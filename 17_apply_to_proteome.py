#!/usr/bin/env python3

import polars as pl
import pandas as pd
from joblib import load
import argparse
import logging
import os
import sys
from io import StringIO

from bird_tool_utils import in_tempdir
import extern
from tqdm import tqdm


# flatten a column that has lists as entries
def flatten_columns(df, cols):
    """Flattens multiple columns in a data frame, cannot specify all columns!"""
    flattened_cols = {}
    for col in cols:
        flattened_cols[col] = pd.DataFrame([(index, value) for (index, values) in df[col].items() for value in values],
                                           columns=['index', col]).set_index('index')
    flattened_df = df.drop(cols, axis=1)
    for col in cols:
        flattened_df = flattened_df.join(flattened_cols[col])
    return flattened_df


# Returns a data frame with 'query' and 'eggNOG_OGs' columns for gene and the predictive cog that provides
def read_annotations(annotations_path, predictive_cogs_series):
    #
    # equivalent of the following R:
    # annotations = train_accessions[,fread(cmd=paste('sed s=^.query=query= eggnog1_results/',accession,".faa.emapper.annotations  |grep -v '^\\#'",sep='')),by=accession]
    # nrow(annotations)
    #
    a1 = pd.read_csv(annotations_path, sep="\t", comment='#', header=None)
    # a1.columns[0] = 'query'
    # a1[:,0]
    expected_columns = [
        "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
        "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
        "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
    ]
    if (len(a1.columns) != len(expected_columns)):
        raise Exception("Unexpected formation (number of columns) in eggnog-mapper .annotations file")
    a1.columns = expected_columns

    a2 = a1.loc[:, list(['query', 'eggNOG_OGs'])]
    a2.loc[:, 'eggNOG_OGs'] = a2['eggNOG_OGs'].str.split(",")

    a3 = flatten_columns(a2, ['eggNOG_OGs'])
    a4 = pd.merge(a3, predictive_cogs_series, on='eggNOG_OGs')

    # dcast - do all at once later
    # a5 = pd.pivot_table(data=a4, columns=predictive_cogs_series, aggfunc=sum, index=None)

    return a4


# Seems to require pandas 1.3.3 (or at least >1.2.2). Dummy first row should be removed first.
def read_multiple_annotations(to_read, predictive_cogs_series):  # read in a number of annotations
    # A dummy is needed so the right otherwise the pivot_table at the end fails with a KeyError
    # when predictive ones aren't found in any genomes.
    collected = pd.DataFrame({"accession": 'dummy', 'eggNOG_OGs': predictive_cogs_series, 'query': 'dummy'})

    for ann in tqdm(to_read, desc="reading annotation files"):
        cogs = read_annotations(ann.strip(), predictive_cogs_series)
        cogs['accession'] = ann
        # if collected is None:
        #     collected = cogs
        # else:
        collected = pd.concat([collected, cogs])
    return collected


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--working-directory', help='working directory', default=None)
    parent_parser.add_argument('--protein-fasta', help='protein fasta file', required=True)

    parent_parser.add_argument('--annotation-table', help='Pre-calculated counts of COG / KO families in one or more genomes, in TSV format')
    # add missing annotations
    parent_parser.add_argument('--add-missing-annotations', help='when columns are missing, add them as 0 [default: croak]', action="store_true")
    
    parent_parser.add_argument('--threads', help='number of threads to use', default=1, type=int)
    
    eggnog_parser = parent_parser.add_mutually_exclusive_group()
    eggnog_parser.add_argument('--eggnog-data-dir',
                               help='eggnog data directory e.g. ~/m/db/eggnog-mapper/2.1.3')
    eggnog_parser.add_argument('--eggnog-annotation-file', help='eggnog .annotation file')

    kofam_parser = parent_parser.add_mutually_exclusive_group()
    kofam_parser.add_argument('--kofam-hmm-path',
                               help='path to kofam hmm e.g. ~/m/db/kofam/2022-01-30/profiles.hmm')
    kofam_parser.add_argument('--kofam-tsv-file', help='hmmsearch output file')

    # whitelist default to 'data/COGs_with_more_than_50_percent_remaining.tsv'
    parent_parser.add_argument('--whitelist', help='whitelist of COGs to use', default='/work/microbiome/bacterial_dating/19_coherent_prediction_workflow/data/COGs_with_more_than_50_percent_remaining.tsv')

    parent_parser.add_argument('--modal-keggs', help='modal keggs', default='/work/microbiome/bacterial_dating/19_coherent_prediction_workflow/data/ModalKEGGs.tsv')
    # ./12_apply_model.py --model the.model -x the.x --training-data-header <(head the.training-data-header)
    parent_parser.add_argument('--models', nargs='+', help='models to use')

    parent_parser.add_argument('--training-data-header', help='header of training data', required=True)
    
    parent_parser.add_argument('--output-predictions', help='output predictions')
    parent_parser.add_argument('--output-annotations', help='output annotations, and exit, do not apply models')
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

    # Read training data header
    eg_data = pl.read_csv(args.training_data_header, separator="\t", has_header=True)
    header = eg_data.select(pl.exclude(['accession', 'false_negative_rate', 'false_positive_rate'])).columns
    # Blacklist these as they aren't in the current ancestral file, not sure why
    header = list([h for h in header if h not in ['COG0411', 'COG0459', 'COG0564', 'COG1344', 'COG4177']])

    if not args.annotation_table:
        if not args.eggnog_data_dir and not args.eggnog_annotation_file:
            raise Exception("Either --annotations, --eggnog-data-dir or --eggnog-annotation-file must be specified")
        if not args.kofam_hmm_path and not args.kofam_tsv_file:
            raise Exception("Either --annotations, --kofam-hmm-path or --kofam-tsv-file must be specified")
    
    if not args.output_predictions and not args.output_annotations:
        raise Exception("Either --output-predictions or --output-annotations must be specified")

    # Gather abs paths so chdir to working directory doesn't break things
    working_directory = None
    if args.working_directory:
        working_directory = os.path.abspath(args.working_directory)
    protein_fasta = os.path.abspath(args.protein_fasta)
    if args.eggnog_data_dir:
        eggnog_data_dir = os.path.abspath(args.eggnog_data_dir)
    if args.eggnog_annotation_file:
        eggnog_annotation_file = os.path.abspath(args.eggnog_annotation_file)
    if args.kofam_hmm_path:
        kofam_hmm_path = os.path.abspath(args.kofam_hmm_path)
    if args.kofam_tsv_file:
        kofam_tsv_file = os.path.abspath(args.kofam_tsv_file)

    modal_keggs = os.path.abspath(args.modal_keggs)

    if args.annotation_table:
        if args.output_annotations:
            logging.error("Cannot specify --annotation-table and --output-annotations")
            sys.exit(1)
        logging.info("Reading annotation table ..")
        raw = pl.read_csv(args.annotation_table, separator='\t')
        if args.add_missing_annotations:
            added_column_count = 0
            for h in header:
                if h not in raw.columns:
                    raw = raw.with_columns(pl.lit(0).alias(h))
                    added_column_count = 0
            if added_column_count > 0:
                logging.info("Added {} extra 0 count columns to annotation table".format(added_column_count))
        d5 = raw.to_pandas()[header]
        accessions = list(raw['accession'])
    else:
        logging.info("Reading or calculating annotations not from a table ..")
        with in_tempdir():
            logging.info("Created temp working directory {}".format(os.getcwd()))
            if args.working_directory:
                logging.info("Changing working directory to {}".format(working_directory))
                if not os.path.exists(working_directory):
                    logging.info("Creating working directory {}".format(working_directory))
                    os.makedirs(working_directory)
                os.chdir(working_directory)

            logging.info("Reading in COG whitelist")
            whitelist1 = open(args.whitelist).read().splitlines()
            logging.info("Found {} COGs in the whitelist".format(len(whitelist1)))

            # Run eggnog-mapper if needed
            predictive_cogs_series = pd.Series(["{}@1|root".format(c) for c in whitelist1], name='eggNOG_OGs')
            if args.eggnog_annotation_file:
                logging.info("Reading EGNOG annotations from supplied file")
                cog_annotations = read_annotations(eggnog_annotation_file, predictive_cogs_series)
            else:
                logging.info("Running eggnog-mapper")
                eggnog_output = 'eggnog_output'
                extern.run(
                    f'EGGNOG_DATA_DIR={eggnog_data_dir} emapper.py -m diamond -i {protein_fasta} --target_orthologs one2one --query_cover 50.0 --evalue 0.0000001 --cpu {args.threads} -o {eggnog_output}'
                )
                cog_annotations = read_multiple_annotations(
                    [eggnog_output+'.emapper.annotations'],
                    predictive_cogs_series)

            # Run hmmsearch in required
            if args.kofam_tsv_file:
                hmmsearch_tblout = kofam_tsv_file
            else:
                # $ ls proteomes/*faa |parallel hmmsearch --tblout kegg_annotations/annotations/{/}.hmmsearch_tblout.csv -o /dev/null --notextw --cpu 1 ~/m/db/kofam/2022-01-30/profiles.hmm {}
                logging.info("Running hmmsearch")
                hmmsearch_tblout = 'hmmsearch_tblout.csv'
                extern.run(
                    f'hmmsearch --tblout {hmmsearch_tblout} -o /dev/null --notextw --cpu {args.threads} {kofam_hmm_path} {protein_fasta}'
                )

            # and then we take the best hit for each protein
            # $ ls kegg_annotations/annotations/*_tblout.csv |parallel -j10 ./7_process_kegg_hmmsearch.py --input-tblout {} --output-csv kegg_annotations/best_hits/{/}.best_hits.csv &>7_process_kegg_hmmsearch.py.log
            hits = extern.run("sed 's/  */\t/g' {} |cut -f1-6 |grep -v '^#'".format(hmmsearch_tblout))
            hmmsearch_best_hits_df = pd.read_csv(StringIO(hits), sep='\t', header=None)
            hmmsearch_best_hits_df = hmmsearch_best_hits_df.loc[:, [0, 2, 4, 5]]
            logging.info("Read in {} annotations".format(len(hmmsearch_best_hits_df)))

            hmmsearch_best_hits_df.columns = ['query', 'target', 'evalue', 'score']
            hmmsearch_best_hits_df2 = hmmsearch_best_hits_df.groupby('query', as_index=False).apply(
                lambda x: x.nsmallest(1, 'evalue')).reset_index().loc[:, ['query', 'target', 'evalue', 'score']]
            logging.info("Found {} unique HMM annotations".format(len(hmmsearch_best_hits_df2)))

            # Process eggNOG annotations
            cog_annotations.loc[:, 'eggNOG_OGs'] = [x.replace('@1|root', '') for x in cog_annotations['eggNOG_OGs']]
            # Remove rows with dummy accession
            # cog_annotations = cog_annotations[cog_annotations['accession'] != 'dummy']
            # fixed_accessions = []
            # r = re.compile(r'.*/(.*)_protein.faa.gz.emapper.annotations')
            # for acc in cog_annotations['accession']:
            #     m = r.match(acc)
            #     if m:
            #         fixed_accessions.append(m.group(1))
            #     else:
            #         raise Exception("Regex {} failed to match accession {}".format(r, acc))
            # cog_annotations.loc[:, 'accession'] = fixed_accessions
            logging.info("Read in {} cog_annotations".format(len(cog_annotations)))

            # Read in modal KEGGs for each COG grouping
            kegg_pairings = pd.read_csv(modal_keggs, sep="\t", header=None, names=['cog', 'kegg'])

            kegg_cog = pd.merge(hmmsearch_best_hits_df, cog_annotations, on='query')
            kegg_cog.loc[:, 'cog=kegg'] = kegg_cog['eggNOG_OGs'] + '=' + kegg_cog['target']
            kegg_pairings.loc[:, 'cog=kegg'] = kegg_pairings['cog'] + '=' + kegg_pairings['kegg']

            medianed = pd.merge(kegg_cog.loc[:, ['cog=kegg']], kegg_pairings, on='cog=kegg')
            medianed.columns = ['cog=kegg', 'cog', 'kegg']
            medianed['accession'] = 'query_genome'

            wide_table = pd.pivot_table(medianed.loc[:, ['accession','cog']],
                                        index='accession',
                                        aggfunc=len,
                                        columns='cog',
                                        values='cog',
                                        fill_value=0)

            whitelist = ['accession'] + whitelist1
            # do intersect because some COGs are in whitelist but not in any genomes
            wide_table2 = wide_table.loc[:, wide_table.columns.intersection(whitelist)]

            # # Read in data
            # # d = pl.read_csv('TableAncestralRoot1.tsv',sep="\t")
            # d = pd.read_csv(args.x, sep="\t")
            # logging.info("Read in input data of shape {}".format(d.shape))

            # # Collapse counts of each COG subfamily
            # d2 = d
            # d2['COG'] = d2['COG'].str.split('_').str[0]
            # d3 = d2.groupby('COG').sum()
            # d4 = d3.transpose()

            wide_table2_melt_plus = pd.concat([wide_table2.melt(), pd.DataFrame({'cog': header, 'value': 0})])
            wide_table2_melt_plus2 = pl.DataFrame(wide_table2_melt_plus).groupby('cog').sum()
            wide_table2_melt_plus2 = wide_table2_melt_plus2.with_columns(pl.lit('query_genome').alias('accession'))
            wide_table2_melt_plus_pivot = wide_table2_melt_plus2.pivot(index='accession', columns='cog', values='value', aggregate_function='sum')

        # Reorder columns to be the same as the training dataset
        d5 = wide_table2_melt_plus_pivot.to_pandas()[header]

        if args.output_annotations:
            d6 = d5
            d6['accession'] = args.protein_fasta
            d6.to_csv(args.output_annotations, sep="\t", index=False)
            logging.info("Wrote annotation table to {}".format(args.output_annotations))

    all_results = []
    for model_path in args.models:
        logging.info("Loading model {}".format(model_path))
        model = load(model_path)
        logging.info("Loaded model {}".format(model_path))

        doing_perceptron = 'Perceptron' in model_path

        preds = model.predict(d5)
        # if doing_perceptron:
        #     probas = pl.lit(-1.0)
        # else:
        #     probas = model.predict_proba(d5)[:,1]

        results = pl.DataFrame({
            'node': accessions if args.annotation_table else args.protein_fasta,
            'preds': preds})
        results = results.select(
            pl.col('node'),
            pl.col('preds').alias('prediction').cast(pl.Int64),
            # pl.col('proba').alias('probability').cast(pl.Float64),
            pl.lit(model_path).alias('model'))
        if doing_perceptron:
            results = results.with_columns(pl.lit(-1.0).alias('probability').cast(pl.Float64))
        else:
            results = results.with_columns(pl.lit(model.predict_proba(d5)[:, 1]).alias('probability').cast(pl.Float64))
        all_results.append(results)

    pl.concat(all_results).write_csv(args.output_predictions, separator="\t")
    logging.info("Wrote predictions to {}".format(args.output_predictions))
