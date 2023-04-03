#!/usr/bin/env python3

import pandas as pd
import polars as pl
import random
import os, re
from tqdm import tqdm
import argparse
import logging

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
    expected_columns = ["query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"]
    if(len(a1.columns) != len(expected_columns)):
        raise Exception("Unexpected formation (number of columns) in eggnog-mapper .annotations file")
    a1.columns = expected_columns

    a2 = a1.loc[:,list(['query','eggNOG_OGs'])]
    a2.loc[:,'eggNOG_OGs'] = a2['eggNOG_OGs'].str.split(",")

    a3 = flatten_columns(a2, ['eggNOG_OGs'])
    a4 = pd.merge(a3, predictive_cogs_series, on='eggNOG_OGs')

    # dcast - do all at once later
    # a5 = pd.pivot_table(data=a4, columns=predictive_cogs_series, aggfunc=sum, index=None)

    return a4

# Seems to require pandas 1.3.3 (or at least >1.2.2). Dummy first row should be removed first.
def read_multiple_annotations(to_read, predictive_cogs_series): # read in a number of annotations
    # A dummy is needed so the right otherwise the pivot_table at the end fails with a KeyError
    # when predictive ones aren't found in any genomes.
    collected = pd.DataFrame({
        "accession": 'dummy',
        'eggNOG_OGs': predictive_cogs_series,
        'query': 'dummy'
    })

    for ann in tqdm(to_read, desc="reading annotation files"):
        cogs = read_annotations(ann.strip(), predictive_cogs_series)
        cogs['accession'] = ann
        # if collected is None:
        #     collected = cogs
        # else:
        collected = pd.concat([collected, cogs])
    return collected

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    # $ ./6_split_families.py -i data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv --training-families data/training_families.txt --testing-families data/testing_families.txt
    parent_parser.add_argument('-i','--input-csvs', nargs='+', required=True)
    parent_parser.add_argument('-o','--output-csv', required=True)
    parent_parser.add_argument('--eggnog-mapper-annotations-key', required=True)
    parent_parser.add_argument('--kegg-best-hits-annotations-key', required=True)
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read all accessions from all input CSVs
    all_annotations = pl.concat([pl.read_csv(i, separator="\t") for i in args.input_csvs])
    all_accessions = all_annotations['accession'].unique()
    logging.info("Found {} accessions in the input CSVs".format(len(all_accessions)))

    # all_accessions = all_accessions[:3] # debug

    # # Read training data annotations
    # training_data_annotations = pd.read_csv('../12_expanded_set_and_cyanos/train_validate_test_sets3.csv', sep="\t")
    # # training_data_annotations = training_data_annotations[~(training_data_annotations['set'].isin(['test','validation']))]
    # training_data_annotations['false_negative_rate'] = 0
    # training_data_annotations['false_positive_rate'] = 0

    whitelist1 = open('data/COGs_with_more_than_50_percent_remaining.tsv').read().splitlines()
    logging.info("Found {} COGs in the whitelist".format(len(whitelist1)))

    # Read eggnog file map
    eggnog_file_map_df = pl.read_csv(args.eggnog_mapper_annotations_key, separator="\t", has_header=False)
    eggnog_file_map = {}
    for acc, path in zip(eggnog_file_map_df['column_1'], eggnog_file_map_df['column_2']):
        eggnog_file_map[acc] = path

    # Read in COG annotations
    cog_annotations = read_multiple_annotations(
        [eggnog_file_map[x] for x in all_accessions],
        pd.Series(["{}@1|root".format(c) for c in whitelist1], name='eggNOG_OGs'))
    cog_annotations.loc[:,'eggNOG_OGs'] = [x.replace('@1|root','') for x in cog_annotations['eggNOG_OGs']]
    # Remove rows with dummy accession
    cog_annotations = cog_annotations[cog_annotations['accession'] != 'dummy']
    fixed_accessions = []
    r = re.compile(r'.*/(.*)_protein.faa.gz.emapper.annotations')
    for acc in cog_annotations['accession']:
        m = r.match(acc)
        if m:
            fixed_accessions.append(m.group(1))
        else:
            raise Exception("Regex {} failed to match accession {}".format(r, acc))
    cog_annotations.loc[:,'accession'] = fixed_accessions
    logging.info("Read in {} cog_annotations".format(len(cog_annotations)))

    # Read in KEGG file map
    kegg_file_map_df = pl.read_csv(args.kegg_best_hits_annotations_key, separator="\t", has_header=False)
    kegg_file_map = {}
    for acc, path in zip(kegg_file_map_df['column_1'], kegg_file_map_df['column_2']):
        kegg_file_map[acc] = path

    # Read in KEGG annotations
    kegg_annotations = None
    for accession in tqdm(all_accessions, desc="reading KEGG annotation files"):
        df = pd.read_csv(kegg_file_map[accession],sep=',')
        df['accession'] = accession
        if kegg_annotations is None:
            kegg_annotations = df
        else:
            kegg_annotations = pd.concat([kegg_annotations, df])

    kegg_annotations['accession_protein'] = kegg_annotations['accession'] + '~' + kegg_annotations['query']
    cog_annotations['accession_protein'] = cog_annotations['accession'] + '~' + cog_annotations['query']

    # Read in modal KEGGs for each COG grouping
    kegg_pairings = pd.read_csv('data/ModalKEGGs.tsv', sep="\t", header=None, names=['cog','kegg'])

    kegg_cog = pd.merge(kegg_annotations, cog_annotations, on='accession_protein')
    kegg_cog.loc[:,'cog=kegg'] = kegg_cog['eggNOG_OGs'] + '=' + kegg_cog['target']
    kegg_pairings.loc[:,'cog=kegg'] = kegg_pairings['cog'] + '=' + kegg_pairings['kegg']

    medianed = pd.merge(kegg_cog.loc[:,['accession_x','cog=kegg']], kegg_pairings, on='cog=kegg')
    medianed.columns = ['accession','cog=kegg','cog','kegg']

    wide_table = pd.pivot_table(medianed.loc[:,['accession','cog']], index='accession', aggfunc=len, columns='cog', values='cog', fill_value=0)

    whitelist = ['accession'] + whitelist1
    # do intersect because some COGs are in whitelist but not in any genomes
    wide_table2 = wide_table.loc[:, wide_table.columns.intersection(whitelist)]

    wide_table2.to_csv(args.output_csv, sep="\t", index=True)
