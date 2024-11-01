#!/usr/bin/env python3
 
import argparse
import logging
import polars as pl
import re
import random

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-csv', required=True)
    parent_parser.add_argument('--output-csv', required=True)
    parent_parser.add_argument('--respiration-genes', required=True)
    parent_parser.add_argument('--set-as-aerobic', action='store_true')
    parent_parser.add_argument('--set-as-new-class', action='store_true')
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read input TSV
    df = pl.read_csv(args.input_csv, separator="\t")
    logging.info("Read in {} rows".format(len(df)))

    # Read in respiration genes
    respiration_genes = pl.read_csv(args.respiration_genes, separator="\t")
    logging.info("Read in {} respiration species' respiration genes".format(len(respiration_genes)))

    # Add a summary column for respiration genes - "anything with a HCO (A,B and B subfamilies, C) or anything with a cytbd (qOR1-4, OR-C, OR-N)"
    # Get the columns we want to check
    # In [2]: respiration_genes.columns
    # Out[2]: 
    # ['Taxonomy',
    # 'Genome',
    # 'A',
    # 'C',
    # 'B1',
    # 'B2',
    # 'B3',
    # 'B4',
    # 'B5',
    # 'B6',
    # 'B7',
    # 'B8',
    # 'B9',
    # 'qOR1',
    # 'qOR2',
    # 'qOR3',
    # 'qOR4a',
    # 'OR-C1a',
    # 'OR-N1',
    # 'OR-N2',
    # 'OR-N3a',
    # 'OR-N4a',
    # 'OR-N5a']
    # So find the sum of all columns except the first two
    # First remove the taxonomy column
    respiration_genes = respiration_genes.select([pl.exclude('Taxonomy')])
    # Melt the dataframe
    respiration_genes2 = respiration_genes.melt(id_vars=['Genome'], value_name='count')
    respiration_genes3 = respiration_genes2.groupby(['Genome']).agg(pl.sum('count').alias('respiration_gene_count'))
    respiration_genes4 = respiration_genes3.select([
        pl.col('Genome').alias('accession'),
        (pl.col('respiration_gene_count') > 0).alias('has_respiration_genes')
    ])
    logging.info("Summary of all respiration gene annotations: {}".format(respiration_genes4.groupby('has_respiration_genes').count()))

    # Merge the respiration gene summary into the main dataframe
    df2 = df.join(respiration_genes4, on='accession', how='left')
    logging.info("Summary of respiration gene annotations in dataset: {}".format(df2.groupby(['oxytolerance','has_respiration_genes']).count().sort('oxytolerance')))

    # Either set as aerobic or exclude those with respiration genes
    if args.set_as_aerobic:
        df3 = df2.select(
            pl.col('accession'),
            pl.when(
                pl.col('has_respiration_genes') == True).
                then('aerobe').
                otherwise(pl.col('oxytolerance')).alias('oxytolerance')
        )
    elif args.set_as_new_class:
        df3 = df2.select(
            pl.col('accession'),
            pl.col('oxytolerance'),
            pl.when(
                pl.col('has_respiration_genes') == True).
                then('anaerobic_with_respiration_genes').
                otherwise(pl.col('oxytolerance')).alias('oxytolerance0')
        )
        df3 = df3.select(
            pl.col('accession'),
            pl.when(
                pl.col('oxytolerance') == 'aerobe').
                then('aerobe').
                otherwise(pl.col('oxytolerance0')).alias('oxytolerance')
        )
    else:
        df3 = df2.filter(
            pl.any([
                pl.col('oxytolerance') == 'aerobe',
                pl.col('has_respiration_genes') == False]))
    logging.info("Summary of respiration gene annotations in dataset after filtering: {}".format(df3.groupby(['oxytolerance']).count().sort('oxytolerance')))

    # Write out the dataframe
    df3.select(['accession','oxytolerance']).write_csv(args.output_csv, separator="\t")

    # # Why do so many anaerobic genomes contain positive calls?
    # df3 = df2.join(respiration_genes.select([
    #     pl.col('Genome').alias('accession'),
    #     pl.col('*')
    # ]), on='accession', how='left')
    # anaerobes = df3.filter(pl.col('oxytolerance') == 'anaerobe')
    # anaerobes2 = anaerobes.filter(pl.col('has_respiration_genes') == True)
    # anaerobes3 = anaerobes2.melt(id_vars=['accession', 'oxytolerance', 'has_respiration_genes','Genome'], value_name='count')
    # anaerobes4 = anaerobes3.filter(pl.col('count') > 0).groupby('variable').count()
    # anaerobes4.write_csv('debug/respiration_genes_in_anaerobes.csv', separator="\t")

    # # Write out all annotations
    # respiration_genes_again = pl.read_csv(args.respiration_genes, separator="\t")
    # df3.join(
    #     respiration_genes_again.select([
    #         pl.col('Genome').alias('accession'),
    #         pl.col('Taxonomy')
    #     ]), on='accession', how='left').write_csv('debug/respiration_genes_all_annotations.csv', separator="\t")

