#!/usr/bin/env python3

import polars as pl
import logging
import argparse

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input', help='from 9_expand_incompletenss_and_contamination4.py', required=True)
    parent_parser.add_argument('--output', help='output predictions', required=True)
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

    # Read in data
    # d = pl.read_csv('TableAncestralRoot1.tsv',sep="\t")
    logging.info("Reading input ..")
    d = pl.read_csv(args.input, separator="\t", has_header=True)
    logging.info("Read in input data of shape {}".format(d.shape))

    # Rename accession to inlcude false_negative_rate and false_positive_rate
    d2 = d.with_columns([
        pl.concat_str(d["accession"], d["false_negative_rate"], d['false_positive_rate'], separator="=").alias('accession2')
    ])

    d3 = d2.select(pl.exclude([
        'accession','false_negative_rate','false_positive_rate'
    ]))

    logging.info("Transposing ..")
    d4 = d3.select([
        pl.exclude('accession2')
    ]).transpose()

    cogs = d3.select(pl.exclude('accession2')).columns
    d5 = d4.with_columns(
        # df.with_column(pl.Series(name="col_list", values=my_list))
        pl.Series('COG', cogs)
    )

    d6 = d5.select([pl.col('COG'), pl.exclude('COG')])
    d6.columns = ['COG'] + list(d3['accession2'])

    d6.write_csv(args.output, separator="\t", has_header=True)
