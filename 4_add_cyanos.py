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
    parent_parser.add_argument('--cyano-species', help='TSV with cyanobacteria species and oxygen tolerance. If not specified, just transform input and do not add any cyanos')
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

    # Read in cyano species
    if args.cyano_species:
        cyanos = pl.read_csv(args.cyano_species, separator="\t")
        # Change aerobic to aerobe
        cyanos2 = cyanos.select([
            pl.col('accession'),
            pl.col('oxytolerance').str.replace('aerobic', 'aerobe'),
        ])
        logging.info("Read in {} cyanos".format(len(cyanos2)))

        to_write = pl.concat(
            [
                df.select([
                    pl.col('accession'),
                    pl.col('Oxygen tolerance').alias('oxytolerance'),
                ]),
                cyanos2
            ],
            how='vertical')
    else:
        to_write = df.select([
            pl.col('accession'),
            pl.col('Oxygen tolerance').alias('oxytolerance'),
        ])

    logging.info("Class counts in input: {}".format(df.groupby('Oxygen tolerance').count()))
    logging.info("Class counts in output: {}".format(to_write.groupby('oxytolerance').count()))

    # Write to TSV
    to_write.write_csv(args.output_csv, separator="\t")


