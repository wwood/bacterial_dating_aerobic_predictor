#!/usr/bin/env python3
 
import argparse
import logging
import polars as pl
import re
import random

"""
Adding cyanobacteria into output. 
"""
def replace_string(match):
    return replacement_mapping[match.group(0)]

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

        replacement_mapping = {
            'thermophilic': 'psychrophilic',
            'mesophilic': 'psychrophilic', 
            'hyperthermophilic': 'psychrophilic'
        }
        
        # pattern = '|'.join(re.escape(key) for key in replacement_mapping.keys())

        # Convert the DataFrame to Pandas to perform the replacement
        cyanos_pd = cyanos.to_pandas()
        cyanos_pd['temperaturerange'] = cyanos_pd['temperaturerange'].replace(replacement_mapping, regex=True)
        
        # Convert the modified Pandas DataFrame back to Polars DataFrame
        cyanos2 = pl.DataFrame(cyanos_pd)
        logging.info("Read in {} cyanos".format(len(cyanos2)))

        to_write = pl.concat(
            [
                df.select([
                    pl.col('accession'),
                    pl.col('Temperature range').alias('temperaturerange'),
                ]),
                cyanos2
            ],
            how='vertical')
    else:
        to_write = df.select([
            pl.col('accession'),
            pl.col('Temperature range').alias('temperaturerange'),
        ])

    logging.info("Class counts in input: {}".format(df.groupby('Temperature range').count()))
    logging.info("Class counts in output: {}".format(to_write.groupby('temperaturerange').count()))

    # Write to TSV
    to_write.write_csv(args.output_csv, separator="\t")


