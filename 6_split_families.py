#!/usr/bin/env python3
 
import argparse
import logging
import os
import polars as pl
import random

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    # $ ./6_split_families.py -i data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv --training-families data/training_families.txt --testing-families data/testing_families.txt
    parent_parser.add_argument('-i','--input-csv', required=True)
    parent_parser.add_argument('--training-families', required=True)
    parent_parser.add_argument('--testing-families', required=True)
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    df = pl.concat([
        pl.read_csv('data/bac120_metadata_r202.tsv', sep="\t"),
        pl.read_csv('data/ar122_metadata_r202.tsv', sep="\t")
    ])

    df = df.filter(pl.col("gtdb_representative") == "t")
    logging.info("Read in {} GTDB reps".format(len(df)))
    # df["phylum"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[1])
    # df["class"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[2])
    # df["order"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[3])
    df = df.with_columns(pl.col("gtdb_taxonomy").str.split(';').arr.get(4).alias("family"))
    # df["genus"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[5])
    # df["sp"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[6])

    # Read input tsv
    data = pl.read_csv(args.input_csv, sep="\t")
    d = data.join(df, on="accession", how="left")

    testing_families = set()
    testing_set_size = 0

    all_families = list(set(d["family"].to_list()))
    logging.info("Found {} families from {} data points".format(len(all_families), len(d)))

    # Randomly order families
    random.shuffle(all_families)

    while testing_set_size < len(d) * 0.2:
        family = all_families.pop()
        testing_families.add(family)
        testing_set_size += len(d.filter(pl.col("family") == family))

    logging.info("Found {} testing families, comprising {} data points".format(len(testing_families), testing_set_size))

    with open(args.training_families, "w") as f:
        f.write("\n".join(all_families))
    with open(args.testing_families, "w") as f:
        f.write("\n".join(testing_families))

    logging.info("Finished")