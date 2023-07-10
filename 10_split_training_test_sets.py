#!/usr/bin/env python3

import polars as pl
from tqdm import tqdm
import argparse
import logging

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    # 10_split_training_test.py --input-file data/all_gene_annotations.added_incompleteness_and_contamination.tsv --training-families data/training_families.txt --testing-families data/testing_families.txt --output-training data/all_gene_annotations.added_incompleteness_and_contamination.training.tsv --output-testing data/all_gene_annotations.added_incompleteness_and_contamination.testing.tsv
    parent_parser.add_argument('-i','--input-file', required=True)
    parent_parser.add_argument('--training-families', required=True)
    parent_parser.add_argument('--testing-families', required=True)
    parent_parser.add_argument('--output-training', required=True)
    parent_parser.add_argument('--output-testing', required=True)
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read GTDB annotations and calc family
    gtdb = pl.concat([
        pl.read_csv('data/bac120_metadata_r202.tsv', separator="\t"),
        pl.read_csv('data/ar122_metadata_r202.tsv', separator="\t")
    ])
    gtdb = gtdb.filter(pl.col("gtdb_representative") == "t")
    logging.info("Read in {} GTDB reps".format(len(gtdb)))
    # df["phylum"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[1])
    # df["class"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[2])
    # df["order"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[3])
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(4).alias("family"))

    # Read testing and training families
    with open(args.testing_families) as f:
        testing_families = set(f.read().splitlines())
    with open(args.training_families) as f:
        training_families = set(f.read().splitlines())

    # Read input file lazily
    logging.info("Reading input data ..")
    df = pl.read_csv(args.input_file, separator="\t")
    logging.info("Read {} rows".format(len(df)))

    # Join with GTDB
    logging.info("Joining with GTDB ..")
    df2 = df.join(gtdb.select(['accession','family']), on="accession", how="left")

    # Prep train df
    train_df = df2.filter(pl.col("family").is_in(training_families))
    logging.info("Training set has {} rows".format(len(train_df)))
    train_df = train_df.drop(["family"])
    train_df.write_csv(args.output_training, separator="\t")

    # Prep test df
    test_df = df2.filter(pl.col("family").is_in(testing_families))
    logging.info("Testing set has {} rows".format(len(test_df)))
    test_df = test_df.drop(["family"])
    test_df.write_csv(args.output_testing, separator="\t")

    logging.info("Done")