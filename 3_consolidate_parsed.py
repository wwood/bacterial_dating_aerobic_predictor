#!/usr/bin/env python3
 
import argparse
import logging
import pandas as pd
import re
import random

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-csv', required=True)
    parent_parser.add_argument('--output-csv', required=True)
    parent_parser.add_argument('--class-mapping-file', required=True)
    parent_parser.add_argument('--multi-class-output-csv', help='output cases where 2+ classes are predicted')
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read in class mapping
    class_mapping_csv = pd.read_csv(args.class_mapping_file)
    class_mapping = {}
    for i, row in class_mapping_csv.iterrows():
        class_mapping[row['member']] = row['class']
    logging.info("Read in {} class mappings: {}".format(len(class_mapping), class_mapping))

    # Read GTDB
    df = pd.concat([
        pd.read_csv('data/bac120_metadata_r202.tsv', sep="\t"),
        pd.read_csv('data/ar122_metadata_r202.tsv', sep="\t")
    ])
    df = df[df["gtdb_representative"] == "t"]
    logging.info("Read in {} GTDB reps".format(len(df)))
    df["phylum"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[1])
    df["class"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[2])
    df["order"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[3])
    df["family"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[4])
    df["genus"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[5])
    df["sp"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[6])

    df["quality"] = df["checkm_completeness"] - (df["checkm_contamination"] * 5)

    df = df[df["ncbi_taxonomy"] != ""]

    df['ncbi_taxonomy2'] = [x for x in df["ncbi_taxonomy"].apply(lambda x:x.split(";")[6][3:])]
    # Choose original type species, if there are more than one.
    r = re.compile('.* .*_[A-Z]+$')
    df['is_original_species'] = [r.match(x) is None for x in df['sp']]
    df_original = df[df['is_original_species'] == True]
    logging.info("Found {} reps of original species".format(len(df_original)))

    # Get best quality rep for each ncbi_taxonomy
    df = df.sort_values(by=["gtdb_type_species_of_genus","quality"], ascending=False)
    df_best = df_original.groupby("ncbi_taxonomy2", as_index=False).first()
    logging.info("Found {} reps of best quality".format(len(df_best)))
    logging.debug("e.g. {}".format(df_best.iloc[0]))
    logging.debug("i.e. {}".format(df_best[['ncbi_taxonomy2','ncbi_taxonomy','gtdb_taxonomy','genus','sp','accession']].iloc[0]))

    bacdive = pd.read_csv(args.input_csv, sep="\t")
    bacdive = bacdive[bacdive["Oxygen tolerance raw"] != ""]
    logging.info("Read in {} rows of bacdive annotation".format(bacdive.shape[0]))
    logging.debug("e.g. {}".format(bacdive.iloc[0]))

    # Remove duplicates from BacDive.
    bacdive = bacdive[[not pd.isna(a) for a in bacdive["Oxygen tolerance raw"]]]
    bacdive.drop_duplicates(inplace=True)
    logging.info("After removing duplicates, have {} bacdive entries".format(bacdive.shape[0]))

    merged = pd.merge(bacdive, df_best, right_on="ncbi_taxonomy2", left_on="Species")
    # print(merged.columns)
    m2 = merged[['Species','gtdb_taxonomy','genus','sp','accession','gtdb_type_species_of_genus','Oxygen tolerance raw']]
    m2.drop_duplicates(inplace=True)

    # Remove some problem cases by hand. These are bacdive duplication data errors somehow? These two are both microaerophile and anaerobe, so not useful for us.
    m3 = m2[~(m2['Species'].isin(['Actinotignum schaalii','Microvirgula aerodenitrificans']))]

    # Ensure there are no duplicate species in m3
    if len(m3.Species.unique()) != len(m3.Species):
        logging.error("Found duplicate species in m3")
        raise Exception("Found duplicate species in m3")

    # Map classes
    def map_classes_def(members):
        classes = set()
        for m in members.split(","):
            # Check there are no class mappings which are not in the mapping
            if m not in class_mapping:
                logging.error("Found class {} not in class mapping".format(c))
                raise Exception("Found class {} not in class mapping".format(c))
            elif class_mapping[m] == 'exclude':
                # Entirely exclude this species, not just this annotation
                logging.debug("Found exclude class for members: {}".format(members))
                return None
            else:
                classes.add(class_mapping[m])
        if len(classes) > 1:
            logging.debug("Found multiple classes for members: {}".format(members))
            return None
        else:
            return list(classes)[0]
    m3["Oxygen tolerance"] = [map_classes_def(c) for c in m3["Oxygen tolerance raw"]]

    # Remove any rows with multiple classes / or that are excluded
    m4 = m3[m3["Oxygen tolerance"].notnull()]
    num_multiple_classes = len(m3) - len(m4)
    logging.info("Found {} rows with multiple classes or were excluded, and {} rows with a single class".format(num_multiple_classes, len(m4)))
    if args.multi_class_output_csv is not None:
        # Write out the multiple class rows
        m3[m3["Oxygen tolerance"].isnull()].to_csv(args.multi_class_output_csv, sep="\t", index=False)

    # Only accept 3 species per GTDB genus. We want type species of the genus, plus at most 2 more
    m4['random_id'] = [random.randint(0, len(m4)) for i in range(len(m4))]
    m4.sort_values(by=["gtdb_type_species_of_genus","random_id"], inplace=True, ascending=False)
    genus_derep = m4.groupby('genus').head(3)

    genus_derep.to_csv(args.output_csv, sep="\t", index=False)
