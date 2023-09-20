#!/usr/bin/env python3
 
import argparse
import logging
import requests
from requests.auth import HTTPBasicAuth
import os
import pandas as pd
import time
import bacdive
import json
from tqdm import tqdm

"""
Script to scrape bacterial information from BacDive. 
"""

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--output-json', required=True)
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    
    from private_credentials import USERNAME, PASSWORD
    client = bacdive.BacdiveClient(USERNAME, PASSWORD)

    df = pd.concat([
        pd.read_csv('data/bac120_metadata_r202.tsv', sep="\t"),
        pd.read_csv('data/ar122_metadata_r202.tsv', sep="\t")
    ])
    df = df[df["gtdb_representative"] == "t"]
    df["phylum"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[1])
    df["class"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[2])
    df["order"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[3])
    df["family"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[4])
    df["genus"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[5])
    df["sp"] = df["gtdb_taxonomy"].apply(lambda x: x.split(";")[6])

    df["quality"] = df["checkm_completeness"] - (df["checkm_contamination"] * 5)

    df = df[df["ncbi_taxonomy"] != ""]
    df = df.sort_values(by=["gtdb_type_species_of_genus","quality"], ascending=False)

    fdf = df.groupby("genus", as_index=False).first()

    sps_list = set([tuple(x.split(" ")) for x in fdf["ncbi_taxonomy"].apply(lambda x:x.split(";")[6][3:])])
    logging.info("Found {} type species".format(len(sps_list)))

    # We only scrape on the genus level, so remove dups
    genus_list = set()
    for sp in sps_list:
        if sp != '':
            genus_list.add(sp[0])

    logging.info("Found {} genera to query".format(len(sps_list)))

    with open(args.output_json, "w") as f:
        for i,genus in enumerate(tqdm(genus_list)):
            logging.debug("Looking for annotations of the genus of i={} {}".format(i, genus))

            query = {'taxonomy': (genus)}
            client.search(**query)
            strains = client.retrieve()
            try:
                l = list([s for s in strains])
                print("\t".join([genus, 'spp.', json.dumps(l)]), file=f)
            except KeyError:
                logging.warning("No strains found for {}".format(genus))
