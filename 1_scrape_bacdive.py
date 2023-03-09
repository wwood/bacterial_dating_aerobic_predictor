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
        pd.read_csv('~/m/db/gtdb/gtdb_release202/bac120_metadata_r202.tsv', sep="\t"),
        pd.read_csv('~/m/db/gtdb/gtdb_release202/ar122_metadata_r202.tsv', sep="\t")
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

    # sps_list = [("Alloiococcus", "otitis")] # 5 strains found, some have oxytolerance
    # sps_list = [('Acetoanaerobium', 'noterae')] # no strains found

    logging.info("Found {} genera to query".format(len(sps_list)))

    with open(args.output_json, "w") as f:
        for i,entry in enumerate(sps_list):
            logging.info("Looking for annotations of the genus of i={} {}".format(i, entry))
            
            # if i < 2057:
            #     continue
            
            if len(entry) != 2:
                continue
            if i % 100 == 0:
                print(str(i))
            genus, species_epithet = entry

            # query = {'taxonomy': (genus, species_epithet)}
            query = {'taxonomy': (genus)}
            client.search(**query)
            strains = client.retrieve()
            try:
                l = list([s for s in strains])
                print("\t".join([genus, species_epithet, json.dumps(l)]), file=f)
            except KeyError:
                logging.warning("No strains found for {}".format(entry))


    # for description in l:
    #     bacdive_id = description['General']['BacDive-ID']
    #     oxygen_tolerance = ''

    #     try:
    #         oxygen_tolerance = description['Physiology and metabolism']['oxygen tolerance']['oxygen tolerance']
    #     except KeyError:
    #         pass

    #     if oxygen_tolerance != '': import IPython; IPython.embed()
    #     print('\t'.join([str(i), genus, species_epithet, str(bacdive_id), oxygen_tolerance]))

    # [s['General']['BacDive-ID'] for s in l]
    # import IPython; IPython.embed()
    
    # response = requests.get('https://api.bacdive.dsmz.de/taxon/%s/%s' % (genus,species_epithet), auth=credentials)
    # time.sleep(0.3)
    # if response.status_code == 200:
    #     results = response.json()
    #     sp2links[entry] = results
        
    # else:
    #     # raise
    #     continue      
# print("Done!")




# sp2raw_info = dict()
# i = 0
# for sp, links in sp2links.items():
    
#     i+=1
#     if i % 50 == 0:
#         print(str(i))
        
#     # if i < 864:
#     #     continue
    
#     url = links["results"][0]["url"]
#     response = requests.get(url, headers=headers, auth=credentials)
#     time.sleep(0.3)
#     if response.status_code == 200:
#         results = response.json()
#         sp2raw_info[sp] = results
#     else:
#         continue
# print("Done!")