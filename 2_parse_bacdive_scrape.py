#!/usr/bin/env python3
 
import argparse
import logging
import json

def parse_from_species_block(input_json):
    #          "General": {
    #   "BacDive-ID": 151694,
      
    #  "Physiology and metabolism": {
    #   "oxygen tolerance": {
    #     "@ref": 56421,
    #     "oxygen tolerance": "microaerophile"
    #   },

    # Collect all data for each species from the genus
    name_to_bacdive_ids = {}
    name_to_oxygen_tolerances = {}
    for j in input_json:
        name = j["Name and taxonomic classification"]["species"]
        if name not in name_to_bacdive_ids:
            name_to_bacdive_ids[name] = []
            name_to_oxygen_tolerances[name] = set()

        if 'oxygen tolerance' in j['Physiology and metabolism']:
            tolerances = []
            ot = j['Physiology and metabolism']['oxygen tolerance']
            if isinstance(ot, dict):
                if 'oxygen tolerance' in ot: # In some instances the oxygen tolerance is not specied
                    for o in ot['oxygen tolerance'].split(','):
                        tolerances.append(o)
            else:
                for entry in ot:
                    if 'oxygen tolerance' in entry: # In some instances the oxygen tolerance is not specified e.g. Butyricimonas faecihominis
                        for o in entry['oxygen tolerance'].split(','):
                            tolerances.append(o)
            if tolerances != []:
                # Only record BacDive IDs which are associated with data
                name_to_bacdive_ids[name].append(str(j['General']['BacDive-ID']))
                for o in tolerances:
                    name_to_oxygen_tolerances[name].add(o.strip())


    # Print per-species
    for name, bacdive_ids in name_to_bacdive_ids.items():
        oxytolerances2 = set()
        for ot in name_to_oxygen_tolerances[name]:
            if ot == 'obligate aerobe':
                oxytolerances2.add('aerobe')
            elif ot == 'obligate anaerobe':
                oxytolerances2.add('anaerobe')
            # elif ot == 'microaerophile':
            #     oxytolerances2.add('aerobe')
            else:
                oxytolerances2.add(ot)
        print('\t'.join([','.join(bacdive_ids), name, ','.join(oxytolerances2), ','.join(name_to_oxygen_tolerances[name])]))

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-json', required=True)
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    print('\t'.join([
        'BacDive-ID',
        'Species',
        'Oxygen tolerance',
        'Oxygen tolerance raw',
    ]))

    with open(args.input_json) as f:
        for l in f:
            (genus,species,json_string) = l.strip().split('\t')
            if l != '[]\n':
                parse_from_species_block(json.loads(json_string))