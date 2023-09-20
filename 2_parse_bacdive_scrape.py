#!/usr/bin/env python3
 
import argparse
import logging
import json
"""
Parsing BacDive information for temperature range. 
"""

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
    name_to_temperature_range = {}
    for j in input_json:
        name = j["Name and taxonomic classification"]["species"]
        if name not in name_to_bacdive_ids:
            name_to_bacdive_ids[name] = []
            name_to_temperature_range[name] = set()

<<<<<<< Updated upstream
        if 'culture temp' in j['Culture and growth conditions']:
            tolerances = []
            tempt = j['Culture and growth conditions']['culture temp']
            if isinstance(tempt, dict):
                if 'range' in tempt: # In some instances the range is not specied
                    for temp in tempt['range'].split(','):
                        tolerances.append(temp)
            else:
                for entry in tempt:
                    if 'range' in entry: # In some instances the Temperature range is not specified
                        for temp in entry['range'].split(','):
                            tolerances.append(temp)
            if tolerances != []:
                # Only record BacDive IDs which are associated with data
                name_to_bacdive_ids[name].append(str(j['General']['BacDive-ID']))
                for temp in tolerances:
                    name_to_temperature_range[name].add(temp.strip())
=======
        if 'temperature range' in j['Culture and growth conditions']:
            tolerances = []
            ot = j['Culture and growth conditions']['temperature range']
            if isinstance(ot, dict):
                if 'temperature range' in ot: # In some instances the temperature range is not specied
                    for o in ot['temperature range'].split(','):
                        tolerances.append(o)
            else:
                for entry in ot:
                    if 'temperature range' in entry: # In some instances the temperature range is not specified e.g. Butyricimonas faecihominis
                        for o in entry['temperature range'].split(','):
                            tolerances.append(o)
            if tolerances != []:
                # Only record BacDive IDs which are associated with data
                name_to_bacdive_ids[name].append(str(j['General']['BacDive-ID']))
                for o in tolerances:
                    name_to_temperature_range[name].add(o.strip())
>>>>>>> Stashed changes


    # Print per-species
    for name, bacdive_ids in name_to_bacdive_ids.items():
<<<<<<< Updated upstream
        tmptolerances = set()
        for ot in name_to_temperature_range[name]:
                tmptolerances.add(ot)
        print('\t'.join([','.join(bacdive_ids), name, ','.join(tmptolerances), ','.join(name_to_temperature_range[name])]))
=======
        temperature_range = set()
        for ot in name_to_temperature_range[name]:
            if ot == 'mesophilic':
                temperature_range.add('mesophilic')
            elif ot == 'thermophilic':
                temperature_range.add('thermophilic')
            elif ot == 'psychrophilic':
                temperature_range.add('psychrophilic')
            elif ot == 'hyperthermophilic':
                temperature_range.add('hyperthermophilic')
            else:
                temperature_range.add(ot)
        print('\t'.join([','.join(bacdive_ids), name, ','.join(temperature_range), ','.join(name_to_temperature_range[name])]))
>>>>>>> Stashed changes

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
<<<<<<< Updated upstream
        'Temperature range',
        'Temperature range raw',
=======
        'Temperature Range',
        'Temperature Range raw',
>>>>>>> Stashed changes
    ]))

    with open(args.input_json) as f:
        for l in f:
            (genus,species,json_string) = l.strip().split('\t')
            if l != '[]\n':
                parse_from_species_block(json.loads(json_string))