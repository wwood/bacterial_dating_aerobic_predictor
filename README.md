A reproducible analysis pipeline for aerobic predictions


Grab JSON data from bacdive's API. This first step untested contemporaneously but script not changed significantly since 20220201.
```
$ ./1_scrape_bacdive.py --output-json data/bacdive_scrape_20220201.json.parsed.csv
```

Parse the JSON into a TSV. Also removes "obligate" from the oxytolerance annotation.
```
$ ./2_parse_bacdive_scrape.py --input-json data/bacdive_scrape_20220201.json >data/bacdive_scrape_20220201.json.parsed.csv
```

Consolidate annotations from each strain together into a single 2-class system - (anaerobic) or (aerobic, microaerophilic, facultative)
```
$ ./3_consolidate_parsed.py --input-csv data/bacdive_scrape_20220201.json.parsed.csv --output-csv data/bacdive_scrape_20220201.json.parsed.anaerobe_vs_rest_consolidated.csv --class-mapping data/anaerobe_vs_rest.mapping --multi-class-output-csv debug/multi_class.csv
```
