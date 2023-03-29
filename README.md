A reproducible analysis pipeline for aerobic predictions

First setup conda env
```
$ conda env create -f env.yml -p env
```

Then activate the env
```
$ conda activate ./env
```

Grab JSON data from bacdive's API. You'll need to create a `private_credentials.py` file in the root directory with the following contents, after registering with bacdive:
```
USERNAME='xxx@xxx'
PASSWORD='xxxx'
```

Then run the scrape
```
$ ./1_scrape_bacdive.py --output-json data/bacdive_scrape_20230315.json
```

Parse the JSON into a TSV. Also removes "obligate" from the oxytolerance annotation.
```
$ ./2_parse_bacdive_scrape.py --input-json data/bacdive_scrape_20230315.json >data/bacdive_scrape_20230315.json.parsed.csv
```

Consolidate annotations from each strain together into a single 2-class system - (anaerobic) or (aerobic, microaerophilic, facultative, ..). NOTE: This step (and the one below) are not reproducible because they involve sorting by randomly generated numbers.
```
$ ./3_consolidate_parsed.py --input-csv data/bacdive_scrape_20230315.json.parsed.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.csv --class-mapping data/anaerobe_vs_rest.mapping
```

Consolidate into a second 2-class system - (anaerobic) or (aerobic) 
```
$ ./3_consolidate_parsed.py --input-csv data/bacdive_scrape_20230315.json.parsed.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.csv --class-mapping data/aerobe_vs_anaerobe.mapping
```

Add cyano annotations, which are not in BacDive
```
$ ./4_add_cyanos.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv --cyano-species data/manual_cyano_annotations3.tsv
$ ./4_add_cyanos.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv --cyano-species data/manual_cyano_annotations3.tsv
```

Integrate manually curated respiration gene info

```
$ ./5_apply_respiration_genes.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.apply_respiration_gene_exclusion.csv --respiration-genes data/aerobic_repiration_by_species_in_gtdb_r202_final.txt

$ ./5_apply_respiration_genes.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.apply_respiration_gene_set_aerobic.csv --respiration-genes data/aerobic_repiration_by_species_in_gtdb_r202_final.txt --set-as-aerobic

$ ./5_apply_respiration_genes.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.apply_respiration_gene_exclusion.csv --respiration-genes data/aerobic_repiration_by_species_in_gtdb_r202_final.txt

$ ./5_apply_respiration_genes.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.apply_respiration_gene_set_aerobic.csv --respiration-genes data/aerobic_repiration_by_species_in_gtdb_r202_final.txt --set-as-aerobic
```

So now we have 6 datasets. To do family-wise GroupKFold validation, we exclude >=20% of the data from training, where no family is both in the training and test sets:
```
$ ./6_split_families.py -i data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv --training-families data/training_families.txt --testing-families data/testing_families.txt
```