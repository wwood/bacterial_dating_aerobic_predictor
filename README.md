Robust machine learning-based prediction of whether a microorganism is aerobic or anaerobic.

This predictor was used in this (currently pre-publication) manuscript, please cite it if appropriate:

Davin, Adrian A., Ben J. Woodcroft, Rochelle M. Soo, Ranjani Murali, Dominik Schrempf, James Clark, Bastien Boussau et al. "An evolutionary timescale for Bacteria calibrated using the Great Oxidation Event." bioRxiv (2023): 2023-08. https://www.biorxiv.org/content/10.1101/2023.08.08.552427v1.full

This repository details how the data was gathered and how the predictor was trained. The predictor itself is available at https://huggingface.co/wwood/aerobicity

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

Then run the scrape (do not include the `$` in the command):
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
and for the predictor not using cyano annotations, just transform the data:
```
$ ./4_add_cyanos.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.no_cyanos.csv
```

Integrate manually curated respiration gene info

```
$ ./5_apply_respiration_genes.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.apply_respiration_gene_exclusion.csv --respiration-genes data/aerobic_repiration_by_species_in_gtdb_r202_final.txt

$ ./5_apply_respiration_genes.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.apply_respiration_gene_set_aerobic.csv --respiration-genes data/aerobic_repiration_by_species_in_gtdb_r202_final.txt --set-as-aerobic

$ ./5_apply_respiration_genes.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.apply_respiration_gene_exclusion.csv --respiration-genes data/aerobic_repiration_by_species_in_gtdb_r202_final.txt

$ ./5_apply_respiration_genes.py --input-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv --output-csv data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.apply_respiration_gene_new_class.csv --respiration-genes data/aerobic_repiration_by_species_in_gtdb_r202_final.txt --set-as-new-class
```

So now we have the datasets. To do family-wise GroupKFold validation, we exclude >=20% of the data from training, where no family is both in the training and test sets:
```
$ ./6_split_families.py -i data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv --training-families data/training_families.txt --testing-families data/testing_families.txt
```

Create a list of dataset files so we can start to do these things in parallel.

TODO: Add no_cyano dataset to this list
```
$ echo \
   bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.csv \
   bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv \
   bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.apply_respiration_gene_exclusion.csv \
   bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.apply_respiration_gene_set_aerobic.csv \
   bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.apply_respiration_gene_exclusion.csv \
   bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.apply_respiration_gene_set_aerobic.csv \
   bacdive_scrape_20230315.json.parsed.anaerobe_vs_rest.with_cyanos.apply_respiration_gene_new_class.csv \
   |sed 's/ /\n/g' \
    >data/dataset_files.txt
$ cat data/dataset_files.txt |parallel ls data/{} -1 # Should not fail
```

Get a list of all genome IDs, maybe not needed for this, but might be useful for generating the eggnog/kegg annotations
```
$ cat data/dataset_files.txt |sed 's/^/data\//g' |parallel tail -n+2 |cut -f1 |sort |uniq >data/all_genome_ids
```

Process KEGG annotations. KEGG annotations are first annotated using e.g.
```
$ ls proteomes/*faa |parallel hmmsearch --tblout kegg_annotations/annotations/{/}.hmmsearch_tblout.csv -o /dev/null --notextw --cpu 1 ~/m/db/kofam/2022-01-30/profiles.hmm {}
```
and then we take the best hit for each protein
```
$ ls kegg_annotations/annotations/*_tblout.csv |parallel -j10 ./7_process_kegg_hmmsearch.py --input-tblout {} --output-csv kegg_annotations/best_hits/{/}.best_hits.csv &>7_process_kegg_hmmsearch.py.log
```

Generate a master list of all genomes in all datasets with their annotations
```
$ ./8_gather_genome_annotations.py -i `cat data/dataset_files.txt |sed s=^=data/=` -o data/all_gene_annotations.tsv --eggnog-mapper-annotations-key data/eggnog_mapper_annotation_paths.tsv --kegg-best-hits-annotations-key data/kegg_best_hits_key
```

Generate incompleteness and contamination for all genomes in all datasets
```
$ ./9_expand_incompletenss_and_contamination4.py --input-file data/all_gene_annotations.tsv --output-file data/all_gene_annotations.added_incompleteness_and_contamination.tsv
```

Generate especially low completeness and contamination datasets
```
$ ./9_expand_incompletenss_and_contamination4.py --low-completeness --input-file data/all_gene_annotations.tsv --output-file data/all_gene_annotations.added_very_low_completeness_and_contamination.tsv
```

Split out training and test datasets based on their families chosen randomly above
```
$ ./10_split_training_test_sets.py --input-file data/all_gene_annotations.added_incompleteness_and_contamination.tsv --training-families data/training_families.txt --testing-families data/testing_families.txt --output-training data/all_gene_annotations.added_incompleteness_and_contamination.training.tsv --output-testing data/all_gene_annotations.added_incompleteness_and_contamination.testing.tsv
```

Do the actual training and testing
```
$ cat data/dataset_files.txt |parallel -j1 ./11_generate_models.py --training-file data/all_gene_annotations.added_incompleteness_and_contamination.training.tsv --testing-file data/all_gene_annotations.added_incompleteness_and_contamination.testing.tsv -y data/{} --model-output-dir data/{}.models --cross-validation-data-output-dir data/{}.cv_data
```

Inspect results of training according to cross validation - run this and then inspect the ipynb file.
```
$ jupyter nbconvert --execute 12_inspect_cross_validation.ipynb --to notebook --inplace
```

Apply the predictor against the ancestral state predictions
```
$ cat data/dataset_files.txt |parallel -j1 ./13_apply_model.py --model data/{}.models/*.model -x data/TableAncestralRoot1.tsv --training-data-header data/all_gene_annotations.added_incompleteness_and_contamination.testing.tsv --output-predictions data/{}.predictions.csv
```

Determine SHAP values for the models
```
jupyter nbconvert --execute 15_inspect_shap_values.ipynb --to notebook --inplace
```
