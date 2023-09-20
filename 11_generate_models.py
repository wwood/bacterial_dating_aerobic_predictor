#!/usr/bin/env python3

from sklearn.model_selection import GroupKFold
from xgboost import XGBClassifier
import polars as pl
import pandas as pd
from joblib import dump
import argparse
import logging
import os

from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier,AdaBoostClassifier,ExtraTreesClassifier
from sklearn.preprocessing import MaxAbsScaler
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression,Perceptron
from sklearn.pipeline import make_pipeline
from sklearn.calibration import CalibratedClassifierCV

"""
Generate and train ML models based on temperature range predictions. 
"""

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    # ./11_generate_models.py --training-file data/all_gene_annotations.added_incompleteness_and_contamination.training.tsv --testing-file data/all_gene_annotations.added_incompleteness_and_contamination.testing.tsv -y data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv --output-dir data/bacdive_scrape_20230315.json.parsed.anaerobe_vs_aerobe.with_cyanos.csv.models
    parent_parser.add_argument('--training-file', help='training file e.g. completeness_removed_contamination_added_training_data.csv', required=True)
    parent_parser.add_argument('--testing-file', help='testing file e.g. completeness_removed_contamination_added_testing_data.csv', required=True)
    parent_parser.add_argument('--model-output-dir', help='output directory', required=True)
    parent_parser.add_argument('--cross-validation-data-output-dir', help='output directory', required=True)
    parent_parser.add_argument('-y', help='annotations to train/predict on', required=True)
    parent_parser.add_argument('--add-calibrated-models', help='add calibrated models', action="store_true")
    parent_parser.add_argument('--pipelines', nargs='+', help='Use these models only')
    # set target column name and levels
    parent_parser.add_argument('--target-column', help='target column name', default='temperaturerange')
    parent_parser.add_argument('--target-levels', nargs='+', help='target levels')
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

    # Ready output directory
    if not os.path.exists(args.model_output_dir):
        os.makedirs(args.model_output_dir)
    if not os.path.exists(args.cross_validation_data_output_dir):
        os.makedirs(args.cross_validation_data_output_dir)

    target_column = args.target_column

    # Read y
    y0 = pl.read_csv(args.y, separator="\t")
    y1 = y0.unique() # There are some duplicates in the cyanos, so dedup
    logging.info("Read y: %s", y1.shape)
    # Log counts of each class
    logging.info("Counts of each class amongst unique accessions: %s", y1.groupby(target_column).agg(pl.count()))

    # Read GTDB
    gtdb = pl.concat([
        pl.read_csv('data/bac120_metadata_r202.tsv', separator="\t"),
        pl.read_csv('data/ar122_metadata_r202.tsv', separator="\t")
    ])
    gtdb = gtdb.filter(pl.col("gtdb_representative") == "t")
    logging.info("Read in {} GTDB reps".format(len(gtdb)))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(1).alias("phylum"))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(2).alias("class"))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(3).alias("order"))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(4).alias("family"))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(5).alias("genus"))

    # Read training data
    d = pl.read_csv(args.training_file,separator="\t")
    logging.info("Read training data: %s", d.shape)

    # Ignore all but training data
    d2 = d.join(gtdb.select(['accession','phylum','class','order','family','genus']), on="accession", how="left")
    d3 = d2.join(y1, on="accession", how="inner") # Inner join because test accessions are in y1 but not in d2
    logging.info("Counts of each class in training/test data: %s", d3.groupby(target_column).agg(pl.count()))

    X = d3.select(pl.exclude(['accession',target_column,'phylum','class','order','family','genus','false_negative_rate','false_positive_rate'])).to_pandas()
    # Map oxytolerance to 0, 1, 2
    if args.target_levels:
        classes_map = {k: i for i, k in enumerate(args.target_levels)}
    else:
        if 'anaerobic_with_respiration_genes' in d3['temperaturerange'].to_list():
            classes_map = {
                'mesophilic': 0,
                'hyperthermophilic': 1,
                'thermophilic': 2,
                'psychrophilic': 3,
            }
        else:
            classes_map = {
                'mesophilic': 0,
                'hyperthermophilic': 1,
                'thermophilic': 2,
                'psychrophilic': 3,
            }

    y = d3.select(pl.col(target_column).apply(lambda x: classes_map[x]).alias(target_column))
    logging.info("Counts of y: %s", y.groupby(target_column).agg(pl.count()))
    y = y.to_pandas()
    
    groups = d3['family'].to_list()

    d_gtdb = d3.to_pandas()

    # Blacklist these as they aren't in the current ancestral file, not sure why
    X = X.drop(['COG0411', 'COG0459', 'COG0564', 'COG1344', 'COG4177'],axis=1)

    n_jobs=64
    classifiers = [
        LogisticRegression(max_iter=1000,n_jobs=n_jobs),
        RandomForestClassifier(n_jobs=n_jobs,n_estimators=1000),
        ExtraTreesClassifier(n_jobs=n_jobs,n_estimators=1000),
        AdaBoostClassifier(),
        GradientBoostingClassifier(learning_rate=0.1),
        Perceptron(),
        GaussianNB()]
    classifiers_names = [
        "LogisticRegression","RandomForest","ExtraTrees",
        "AdaBoostClassifier","GradientBoosting"
        ,"Perceptron"
        ,"GaussianNB"]
    pipes = [
        make_pipeline(
        #StandardScaler(),
        MaxAbsScaler(),    
        classifier
        ) for classifier in classifiers]

    # Add in xgboost
    classifiers = [XGBClassifier(n_jobs=n_jobs, use_label_encoder=False)] + classifiers
    classifiers_names = ["XGBoost"] + classifiers_names

    pipes = [
        make_pipeline(
        #StandardScaler(),
        MaxAbsScaler(),
        classifier
        ) for classifier in classifiers]
    
    # Add calibrated models
    names_and_pipes = []
    for (name, pipe) in zip(classifiers_names, pipes):
        if args.pipelines is None or name in args.pipelines:
            names_and_pipes.append((name, pipe))
        if args.add_calibrated_models:
            names_and_pipes.append((name + "_Isotonic", CalibratedClassifierCV(pipe, cv=5, method='isotonic')))
            names_and_pipes.append((name + "_Sigmoid", CalibratedClassifierCV(pipe, cv=5, method='sigmoid')))
    logging.info("Using these pipelines: {}".format([np[0] for np in names_and_pipes]))

    gkf = GroupKFold(n_splits=5)
    for i, (train, test) in enumerate(gkf.split(X, y, groups=groups)):
        for (model_name, model) in names_and_pipes:
            logging.info("Fold %i, Training model %s .." % (i, model_name))
            model.fit(X.iloc[train], y.iloc[train, 0])
            y_pred = model.predict(X.iloc[test])
            
            y_actual = y.iloc[test, 0].values

            if "Perceptron" in model_name:
                # Perception doesn't have predict_proba
                df1 = pd.DataFrame(y_pred, columns=['prediction'])
            else:
                pp = model.predict_proba(X.iloc[test])
                df1 = pd.DataFrame(
                    pp,
                    columns=[f"probability_{classes_map[k]}" for k in classes_map.keys()],
                )
                df1['prediction'] = y_pred

            df1['accession'] = d_gtdb.loc[test, 'accession'].values
            df1['y_actual'] = y_actual
            df1['false_negative_rate'] = d_gtdb.loc[test, 'false_negative_rate'].values
            df1['false_positive_rate'] = d_gtdb.loc[test, 'false_positive_rate'].values
            df1['predictor'] = model_name

            print("\t".join([
                str(i),
                model_name,
                str(accuracy_score(y_actual, y_pred)),
            ]))

            df1.to_csv('{}/prediction_probabilities_{}_{}.csv'.format(
                args.cross_validation_data_output_dir, model_name, i), index=False, sep="\t", header=True)



    # Generate final predictors that include no cross-validation removal of samples
    logging.info("Creating final predictor")
    for (model_name, model) in names_and_pipes:
        logging.info("Creating non-cross-validation predictor for {}".format(model_name))
        model.fit(X, y.iloc[:, 0])

        model_filename = os.path.join(args.model_output_dir, '{}.model'.format(model_name))
        dump(model, model_filename)

    logging.info("Done")