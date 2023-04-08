#!/usr/bin/env python3

from sklearn.model_selection import GroupKFold, cross_val_score
from xgboost import XGBClassifier
import polars as pl
import pandas as pd
import numpy as np
from joblib import dump
import argparse
import logging
import os

from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier,AdaBoostClassifier,ExtraTreesClassifier
from sklearn.preprocessing import StandardScaler,MaxAbsScaler
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression,Perceptron
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.pipeline import make_pipeline

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
    parent_parser.add_argument('-y', help='oxytolerance annotations', required=True)
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

    # Read y
    y1 = pl.read_csv(args.y, separator="\t")
    logging.info("Read y: %s", y1.shape)
    # Log counts of each class
    logging.info("Counts of each class amongst unique accessions: %s", y1.groupby("oxytolerance").agg(pl.count()))

    # Read GTDB
    gtdb = pl.concat([
        pl.read_csv('data/bac120_metadata_r202.tsv', separator="\t"),
        pl.read_csv('data/ar122_metadata_r202.tsv', separator="\t")
    ])
    gtdb = gtdb.filter(pl.col("gtdb_representative") == "t")
    logging.info("Read in {} GTDB reps".format(len(gtdb)))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').arr.get(1).alias("phylum"))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').arr.get(2).alias("class"))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').arr.get(3).alias("order"))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').arr.get(4).alias("family"))
    gtdb = gtdb.with_columns(pl.col("gtdb_taxonomy").str.split(';').arr.get(5).alias("genus"))

    # Read training data
    d = pl.read_csv(args.training_file,separator="\t")
    logging.info("Read training data: %s", d.shape)

    # Ignore all but training data
    d2 = d.join(gtdb.select(['accession','phylum','class','order','family','genus']), on="accession", how="left")
    d3 = d2.join(y1, on="accession", how="inner") # Inner join because test accessions are in y1 but not in d2
    logging.info("Counts of each class in training/test data: %s", d3.groupby("oxytolerance").agg(pl.count()))

    X = d3.select(pl.exclude(['accession','oxytolerance','phylum','class','order','family','genus','false_negative_rate','false_positive_rate'])).to_pandas()
    y = d3.select(pl.when(pl.col('oxytolerance')== 'aerobe').then(1).otherwise(0)).to_pandas()
    groups = d3['family'].to_list()

    d_gtdb = d3.to_pandas()

    # Blacklist these as they aren't in the current ancestral file, not sure why
    X = X.drop(['COG0411', 'COG0459', 'COG0564', 'COG1344', 'COG4177'],axis=1)

    #bunch of semi random classifiers from sklearn 
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


    gkf = GroupKFold(n_splits=5)
    for i, (train, test) in enumerate(gkf.split(X, y, groups=groups)):
        for model, model_name in zip(pipes, classifiers_names):
            logging.info("Fold %i, Training model %s .." % (i, model_name))
            model.fit(X.iloc[train], y.iloc[train, 0])
            y_pred = model.predict(X.iloc[test])
            
            y_actual = y.iloc[test, 0]

            if model_name == "Perceptron":
                # Perception doesn't have predict_proba
                df1 = pd.DataFrame(y_pred, columns=['prediction'])
            else:
                pp = model.predict_proba(X.iloc[test])
                df1 = pd.DataFrame(
                    pp,
                    columns=['prob_anaerobe', 'prob_aerobe'],
                )
                df1['prediction'] = y_pred

            df1['accession'] = d_gtdb.loc[test, 'accession'].values
            # df1['y_actual'] = y_actual.array
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
    # results = {}
    for model, model_name in zip(pipes, classifiers_names):
        logging.info("Creating non-cross-validation predictor for {}".format(model_name))
        model.fit(X, y.iloc[:, 0])

        model_filename = os.path.join(args.model_output_dir, '{}.model'.format(model_name))
        dump(model, model_filename)

    logging.info("Done")