#!/usr/bin/env python3

from sklearn.model_selection import GroupKFold, cross_val_score
from xgboost import XGBClassifier
import polars as pl
import pandas as pd
import numpy as np
from joblib import dump, load
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

    # ./12_apply_model.py --model the.model -x the.x --training-data-header <(head the.training-data-header)
    parent_parser.add_argument('--models', nargs='+', help='models to use', required=True)
    parent_parser.add_argument('-x', help='table of inputs to use', required=True)
    parent_parser.add_argument('--training-data-header', help='header of training data', required=True)
    parent_parser.add_argument('--output-predictions', help='output predictions', required=True)
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')


    # Read in data
    # d = pl.read_csv('TableAncestralRoot1.tsv',sep="\t")
    d = pd.read_csv(args.x, sep="\t")
    logging.info("Read in input data of shape {}".format(d.shape))

    # Collapse counts of each COG subfamily
    d2 = d
    d2['COG'] = d2['COG'].str.split('_').str[0]
    d3 = d2.groupby('COG').sum()
    d4 = d3.transpose()

    # Read training data header
    eg_data = pl.read_csv(args.training_data_header, separator="\t", has_header=True)
    header = eg_data.select(pl.exclude([
        'accession','false_negative_rate','false_positive_rate'
    ])).columns
    # Blacklist these as they aren't in the current ancestral file, not sure why
    header = list([h for h in header if h not in ['COG0411', 'COG0459', 'COG0564', 'COG1344', 'COG4177']])
    
    # Reorder columns to be the same as the training dataset
    d5 = d4[header]

    all_results = []
    for model_path in args.models:
        logging.info("Loading model {}".format(model_path))
        model = load(model_path)
        logging.info("Loaded model {}".format(model_path))

        doing_perceptron = 'Perceptron' in model_path

        preds = model.predict(d5)
        if doing_perceptron:
            probas = pl.lit(-1.0)
        else:
            probas = model.predict_proba(d5)[:,1]

        results = pl.DataFrame({
            'node': list(d5.index.values),
            'preds': preds,
            'proba': -1
        })
        results = results.select(
            pl.col('node'),
            pl.col('preds').alias('prediction').cast(pl.Int64),
            pl.col('proba').alias('probability').cast(pl.Float64),
            pl.lit(model_path).alias('model')
        )
        all_results.append(results)

    pl.concat(all_results).write_csv(args.output_predictions, separator="\t")
    logging.info("Wrote predictions to {}".format(args.output_predictions))
    