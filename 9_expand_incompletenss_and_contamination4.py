#!/usr/bin/env python3

import pandas as pd
import random
import os
from tqdm import tqdm
import argparse
import logging


def remove_completeness_and_add_contamination(row, fraction_to_remove, fraction_to_add, cog_weights,
                                              cog_weight_indices):
    # Make an array of indices
    indices = [[i] * count for i, count in enumerate(row)]

    # Flatten the list
    indices = list([item for sublist in indices for item in sublist])

    # Choose a random fraction_to_remove of these
    # Return a new row with these removed
    num_to_remove = int(fraction_to_remove * len(indices))
    row_remaining = row.copy()
    for i in random.sample(indices, num_to_remove):
        row_remaining[i] -= 1

    # false positive rate is independent of genome size
    num_to_add = int(fraction_to_add * len(cog_weights))
    if num_to_add > 0:
        chosens = random.choices(cog_weight_indices, weights=cog_weights, k=num_to_add)
        for chosen in chosens:
            row_remaining[chosen] += 1

    return row_remaining


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    # parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--low-completeness', help='generate very incomplete data', action="store_true")
    parent_parser.add_argument('--input-file',
                               help='input file e.g. x_and_y_training_data_annotations.csv',
                               required=True)
    parent_parser.add_argument('--output-file',
                               help='output file e.g. completeness_removed_contamination_added_training_data.csv',
                               required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    x_and_y = pd.read_csv(args.input_file, sep="\t")

    # Create a new dataset with varying amounts of completeness
    # For each row in the training set, add the same datapoint with 5 levels of varying completeness
    # d = x_and_y.iloc[:2,:].copy()
    fractions_to_add = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    if args.low_completeness:
        fractions_to_remove = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    else:
        fractions_to_remove = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    logging.info("Using fractions to add {} and fractions to remove {}".format(fractions_to_add, fractions_to_remove))

    # Create a single vector containing the number of counts of each COG across the entire dataset, for false positive generation
    cog_weights = list(x_and_y.iloc[:, 1:].sum(axis=0).values)
    cog_weight_indices = list([i for i, _ in enumerate(cog_weights)])

    num_genomes = len(x_and_y)

    output_file = args.output_file
    if os.path.exists(output_file):
        os.remove(output_file)

    first_genome = True
    for i in tqdm(range(num_genomes)):
        this_genomes_data = None
        for fraction_to_remove in fractions_to_remove:
            for fraction_to_add in fractions_to_add:
                to_add = pd.concat([
                    x_and_y.iloc[i, :1],
                    pd.Series({'false_negative_rate': fraction_to_remove}),
                    pd.Series({'false_positive_rate': fraction_to_add}),
                    remove_completeness_and_add_contamination(x_and_y.iloc[i, 1:], fraction_to_remove, fraction_to_add,
                                                              cog_weights, cog_weight_indices)
                ])
                to_add_df = pd.DataFrame(to_add).transpose()
                if this_genomes_data is None:
                    this_genomes_data = to_add_df
                else:
                    this_genomes_data = pd.concat([this_genomes_data, to_add_df])

        this_genomes_data.to_csv(output_file, sep="\t", index=False, mode='a', header=first_genome)
        first_genome = False