#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2021 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2021"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os

import pandas as pd
import extern
from io import StringIO

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-tblout', help='path to hmmsearch --tblout', required=True)
    parent_parser.add_argument('--output-csv', help='output one KEGG annotation per gene here', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info("Reading annotations from {}".format(args.input_tblout))
    tbl = args.input_tblout
    hits = extern.run("sed 's/  */\t/g' {} |cut -f1-6 |grep -v '^#'".format(tbl))
    df = pd.read_csv(StringIO(hits), sep='\t', header=None)
    df = df.loc[:,[0,2,4,5]]
    logging.info("Read in {} annotations".format(len(df)))

    df.columns = ['query', 'target', 'evalue', 'score']
    df2 = df.groupby('query', as_index=False).apply(lambda x: x.nsmallest(1, 'evalue')).reset_index().loc[:, ['query', 'target', 'evalue', 'score']]
    logging.info("Found {} unique annotations".format(len(df2)))

    df2.to_csv(args.output_csv, index=False)