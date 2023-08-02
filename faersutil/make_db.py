# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

CLI preprocesser

@author: tadahaya
"""
import sys
import os
import datetime
import argparse
import numpy as np
import pandas as pd
import glob

from tqdm.auto import trange, tqdm

# original packages in src
from .src import db_handler as dh

# setup
if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

parser = argparse.ArgumentParser(description='CLI preprocess')
parser.add_argument('--note', type=str, help='preprocessing FAERS raw data (xml, sgml)')
parser.add_argument(
    'workdir',
    type=str,
    help='working directory that contains unzipped FAERS raw directories'
    )
parser.add_argument(
    '-s', '--sql_only', action='store_true',
    help='whether sql only or not'
    )
args = parser.parse_args()

### main ###
def main():
    """ main function """
    if args.sql_only:
        print("=== sql database preparation ===")
        print("> db preparation only")
        prep_database()
        print("> DONE")
    else:
        print("=== sql database preparation ===")
        prep_database()
        print("> DONE")   


### prepare database ###

def prep_database():
    """ prepare database from clean data """
    raise NotImplementedError

    # init
    now = datetime.datetime.now().strftime('%Y%m%d')
    fileout = args.workdir + SEP + f"sqlite_{now}.db"
    dat = DBhandler()
    dat.set_path(fileout)
    # prepare rxn_table
    tmp_filein = glob.glob(__file__.replace("cli_preprocess.py", "reaction") + SEP + "*.txt")[0]
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0)
    dat.make_rxn_table(df)
    del df
    print("> rxn_table is ready")
    # prepare 
    tmp_filein = glob.glob(args.workdir + SEP + "clean_*.txt")[0]
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0)
    dat.make_case_table(df)
    del df
    print("> case_table is ready")


if __name__ == '__main__':
    main()     
