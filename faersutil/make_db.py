# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

make sqlite3 database by integrating FAERS and other data

< How to use >
- 

< Note >
- 

@author: tadahaya
"""
import os
import datetime
import argparse
import numpy as np
import pandas as pd
import glob
from itertools import chain

from tqdm.auto import trange, tqdm

# original packages in src
from .src import synodict as dh


### setup ###
if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

parser = argparse.ArgumentParser(description='preprocessing of FAERS raw data')
parser.add_argument('--note', type=str, help='preprocessing FAERS raw data (xml, sgml)')
parser.add_argument(
    'workdir',
    type=str,
    help='working directory that contains unzipped FAERS raw directories'
    )
parser.add_argument(
    '-c', '--clean_only', action='store_true',
    help='whether cleansing only or not'
    )
parser.add_argument(
    '-p', '--parse_only', action='store_true',
    help='whether parse only or not'
    )
parser.add_argument(
    '-d', '--drug_only', action='store_true',
    help='whether drug curation only or not'
    )
args = parser.parse_args()


### main ###
def main():
    """ main function """
    if args.sql_only:
        print("=== sql database preparation ===")
        print("> db preparation only")
        init_database()
        print("> DONE")
    else:
        print("=== sql database preparation ===")
        init_database()
        print("> DONE")   



### prepare database ###

def init_database():
    """
    initialize database from FAERS clean data
    
    """
    # init
    now = datetime.datetime.now().strftime('%Y%m%d')
    fileout = args.workdir + SEP + f"sqlite_{now}.db"
    dat = dh.DBhandler()
    dat.set_path(fileout)
    # prepare rxn_table
    tmp_filein = glob.glob(__file__.replace("make_db.py", f"data{SEP}reaction{SEP}*.txt"))[0]
    if len(tmp_filein)==0:
        raise ValueError("!! No MedDRA data: check faersutil/data/reaction !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0)
    dat.make_rxn_table(df)
    del df
    print("> rxn_table is ready")
    # prepare 
    tmp_filein = glob.glob(args.workdir + SEP + "clean_*.txt")
    if len(tmp_filein)==0:
        raise ValueError("!! No clean FAERS data: use 'preprocess' before this !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0)
    dat.make_case_table(df)
    del df
    print("> case_table is ready")


def integrate():
    """ integrate FAERS and OHDSI data """
    # url setting
    path_faers = glob.glob(args.workdir + SEP + "clean_*.txt")
    path_ohdsi = glob.glob(args.workdir + SEP + "Drug_dict_*.txt")
    if len(path_faers)==0:
        raise ValueError("!! No clean FAERS data: use 'preprocess' before this !!")
    else:
        path_faers = sorted(path_faers, reverse=True)[0]
    if len(path_ohdsi)==0:
        raise ValueError("!! No Drug_dict from OHDSI data: use 'preprocess' before this !!")
    else:
        path_ohdsi = sorted(path_ohdsi, reverse=True)[0]
    # prep base dict
    ohdsi = pd.read_csv(path_ohdsi, sep="=t", index_col=0)
    base_dic = dict(zip(list(ohdsi["key"]), list(ohdsi["value"])))
    ohdsi = ohdsi[ohdsi["representative"]==1]
    # load FAERS data
    faers = pd.read_csv(path_faers, sep="=t", index_col=0)
    whole = set(faers["reactions"].map(lambda x: set(x.splite("///"))))
    whole = set(chain.from_iterable(faers.values.tolist()))
    ## use chain.from_iterable is much faster than for loop
    whole = sorted(list(whole))






if __name__ == '__main__':
    main()     
