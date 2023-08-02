# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

preprocessor of FAERS raw data (xml)

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
from .src import xml_loader as xm
from .src import cleanser as cl

# setup
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
args = parser.parse_args()

### main ###

def main():
    """ main function """
    if args.parse_only:
        print("=== parse xml files ===")
        print("> parse only")
        parse_xml()
        print("> DONE")
        # print("=== parse sgml files ===")
        # parse_sgml()
        # print("> DONE")
    elif args.clean_only:
        print("=== clean and merge files ===")
        print("> cleansing only")
        clean_and_merge()
        print("> DONE")
    else:
        print("=== parse xml files ===")
        parse_xml()
        print("> DONE")
        print("=== clean and merge files ===")
        clean_and_merge()
        print("> DONE")   

### xml to tsv ###

def parse_xml():
    """
    parse xml files
    Note files before 2014Q2 and after it are different format
    and need different modules
    
    """
    # url setting
    path_list = glob.glob(args.workdir + SEP + "faers_xml*")
    keys = []
    ## both q and Q exist to indicate quarter
    for p in path_list:
        tmp = p.split("_")[-1].replace("q", ".").replace("Q", ".")
        keys.append(float(tmp))
    outdir = args.workdir + SEP + "parsed"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # parse
    for p_dir, key in tqdm(zip(path_list, keys)):
        path_xmls = glob.glob(p_dir + f"{SEP}xml{SEP}*.xml")
        if len(path_xmls)==0:
            # 2020~
            path_xmls = glob.glob(p_dir + f"{SEP}XML{SEP}*.xml")
        if key > 2014.2:
            # from 2014Q3
            for i, path in enumerate(path_xmls):
                dat = xm.Loader2014Q3_()
                key = str(key).replace(".", "q")
                fileout = outdir + SEP + f"parsed_{key}_{i}.txt"
                try:
                    df = dat.load_xml(path, fileout=fileout, to_csv=True, sep="\t")
                except Exception as e:
                    print(f"!! An Error Ocurred: {e} !!")
                    print(path)
        else:
            pass
            # note: skip before 2014Q3 since the current program became inappropriate
            # from 2012Q4 to 2014Q2
            # for i, path in enumerate(path_xmls):
            #     dat = xm.Loader2012Q4_2014Q2()
            #     key = str(key).replace(".", "q")
            #     fileout = outdir + SEP + f"parsed_{key}_{i}.txt"
            #     try:
            #         df = dat.load_xml(path, fileout=fileout, to_csv=True, sep="\t")
            #     except:
            #         print("!! An Error Ocurred !!")
            #         print(path)


# def parse_sgml():
#     """
#     parse sgml files, old format but can be read like xml
    
#     """
#     # url setting
#     path_list = glob.glob(args.workdir + SEP + "aers_sgml*")
#     ## both q and Q exist to indicate quarter
#     keys = []
#     for p in tqdm(path_list):
#         tmp = p.split("_")[-1].replace("q", ".").replace("Q", ".")
#         keys.append(float(tmp))
#     outdir = args.workdir + SEP + "parsed"
#     if not os.path.exists(outdir):
#         os.makedirs(outdir)
#     # parse
#     for p_dir, key in zip(path_list, keys):
#         path_sgmls = glob.glob(p_dir + f"{SEP}sgml{SEP}*.SGM")
#         # from 2012Q4 to 2014Q2
#         for i, path in enumerate(path_sgmls):
#             dat = xm.Loader2012Q4_2014Q2()
#             key = str(key).replace(".", "q")
#             fileout = outdir + SEP + f"parsed_{key}_{i}.txt"
#             try:
#                 df = dat.load_xml(path, fileout=fileout, to_csv=True, sep="\t")
#             except:
#                 print("!! An Error Ocurred !!")
#                 print(path)

### cleansing ###

def clean_and_merge():
    """ cleansing parsed files """
    path_list = glob.glob(args.workdir + SEP + "parsed" + SEP + "*.txt")
    now = datetime.datetime.now().strftime('%Y%m%d')
    fileout = args.workdir + SEP + f"clean_{now}.txt"
    results = []
    for path in tqdm(path_list):
        fname = path.split("parsed_")[-1]
        key = fname.split("_")[0].replace("q", ".")
        df = pd.read_csv(path, sep="\t").dropna(subset=["Active Substances"]) # Some records miss it
        dat = cl.Cleanser()
        dat.set_data(data=df)
        dat.data_cleansing()
        dat.set_exception()
        dat.exclude_exception()
        res = dat.get_data()
        res.loc[:, "Stored Year"] = [float(key)] * res.shape[0] # add store date
        results.append(res)
    results = pd.concat(results, axis=0, join="inner")
    # set to str
    results.loc[:,"Active Substances"] = results.loc[:,"Active Substances"].map(list).apply(_concat)
    results.loc[:,"Reactions"] = results.loc[:,"Reactions"].map(list).apply(_concat)
    # care after concatenation
    ## check duplicates
    results = results.drop_duplicates(subset=["Case ID"], keep='first')
    results = results.reset_index(drop=True)    
    # export
    col = [v.replace(" ", "_").lower() for v in list(results.columns)]
    results.columns = col
    results.to_csv(fileout, sep="\t")


def _concat(lst:list, gap:str="///"):
    """ concat the contents of a list with the indicated gap """
    if len(lst) > 1:  
        res = lst[0]
        for v in lst[1:]:
            res += gap + v
        return res
    elif len(lst)==1:
        return lst[0]
    else:
        return np.nan


if __name__ == '__main__':
    main()