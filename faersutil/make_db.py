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
import itertools

from tqdm.auto import trange, tqdm

# original packages in src
from .src import synodict as dh
from .src import chem_editor as ce


### setup ###
if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

parser = argparse.ArgumentParser(description='convert FAERS data into a sqlite database')
parser.add_argument('--note', type=str, help='convert FAERS data into a sqlite database')
parser.add_argument(
    'workdir',
    type=str,
    help='working directory that contains outputs of preprocess such as cleansed FAERS data'
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


def update_drugdict():
    """ update drug-dict using name identification with chem_editor """
    now = datetime.datetime.now().strftime('%Y%m%d')
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
    ohdsi = pd.read_csv(path_ohdsi, sep="\t", index_col=0)
    base_dic = dict(zip(list(ohdsi["key"]), list(ohdsi["value"])))
    rep = list(ohdsi[ohdsi["representative"]==1]["key"])
    # load FAERS data
    faers = pd.read_csv(path_faers, sep="\t", index_col=0)
    faers = faers["active_substances"].map(lambda x: set(x.split("///")))
    faers = set(itertools.chain.from_iterable(faers.values.tolist()))
    # note: chain is much faster than loop
    # get remaining
    whole = sorted(list(faers))
    remain = []
    for w in tqdm(whole):
        try:
            tmp = base_dic[w]
        except KeyError:
            remain.append(w)
    # edit remaining
    conved, summary = ce.main(remain)
    summary.to_csv(args.workdir + SEP + f"whole_drug_conversion_{now}.txt", sep="\t")
    new_dic = dict()
    for r, c in zip(remain, conved):
        if r!=c:
            try:
                tmp = base_dic[c]
                new_dic[r] = tmp
            except KeyError:
                pass
    base_dic.update(new_dic)
    # export
    res = pd.DataFrame({"key":base_dic.keys(), "value":base_dic.values()})
    res.loc[:, "representative"] = 0
    res.loc[res["key"].isin(rep), "representative"] = 1
    res.to_csv(args.workdir + SEP + f"Drug_dict_updated_{now}.txt", sep="\t")


def prep_drug_rxn():
    """
    prepare drug-rxn table
    Due to memory restriction, data is exported by each stored_year
    Note time-consuming
    
    """
    # init
    now = datetime.datetime.now().strftime('%Y%m%d')
    outdir = args.workdir + SEP + "drug_rxn"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # prep reaction dict
    path_rxn = glob.glob(__file__.replace("make_db.py", f"data{SEP}reaction{SEP}*.txt"))
    path_rxn = sorted(path_rxn, reverse=True)[0]
    rxns = pd.read_csv(path_rxn, sep="\t", index_col=0)
    for c in rxns.columns:
        rxns.loc[:, c] = rxns.loc[:, c].map(lambda x: x.lower())
    rxns.loc[:, "rxn_id"] = list(range(rxns.shape[0]))
    dic_rxn = dict(zip(list(rxns["PT"]), list(rxns["rxn_id"])))
    # prep drug dict
    path_ohdsi = glob.glob(args.workdir + SEP + "Drug_dict_updated_*.txt")
    if len(path_ohdsi)==0:
        raise ValueError("!! No Drug_dict_updated data: use 'update_drugdict' before this !!")
    else:
        path_ohdsi = sorted(path_ohdsi, reverse=True)[0]
    ohdsi = pd.read_csv(path_ohdsi, sep="\t", index_col=0)
    dic_drug = dict(zip(list(ohdsi["key"]), list(ohdsi["value"])))
    # load FAERS data
    path_faers = glob.glob(args.workdir + SEP + "clean_*.txt")
    if len(path_faers)==0:
        raise ValueError("!! No clean FAERS data: use 'preprocess' before this !!")
    else:
        path_faers = sorted(path_faers, reverse=True)[0]
    faers = pd.read_csv(path_faers, sep="\t", index_col=0)[
        ["case_id", "active_substances", "reactions", "stored_year"]
        ]
    # expand and convert one by one
    stored_year = list(set(list(faers["stored_year"])))
    for s in stored_year:
        print(s)
        arrays = faers[faers["stored_year"]==s].values
        drug_rxn = []
        ids = []
        for i in trange(arrays.shape[0]):
            cid = arrays[i, 0] # case_id
            drugs = arrays[i, 1].split("///") # active_substances
            rxns = arrays[i, 2].split("///") # reactions
            drugs = [dic_drug.get(k, np.nan) for k in drugs]
            rxns = [dic_rxn.get(k, np.nan) for k in rxns]
            drug_rxn += list(itertools.product(drugs, rxns)) # cartesian product
            ids += [cid] * len(drugs) * len(rxns)
        df_tmp = pd.DataFrame(drug_rxn, columns=["active_substances", "reactions"])
        df_tmp.loc[:, "case_id"] = ids
        df_tmp = df_tmp.dropna() # delete not found keys
        df_tmp.loc[:, "stored_year"] = s
        df_tmp = df_tmp.reset_index(drop=True)
        df_tmp.to_csv(outdir + SEP + f"drug_rxn_{s}_{now}.txt", sep="\t")


if __name__ == '__main__':
    main()     
