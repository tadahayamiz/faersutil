# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

make sqlite3 database by integrating FAERS and other data

< database content >
- case table
    - case_id (unique, int)
    - case_name (str)
    - sex (str)
    - event_date (int)
    - event_country (str)
    - patient_age (int)
    - qualification (int)
    - stored_year (float)
- quality_table
    - qual_id (unique, int)
    - qual_name (str)
- drug table
    - whole_drug_id (unique, int)
    - drug_id (int)
    - drug_name (str)
    - representative (bool)
- rxn table
    - rxn_id (unique, int)
    - pt_name (str)
    - hlt_name (str)
    - hlgt_name (str)
    - soc_name (str)
- drug-rxn table
    - relation_id (unique, int)
    - drug_id (int)
    - rxn_id (int)
    - case_id (int)

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
SEP = os.sep
parser = argparse.ArgumentParser(description='convert FAERS data into a sqlite database')
parser.add_argument('--note', type=str, help='convert FAERS data into a sqlite database')
parser.add_argument(
    'workdir',
    type=str,
    help='working directory that contains outputs of preprocess such as cleansed FAERS data'
    )
args = parser.parse_args()

### main ###
def main():
    """ main function """
    print("=== update drug-dict ===")
    update_drugdict()
    print("> completed")
    print("=== prepare drug-rxn relationship ===")
    prep_drug_rxn()
    print("> completed")
    print("=== initialize database ===")
    init_database()
    print("> completed")
    print("=== update database ===")
    # update dabase
    print("> completed")


### prepare database ###

def update_drugdict():
    """ update drug-dict using name identification with chem_editor """
    print("prepare data", end="...")
    now = datetime.datetime.now().strftime('%Y%m%d')
    path_faers = glob.glob(args.workdir + SEP + "cleansed" + SEP + "clean_*.txt")
    path_ohdsi = glob.glob(args.workdir + SEP + "ohdsi" + SEP + "Drug_dict_*.txt")
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
    print("DONE")
    # get remaining
    print("editing...")
    whole = sorted(list(faers))
    remain = []
    for w in whole:
        try:
            tmp = base_dic[w]
        except KeyError:
            remain.append(w)
    # edit remaining
    conved, summary = ce.main(remain)
    summary.to_csv(
        args.workdir + SEP + "ohdsi" + SEP + f"whole_drug_conversion_{now}.txt", sep="\t"
        )
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
    res.to_csv(args.workdir + SEP + "ohdsi" + SEP + f"Drug_dict_updated_{now}.txt", sep="\t")
    print("DONE")


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
    print("prepare reaction dict", end="...")
    path_rxn = glob.glob(__file__.replace("make_db.py", f"data{SEP}reaction{SEP}*.txt"))
    path_rxn = sorted(path_rxn, reverse=True)[0]
    rxns = pd.read_csv(path_rxn, sep="\t", index_col=0)
    for c in rxns.columns:
        rxns.loc[:, c] = rxns.loc[:, c].map(lambda x: x.lower())
    rxns.loc[:, "rxn_id"] = list(range(rxns.shape[0]))
    dic_rxn = dict(zip(list(rxns["PT"]), list(rxns["rxn_id"])))
    print("DONE")
    # prep drug dict
    print("prepare drug dict", end="...")
    path_ohdsi = glob.glob(args.workdir + SEP + "ohdsi" + SEP + "Drug_dict_updated_*.txt")
    if len(path_ohdsi)==0:
        raise ValueError("!! No Drug_dict_updated data: use 'update_drugdict' before this !!")
    else:
        path_ohdsi = sorted(path_ohdsi, reverse=True)[0]
    ohdsi = pd.read_csv(path_ohdsi, sep="\t", index_col=0)
    dic_drug = dict(zip(list(ohdsi["key"]), list(ohdsi["value"])))
    print("DONE")
    # load FAERS data
    print("prepare FAERS data", end="...")
    path_faers = glob.glob(args.workdir + SEP + "cleansed" + SEP + "clean_*.txt")
    if len(path_faers)==0:
        raise ValueError("!! No clean FAERS data: use 'preprocess' before this !!")
    else:
        path_faers = sorted(path_faers, reverse=True)[0]
    faers = pd.read_csv(path_faers, sep="\t", index_col=0)[
        ["case_id", "active_substances", "reactions", "stored_year"]
        ]
    print("DONE")
    # convert and expand records
    stored_year = list(set(list(faers["stored_year"])))
    print("convert and expand records ...")
    record0 = 0
    record1 = 0
    relation0 = 0
    relation1 = 0
    for s in stored_year:
        print(s)
        arrays = faers[faers["stored_year"]==s].values
        record0 += arrays.shape[0] # num records before NaN treatment
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
        del arrays
        df_tmp = pd.DataFrame(drug_rxn, columns=["active_substances", "reactions"])
        df_tmp.loc[:, "case_id"] = ids
        relation0 += df_tmp.shape[0] # num relation before NaN treatment
        df_tmp = df_tmp.dropna() # delete not found keys
        record1 += df_tmp["case_id"].nunique()
        relation1 += df_tmp.shape[0] # num relation after NaN treatment
        df_tmp.loc[:, "stored_year"] = s
        df_tmp = df_tmp.reset_index(drop=True)
        df_tmp.to_csv(outdir + SEP + f"drug_rxn_{s}_{now}.txt", sep="\t")
        del df_tmp
    summary = pd.DataFrame(
        {"count":[record0, record1, relation0, relation1]},
        index=["n_record_before", "n_record_after", "n_relation_before", "n_relation_after"]
        )
    summary.to_csv(outdir + SEP + f"count_drug_rxn_{now}.txt", sep="\t")
    print("DONE")


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
    print("prepare reaction table", end="...")
    tmp_filein = glob.glob(__file__.replace("make_db.py", f"data{SEP}reaction{SEP}*.txt"))[0]
    if len(tmp_filein)==0:
        raise ValueError("!! No MedDRA data: check faersutil/data/reaction !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0)
    dat.make_rxn_table(df)
    del df
    print("DONE")
    # prepare 
    print("prepare case table", end="...")
    tmp_filein = glob.glob(args.workdir + SEP + "cleansed" + SEP + "clean_*.txt")
    if len(tmp_filein)==0:
        raise ValueError("!! No clean FAERS data: use 'preprocess' before this !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0)
    dat.make_case_table(df)
    del df
    print("DONE")


if __name__ == '__main__':
    main()     
