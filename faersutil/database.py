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
import time
import argparse
import numpy as np
import pandas as pd
import glob
import itertools

from tqdm.auto import trange, tqdm

# original packages in src
from .src import chem_editor as ce
from .src import db_handler as dh

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
    # update drug dict
    update_drugdict()
    # prepare drug-rxn data
    prep_drug_rxn()
    # generate database
    make_database()


### prepare database ###

def update_drugdict():
    """
    update drug-dict using name identification with chem_editor
    
    Returns
    -------
    tsv, the updated drug dict
    
    """
    print("prepare data", end="...")
    now = datetime.datetime.now().strftime('%Y%m%d')
    path_faers = glob.glob(args.workdir + SEP + "clean" + SEP + "clean_*.txt")
    outdir = args.workdir + SEP + "curated"
    path_ohdsi = glob.glob(outdir + SEP + "Drug_dict_*.txt")
    path_ohdsi = [v for v in path_ohdsi if "Drug_dict_updated" not in v]
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
    summary.to_csv(outdir + SEP + f"whole_drug_conversion_{now}.txt", sep="\t")
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
    res.loc[:, "drug_dict_id"] = list(range(res.shape[0])) # add index for DB construction
    res.to_csv(outdir + SEP + f"Drug_dict_updated_{now}.txt", sep="\t")
    print("> completed")


def prep_drug_rxn():
    """
    prepare drug-rxn table
    Due to memory restriction, data is exported by each stored_year
    Note time-consuming

    Returns
    -------
    - tsv files, drug-rxn tables for each stored_year, pached in drug_rxn dir
    - tsv, summarizing the number of records before and after conversion
    
    """
    # init
    now = datetime.datetime.now().strftime('%Y%m%d')
    outdir = args.workdir + SEP + "drug_rxn"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # prep reaction dict
    print("prepare reaction dict", end="...")
    path_rxn = glob.glob(__file__.replace("database.py", f"data{SEP}reaction{SEP}*.txt"))
    path_rxn = sorted(path_rxn, reverse=True)[0]
    rxns = pd.read_csv(path_rxn, sep="\t", index_col=0)
    for c in rxns.columns:
        rxns.loc[:, c] = rxns.loc[:, c].map(lambda x: x.lower())
    ## prep rxn_id that corresponds to pt
    pt = list(set(list(rxns["PT"])))
    dic_rxn = dict(zip(pt, list(range(len(pt)))))
    rxns.loc[:, "rxn_id"] = [dic_rxn[i] for i in list(rxns["PT"])]
    rxns.loc[:, "unique_rxn_id"] = list(range(rxns.shape[0]))
    rxns.to_csv(args.workdir + SEP + "curated" + f"/rxn_table_{now}.txt", sep="\t")
    ## save w/ ID version
    print("DONE")
    # prep drug dict
    print("prepare drug dict", end="...")
    path_ohdsi = glob.glob(args.workdir + SEP + "curated" + SEP + "Drug_dict_updated_*.txt")
    if len(path_ohdsi)==0:
        raise ValueError("!! No Drug_dict_updated data: use 'update_drugdict' before this !!")
    else:
        path_ohdsi = sorted(path_ohdsi, reverse=True)[0]
    ohdsi = pd.read_csv(path_ohdsi, sep="\t", index_col=0)
    dic_drug = dict(zip(list(ohdsi["key"]), list(ohdsi["value"])))
    print("DONE")
    # load FAERS data
    print("prepare FAERS data", end="...")
    path_faers = glob.glob(args.workdir + SEP + "clean_*.txt")
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
    print("convert and expand records")
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
        df_tmp = pd.DataFrame(drug_rxn, columns=["active_substances", "reactions"])
        df_tmp.loc[:, "case_id"] = ids
        relation0 += df_tmp.shape[0] # num relation before NaN treatment
        df_tmp = df_tmp.dropna() # delete not found keys
        record1 += df_tmp["case_id"].nunique()
        relation1 += df_tmp.shape[0] # num relation after NaN treatment
        df_tmp.loc[:, "stored_year"] = s
        df_tmp = df_tmp.reset_index(drop=True)
        df_tmp.to_csv(outdir + SEP + f"drug_rxn_{s}_{now}.txt", sep="\t")
    summary = pd.DataFrame(
        {"count":[record0, record1, relation0, relation1]},
        index=["n_record_before", "n_record_after", "n_relation_before", "n_relation_after"]
        )
    summary.to_csv(outdir + SEP + f"count_drug_rxn_{now}.txt", sep="\t")
    print("> completed")


def make_database():
    """
    generate database from FAERS clean data
    
    """
    # init
    start = time.time()
    now = datetime.datetime.now().strftime('%Y%m%d')
    fileout = args.workdir + SEP + f"sqlite_{now}.db"
    dat = dh.DBhandler()
    dat.set_path(fileout)
    # qualification table
    print("prepare qualification table")
    tmp_filein = glob.glob(args.workdir + SEP + "clean" + SEP + "qualification_*.txt")
    if len(tmp_filein)==0:
        raise ValueError("!! No qualification table: use 'preprocess' before this !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    dtypes = {"qual_id":int, "qual_name":str}
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0, dtype=dtypes)
    dat.make_qualification_table(df)
    del df
    print("> DONE")
    # rxn_table
    print("prepare rxn table")
    tmp_filein = glob.glob(args.workdir + SEP + "curated" + SEP + f"rxn_table_*.txt")
    if len(tmp_filein)==0:
        raise ValueError("!! No rxn_table: use 'prep_drug_rxn' before this !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    dtypes = {"rxn_id":int}
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0)
    dat.make_rxn_table(df)
    del df
    print("> DONE")
    # drug table
    print("prepare drug table")
    tmp_filein = glob.glob(args.workdir + SEP + "curated" + SEP + "Drug_curated_*.txt")
    if len(tmp_filein)==0:
        raise ValueError("!! No curated drug information: use 'preprocess' before this !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    dtypes = {
        "concept_name":str, "concept_id":int, "CID":int, "CanonicalSMILES":str,
        "IUPACName":str, "MolecularFormula":str, "MolecularWeight":float, 
        "TPSA":float, "XLogP":float, "category":int,
        }
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0, dtype=dtypes)
    dat.make_drug_table(df)
    del df
    print("> DONE")
    # drug_dict
    print("prepare drug dict")
    tmp_filein = glob.glob(args.workdir + SEP + "curated" + SEP + "Drug_dict_updated_*.txt")
    if len(tmp_filein)==0:
        raise ValueError("!! No updated drug dict information: use 'preprocess' before this !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    dtypes = {
        "drug_dict_id":int, "key":str, "value":int, "representative":int,
        }
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0, dtype=dtypes)
    dat.make_drug_dict(df)
    del df
    print("> DONE")
    # case table
    print("prepare case table")
    print("<< time-consuming step >>")
    tmp_filein = glob.glob(args.workdir + SEP + "clean" + SEP + "clean_*.txt")
    if len(tmp_filein)==0:
        raise ValueError("!! No clean FAERS data: use 'preprocess' before this !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=True)[0]
    dtypes = {
        "case_id":int, "active_substances":str, "reactions":str, "sex":str,
        "event_date":int, "event_country":str, "patient_age":int, 
        "qualification":int, "stored_year":float
        }
    df = pd.read_csv(tmp_filein, sep="\t", index_col=0, dtype=dtypes)
    dat.make_case_table(df)
    del df
    print("> DONE")
    # drug_rxn_table
    print("prepare drug-rxn table")
    print("<< time-consuming step >>")
    tmp_filein = glob.glob(args.workdir + SEP + "drug_rxn" + SEP + "drug_rxn_*.txt")
    if len(tmp_filein)==0:
        raise ValueError("!! No updated drug-rxn relationship: use 'preprocess' before this !!")
    else:
        tmp_filein = sorted(tmp_filein, reverse=False)
    dtypes = {
        "drug_rxn_id":int, "case_id":int, "drug_id":int, "rxn_id":int,
        }
    for t in tqdm(tmp_filein):
        df = pd.read_csv(t, sep="\t", index_col=0, dtype=dtypes)
        dat.make_drug_rxn_table(df, if_exists="append")
        del df
    # logging
    dat._to_history(target_table="drug_rxn_table", description="newly create")
    dat.head("drug_rxn_table")
    print("> DONE")
    # summary
    print("> completed")
    elapsed = time.time() - start
    h, rem = divmod(elapsed, 3600)
    m, s = divmod(rem, 60)
    print(f"elapsed time: {int(h)} hr {int(m)} min {s:.1f} sec")

if __name__ == '__main__':
    main()     
