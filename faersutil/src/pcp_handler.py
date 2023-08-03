# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

obtain chemical information with pubchempy

requirement:
- pubchempy
- tqdm

@author: tadahaya
"""
import os
import pandas as pd
import numpy as np

from tqdm.auto import tqdm, trange

import pubchempy as pcp

EXCLUDE_KEYS = {""," "," Certified Reference Material"}
KEYS = [
    "CID",
    "CanonicalSMILES",
    "IUPACName",
    "MolecularFormula",
    "MolecularWeight",
    "Synonym",
    "TPSA",
    "XLogP"
    ]

if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

def _del_key(lst:list):
    """ exclude elements indicated by the excluded set """
    return {v for v in lst if v not in EXCLUDE_KEYS}

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

def pull_info(
        chem_list:list, key:str="concept_name", gap:str="///", namespace:str="name"
        ):
    """
    obtain chemical information using pubchempy

    Parameters
    ----------
    chem_list: list
        a list of chemicals of analyzed
    
    key: str
        indicates the column name for the given chemical list in the output
    
    gap: str
        the gap that separates each element in Synonym

    namespace: str
        indicates how to indicate the chemicals in the given list
        name or cid can be indicated
    
    """
    assert (namespace=="name") | (namespace=="cid")
    n = len(chem_list)
    prop = ['iupacname', 'molecularformula', 'molecularweight', 'xlogp', 'tpsa', 'canonicalsmiles']
    res = []
    posi = []
    nega = []
    not_yet = []
    # pubchempy
    print("=== Search with PubChempy ===")
    for w in tqdm(chem_list):
        try:
            x = pcp.get_properties(prop, w, namespace)
            if len(x) > 0:
                tmp = x[0] # [dict(hit1),dict(hit2),...]
                z = pcp.get_synonyms(tmp["CID"], namespace="cid")
                if len(z)==0:
                    syno = set()
                else:
                    syno = z[0]["Synonym"] # list
                    syno = _del_key(syno) # returned as set
                    syno = {w} | syno
                    syno = {s.lower() for s in syno}
                tmp["Synonym"] = syno
                if "XLogp" not in set(tmp.keys()):
                    tmp["XLogP"] = np.nan
                res.append(tmp) # add a dict
                posi.append(w)
            else:
                nega.append(w)
        except KeyboardInterrupt:
            print("!! Forced termination !!")
            print("> save the current records")
            not_yet = [v for v in chem_list if v not in posi + nega]
            break
        except Exception as e:
            print("!! an ERROR: {} !!".format(e.args))
            print("> save the current records")
            not_yet = [v for v in chem_list if v not in posi + nega]
            break
    print("> extracted {0}/{1} compounds".format(len(posi), n))
    if len(res)==0:                
        return pd.DataFrame(columns=[key] + KEYS), nega, not_yet 
    else:
        # remove overlapping synonyms
        print("=== Remove overlapping Synonyms ===")
        num_posi = len(posi)
        for i in trange(num_posi):
            for j in range(i + 1, num_posi):
                inter = res[i]["Synonym"] & res[j]["Synonym"]
                if len(inter) > 0:
                    res[i]["Synonym"] -= inter
                    res[j]["Synonym"] -= inter
        data = pd.DataFrame(res)
        data.loc[:, key] = posi
        data.reset_index(drop=True)
        data = data[[key] + KEYS]
        # set to str
        data.loc[:, "Synonym"] = data.loc[:, "Synonym"].map(list).apply(_concat, gap=gap)
        return data, nega, not_yet


def main(
        chem_list:list, key:str="concept_name", gap:str="///",
        fileout:str="", sep:str="\t", namespace:str="name"
        ):
    """
    repeat pull_info until all chemicals in the given list are searched

    Parameters
    ----------
    chem_list: list
        a list of chemicals of analyzed

    key: str
        indicates the column name for the given chemical list in the output
    
    gap: str
        the gap that separates each element in Synonym

    namespace: str
        indicates how to indicate the chemicals in the given list
        name or cid can be indicated
    
    """
    if len(chem_list) != len(set(chem_list)):
        raise ValueError(
            "!! chem_list has duplicated names: remove the duplicates !!"
            )
    not_yet = chem_list.copy()
    results = []
    i = 0
    while len(not_yet) > 0:
        prev = len(not_yet)
        res, nega, not_yet = pull_info(
            not_yet, key=key, gap=gap, namespace=namespace
            )
        if prev==len(not_yet):
            print("!! NO PROGRESS: Check Bug !!")
            break
        i += 1
        results.append(res)
    results = pd.concat(results, axis=0, join="inner")
    if len(fileout) > 0:
        results.to_csv(fileout, sep=sep)
    return results, nega