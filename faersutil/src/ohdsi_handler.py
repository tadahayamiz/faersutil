# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

handling CONCEPT.csv from OHDSI
https://athena.ohdsi.org/search-terms/start

< Note >
- concept_name of Ingredient has duplicates while concept_id is unique

@author: tadahaya
"""
import pandas as pd
import numpy as np

from . import pcp_handler as ph

class OHDSIhandler():
    def __init__(self):
        self.url = ""
        self.df = None
        self.pubchem = None


    def set_path(self, url:str=""):
        """ set the path for CONCEPT.csv """
        if len(url)==0:
            raise ValueError("!! Provide url for CONCEPT.csv !!")
        self.url = url


    def load_df(self, df:pd.DataFrame=None):
        """ load CONCEPT.csv as a df """
        if df is None:
            self.df = pd.read_csv(self.url, sep="\t", index_col=0)
        else:
            self.df = df


    def load_pubchem(self, df:pd.DataFrame=None):
        """ load PubChem search result a df """
        if df is None:
            raise ValueError("!! Provide PubChem result as a df !!")
        else:
            self.df = df


    def extract_ingredient(self, fileout:str=""):
        """ extract ingredient from raw CONCEPT.csv """
        self.df = self.df[self.df["domain_id"]=="Drug"] # restrict Drug domain
        self.df = self.df[self.df["concept_class_id"]=="Ingredient"]
        self.df = self.df.drop_duplicates(subset=["concept_name"], keep="first")
        self.df = self.df.reset_index(drop=True)
        if len(fileout) > 0:
            self.df.to_csv(fileout, sep="\t")

    
    def get_pubchem(
            self, chem_list:list=[], key:str="concept_name", gap:str="///",
            fileout:str="", sep:str="\t", namespace:str="name"
            ):
        """ get pubchem information """
        if len(chem_list)==0:
            chem_list = list(self.df["concept_name"])
        res, nega = ph.main(
            chem_list=chem_list, key=key, gap=gap, 
            fileout=fileout, sep=sep, namespace=namespace
            )
        self.pubchem = res
        if len(fileout) > 0:
            self.pubchem.to_csv(fileout, sep="\t")


    def integrate_pubchem(
            self, pubchem:pd.DataFrame=None, fileout:str="",
            def_small:dict={
                "min_len_SMILES":4, "min_len_MolecularFormula":3,
                "min_MolecularWeight":40, "max_MolecularWeight":800,
                }
            ):
        """
        search ingredients in PubChem
        and integrate the results with concept_ingredient
        note we categorize the compounds as follows:
        - 0: all compounds (OHDSI ingredients)
        - 1: PubChem (positive in PubChem search, having SMILES)
        - 2: small molecules (based on MW etc.)

        """
        # preparation
        chem_list = list(self.df["concept_name"])
        ohdsi_id = list(self.df["concept_id"])
        dic = dict(zip(chem_list, ohdsi_id))
        if self.pubchem is None:
            if pubchem is None:
                raise ValueError(
                    "!! Provide pubchem result as an argument or run get_pubchem before this !!"
                    )
        # integration
        self.pubchem.loc[:, "concept_id"] = self.pubchem.loc[:, "concept_name"].map(lambda x: dic[x])
        remains = [v for v in chem_list if v not in list(self.pubchem["concept_name"])]
        remains_id = [dic[v] for v in remains]
        tmp = pd.DataFrame({"concept_name":remains, "concept_id":remains_id})
        self.df = pd.concat([self.pubchem, tmp], axis=0, join="outer")
        self.df = self.df.reset_index(drop=True)
        self.df = self.df.fillna(
            {
                "CID":0, "CanonicalSMILES":"N", "IUPACName":"N",
                "MolecularFormula":"N", "Synonym":"N"
                }
            )
        # align dtypes
        self.df.loc[:, "concept_name"] = self.df.loc[:, "concept_name"].astype(str)
        self.df.loc[:, "CanonicalSMILES"] = self.df.loc[:, "CanonicalSMILES"].astype(str)
        self.df.loc[:, "IUPACName"] = self.df.loc[:, "IUPACName"].astype(str)
        self.df.loc[:, "MolecularFormula"] = self.df.loc[:, "MolecularFormula"].astype(str)
        self.df.loc[:, "Synonym"] = self.df.loc[:, "Synonym"].astype(str)
        self.df.loc[:, "concept_id"] = self.df.loc[:, "concept_id"].astype(int)
        self.df.loc[:, "CID"] = self.df.loc[:, "CID"].astype(int)
        self.df.loc[:, "MolecularWeight"] = self.df.loc[:, "MolecularWeight"].astype(float)
        self.df.loc[:, "TPSA"] = self.df.loc[:, "TPSA"].astype(float)
        self.df.loc[:, "XLogP"] = self.df.loc[:, "MolecularWeight"].astype(float)
        # add category
        ## all compounds
        self.df.loc[:, "category"] = 0
        ## PubChem
        self.df.loc[self.df["CID"] > 0, "category"] = 1
        ## smal molecules
        tmp = self.df.copy()
        tmp.loc[:, "len_SMILES"] = tmp["CanonicalSMILES"].map(lambda x: len(x))
        tmp.loc[:, "len_MF"] = tmp["MolecularFormula"].map(lambda x: len(x))

        for t in list(tmp["MolecularWeight"]):
            print(type(t), t)

        tmp = tmp[tmp["MolecularWeight"] >= def_small["min_MolecularWeight"]]
        tmp = tmp[tmp["MolecularWeight"] <= def_small["max_MolecularWeight"]]
        tmp = tmp[tmp["len_SMILES"] > def_small["min_len_SMILES"]]
        tmp = tmp[tmp["len_MF"] > def_small["min_len_MolecularFormula"]]
        self.df.loc[tmp.index, "category"] = 2
        del tmp
        if len(fileout) > 0:
            self.df.to_csv(fileout, sep="\t")