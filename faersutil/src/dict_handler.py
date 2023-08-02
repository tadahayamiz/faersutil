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

if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")


class SynonymDict():
    def __init__(self):
        self.encoder = dict() # encode many to one
        self.decoder = dict() # decode to the representatives
    
    def get_encoder(self):
        return self.encoder

    def get_decoder(self):
        return self.decoder

    def from_df(
        self, df:pd.DataFrame, key:str="", value:str="", synonym:str="", sep:str="///"
        ):
        """
        prepare SynonymDict from dataframe with the following columns:
        names, which indicate the representative keys
        values, which indicate the IDs corresponding to the keys
        synonyms, which indicate the Synonyms of the representatives

        Parameters
        ----------
        df: pd.DataFrame
        
        key: str
            the column name for the representative keys

        value: str
            the column name for the ID values corresponding to the keys
            if it's not given, integers from 0 will be employed

        synonym: str
            the column name for the synonyms
            each element is str bridged with the indicated separator as sep
        
        """
        if (len(key)==0) | (len(synonym)==0):
            raise ValueError("!! Provide keys or synonyms !!")
        assert (key in df.columns) & (synonym in df.columns)
        if len(value)==0:
            self.encoder = dict(zip(list(df[key]), list(range(df.shape[0]))))
        else:
            self.encoder = dict(zip(list(df[key]), list(df[value])))
        self.decoder = dict(zip(self.encoder.values(), self.encoder.keys()))
        # add synonyms
        df = df[[key, value, synonym]] # sort
        for v in df.values:
            tmp = v[2].split(sep)
            tmp_dic = dict(zip(tmp, [v[1]] * len(tmp)))
            self.encoder.update(tmp_dic)

    def from_lists(
        self, keys:list=[], values:list=[], synonyms:list=[], sep:str="///"
        ):
        """
        prepare SynonymDict from dataframe with the following columns:
        names, which indicate the representative keys
        values, which indicate the IDs corresponding to the keys
        synonyms, which indicate the Synonyms of the representatives

        Parameters
        ----------        
        keys: list
            contain the representative keys

        values: list
            contain the ID values corresponding to the keys

        synonyms: list
            contain the snynonyms of the representative keys
            each element is str bridged with the indicated separator as sep
        
        """
        assert len(keys)==len(synonyms)
        assert (type(keys)==list) & (type(synonyms)==list)
        if len(values)==0:
            self.encoder = dict(zip(keys, list(range(len(keys)))))
        else:
            self.encoder = dict(zip(keys, values))
        self.decoder = dict(zip(self.encoder.values(), self.encoder.keys()))
        # add synonyms
        for k, v, s in zip(keys, values, synonyms):
            tmp = s.split(sep)
            tmp_dic = dict(zip(tmp, [v] * len(tmp)))
            self.encoder.update(tmp_dic)


class OHDSIhandler():
    def __init__(self):
        pass



import pandas as pd

df = pd.read_csv(r"D:\OHDSI\CONCEPT.csv", sep="\t") # tsvだったので注意
print(df.shape)
df.head()

df2 = df[df["domain_id"]=="Drug"]
print(df2.shape, df.shape)
df2.head()

df3 = df2[df2["concept_class_id"]=="Ingredient"]
print(df3.shape, df2.shape)
df3.head()

df3 = df3.reset_index(drop=True)
df3.to_csv(r"D:\OHDSI\CONCEPT_drug_ingredient_230731.txt", sep='\t')