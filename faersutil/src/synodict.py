# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

dict handling Synonyms

@author: tadahaya
"""
import pandas as pd

class SynoDict():
    def __init__(self):
        self.encoder = dict() # encode many to one value
        self.decoder = dict() # decode to the representatives
    

    def set_encoder(self, encoder:dict):
        self.encoder = encoder


    def set_decoder(self, decoder:dict):
        self.decoder = decoder


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
        # prep
        if len(value)==0:
            encoder = dict(zip(list(df[key]), list(range(df.shape[0]))))
        else:
            encoder = dict(zip(list(df[key]), list(df[value])))
        decoder = dict(zip(encoder.values(), encoder.keys()))
        # update
        self.encoder.update(encoder)
        self.decoder.update(decoder)
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
        # prep
        if len(values)==0:
            encoder = dict(zip(keys, list(range(len(keys)))))
        else:
            encoder = dict(zip(keys, values))
        decoder = dict(zip(encoder.values(), encoder.keys()))
        # update
        self.encoder.update(encoder)
        self.decoder.update(decoder)
        # add synonyms
        for k, v, s in zip(keys, values, synonyms):
            tmp = s.split(sep)
            tmp_dic = dict(zip(tmp, [v] * len(tmp)))
            self.encoder.update(tmp_dic)
            # ToDo: add a system to check duplicates when updated


    def to_csv(self, fileout:str="", sep:str="\t"):
        """ export encoder as a df """
        assert len(fileout) > 0
        df = pd.DataFrame({
            "key":list(self.encoder.keys()),
            "value":list(self.encoder.values())
            })
        # add representative to judge representative or synonyms
        df.loc[:, "representative"] = 0
        df.loc[list(self.decoder.values()), "representative"] = 1
        # dtypes
        df.loc[:, "key"] = df.loc[:, "key"].astype(str)
        df.loc[:, "value"] = df.loc[:, "value"].astype(int)
        df.loc[:, "representative"] = df.loc[:, "representative"].astype(int)
        # export
        df.to_csv(fileout, sep=sep)