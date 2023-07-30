# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

utilities for chemical data handling

@author: tadahaya
"""
import json
import pandas as pd
import numpy as np
import re
import scipy.stats as st
from itertools import chain
import pubchempy as pcp
from pathlib import Path
import os
import csv

from .converter import SynoDict
from .info_handler.chem_info import Info

if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

class Chem():
    def __init__(self):
        self.info = Info()
        self.ref = None # an SynoDict object
        self.encoder = dict() # a dict for key to value conversion
        self.decoder = dict() # a dict for value to key conversion
        self.not_found = dict() # a dict of chemicals not found in PubChem

    ### info
    def set_info(self,data=None,url="",col="Synonym",sep=";",use_stored=True):
        """
        obtains compound information from a local file like "chemical-info.txt"
        
        Parameters
        ----------
        data: dataframe
            chemical-info_XXXX.txt

        url: str
            a path of chemical-info_XXXX.txt

        col: str
            indicates the columns of synonyms

        sep: str
            indicates the separator keyword for set conversion

        use_stored: boolean
            whether the stored chemical-info_whole data is employed or not

        """
        if use_stored:        
            paths = self.info.get_stored_path()
            url0 = sorted([v for v in paths if "chemical-info_whole" in v])[-1]
            self.info.load_data(url=url0)
            if (data is not None) or (len(url) > 0):
                temp = Info()
                temp.load_data(data,url)
                self.info.combine(data_list=[temp.get_data()])
        else:
            self.info.load_data(data,url)
        self.info.update_info(replace=True)


    def get_info(self,chem_list,namespace="name",delete=[],replace=[]):
        """
        get pubchem information of the chemicals in the given list

        Parameters
        ----------
        chem_list: list
            a list of the chemicals of interest

        namespace: str
            indicates the type of chem_list element, "name" or "cid"

        delete,replace: list
            indicate additional elements to be treated in deletion or replacement modifications

        """
        return self.info.get_info(chem_list,namespace,delete,replace)


    def update_info(self,replace=False,whole_update=False,col="Synonym",sep=";"):
        """
        update chemical-info

        Parameters
        ----------
        replace: boolean
            indicates whether data is replaced with temp data

        whole_update: boolean
            indicates whether chemical-info_whole.txt is updated

        """
        self.info.update_info(replace=replace)
        if whole_update:
            paths = self.info.get_stored_path()
            url = sorted([v for v in paths if "chemical-info_whole" in v])[-1]
            temp = Info()
            temp.load_data(url=url)
            self.info.combine(data_list=[temp.temp])
            self.info.update_info()
            self.info.export(url.split("chemical-info_whole")[0] + "chemical-info_whole.txt",col,sep)


    ### name integration
    def set_ref(self,ref=None,value_col="CID",synonym_col="Synonym"):
        """
        set reference data
        note reference data should be a SynoDict object

        Parameters
        ----------
        ref: SynoDict
            a SynoDict object derived from a chemical-info file

        """
        if ref is None:
            raise ValueError("!! No prepared reference: prep_ref() before this process !!")
        else:
            if type(ref)==type(SynoDict()):
                self.ref = ref
            else:
                raise TypeError("!! Wrong type: ref should be SynoDict !!")


    def prep_ref(self,info=None,value_col="CID",synonym_col="Synonym"):
        """
        prepare reference data

        Parameters
        ----------
        info: dataframe
            a dataframe of info file
            should have been processed by str2set()

        """
        if info is None:
            info = self.info.get_data()
        else:
            col = list(info.columns)
            if (value_col not in col) or (synonym_col not in col):
                raise KeyError("!! Wrong column names: check column names of info and each col argument !!")
            if type(info[synonym_col][0])==str:
                temp = Info()
                temp.load_data(data=info)
                info = temp.get_data()
        self.ref = SynoDict(keys=list(info.index),values=list(info[value_col]),synonyms=list(info[synonym_col]),processing=True)


    def prep_identifier(self,chem_list,ref=None):
        """
        prepare identifier that achieves name identification based on a given list

        Parameters
        ----------
        chem_list: list
            a list of the chemicals of interest

        ref: SynoDict
            a SynoDict object derived from a chemical-info file

        """
        if ref is None:
            if self.ref is None:
                raise ValueError("!! Give a reference (SynoDict) or prep_ref() before this process !!")
        else:
            self.ref = ref
        self.encoder = self.ref.fix(chem_list)
        key = list(self.encoder.keys())
        val = list(self.encoder.values())
        val2 = sorted(list(set(val)))
        key2 = []
        for v in val2:
            temp = [i for i,w in enumerate(val) if w==v]
            temp2 = [key[i] for i in temp]
            temp3 = [len(k) for k in temp2]
            key2.append(temp2[np.argmin(temp3)])
        self.decoder = dict(zip(val2,key2)) # the shortest key is selected
        self.not_found = self.ref.get_not_found()

        # export
        temp1 = pd.DataFrame({"name":key,"CID":val})
        url1 = __file__.replace("name_identifier.py",f"dwh{SEP}identifier.txt")
        temp1.to_csv(url1,sep="\t")
        temp2 = pd.DataFrame({"not_found":list(self.not_found.keys()),"CID":list(self.not_found.values())})
        url2 = __file__.replace(f"identifier{SEP}name_identifier.py","exception_list{SEP}not_found.txt")
        temp2.to_csv(url2,sep="\t")


    def set_identifier(self,data=None,url=""):
        """ set identifier dict """
        if data is None:
            if len(url)==0:
                raise ValueError("!! Give data or its url !!")
            else:
                data = pd.read_csv(url,sep="\t",index_col=0)
        key = list(data["name"])
        val = list(data["CID"])
        self.encoder = dict(zip(key,val))
        self.decoder = dict(zip(val,key))


    def get_not_found(self):
        """ get not found dictionary """
        return self.not_found


    def encode(self,chem_list):
        """
        name identification by converting names to IDs

        Parameters
        ----------
        chem_list: list
            a list of the chemicals of interest

        """
        return [self.encoder[v] for v in chem_list]


    def encode_set(self,chem_set):
        """
        name identification by converting names to IDs

        Parameters
        ----------
        chem_set: set
            a set of the chemicals of interest

        """
        return {self.encoder[v] for v in chem_set}


    def decode(self,chem_list):
        """
        name identification by converting names to IDs

        Parameters
        ----------
        chem_list: list
            a list of the chemicals of interest

        """
        return [self.decoder[v] for v in chem_list]


    def decode_set(self,chem_set):
        """
        name identification by converting names to IDs

        Parameters
        ----------
        chem_set: set
            a set of the chemicals of interest

        """
        return {self.decoder[v] for v in chem_set}
