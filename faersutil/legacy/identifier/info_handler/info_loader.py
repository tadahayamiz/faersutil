# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

Chem o-- Info o-- InfoLoader

@author: tadahaya
"""
import json
import pandas as pd
import numpy as np
import re
from itertools import chain
from collections import defaultdict
import pubchempy as pcp
from pathlib import Path
import csv
import os

if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

class InfoLoader():
    """ data loader for Info """
    def __init__(self,url="",extension="txt",key="chemical-info"):
        self.__data = pd.DataFrame()
        if len(url)==0:
            url = __file__.replace(f"info_handler{SEP}info_loader.py","dwh")
        self.stored_path = []
        self.set_stored_path(url,extension,key)
        if len(self.stored_path)==0:
            print("!! Notice no stored chemical-info files !!")


    def set_stored_path(self,url="",extension="txt",key="chemical-info"):
        """ set stored file paths """
        p = Path(url)
        self.stored_path = [v for v in list(map(lambda x: x.as_posix(),list(p.glob("*.{}".format(extension))))) if key in v]


    def get_stored_path(self):
        """ return stored file paths """
        return self.stored_path


    def get_data(self):
        return self.__data


    def load_stored(self,include=[],exclude=[],col="Synonym",sep=";"):
        """ load stored data """
        lst = []
        for v in self.stored_path:
            flag = False
            for i in include:
                if i in v:
                    flag = True
                    break
            if flag:
                for e in exclude:
                    if e in v:
                        flag = False
                        break
            if flag:
                lst.append(v)
        temp = [pd.read_csv(v,sep="\t",index_col=0) for v in list(set(lst))]
        if len(temp)!=0:
            self.__data = pd.concat(temp,join="inner",axis=0)
            self.__data = self.aggregate(self.__data,col_synonym=col,sep=sep)
        

    def load_data(self,data=None,url="",col="Synonym",sep=";"):
        """
        obtains compound information from a local file like "chemical-info.txt"
        
        Parameters
        ----------
        data: dataframe
            chemical-info_XXXX.txt

        url: str
            a path of chemical-info_XXXX.txt

        """
        self.__data = data
        if self.__data is None:
            if len(url) > 0:
                self.__data = pd.read_csv(url,sep="\t",index_col=0).fillna("")
            else:
                raise ValueError("!! Indicate data or url for loading !!")
        self.__data = self.__data.fillna("")
        self.__data = self.aggregate(self.__data)


    def str2set(self,data,col="Synonym",sep=";"):
        """ convert string separated by sep into set """
        data[col] = data[col].map(lambda x: set(x.split(sep)))
        return data


    def set2str(self,data,col="Synonym",sep=";"):
        """ convert string separated by sep into set """
        
        def func(temp):
            temp = list(temp)
            res = temp[0]
            for v in temp[1:]:
                res += sep + v
            return res

        data[col] = data[col].map(func)
        return data


    def aggregate(self,data,col_cid="CID",col_synonym="Synonym",sep=";"):
        """
        aggregate records with the same CID
        Note that synonyms of the outcomes are set
        
        """
        dat = data.copy()
        if type(dat[col_synonym][0])==str:
            dat[col_synonym] = dat[col_synonym].map(lambda x: x.lower())
            dat = self.str2set(data=dat,col=col_synonym,sep=sep)
        cid = list(dat[col_cid])
        dup = [x for x in set(cid) if cid.count(x) > 1]
        if len(dup)==0:
            return dat
        else:
            remain = []
            ap = remain.append
            for d in dup:
                temp = dat[dat[col_cid]==d]
                syno = set(chain.from_iterable(temp[col_synonym]))
                idx = list(temp.index)
                n_idx = [len(v) for v in idx]
                ap(temp.iloc[np.argmin(n_idx),:])
            fxn = lambda x: x not in dup
            remain = pd.concat(remain,join="inner",axis=1).T
            single = dat[dat[col_cid].map(fxn)]
            return pd.concat([single,remain],join="inner",axis=0)
