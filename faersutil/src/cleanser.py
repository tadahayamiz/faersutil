# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

FAERSUtil o-- Cleanser

@author: tadahaya
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from scipy import stats
import datetime
import matplotlib.pyplot as plt
from itertools import chain
import statsmodels.stats.multitest as multitest

if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

class Cleanser():
    def __init__(self):
        self.data = pd.DataFrame()
        self.n_record = tuple()
        self.exception = dict()


    def get_data(self):
        return self.data

    
    def get_record(self):
        return self.n_record


    def get_exception(self):
        return self.exception


    def set_data(self,data=None,url=""):
        """ set cleansed data """
        self.data = data
        if self.data is None:
            if len(url)==0:
                raise ValueError("!! Give data or its url !!")
            else:
                import pickle
                with open(url,"rb") as f:
                    self.data = pickle.load(f)


    def set_xml(self, data=[], url="", from_pickle=False, from_csv=False):
        """
        code for making list of data files derived from FAERS xml (output of load_xml)
            
        Parameters
        ----------
        url: str
            path of folder which contains pickle or csv files of data derived from FAERS xml
            
        from_pickle: boolean
            whether result is exported as pickle
            
        from_csv: boolean
            whether result is exported as csv
                
        """
        if len(data)==0:
            if len(url)==0:
                raise ValueError("!! No data: give a data list or a url for the directory containing data !!")
            else:
                p = Path(url)
                if from_pickle:
                    filelist = list(map(lambda x: x.as_posix(),list(p.glob("*.pkl"))))
                    if len(filelist)==0:
                        raise ValueError("!! No pickle files: check the inidicated directory !!")
                    data = [pd.read_pickle(v) for v in filelist]
                elif from_csv:
                    filelist = list(map(lambda x: x.as_posix(),list(p.glob("*.csv"))))
                    if len(filelist)==0:
                        raise ValueError("!! No csv files: check the inidicated directory !!")
                    data = [pd.read_csv(v,index_col=0,low_memory=True,dtype=str) for v in filelist]
                else:
                    raise ValueError("!! Turn on check from_pickle or from_csv !!")
        self.data = pd.concat(data)
    

    def data_cleansing(self, data=None):
        """
        data cleaning of concatenated data derived from FAERS xml files
        inspired by Scientific Data volume 3, Article number: 160026 (2016), Banda J, et al.
        
        Parameters
        ----------
        data: dataframe
            raw FAERS data after load_xml()

        """
        if data is not None:
            self.data = data
        n_all = self.data.shape[0]

        # check 'Event Date', etc.
        self.data = self.data.drop_duplicates(keep='first')
        
        # check 'case ID'
        self.data = self.data.drop_duplicates(subset=["Case ID"], keep='first')
        self.data = self.data.reset_index(drop=True)
        self.data.loc[:, "Active Substances"] = self.data.loc[:, "Active Substances"].map(lambda x: set(x.split(";")))
        self.data.loc[:, "Reactions"] = self.data.loc[:, "Reactions"].map(lambda x: set(x.split(";")))
        
        # treat Active Substances with hard coding ("or")
        self.data.loc[:, "Active Substances"] = self.data.loc[:, "Active Substances"].map(treat_or) # treat only "or"
        
        # align
        empty_set = {""}
        self.data.loc[:, "Active Substances"] = self.data.loc[:, "Active Substances"].map(lambda x: x - empty_set)
        self.data = self.data[self.data["Active Substances"].map(lambda x: len(x) > 0)]
        self.n_record = (n_all,self.data.shape[0])

        # arrange type, mainly due to np.nan
        dic = {"Female":"F", "Male":"M", np.nan:"N"}
        self.data.loc[:, "Sex"] = self.data.loc[:, "Sex"].map(lambda x: dic[x])
        self.data = self.data.fillna({"Event Country":"N", "Patient Age":0, "Qualification":"N"})
        self.data.loc[:, "Patient Age"] = self.data.loc[:, "Patient Age"].astype(int)

        # correct date
        self._correct_date()


    def set_exception(self):
        """ load exception lists """
        if len(self.n_record)==0:
            self.n_record = (self.data.shape[0],self.data.shape[0])
        p = Path(__file__.replace(f"cleanser.py", "exception_list"))
        paths = list(map(lambda x: x.as_posix(),list(p.glob("*.txt"))))
        if len(paths) > 0:
            for v in paths:
                temp = pd.read_csv(v,sep="\t")
                fname = v.split(SEP)[-1].replace(".txt","")
                self.exception[fname] = set(temp[fname])


    def exclude_exception(self, chem_only=True):
        """
        exclude exception based on not found drugs in PubChem
        keep chemicals only, or chemicals and proteins
        
        Parameters
        ----------
        chem_only: boolean
            exclude all records that are not chemicals
        
        """
        if len(self.exception)==0:
            print("!! No exception: check exception or set_exception() before this process !!")
        else:
            exceptions = self.exception["not_found"]
            if chem_only:
                pass
            else:
                for k,v in self.exception.items():
                    if k!="not_found":
                        exceptions -= v
            fxn = lambda x: x - exceptions
            self.data.loc[:, "Active Substances"] = self.data.loc[:, "Active Substances"].map(fxn)
            empty_set = {""}
            self.data.loc[:, "Active Substances"] = self.data.loc[:, "Active Substances"].map(lambda x: x - empty_set)
            self.data = self.data[self.data["Active Substances"].map(lambda x: len(x) > 0)]
            self.n_record = (self.n_record[0], self.data.shape[0])


    def _correct_date(self, key_date:str="Event Date"):
        """
        correct date
        Event Date is in {NaN, YYYY, YYYYMM, YYYYMMDD}
        Align all date to int YYYYMMDD format or 0 (NaN)
        Convert YYYY to YYYY0701, YYYYMM to YYYYMM15, NaN to 0

        """
        # replace NaN and set int
        self.data.loc[:, key_date] = self.data.loc[:, key_date].replace(np.nan, 0).astype(int)
        # add Order for ease
        self.data.loc[:, "Order"] = self.data[key_date].map(lambda x: len(str(x)))
        # convert YYYY
        self.data.loc[self.data["Order"]==4, key_date] = self.data[key_date].map(lambda x: x * 10000 + 701)
        # convert YYYYMM
        self.data.loc[self.data["Order"]==6, key_date] = self.data[key_date].map(lambda x: x * 100 + 15)
        del self.data["Order"]


def treat_or(x):
    """ processing drugs containing ' or ' """
    return set(chain.from_iterable(v.split(" or ") for v in x))

def treat_and(x):
    """ processing drugs containing ' and ' """
    return set(chain.from_iterable(v.split(" and ") for v in x))

def treat_and_or(x):
    """ processing drugs containing ' and or ' """
    return set(chain.from_iterable(v.split(" and") for v in x))