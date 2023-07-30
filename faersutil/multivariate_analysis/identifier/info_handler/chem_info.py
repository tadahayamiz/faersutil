# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

Chem o-- Info

Note: keys obtained with pubchempy were changed ~2020.

@author: tadahaya
"""
import json
import pandas as pd
import numpy as np
import scipy.stats as st
from itertools import chain
import pubchempy as pcp
from pathlib import Path
import csv
import re
from tqdm import tqdm
import sys

from .info_loader import InfoLoader
from .info_editor import Editor

EXCLUDE_KEYS = {""," "," Certified Reference Material"}
KEYS = ["CID","CanonicalSMILES","IUPACName","MolecularFormula","MolecularWeight","Synonym","TPSA","XLogP"]

class Info():
    def __init__(self,url="",extension="txt",key="chemical-info"):
        self.__iloader = InfoLoader(url,extension,key)
        self.__editor = Editor()
        self.temp = pd.DataFrame(columns=KEYS)
        self.nega = []
        self.data = pd.DataFrame(columns=KEYS)


    def get_data(self):
        """ return data """
        return self.data


    def get_nega(self):
        """ return nega """
        return self.nega


    def combine(self,data_list=[],exclude_temp=False):
        """ combine info data """
        if len(data_list)==0:
            raise ValueError("!! Give a list of information dataframes !!")
        if not exclude_temp:
            data_list += [self.temp]
        merge = pd.concat(data_list,axis=0,join="inner",sort=False)
        self.temp = self.__iloader.aggregate(merge)


    def export(self,fileout="",col="Synonym",sep=";"):
        """ export data """
        if self.data.shape[0]==0:
            raise ValueError("!! No data is stored: update_info before exporting data !!")
        if len(fileout)==0:
            raise ValueError("!! No fileout: give a path for fileout !!")
        temp = self.__iloader.set2str(self.data,col,sep)
        temp.to_csv(fileout,sep="\t")


    def get_stored_path(self):
        """ return stored file paths """
        return self.__iloader.get_stored_path()


    def load_stored(self,include=[],exclude=[],col="Synonym",sep=";"):
        """
        load stored info data
        
        Parameters
        ----------
        include,exclude: list
            keywords to be included and excluded from the stored data, respectively

        col: str
            indicates the columns of synonyms

        sep: str
            indicates the separator keyword for set conversion

        """
        if (len(include) > 0) or (len(exclude)) > 0:
            self.__iloader.load_stored(include,exclude)
        self.temp = self.__iloader.get_data()


    def load_data(self,data=None,url="",col="Synonym",sep=";"):
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

        """
        self.__iloader.load_data(data,url)
        self.temp = self.__iloader.get_data()


    def update_info(self,replace=False):
        """
        update chemical-info
        
        """
        if self.temp.shape[1]==0:
            raise ValueError("!! temp is empty !!")
        if self.data.shape[0]==0:
            self.data = self.temp
        else:
            if replace:
                self.data = self.temp
            else:
                merge = pd.concat([self.data,self.temp],axis=0,join="inner",sort=False)
                self.data = self.__iloader.aggregate(merge)


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
        # 1st raw
        print("------ 1st run ------")
        dat1,self.nega = self.pull_info(chem_list,namespace)
        dat2 = pd.DataFrame(columns=dat1.columns)
        dat3 = pd.DataFrame(columns=dat1.columns)

        # 2nd delete
        if len(self.nega) > 0:
            print("------ 2nd run ------")
            self.__editor.to_delete()
            if len(delete) > 0:
                self.__editor.add(delete)
            dat2,self.nega = self.pull_info(self.nega,namespace)

            # 3rd replacement
            if len(self.nega) > 0:        
                print("------ 3rd run ------")
                self.__editor.to_replace()
                if len(replace) > 0:
                    self.__editor.add(replace)
                dat3,self.nega = self.pull_info(self.nega,namespace)

        self.temp = pd.concat([dat1,dat2,dat3],join="inner",axis=0)
        self.temp = self.__iloader.aggregate(self.temp)
        if self.temp.empty:
            self.temp = pd.DataFrame(columns=KEYS)
        return self.temp,self.nega


    def pull_info(self,chem_list,namespace="name"):
        """ obtain chemical information """
        n = len(chem_list)
        prop = ['iupacname', 'molecularformula', 'molecularweight', 'xlogp', 'tpsa', 'canonicalsmiles']
        res = []
        posi = []
        nega = []
        ap = res.append
        ap2 = posi.append
        ap3 = nega.append
        if (namespace=="name") or (namespace=="cid"):
            for v in tqdm(chem_list):
                edited = self.__editor.edit(v) # word editing
                if len(edited) > 0: # edit
                    for w in edited:
                        try:
                            x = pcp.get_properties(prop,w,namespace)
                            if len(x) > 0:
                                temp = x[0] # [dict(hit1),dict(hit2),...]
                                z = pcp.get_synonyms(temp["CID"],namespace="cid")
                                syno0 = z[0]["Synonym"] # list
                                syno0 = self.del_key(syno0)
                                syno = {v} | set(syno0)
                                syno = {s.lower() for s in syno}
                                temp["Synonym"] = syno
                                if "XLogP" not in list(temp.keys()):
                                    temp["XLogP"] = np.nan
                                ap([temp[k] for k in KEYS])
                                ap2(v)
                                break
                        except KeyboardInterrupt:
                            print("!! Forced termination !!")
                            sys.exit()
                        except Exception as e:
                            print("!! an ERROR: {} !!".format(e.args))
                            break
                    else:
                        ap3(v)
                else:
                    ap3(v)
        else:
            raise ValueError("!! Wrong namespace: choose 'name' or 'cid' !!")
        print("> extracted {0}/{1} compounds".format(len(posi),n))
        if len(res)==0:                
            return pd.DataFrame(columns=KEYS),nega
        else:
            data = pd.DataFrame(res,columns=KEYS,index=posi)
            return data,nega


    def del_key(self,lst):
        """ exclude elements indicated by the excluded set """
        return [v for v in lst if v not in EXCLUDE_KEYS]