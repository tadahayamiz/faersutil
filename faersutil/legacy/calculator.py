# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

FAERSUtil o-- Calculator

@author: tadahaya
"""

import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
import os
import time
from pathlib import Path
from scipy import stats
import datetime
import matplotlib.pyplot as plt
from itertools import chain
import statsmodels.stats.multitest as multitest
import pickle
from tqdm import tqdm

MEDDRA_TERM = ["SOC","HLGT","HLT","PT"]

class Calculator():
    def __init__(self):
        self.data = pd.DataFrame()
        self.rxn = []
        self.drug = []
        self.meddra = pd.DataFrame()
        self.focused_pt = set()
        self.res = pd.DataFrame()
        self.__stats = CommonFxn()


    def set_data(self,data=None,url=""):
        """ setter """
        self.data = data
        if self.data is None:
            if len(url)==0:
                raise ValueError("!! Give data or its url !!")
            else:
                with open(url,"rb") as f:
                    self.data = pickle.load(f)
        self.drug = sorted(set(chain.from_iterable(self.data["Active Substances"])))
        self.rxn = sorted(set(chain.from_iterable(self.data["Reactions"])))


    def set_meddra(self,data=None,url="",ignore=True):
        """
        setter for MedDRA data
        
        Parameters
        ----------
        ignore: boolean
            whether the terms not related to drug induced events are removed or not
        
        """
        self.meddra = data
        if self.meddra is None:
            if len(url)==0:
                temp = __file__.replace("calculator.py","")
                p = Path(temp)
                fpath = [v.as_posix() for v in p.glob("*.txt") if "MedDRA" in v.as_posix()]
                if len(fpath) > 0:
                    url = fpath[-1]
                else:
                    raise ValueError("!! Give data or url for MedDRA data !!")
            self.meddra = pd.read_csv(url) # edit on 210105
    
        col = list(self.meddra.columns)
        for v in col:
            self.meddra[v] = self.meddra[v].map(lambda x: x.lower())
        if ignore:
            # ignore the following:
                # "SOC: congenital, familial and genetic disorders"
                # "SOC: product issues"
                # all terms containing "congenital"
            med = self.meddra
            med = med[med["SOC"]!="congenital, familial and genetic disorders"]
            med = med[med["SOC"]!="product issues"]
            for v in list(med.columns):
                med = med[~med[v].str.contains("congenital")]
            self.meddra = med


    def get_data(self):
        """ getter """
        return self.data


    def check_drug(self,interest=[]):
        """ check whether drugs related to the indicated keywords """
        if len(self.drug)==0:
            raise ValueError("!! Set data before this process !!")
        if type(interest)==str:
            temp = [v for v in self.drug if interest in v]
            for v in temp:
                print(v)
            print(">> {} candidates".format(len(temp)))
        else:
            for i in interest:
                temp = [v for v in self.drug if i in v]
                for v in temp:
                    print(v)
                print(">> {0} candidates for {1}".format(len(temp),i))
                print("")


    def check_rxn(self,interest=[],layer=["SOC","HLGT","HLT"]):
        """
        check whether adverse events related to the indicated keywords
                
        """
        if self.meddra.empty:
            raise ValueError("!! Set MedDRA data before this process !!")
        if len(interest)==0:
            raise ValueError("!! Give adverse events of interest as a list !!")
        if len(layer)==0:
            layer = MEDDRA_TERM
        else:
            layer = [v for v in layer if v in MEDDRA_TERM]
        for l in layer:
            print("------------ {} ------------".format(l))
            temp = list(self.meddra[l])
            if type(interest)==str:
                hit = [v for v in temp if interest in v]
            else:
                hit = []               
                for i in interest:
                    hit += [v for v in temp if i in v]
            hit = set(hit)
            for v in hit:
                print(v)
            print("")

    # revise on 210217
    def set_rxn(self,interest=[],layer="SOC"):
        """
        set PT (adverse events) based on the indicated SOC, HLGT, or HLT

        Parameters
        ----------
        interest: str
            indicates the adverse events of interest

        """
        if layer not in MEDDRA_TERM:
            raise ValueError("Wrong layer: check MedDRA layer")
        temp = self.meddra.copy()
        res = set()
        for i in interest:
            #temp = temp[temp[layer]==i]["PT"] # revise
            test = temp[temp[layer]==i]["PT"]
            if test.empty:
                print("CAUTION: `{}` is wrong key (skipped)".format(i))
            else:
                res = res | set(list(test))
        self.focused_pt = res
        print("{} terms were selected".format(len(self.focused_pt)))


    def generate_table(self):
        """
        generate contingency table
        take a long time

        ### contingency table ###
                            focusing adverse effect  not focusing AE     sum
        
        focusing drug                 n11                n12             n1p
        
        not focusing drug             n21                n22             n2p
        
        sum                           np1                np2             npp

        """
        if len(self.focused_pt)==0:
            raise ValueError("!! Run set_rxn() before this process !!")

        # total records
        npp = self.data.shape[0]

        # np1,np2
        fxn = lambda x: len(x & self.focused_pt) > 0
        temp = self.data[self.data["Reactions"].map(fxn)]
        idx_np1 = set(temp.index)
        np1 = len(idx_np1)
        np2 = npp - np1

        # others
        table = []
        ap = table.append
        for d in tqdm(self.drug):
            fxn = lambda x: len(x & set([d])) > 0
            temp = self.data[self.data["Active Substances"].map(fxn)]
            idx_n1p = set(temp.index)
            n1p = len(idx_n1p)
            n2p = npp - n1p
            idx_n11 = idx_np1 & idx_n1p
            n11 = len(idx_n11)
            n12 = n1p - n11
            n21 = np1 - n11
            n22 = n2p - n21
            ap([d,n11,n12,n21,n22,n1p,n2p])
        table = pd.DataFrame(table,columns=["drug","n11","n12","n21","n22","n1p","n2p"])
        table = table.set_index("drug",drop=True)
        n_drug = len(self.drug)
        table["np1"] = np.full(n_drug,np1,dtype=int)
        table["np2"] = np.full(n_drug,np2,dtype=int)
        table["npp"] = np.full(n_drug,npp,dtype=int)
        self.res = table


    def calc_stat(self,table=None):
        """
        calculate statistics from the contingency table
        
        """
        if table is None:
            if self.res.empty:
                raise ValueError("!! Give table or run generate_table() before this process !!")
            else:
                table = self.res.copy()
        self.drug = sorted(list(table.index))
        n11,n12,n21,n22 = table["n11"].values,table["n12"].values,table["n21"].values,table["n22"].values
        ror = self.__stats.calc_ror(n11,n12,n21,n22)
        lci,uci = self.__stats.calc_ci(n11,n12,n21,n22)
        #chi2,p,q_bh,q_holm = self.__stats.calc_chi2(n11,n12,n21,n22)
        table["lower_CI"] = lci
        table["upper_CI"] = uci
        table["ROR"] = ror
        #table["chi-square"] = chi2
        #table["adjusted_p_value_BH"] = q_bh
        #table["adjusted_p_value_Holm"] = q_holm
        #table["p_value"] = p
        table = table.sort_values("lower_CI",ascending=False)
        self.res = table


    def integrate_drug(self,name="",drugs=[]):
        """
        merge counts of the indicated drugs
        
        Parameters
        ----------
        name: str
            indicates the drug name after integration

        drugs: list
            indicates the drugs to be integrated
        
        """
        if len(name)==0:
            name = "integrated_drug"
        merge = [v for v in drugs if v in self.drug]
        if len(merge) > 1:
            rest = [v for v in self.drug if v not in merge]
            temp = self.res.loc[merge,:]
            npp = temp["npp"].values[0]
            np1 = temp["np1"].values[0]
            np2 = temp["np2"].values[0]
            n1p = np.sum(temp["n1p"])
            n2p = npp - n1p
            n11 = np.sum(temp["n11"])
            n12 = np.sum(temp["n12"])
            n21 = np1 - n11
            n22 = np2 - n12
            col = ["n11","n12","n21","n22","n1p","n2p","np1","np2","npp"]
            rest = self.res.loc[rest,col]
            df = pd.DataFrame([n11,n12,n21,n22,n1p,n2p,np1,np2,npp],index=col,columns=[name]).T
            rest = pd.concat([rest,df],join="inner",axis=0)
            self.calc_stat(rest)
        else:
            print("CAUTION: no integration because the indicated drugs are not enough")


class CommonFxn():
    def __init__(self):
        self.ror = np.array([])


    def calc_ror(self,n11,n12,n21,n22):
        """ calculate ROR """
        self.ror = (n11*n22)/(n12*n21)
        return self.ror


    def calc_ci(self,n11,n12,n21,n22):
        """ calculate confidence interval """
        if len(self.ror)==0:
            ror = self.calc_ror(n11,n12,n21,n22)
        temp = 1.96*np.sqrt(1/n11 + 1/n12 + 1/n21 + 1/n22)
        upperCI = np.exp(np.log(self.ror) + temp)
        lowerCI = np.exp(np.log(self.ror) - temp)
        return lowerCI,upperCI


    def calc_chi2(self,n11,n12,n21,n22):
        """ calculate chi-squared values """
        contingency = np.stack([n11,n12,n21,n22],axis=1)
        lst_chi = [stats.chi2_contingency(np.reshape(v,(2,2))) for v in contingency]
        chi = np.array([v[0] for v in lst_chi])
        p = np.array([v[1] for v in lst_chi])
        q_bh = multitest.multipletests(p,alpha=0.05,method="fdr_bh")[1]
        q_holm = multitest.multipletests(p,alpha=0.05,method="holm")[1]
        return chi,p,q_bh,q_holm


"""
*** Note ***
len([condition]) was faster than sum(condition) and len(list(filter))

"""