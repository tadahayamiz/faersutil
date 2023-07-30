#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 00:53:19 2021

@author: docker
"""
import pandas as pd
from itertools import chain
import numpy as np
from tqdm import tqdm
import os

if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")

class DateNormalizer():
    def __init__(self):
        self.raw_whole = None
    
    def set_info(self,raw_path:str,datenorm_path:str,identifier_path:str,drugbank_info):
        """
        raw_path : '/path/to/narrawed/and/cleaning/data/'
            contains e.g. ADR16Q3.pkl under the path
        datenorm_path : '/path/to/date/normalized/narrawed/and/cleaning/data/'
            contains e.g. ADR16Q3.pkl under the path
        identifier_path : '/mnt/FAERS/github/faersutil/faersutil/identifier/dwh/identifier.txt'
            processed with identifier_processor.py
        drugbank_info : DataFrame
            columns = ["ID","name","cas","dates","date"]
        """
        self.raw_path = raw_path
        self.datenorm_path = datenorm_path
        self.ident_df = pd.read_table(identifier_path,index_col=0)
        self.db_info = drugbank_info
    
    def set_whole(self,raw_whole,datenorm_whole):
        """
        note that raw_whole and datenorm_whle have the same index
        """
        self.raw_whole = raw_whole
        self.datenorm_whole = datenorm_whole
    
    def reflect_drugbank(self):
        """
        1. defne drug marketing date with DrugBanks
        2. drop unreliable data
        """
        total_comp = sorted(set(chain.from_iterable(self.raw_whole["Active Substances"]))) # 4304
        
        db_overlap = [] # the coumpound IDs registered in drugbank
        db_date = []
        for c in total_comp:
            names = self.ident_df[self.ident_df["CID"]==c]["name"].tolist()
            drug_tmp = self.db_info[self.db_info["name"].isin(names)].sort_values("dates")
            if len(drug_tmp) > 0:
                d = str(drug_tmp["date"].tolist()[0]).split(" ")[0] # if there is more than one, select the oldest one
                tmp = d.split("-")
                norm_d = tmp[0]+tmp[1]+tmp[2]
                db_overlap.append(c)
                db_date.append(norm_d)
            else:
                pass
        self.db_dic = dict(zip(db_overlap,db_date))
        
        """drop unreliable records"""
        self.unreliable_idx = []
        # 1. registered date is inappropriate (e.g. '10180302','10190205')
        reliable = self.datenorm_whole.dropna(subset=["Event Date"])
        diff_idx = list(set(self.datenorm_whole.index)-set(reliable.index))
        self.unreliable_idx.extend(diff_idx)
        
        # 2. records reported prior to drugbank marketing date
        # multi loop (1.58 it/s) takes a few minutes
        for c in tqdm(list(self.db_dic.keys())):
            fxn = lambda x : len(x & set([c]))>0
            tmp_df = reliable[reliable["Active Substances"].map(fxn)]
            f = self.db_dic.get(c)
            unreliable = tmp_df[tmp_df['Event Date']<f] # reported before its marketing date
            self.unreliable_idx.extend(unreliable.index.tolist())
        # narraw down the relaibale records
        reliable_idx = set(self.datenorm_whole.index.tolist())-set(self.unreliable_idx)
        self.reliable = self.datenorm_whole.loc[list(reliable_idx)] # 3319935
        print("DrugBank information reflection completed")
        print(len(self.raw_whole),"-->",len(self.reliable))
        print("")
        
    def firstdate_definition(self):
        """
        1. define drug marketng date with FAERS reported date
        2. merge drugbank info and firstdate info
        3. create the drug and its first marketing date dict
        """
        self.reliable.index = [i for i in range(len(self.reliable))] # reindex
        total_com = sorted(set(chain.from_iterable(self.reliable["Active Substances"])))
        com = self.reliable["Active Substances"].tolist()
        first_date = []
        target_size = []
        reported_size = []
        # time consuming
        for c in tqdm(total_com):
            idxs = [] # records index with focusing compound
            for i in range(len(com)):
                if c in com[i]:
                    idxs.append(i)
                else:
                    pass
            target_size.append(len(idxs))
            tmp_df = self.reliable.loc[idxs]
            #tmp_df['Event Date'] = tmp_df['Event Date'].replace('',np.nan)
            #tmp_df = tmp_df.dropna() # not missing the value
            reported_size.append(len(tmp_df))
            first = tmp_df.sort_values("Event Date",ascending=False)["Event Date"].tolist()[-1]
            first_date.append(first)
        total_first_dic = dict(zip(total_com,first_date))
        #pd.to_pickle(first_dic,'/mnt/FAERS/workspace/211127_renew/211128_whole_def/result/total_first_dic.pkl')
        
        use_first = sorted(list(set(total_com)-set(list(self.db_dic.keys()))))
        final_k = []
        final_v = []
        for t in use_first:
            final_k.append(t)
            final_v.append(total_first_dic.get(t))
        for i,k in enumerate(self.db_dic):
            final_k.append(k)
            final_v.append(self.db_dic.get(k))
        self.final_dic = dict(zip(final_k,final_v))
    
    def dropduplicate(self):
        # focus on the records age and sex are not missing
        df1 = self.reliable[~self.reliable["Patient Age"].isin([""])]
        df2 = df1[~df1["Sex"].isin([""])]
        print(len(self.reliable),"; original size")
        print(len(df1),"; age registered")
        print(len(df2),"; age+sex registered")
        # drop duplication
        target_idx = df2.index.tolist()
        other_idx = sorted(list(set(self.reliable.index.tolist())-set(target_idx)))
        
        fxn = lambda x : str(x)
        df2['Active Substances'] = df2['Active Substances'].map(fxn)
        df2['Reactions'] = df2['Reactions'].map(fxn)
        
        df3 = df2.drop_duplicates(subset=['Active Substances','Reactions','Sex','Event Date','Patient Age'],keep='first')
        target_idx2 = sorted(df3.index.tolist())
        print(len(target_idx2),"; after drop duplication")
        print(len(df2),"-->",len(df3))
        target_idx2.extend(other_idx)
        after = self.reliable.loc[target_idx2]
        
        after = after.sort_values("Event Date",ascending=False)
        after.index = [i for i in range(len(after))]
        self.final_whole = after
    
    