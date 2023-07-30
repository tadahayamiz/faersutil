#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 20 07:14:42 2021

Focus on records registered with one and two drugs, and infer the effect of drug interaction from the changes in the ratio of SOC for each

@author: docker
"""

import pandas as pd
from itertools import chain
import copy
import collections
from tqdm import tqdm
import matplotlib.pyplot as plt
from tqdm import tqdm

class MixtureComparison():
    def __init__(self):
        self.whole = pd.DataFrame()
        
    def set_reference(self,whole,syn_dic,soc_dic=None,exception=["Product issues","Congenital, familial and genetic disorders"]):
        """
        whole : whole records before converting
        syn_dic : dictionary about compound name and its ID
        """
        self.whole = whole
        self.syn_dic = syn_dic
      
        whole_comp = sorted(chain.from_iterable(self.whole["Active Substances"]))
        self.original_c = dict(collections.Counter(whole_comp))
        
        # SOC dict process
        use_k = []
        use_v = []
        for i,k in enumerate(soc_dic):
            if k in exception:
                pass
            else:
                use_k.append(k)
                use_v.append(soc_dic.get(k))
        self.soc_dic = dict(zip(use_k,use_v))
    
    def set_target(self,total_target,base_target):
        self.total_target = total_target
        self.base_target = base_target
    
    def target_definition(self,total_th=10000,single_th=5000):
        self.total_target = []
        self.base_target = []
        for i,k in enumerate(self.original_c):
            if self.original_c.get(k)>total_th:
                self.total_target.append(k)
            else:
                pass
        print("compounds registered over %d records ;"%total_th,len(self.total_target))
        fxn_single = lambda x : len(x)==1
        single = self.whole[self.whole["Active Substances"].map(fxn_single)]
        for t in tqdm(self.total_target):
            fxn = lambda x : len(x & set([t]))==1 # registered as single 
            tmp_df = single[single['Active Substances'].map(fxn)]
            if len(tmp_df) > single_th:
                self.base_target.append(t)
            else:
                pass

    def base_soc(self,soc_dic,exception=["Product issues","Congenital, familial and genetic disorders"],single_id=1983):
        """
        focus on single agent
        plot the SOC distribution
        """
        self.single_id = single_id
        self.name = self.syn_dic.get(single_id)
        if self.soc_dic is None:
            use_k = []
            use_v = []
            for i,k in enumerate(soc_dic):
                if k in exception:
                    pass
                else:
                    use_k.append(k)
                    use_v.append(soc_dic.get(k))
            self.soc_dic = dict(zip(use_k,use_v))
        else:
            pass
        
        # focus on the single agent
        fxn1 = lambda x : len(x)==1
        self.single = self.whole[self.whole["Active Substances"].map(fxn1)]
        fxn2 = lambda x : len(x & set([single_id]))==1
        self.base_single = self.single[self.single["Active Substances"].map(fxn2)]
        
        self.base_size = soc_assign(self.base_single,self.soc_dic)
        
        # plot the SOC distribution
        fig,ax = plt.subplots(figsize=(15,10))
        x = [i for i in range(len(self.soc_dic))]
        plt.bar(x,self.base_size,alpha=0.5,label="single agent (%d)"%len(self.base_single))
        #plt.bar(x,comb,alpha=0.5,label="combination (%d) ; "%len(pair_df)+p_name)
        plt.xticks(x,list(self.soc_dic.keys()),rotation=90)
        plt.ylabel("ratio")
        plt.title(self.name)
        plt.legend(loc="best",fontsize=12)
        plt.show()
    
    def comb_comparison(self,threshold=100,max_change=2,min_change=2):
        """
        after performing base_soc()
        focus on the drug combination registered as two agents
        plot SOC distribution 1.combination, 2.base, 3.pair
        """
        fxn_double = lambda x : len(x)==2
        double = self.whole[self.whole["Active Substances"].map(fxn_double)]
        
        interesting_pair = []
        interesting_soc = []
        for t in tqdm(self.total_target):
            p_name = self.syn_dic.get(t) 
            fxn_d = lambda x : len(x & set([self.single_id,t]))==2
            common_df = double[double["Active Substances"].map(fxn_d)] # base & pair records
            if len(common_df) > threshold:
                common_size = soc_assign(common_df,self.soc_dic)
                fxn3 = lambda x : len(x & set([t]))==1
                t_single = self.single[self.single["Active Substances"].map(fxn3)] # pair single records
                t_size = soc_assign(t_single,self.soc_dic)
                
                change_idx = change_detector(self.base_size,t_size,common_size,max_change=max_change,min_change=min_change)
                if len(change_idx)>0:
                    print(p_name)
                    interesting_pair.append(p_name)
                    interesting_soc.append([list(self.soc_dic.keys())[j] for j in change_idx])
                    # plot the SOC distribution
                    fig,ax = plt.subplots(figsize=(12,8))
                    x = [i*2 for i in range(len(self.soc_dic))]
                    plt.bar([t-0.4 for t in x],common_size,alpha=0.9,label="combination (%d)"%len(common_df))
                    plt.bar([t+0.4 for t in x],self.base_size,alpha=0.5,label="base single (%d) ; "%len(self.base_single)+self.name)
                    plt.bar([t+0.4 for t in x],t_size,alpha=0.5,label="pair single (%d) ; "%len(t_single)+p_name)
                    plt.xticks(x,list(self.soc_dic.keys()),rotation=90)
                    plt.ylabel("ratio")
                    plt.legend(loc="best",fontsize=12)
                    plt.title(self.name+"+"+p_name)
                    plt.show()
                else:
                    pass
            else:
                pass
        self.interest_dic = dict(zip(interesting_pair,interesting_soc))
    
    def focus_ab(self,drug1,drug2,name_id=False,max_change=2,min_change=2):
        """
        identify the effects of two specific drugs
        1. drug1 and drug2
        2. drug1 single
        3. drug2 single
        """
        if name_id:
            pass
        else:
            rev_dic = dict(zip(list(self.syn_dic.values()),list(self.syn_dic.keys())))
            drug1 = rev_dic.get(drug1)
            drug2 = rev_dic.get(drug2)
        
        name = self.syn_dic.get(drug1) # drug1 name
        p_name = self.syn_dic.get(drug2) # drug2 name
        
        fxn_double = lambda x : len(x)==2
        fxn1 = lambda x : x == set([drug1])
        fxn2 = lambda x : x == set([drug2])
        double = self.whole[self.whole["Active Substances"].map(fxn_double)]
        self.single1 = self.whole[self.whole["Active Substances"].map(fxn1)]
        self.single2 = self.whole[self.whole["Active Substances"].map(fxn2)]
        fxn_d = lambda x : len(x & set([drug1,drug2]))==2
        self.common_df = double[double["Active Substances"].map(fxn_d)] # registered as just two agents
        
        # obtain SOC freq
        common_res = soc_assign(self.common_df,self.soc_dic)
        res1 = soc_assign(self.single1,self.soc_dic)
        res2 = soc_assign(self.single2,self.soc_dic)
        
        # detect the freq change signal
        change_idx = change_detector(res1,res2,common_res,max_change=max_change,min_change=min_change)
        self.signal_res = []
        for j in change_idx:
            tmp_dic = {list(self.soc_dic.keys())[j]:(common_res[j],res1[j],res2[j])}
            self.signal_res.append(tmp_dic)
        
        # plot the SOC distribution
        fig,ax = plt.subplots(figsize=(12,8))
        x = [i*2 for i in range(len(self.soc_dic))]
        plt.bar([t-0.4 for t in x],common_res,alpha=0.9,label="combination (%d)"%len(self.common_df))
        plt.bar([t+0.4 for t in x],res1,alpha=0.5,label="base single (%d) ; "%len(self.single1)+name)
        plt.bar([t+0.4 for t in x],res2,alpha=0.5,label="pair single (%d) ; "%len(self.single2)+p_name)
        plt.xticks(x,list(self.soc_dic.keys()),rotation=90)
        plt.ylabel("ratio")
        plt.legend(loc="best",fontsize=12)
        plt.title(name+"+"+p_name)
        plt.show()
    
    def focus_multiab(self,drug1,drug2,use_id=False,max_change=2,min_change=2):
        """
        1. records contain drug1 and drug2
        2. records contain drug1
        3. records contain drug2
        """
        if use_id:
            pass
        else:
            rev_dic = dict(zip(list(self.syn_dic.values()),list(self.syn_dic.keys())))
            drug1 = rev_dic.get(drug1)
            drug2 = rev_dic.get(drug2)
        
        name = self.syn_dic.get(drug1) # drug1 name
        p_name = self.syn_dic.get(drug2) # drug2 name
        
        fxn_common = lambda x : len(x & set([drug1,drug2]))==2
        self.common_df = self.whole[self.whole["Active Substances"].map(fxn_common)] # contain both drugs
        other_idx = list(set(self.whole.index.tolist()) - set(self.common_df.index.tolist()))
        other_df = self.whole.loc[other_idx]
        
        fxn1 = lambda x : len(x & set([drug1]))==1
        fxn2 = lambda x : len(x & set([drug2]))==1
        self.df1 = other_df[other_df["Active Substances"].map(fxn1)]
        self.df2 = other_df[other_df["Active Substances"].map(fxn2)]
        
        
        # obtain SOC freq
        common_res = soc_assign(self.common_df,self.soc_dic)
        res1 = soc_assign(self.df1,self.soc_dic)
        res2 = soc_assign(self.df2,self.soc_dic)
        
        # detect the freq change signal
        change_idx = change_detector(res1,res2,common_res,max_change=max_change,min_change=min_change)
        self.signal_res = []
        for j in change_idx:
            tmp_dic = {list(self.soc_dic.keys())[j]:(common_res[j],res1[j],res2[j])}
            self.signal_res.append(tmp_dic)
        
        # plot the SOC distribution
        fig,ax = plt.subplots(figsize=(12,8))
        x = [i*2 for i in range(len(self.soc_dic))]
        plt.bar([t-0.4 for t in x],common_res,alpha=0.9,label="contain both (%d)"%len(self.common_df))
        plt.bar([t+0.4 for t in x],res1,alpha=0.5,label="contain only(%d) ; "%len(self.df1)+name)
        plt.bar([t+0.4 for t in x],res2,alpha=0.5,label="contain only (%d) ; "%len(self.df2)+p_name)
        plt.xticks(x,list(self.soc_dic.keys()),rotation=90)
        plt.ylabel("ratio")
        plt.legend(loc="best",fontsize=12)
        plt.title(name+"+"+p_name)
        plt.show()
        
        
    def narraw_record(self,soc_dic,exception=["Product issues","Congenital, familial and genetic disorders"]):
        self.soc_dic = soc_dic
        use_k = []
        use_v = []
        for i,k in enumerate(self.soc_dic):
            if k in exception:
                pass
            else:
                use_k.append(k)
                use_v.append(self.soc_dic.get(k))
        use_dic = dict(zip(use_k,use_v))
        
        focus_pt = set()
        for t in use_v:
            focus_pt = focus_pt.union(t)
        
        # focus on the reactions in interest
        fxn = lambda x : len(x & set(focus_pt))>0
        self.df = self.whole[self.whole["Reactions"].map(fxn)]

def soc_assign(df,soc_dic:dict):
    soc_size = []
    for i,k in enumerate(soc_dic):
        tmp_fxn = lambda x : len(x & soc_dic.get(k))>0
        tmp_df = df[df["Reactions"].map(tmp_fxn)]
        soc_size.append(len(tmp_df)/len(df))
    return soc_size

def change_detector(x,y,z,signal_threshold=0.2,max_change=2,min_change=2):
    """
    x : single agent 1
    y : single agent 2
    z : combination as agent1 & agent2
    """
    base = [(s,t) for s,t in zip(x,y)]
    change_idx = [] # idx show the significant change
    for i in range(len(z)):
        comb = z[i]
        if comb > max(base[i])*max_change:
            if comb > signal_threshold:
                change_idx.append(i)
            else:
                pass # the combination signal is poor
        elif comb < min(base[i])/min_change:
            if min(base[i])>signal_threshold:
                change_idx.append(i)
            else:
                pass # the single agent signal is poor
        else:
            pass
    return change_idx
        
        
                
    