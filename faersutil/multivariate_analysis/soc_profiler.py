#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 01:05:22 2021

Generate and analyze profiles from the ratio of SOC reported for each drug

@author: docker
"""
import pandas as pd
import matplotlib.pyplot as plt
from itertools import chain,combinations
import collections
from tqdm import tqdm
import seaborn as sns
import numpy as np

class SocProfile():
    def __init__(self):
        self.whole = None
        self.size_dic = None
    
    def set_reference(self,whole,syn_dic,soc_dic=None,exception=["Congenital, familial and genetic disorders","General disorders and administration site conditions","Infections and infestations","Investigations","Product issues","Social circumstances","Surgical and medical procedures","without MedDRA"]):
        """
        whole : whole records before converting
        syn_dic : dictionary about compound name and its ID
        """
        self.whole = whole
        self.syn_dic = syn_dic
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
    
    def narrow_records(self):
        """
        drop records which contain only exceptional PT
        """
        focused_pt = set()
        for i,k in enumerate(self.soc_dic):
            focused_pt = focused_pt.union(self.soc_dic.get(k)) # update fucused_pt
        fxn = lambda x : len(x & focused_pt) > 0
        self.target = self.whole[self.whole["Reactions"].map(fxn)] # contain any focusing pt
        print(len(self.whole),"-->",len(self.target))
        
        # plot the frequency
        self.soc_list = list(self.soc_dic.keys())
        soc_size = []
        for s in self.soc_list:
            s_pt = self.soc_dic.get(s)
            fxn = lambda x : len(x & s_pt) > 0
            target = self.target[self.target["Reactions"].map(fxn)]
            soc_size.append(len(target))
        plt.bar([i for i in range(len(soc_size))],soc_size)
        plt.xticks([i for i in range(len(soc_size))],self.soc_list,rotation=90)
        plt.ylabel("frequency")
        plt.show()
    
    def focus_a(self,name1):
        rev_dic = dict(zip(self.syn_dic.values(), self.syn_dic.keys()))
        drug_id = rev_dic.get(name1)
        fxn = lambda x : len(x & set([drug_id]))>0
        df = self.target[self.target["Active Substances"].map(fxn)]
        res = soc_assign(df, self.soc_dic)
        
        # plot
        # plot the SOC distribution
        fig,ax = plt.subplots(figsize=(12,8))
        x = [i for i in range(len(self.soc_dic))]
        plt.bar([t for t in x],res,alpha=0.9,label="contain (%d) ;"%len(df)+name1)
        plt.xticks(x,list(self.soc_dic.keys()),rotation=90)
        plt.ylabel("ratio")
        plt.legend(loc="best",fontsize=12)
        plt.title(name1)
        plt.show()
    
    def focus_ab(self,name1,name2):
        rev_dic = dict(zip(self.syn_dic.values(), self.syn_dic.keys()))
        drug1 = rev_dic.get(name1)
        drug2 = rev_dic.get(name2)
        
        both_fxn = lambda x : len(x & set([drug1,drug2]))==2
        self.both = self.target[self.target["Active Substances"].map(both_fxn)]
        
        other_idx = list(set(self.target.index.tolist()) - set(self.both.index.tolist()))
        other = self.target.loc[other_idx]
        
        fxn1 = lambda x : len(x & set([drug1]))==1 # contain drug1
        fxn2 = lambda x : len(x & set([drug2]))==1 # contain drug2
        self.df1 = other[other["Active Substances"].map(fxn1)]
        self.df2 = other[other["Active Substances"].map(fxn2)]
        
        common_reaction = sorted(chain.from_iterable(self.both["Reactions"]))
        c = dict(collections.Counter(common_reaction))
        sc = sorted(c.items(),key=lambda x : x[1],reverse=True)
        counter=0
        print("---frequenly reported PT---")
        for i in range(10):
            if counter == 10:
                break
            else:
                print(sc[i][0],sc[i][1],"(%f2)"%(sc[i][1]/len(self.both)))
                counter+=1
                
    
    def focus_multi(self,id_list:list,soc_selection=False,target_soc=["Immune system disorders","Vascular disorders", "Blood and lymphatic system disorders"]):
        for drug_id in tqdm(id_list):
            name = self.syn_dic.get(drug_id)
            fxn = lambda x : len(x & set([drug_id]))>0
            df = self.target[self.target["Active Substances"].map(fxn)]
            if soc_selection:
                tmp_v = [self.soc_dic.get(k) for k in target_soc]
                tmp_soc_dic = dict(zip(target_soc,tmp_v))
                res = soc_assign(df, tmp_soc_dic)
                
                # plot
                # plot the SOC distribution
                fig,ax = plt.subplots(figsize=(12,8))
                x = [i for i in range(len(tmp_soc_dic))]
                plt.bar([t for t in x],res,alpha=0.9,label="contain (%d) ;"%len(df)+name)
                plt.xticks(x,list(tmp_soc_dic.keys()),rotation=90)
                plt.ylabel("ratio")
                plt.legend(loc="best",fontsize=12)
                plt.title(name)
                plt.show()
            else:
                res = soc_assign(df, self.soc_dic)
                
                # plot
                # plot the SOC distribution
                fig,ax = plt.subplots(figsize=(12,8))
                x = [i for i in range(len(self.soc_dic))]
                plt.bar([t for t in x],res,alpha=0.9,label="contain (%d) ;"%len(df)+name)
                plt.xticks(x,list(self.soc_dic.keys()),rotation=90)
                plt.ylabel("ratio")
                plt.legend(loc="best",fontsize=12)
                plt.title(name)
                plt.show()
            
    def comb_comparison(self,id_list,use_target=True,target_soc=["Immune system disorders","Vascular disorders", "Blood and lymphatic system disorders"],bothsize_th=10):
        pair = combinations(id_list,2) # total combination
        self.total_signal = []
        for p in tqdm(pair):
            drug1 = p[0]
            name1 = self.syn_dic.get(drug1)
            drug2 = p[1]
            name2 = self.syn_dic.get(drug2)
            
            fxn_both = lambda x : len(x & set([drug1,drug2]))==2
            both = self.target[self.target["Active Substances"].map(fxn_both)] # contain both compounds
            
            # stop if the number of records contain both compounds is small
            if len(both)<bothsize_th: 
                continue
            
            other_idx = list(set(self.target.index.tolist()) - set(both.index.tolist()))
            other = self.target.loc[other_idx]
            fxn1 = lambda x : len(x & set([drug1]))==1 # contain drug1
            fxn2 = lambda x : len(x & set([drug2]))==1 # contain drug2
            df1 = other[other["Active Substances"].map(fxn1)]
            df2 = other[other["Active Substances"].map(fxn2)]
            
            both_res = soc_assign(both,self.soc_dic)
            res1 = soc_assign(df1,self.soc_dic)
            res2 = soc_assign(df2,self.soc_dic)
            
            # detect the freq change signal
            change_idx = change_detector(res1,res2,both_res,signal_threshold=0.2,max_change=1.5,min_change=1.5)
            signal_res = []
            signal_key = [list(self.soc_dic.keys())[l] for l in change_idx]
            for j in change_idx:
                tmp_dic = {list(self.soc_dic.keys())[j]:(drug1,drug2,both_res[j],res1[j],res2[j])}
                signal_res.append(tmp_dic)
            if len(signal_res)>0: # detect signal about any SOC
                self.total_signal.append(signal_res)
            else:
                pass
            if use_target:
                ident = set(signal_key) & set(target_soc)
                if len(ident)>0:
                    # plot the SOC distribution
                    fig,ax = plt.subplots(figsize=(12,8))
                    x = [i*2 for i in range(len(self.soc_dic))]
                    plt.bar([t-0.4 for t in x],both_res,alpha=0.9,label="contain both (%d)"%len(both))
                    plt.bar([t+0.4 for t in x],res1,alpha=0.5,label="contain only(%d) ; "%len(df1)+name1)
                    plt.bar([t+0.4 for t in x],res2,alpha=0.5,label="contain only (%d) ; "%len(df2)+name2)
                    plt.xticks(x,list(self.soc_dic.keys()),rotation=90)
                    plt.ylabel("ratio")
                    plt.legend(loc="best",fontsize=12)
                    plt.title(name1+"+"+name2)
                    ax.set_axisbelow(True)
                    ax.grid(c='gainsboro')
                    plt.show()
                else:
                    pass
            else:
                if len(signal_res)>0:
                    # plot the SOC distribution
                    fig,ax = plt.subplots(figsize=(12,8))
                    x = [i*2 for i in range(len(self.soc_dic))]
                    plt.bar([t-0.4 for t in x],both_res,alpha=0.9,label="contain both (%d)"%len(both))
                    plt.bar([t+0.4 for t in x],res1,alpha=0.5,label="contain only(%d) ; "%len(df1)+name1)
                    plt.bar([t+0.4 for t in x],res2,alpha=0.5,label="contain only (%d) ; "%len(df2)+name2)
                    plt.xticks(x,list(self.soc_dic.keys()),rotation=90)
                    plt.ylabel("ratio")
                    plt.legend(loc="best",fontsize=12)
                    plt.title(name1+"+"+name2)
                    ax.set_axisbelow(True)
                    ax.grid(c='gainsboro')
                    plt.show()
                else:
                    pass

    def create_freqmatrix(self,threshold=100):
        self.total_comp = sorted(set(chain.from_iterable(self.target["Active Substances"]))) # unique total compounds
        print(len(self.total_comp),"compounds are registered in target records")
        c_dic = dict(collections.Counter(sorted(chain.from_iterable(self.target["Active Substances"]))))
        sc = sorted(c_dic.items(),key=lambda x : x[1],reverse=True)
        self.analysis_comp = []
        for t in sc:
            if t[1] > threshold:
                self.analysis_comp.append(t[0])
            else:
                break
        print("componds",len(self.total_comp),"-->",len(self.analysis_comp))
        
        self.pt_list = [self.soc_dic.get(x) for x in self.soc_list]
        self.matrix = []
        size = []
        for c in tqdm(self.analysis_comp):
            fxn = lambda x : len(x & set([c]))>0
            df = self.target[self.target["Active Substances"].map(fxn)] # contain the compound
            size.append(len(df))
            df_pts = df["Reactions"].tolist()
            count = [0]*len(self.soc_dic)
            for pt in df_pts:
                for i in range(len(self.pt_list)):
                    if len(pt & self.pt_list[i]) > 0:
                        count[i] += 1
                    else:
                        pass
            self.matrix.append(count)
        
        self.size_dic = dict(zip(self.analysis_comp,size))
                
        # conver to dataframe
        self.count_df= pd.DataFrame(self.matrix)
        self.count_df.columns = self.soc_list
        self.count_df.index = self.analysis_comp
        
    def ratio_norm(self,df=None,size_dic=None):
        """
        column : SOC list
        index : compounds list
        """
        # set target data
        if df is None:
            df = self.count_df
        else:
            pass
        # registeration number for each compound
        if size_dic is None:
            size = []
            for c in tqdm(df.index.tolist()):
                fxn = lambda x : len(x & set([c]))>0
                tmp_df = self.target[self.target["Active Substances"].map(fxn)]
                size.append(len(tmp_df))
            self.size_dic = dict(zip(df.index.tolist(),size))
        else:
            self.size_dic = size_dic
        print("coumpounds freq dict prepared")
        
        # normalize with records size
        self.norm_matrix = []
        for t in tqdm(df.index.tolist()):
            tmp = df.loc[t].tolist()
            self.size_dic.get(t)
            new_v = [x/self.size_dic.get(t) for x in tmp]
            self.norm_matrix.append(new_v)
        
        self.norm_matdf = pd.DataFrame(self.norm_matrix)
        self.norm_matdf.columns = self.soc_list
        self.norm_matdf.index = df.index.tolist()
        # plot
        sns.clustermap(self.norm_matdf,col_cluster=False)
        plt.show()
    
    def center_norm(self,df=None,threshold=100):
        # set target data
        if df is None:
            df = self.norm_matdf
        else:
            pass
        sc = sorted(self.size_dic.items(),key=lambda x : x[1],reverse=True)
        target_comp = []
        for t in sc:
            if t[1]>threshold:
                target_comp.append(t[0])
            else:
                break
        narrow = df.loc[target_comp]
        print(len(df),"-->",len(target_comp))
        # centralization
        for i in range(len(narrow.T)):
            cent_fxn = lambda x : x-np.mean(narrow.iloc[:,i])
            narrow.iloc[:,i] = narrow.iloc[:,i].apply(cent_fxn)
        self.norm_matdf2 = narrow
        # plot
        sns.clustermap(self.norm_matdf2,col_cluster=False)
        plt.show()
    
    def freq_norm(self,df=None,foldchange=True):
        # set target data
        if df is None:
            df = self.norm_matdf
        else:
            pass
        self.whole_freq = soc_assign(self.target,self.soc_dic)
        if foldchange: # cal foldchange
            self.norm_matrix2 = []
            for t in tqdm(df.index.tolist()):
                tmp = df.loc[t].tolist()
                new_v = []
                for i in range(len(tmp)):
                    new_v.append(tmp[i]/self.whole_freq[i])
                self.norm_matrix2.append(new_v)
            
            self.norm_matdf2 = pd.DataFrame(self.norm_matrix2)
            self.norm_matdf2.columns = self.soc_list
            self.norm_matdf2.index = df.index.tolist()
        else: # calc difference
            self.norm_matrix2 = []
            for t in tqdm(df.index.tolist()):
                tmp = df.loc[t].tolist()
                new_v = []
                for i in range(len(tmp)):
                    new_v.append(tmp[i] - self.whole_freq[i])
                self.norm_matrix2.append(new_v)
            
            self.norm_matdf2 = pd.DataFrame(self.norm_matrix2)
            self.norm_matdf2.columns = self.soc_list
            self.norm_matdf2.index = df.index.tolist()
        # plot
        sns.clustermap(self.norm_matdf2,col_cluster=False)
        plt.show()
    
    def narrow_profile(self,threshold=100,target_soc=["Hepatobiliary disorders","Gastrointestinal disorders"],SOC_limitation=False,df=None):
        if df is None:
            df = self.norm_matdf2
        else:
            pass
        sc = sorted(self.size_dic.items(),key=lambda x : x[1],reverse=True)
        target_comp = []
        for t in sc:
            if t[1]>threshold:
                target_comp.append(t[0])
            else:
                break
        print(len(df),"-->",len(target_comp))
        if SOC_limitation:
            df = df[target_soc]
        else:
            pass
        self.narrow = df.loc[target_comp]
        sns.clustermap(self.narrow,col_cluster=False)
        plt.show()
    

def soc_assign(df,soc_dic:dict):
    if len(df)==0:
        soc_size = [0]*len(soc_dic)
    else:
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