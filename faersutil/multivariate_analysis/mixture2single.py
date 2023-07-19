#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 05:39:49 2021

To understand the impact of compounds that are frequently used together as a combination drug as a single drug, change from a combination drug to a single drug and reflect it in the record

@author: docker
"""
import pandas as pd
from itertools import chain
import copy
import collections

class Mixture2Single():
    def __init__(self):
        self.whole = pd.DataFrame()
        self.mixture = None
    
    def set_reference(self,whole,syn_dic):
        """
        whole : whole records before converting
        syn_dic : dictionary about compound name and its ID
        """
        self.whole = whole
        self.syn_dic = syn_dic
        
        whole_comp = sorted(chain.from_iterable(self.whole["Active Substances"]))
        self.original_c = dict(collections.Counter(whole_comp))
    
    def set_mixture(self,data=[]):
        if len(data)==0:
            self.mixture = [("trimethoprim","sulfamethoxazole"),("tazobactam sodium","piperacillin sodium"),("clavulanate potassium","amoxicillin"),("sacubitril","valsartan"),("salmeterol xinafoate","fluticasone propionate"),("carbidopa","levodopa"),("hydrocodone bitartrate","acetaminophen"),("ledipasvir","sofosbuvir"),("sodium chloride","sodium lactate","calcium chloride")]
        else:
            self.mixture = data
        print("--- mixture compounds set ---")
        print(self.mixture)
        print("------")
    
    def set_synonym(self,data=[]):
        if len(data)==0:
            self.synonym = [("zolpidem","zolpidem tartrate"),("mycophenolate mofetil hydrochloride","mycophenolate mofetil")]
        else:
            self.synonym = data
        print("--- synonym compounds set ---")
        print(self.synonym)
        print("------")
    
    def mix2single(self):
        """
        A-B --> C
        """
        # convert compounds name to ID
        rev_dic = dict(zip(list(self.syn_dic.values()),list(self.syn_dic.keys())))
        self.mix_id_list = []
        for t in self.mixture:
            tmp = []
            for i in range(len(t)):
                tmp.append(rev_dic.get(t[i]))
            self.mix_id_list.append(tuple(tmp))
        # copy the original records and process
        df = copy.deepcopy(self.whole)
        
        #max_id = max(set(chain.from_iterable(df["Active Substances"])))
        max_id = max(list(self.syn_dic.keys()))
        print("current largest ID",max_id)
        
        aftername = []
        afterid = []
        for i in range(len(self.mix_id_list)):
            pair = self.mixture[i]
            # register new name
            tmp_name = pair[0]
            for j in range(1,len(pair)):
                tmp_name = tmp_name+" --- "+pair[j]
            aftername.append(tmp_name)
            add_id = max_id+1+i
            # register new id
            afterid.append(add_id)
            
            interest = set(self.mix_id_list[i])
            def converter(x):
                """
                interest : set([drug1 ID, drug2 ID])  e.g. {6047, 34359}
                """
                if len(x & interest)==len(pair):
                    new_x = (x-interest).union(set([add_id]))
                    return new_x
                else:
                    return x
            # apply
            df["Active Substances"] = df["Active Substances"].apply(converter)
        
        self.add_dic = dict(zip(afterid,aftername))
        self.result = df
        self.after_c = dict(collections.Counter(sorted(chain.from_iterable(self.result["Active Substances"]))))
    
    def mixexe_check(self):
        for t in self.mix_id_list:
            for s in t:
                print(self.syn_dic.get(s),end=" : ")
                print(self.original_c.get(s),"-->",self.after_c.get(s))
    
    def syn2single(self):
        """
        after mixture converting
        A-B --> A convert to majority
        """
        # convert compounds name to ID
        rev_dic = dict(zip(list(self.syn_dic.values()),list(self.syn_dic.keys())))
        self.syn_id_list = []
        for t in self.synonym:
            tmp = []
            for i in range(len(t)):
                tmp.append(rev_dic.get(t[i]))
            self.syn_id_list.append(tuple(tmp))
        print(self.syn_id_list)
        # copy the original records and process
        df2 = copy.deepcopy(self.result) # after mixture processing
        
        for i in range(len(self.syn_id_list)):
            pair = self.syn_id_list[i]
            m = 0
            max_id = None # pick up the majority name
            for j in range(len(pair)):
                if self.after_c.get(pair[j]) > m:
                    m = self.after_c.get(pair[j])
                    max_id = pair[j]
                else:
                    pass

            interest = set(pair)-set([max_id]) # others
            def converter2(x):
                """
                interest : set([drug1 ID, drug2 ID])  e.g. {6047, 34359}
                """
                if len(x & interest)>0:
                    common = x & interest
                    new_x = (x-common).union(set([max_id])) # convert to the majority one
                    return new_x
                else:
                    return x
            # apply
            df2["Active Substances"] = df2["Active Substances"].apply(converter2)
        
        self.result2 = df2
        self.after_c2 = dict(collections.Counter(sorted(chain.from_iterable(self.result2["Active Substances"]))))
        
    def synexe_check(self):
        for t in self.syn_id_list:
            for s in t:
                print(self.syn_dic.get(s),end=" : ")
                print(self.original_c.get(s),"-->",self.after_c2.get(s))
    
    def merge_dict(self):
        self.syn_dic.update(self.add_dic)
        return self.syn_dic
    


