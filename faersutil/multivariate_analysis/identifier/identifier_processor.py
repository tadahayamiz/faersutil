#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 10:31:07 2021

reflect additional condition and process identifer.txt

@author: docker
"""
import copy
import pandas as pd
import collections

DEL_COMP = ["acetate","anhydrous",
            "benzonate","benzoate",
            "carbonate","chloride","calcium",
            "dihydrate","dihydrochloride","dimer",
            "fumarate",
            "hemifumarate","hydrobromide","hydrochloride",
            "maleate","magnesium","methanesulfonate","mesylates","mesylate",
            "pottasium",
            "sodium","sulfate","saccharate",
            "trifluoroacetate","trifluoroacetic acid","trifydrate"]

class IdentProcessor():
    def __init__(self):
        self.legacy = None
    
    def set_legacy(self,data_path):
        """
        legacy : path/to/processing/target/legacy/identifier.txt
        """
        self.legacy = pd.read_table(data_path,index_col=0)
        self.legacy.index = [i for i in range(len(self.legacy))]
        
    
    def multi_update(self,keywords=DEL_COMP):
        """update the identifier.txt"""
        for w in keywords:
            df = self.legacy
            self.legacy = update(df,w)
    
    def final_process(self):
        c = sorted(dict(collections.Counter(self.legacy["name"].tolist())).items(),key=lambda x : x[1],reverse=True)
        dup_pro = pd.DataFrame()
        for x in c:
            if x[1] > 1:
                name = x[0]
                tmp_df = self.legacy[self.legacy["name"]==name]
                cids = tmp_df["CID"].tolist()
                target = self.legacy[self.legacy["CID"].isin(cids)]
                target["CID"] = [min(cids)]*len(target)
                dup_pro = pd.concat([dup_pro,target])
            else:
                break
        dup_pro = dup_pro.drop_duplicates()
        other = self.legacy[~self.legacy["name"].isin(dup_pro["name"].tolist())]
        self.final = pd.concat([other,dup_pro])

        # remove distinction of head " "
        print("--- remove distinction of head space ---")
        total_comp = self.final["name"].tolist()
        
        change = []
        new_id = []
        for c in total_comp:
            try:
                if c[0]==" ":
                    if c[1::] in total_comp:
                        change.append(c)
                        new_id.append(self.final[self.final["name"]==c[1::]]["CID"].tolist()[0])
                    else:
                        pass
                else:
                    pass
            except:
                pass
        # update
        change_df = pd.DataFrame({"name":change,"CID":new_id})
        other = self.final[~self.final["name"].isin(change)]
        self.final = pd.concat([change_df,other])
        self.final = self.final.sort_values("CID",ascending=False)
            
    
    def create_iddic(self):
        unique_ids = self.final.dropna()["CID"].unique().tolist()
        representative = []
        for cid in unique_ids:
            tmp_df = self.final[self.final["CID"]==cid]
            names = tmp_df["name"].tolist()
            min_size = 1000
            rep = None
            for name in names:
                if len(name.split(" ")) < min_size: # select smaller size one
                    min_size = len(name.split(" "))
                    rep = name
                else:
                    pass
            representative.append(rep)
        self.id_dic = dict(zip(unique_ids,representative))

def update(legacy,w):
    # total compounds
    print("process start :",w)
    total = legacy["name"].dropna().tolist()
    contain = [] # contain the key words
    not_contain = [] # not contain (can be base name)
    for name in total:
        if w in name.split(" "):
            contain.append(name)
        else:
            not_contain.append(name)
    print(len(contain),"compounds contain the target keywords")
    if len(contain)==0:
        print("completed")
        print("")
        return legacy
    else:
        base = []
        pair = []
        no_overlap = []
        for name in contain:
            common = set(set(name.split(" "+w)).union(set(name.split(w+" ")))) & set(not_contain)
            if len(common) > 0:
                base.append(list(common)[0])
                pair.append(name)
            else:
                no_overlap.append(name)
        print(len(pair),"/",len(contain),"compounds assined to its base compound")
        if len(pair)==0:
            print("completed")
            print("")
            return legacy
        else:
            update_count = 0
            update_df = pd.DataFrame()
            for i in range(len(pair)):
                b = base[i]
                p = pair[i]
                tmp_legacy = legacy[legacy["name"].isin([b,p])]
                if len(tmp_legacy["CID"].unique().tolist())==1: # already regarded as the same CID
                    pass
                else:
                    p_id = legacy[legacy["name"]==p]["CID"].tolist()[0]
                    same_pid = legacy[legacy["CID"]==p_id] # the target already hold other synonym compounds
            
                    new_cid = legacy[legacy["name"]==b]["CID"].tolist()[0] # base compound CID
                    new_cids = [new_cid]*(len(same_pid)+1)
                    new_df = pd.concat([legacy[legacy["name"]==b],same_pid]) 
                    new_df["CID"]=new_cids # convert CID to base one
                    update_df = pd.concat([update_df,new_df])
                    update_count += 1
            
            print(update_count,"update its CID")
            #update_df = merge_dup(update_df) # drop duplication
            update_df = update_df.drop_duplicates() # drop duplication
            updated = update_df["name"].tolist()
            other_original = legacy[~legacy["name"].isin(updated)]
            final = pd.concat([other_original,update_df])
            final = final.sort_values("CID",ascending=False)
            final.index = [i for i in range(len(final))]
            print("completed")
            print("")
            return final
    