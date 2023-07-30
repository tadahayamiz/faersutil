#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 04:43:29 2021

@author: docker
"""
import glob
import numpy as np
import pandas as pd
from itertools import chain
from tqdm import tqdm


class AfterProcessor():
    def __init__(self):
        self.raw_whole = None
    
    def date_norm(self,pkl_path:str,fileout:str):
        """
        pkl_path: e.g. /mnt/FAERS/workspace/211127_renew/211128_whole_def/several_res/*.pkl
            the path that stores the cleansed and procesed date for each quater
        fileout: e.g. '/mnt/FAERS/workspace/211127_renew/211128_whole_def/several_res_datenorm/'
        
        --> datenorm_whole = merge(each pickle file if fileout)
        
        """
        l = glob.glob(pkl_path)
        for path in tqdm(l):
            name = path.split('/')[-1]
            test = pd.read_pickle(path)
            #test = test.fillna(0)
            test['Event Date'] = test['Event Date'].replace("",0)
            test['Event Date'] = test['Event Date'].apply(str)
            test_date = test['Event Date'].tolist()
            date = date_creater(path)
        
            final_dates = []
            ap = final_dates.append
            for i in range(len(test)):
                ap(date_processor(test_date[i],date))
            test['Event Date'] = final_dates # update the  column
            test2 = test.dropna(subset=['Event Date'],axis=0)
        
            pd.to_pickle(test2,fileout+name)
        
    def set_raw(self,raw_whole,datenorm_whole):
        """
        l = glob.glob(xml path)
        for i in range(len(l)):
            # cleansing
            dat = faersutil.FAERS()
            dat.data_cleansing(url=l[i])
            raw = dat.get_data()
            # narrow 1
            dat = faersutil.FAERS()
            dat.set_data(raw)
            dat.narrow_record() # chemical only
            nar = dat.get_data().sort_values('Qualification')
            # narrow 2
            target_list = ['Physician','Pharmacist','Other Health Professional','Lawyer']
            df = nar[nar['Qualification'].isin(target_list)]
            # ID transformation
            dat = faersutil.FAERS()
            dat.set_data(df)
            dat.set_identifier(url='/mnt/FAERS/github/faersutil/faersutil/identifier/dwh/identifier.txt')
            dat.identify_name()
            res = dat.get_data()
            # save
            pd.to_pickle(res,outpath+str(l[i].split('/')[-1])+str('.pkl'))
        
        --> raw_whole = merge(each pickle file)
        
        """
        self.raw_whole = raw_whole
        self.raw_whole.index = [i for i in range(len(self.raw_whole))]
        self.datenorm_whole = datenorm_whole
        self.datenorm_whole.index = [i for i in range(len(self.datenorm_whole))]
        
    def prep_idconvert(self,identifier_txt):
        """
        create CID and its representative compound name dict
        """
        df = identifier_txt
        unique_id = df["CID"].unique().tolist()
        name = []
        for i in unique_id:
            tmp = sorted(df[df["CID"]==i]["name"].tolist())
            tmp2 = [x.split(" ") for x in tmp]
            min_size=100
            rep = None
            for j in range(len(tmp2)):
                if len(tmp2[j]) < min_size:
                    rep = tmp[j]
                    min_size = len(tmp2[j])
                else:
                    pass
            name.append(rep)
        self.id_dic = dict(zip(unique_id,name))
        
    def prep_dateinfo(self,drugbank_df,id_dic=None):
        self.id_dic = id_dic
        """1. focus on compounds registered in DrugBank"""
        total_com = sorted(set(chain.from_iterable(self.raw_whole["Active Substances"])))
        total_name = [self.id_dic.get(x) for x in total_com] # 4309
        # compounds registered in drugbank
        drugbank_overlap = set(drugbank_df["name"]) & set(total_name) # 908
        
        target = drugbank_df[drugbank_df['name'].isin(drugbank_overlap)]
        date = target['date'].tolist() # extract overlap compounds
        norm_date = []
        for t in date:
            x = str(t).split(' ')[0]
            tmp = x.split('-')
            p = tmp[0]+tmp[1]+tmp[2]
            norm_date.append(p)
        idx = target['name'].tolist()
        
        enc_dic = dict(zip(list(self.id_dic.values()),list(self.id_dic.keys())))
        idx_id = [enc_dic.get(x) for x in idx]
        self.db_dic = dict(zip(idx_id,norm_date))
        
        """2. trimm unreliable records"""
        com = self.datenorm_whole["Active Substances"].tolist()
        # multi loop (1.80 it/s) faster
        unreliable_idx = []
        for c in tqdm(list(self.db_dic.keys())):
            idxs = []
            for i in range(len(com)):
                if c in com[i]:
                    idxs.append(i)
                else:
                    pass
            tmp_df = self.datenorm_whole.loc[idxs]
            f = self.db_dic.get(c)
            unreliable = tmp_df[tmp_df['Event Date']<f]
            unreliable_idx.extend(unreliable.index.tolist())
        reliable_idx = set(self.datenorm_whole.index.tolist())-set(unreliable_idx)
        self.reliable = self.datenorm_whole.loc[list(reliable_idx)] # 3319935
        
        """3. focus on comounds not registered in DB"""
        total_com = sorted(set(chain.from_iterable(self.reliable["Active Substances"])))
        com = self.reliable["Active Substances"].tolist()
        first_date = []
        target_size = []
        reported_size = []
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
        self.total_first_dic = dict(zip(total_com,first_date))
        
        """4. concatenate DrugBank and registered Event Date"""
        use_first = sorted(list(set(total_com)-set(list(self.db_dic.keys()))))
        final_k = []
        final_v = []
        for t in use_first:
            final_k.append(t)
            final_v.append(self.total_first_dic.get(t))
        for i,k in enumerate(self.db_dic):
            final_k.append(k)
            final_v.append(self.db_dic.get(k))
        self.final_dic = dict(zip(final_k,final_v))
        
def date_creater(path:str):
    """
    path : /mnt/FAERS/workspace/211010_renew/211010_data_processing/data/faers_data/ADR10Q2.pkl
    """
    tmp = path.split('ADR')[-1].split('.pkl')[0]
    y = '20'+tmp.split('Q')[0]
    q = tmp.split('Q')[1]
    if len(y)==3:
        y = '200'+y
    else:
        pass
    if q == '1':
        q = '0401'
    elif q == '2':
        q = '0701'
    elif q == '3':
        q = '1001'
    elif q == '4':
        q = '0101'
    return y+q

def date_processor(eventdate:str,date:str,remove=['02050409','02190528','10100610','10180302','10190205','15001120','16000223','17001130'],drop=['0','1683','21040724']):
    if eventdate in drop:
        return date
    #elif eventdate in remove:
        #return np.nan
    else:
        if len(eventdate)==8:
            if eventdate < '19000101':
                return np.nan
            else:
                return eventdate
        elif len(eventdate)==6:
            return eventdate+'01'
        elif len(eventdate)==4:
            return eventdate+'0101'
        else:
            print(eventdate)
        