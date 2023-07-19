#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 06:32:26 2021

set pure liver injury inducer drug as covariate
perfomr logstic regression and predict the adverse events emergence
focus on the records with high predicted value but no actual adverse events

@author: docker
"""
import pandas as pd
import random
import numpy as np
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
from tqdm import tqdm
import math

class LogisticWeaken():
    def __init__(self):
        self.base_records = pd.DataFrame()
        self.y_prob_df = pd.DataFrame()
        self.final_table = pd.DataFrame()
        
    def set_reference(self,base_records,pure_liid:list):
        self.base_records = base_records
        self.pure_liid = pure_liid
    
    def regression_preparation(self,on_on_df,seed=123):
        # select event off records as to the same size to the event on records
        random.seed(seed)
        r_s = random.sample(self.base_records.index.tolist(),len(on_on_df))
        off_any = self.base_records.loc[r_s]
        
        liid_list = sorted(list(self.pure_liid))
        table = np.zeros((len(off_any),len(liid_list)),dtype=int)
        # time comsuming
        off_any_sub = off_any['Active Substances'].tolist()
        for i in range(len(off_any_sub)):
            tmp = off_any_sub[i]
            for x in tmp:
                if x in liid_list:
                    ix = liid_list.index(x)
                    table[i][ix] = 1
                else:
                    pass
        off_any_df = pd.DataFrame(table,columns=liid_list,index=off_any.index.tolist())
        off_any_df['Event'] = [0]*len(off_any_df)
        
        # concat and create regression target data
        self.reg_data = pd.concat([on_on_df,off_any_df])
    
    def perform_regression(self,do_plot=False,method='lasso'):
        x = self.reg_data.iloc[:,0:-1]
        y = self.reg_data['Event']
        # prepare model
        if method == 'lasso':
            print("-- LASSO regression --")
            lr = LogisticRegression(penalty="l1",solver='liblinear')
            lr.fit(x,y)
        elif method == 'ridge':
            print("-- Ridge regression --")
            lr = LogisticRegression(penalty="l2",solver='liblinear',C=0.01)
            lr.fit(x,y)
        else:
            raise ValueError("!! inappropriate method !!")
        #coef = lr.coef_[0]
        #intercept = lr.intercept_
        # summarize
        #y_pred = lr.predict(x)
        y_prob = lr.predict_proba(x)
        self.y_prob_df = pd.DataFrame(y_prob)
        self.y_prob_df.index = self.reg_data.index.tolist()
        if do_plot:
            # plot the predicted value distribution
            plt.plot(sorted(self.y_prob_df[1].tolist(),reverse=True))
            plt.xlabel('predicted value ranking')
            plt.ylabel('predicted value')
            plt.show()
        else:
            pass
    
    def records_selection(self,target_comp:list):
        """
        target_comp = pd.read_pickle('/mnt/FAERS/workspace/211010_renew/211022_propensity_score_approach/211030_iteration_trial/result/827_target_drug.pkl') : contains the analytical target compounds
        """
        interest_idx = self.y_prob_df[(self.y_prob_df[1]>0.8)].index.tolist()
        interest_df = self.reg_data.loc[interest_idx]
        # no adverse effect
        self.interest_record = self.base_records.loc[interest_df[interest_df['Event']==0].index.tolist()]
        
        self.interest_comp = self.interest_record['Active Substances'].tolist()
        comp_candi = []
        for t in self.interest_record['Active Substances'].tolist():
            comp_candi.extend(t)
        comp_candi = sorted(list(set(comp_candi)-set(self.pure_liid)))
        # reflect target comp
        self.final_candi = sorted(list(set(target_comp)&set(comp_candi)))
        print("focus on the interesting records")
        print(len(interest_df),"records with high predicted value")
        print(len(self.interest_record),"records with high predicted value but no adverse event")
        print("")
    
    def weaken_drug_table(self,freq_dic):
        """
        generate contingency table

        ### contingency table ###
                                    focusing adverse effect  not focusing AE     sum
        
        drug in interesting records           n11                n12             n1p
        
        drug in base whole records            n21                n22             n2p
        
        sum                                   np1                np2             npp
        
        Note that both records are limited to those that dont have the target adverse events

        """
        # calc table
        #freq_dic = pd.read_pickle('/mnt/FAERS/workspace/211010_renew/211022_propensity_score_approach/211030_iteration_trial/result/other_5025_freq_dic.pkl')
        table = []
        ap = table.append
        for c in self.final_candi:
            # count n11
            n11 = 0
            for t in self.interest_comp:
                if c in t:
                    n11 +=1
                else:
                    pass
            n21 = freq_dic.get(c)
            n1p = len(self.interest_record)
            n2p = len(self.base_records)
            n12 = n1p-n11
            n22 = n2p-n21
            np1 = n11+n21
            np2 = n12+n22
            npp = n1p+n2p
            ap([c,n11,n12,n21,n22,n1p,n2p,np1,np2,npp])
        
        table = pd.DataFrame(table,columns=["drug","n11","n12","n21","n22","n1p","n2p","np1","np2","npp"])
        table = table.set_index("drug",drop=True)
        
        # calculate odds ratio
        n11,n12,n21,n22 = table["n11"].values,table["n12"].values,table["n21"].values,table["n22"].values
        ror,lci,uci = calc_odds(n11,n12,n21,n22)
        table["lower_CI"] = lci
        table["upper_CI"] = uci
        table["ROR"] = ror
        
        # calculate IC lower CI
        l1 = table['n11'].tolist()
        l2 = table['n12'].tolist()
        l3 = table['n21'].tolist()
        l4 = table['n22'].tolist()
        
        stat_res = []
        for i in range(len(table)):
            res = calc_ic(l1[i],l2[i],l3[i],l4[i])
            stat_res.append(res)
        table['IC_lowerCI']=stat_res
        table = table.sort_values('IC_lowerCI',ascending=False)
        self.final_table = table


def calc_odds(n11,n12,n21,n22):
    """ calculate confidence interval """
    ror = (n11*n22)/(n12*n21)
    temp = 1.96*np.sqrt(1/n11 + 1/n12 + 1/n21 + 1/n22)
    upperCI = np.exp(np.log(ror) + temp)
    lowerCI = np.exp(np.log(ror) - temp)
    return ror,lowerCI,upperCI

# calculate IC lower CI
def calc_ic(n11,n12,n21,n22):
    """
    return the IC 95% lower confidence interval signal
    if the signal > 0, the signal is statistically significant
    """
    n1p = n11 + n12
    n2p = n21 + n22
    np1 = n11 + n21
    #np2 = n12 + n22
    npp = n1p + n2p
    a1 = 1
    b1 = 1
    a = 2
    b = 2
    r11 = 1
    #r12 = 1
    #r21 = 1
    #r22 = 1
    r = r11*(npp+a)*(npp+b) / ((n1p+a1)*(n1p+b1))
    # calc the IC11 expected value
    e_deno = (n11+r11)*(npp+a)*(npp+b)
    e_nume = (npp+r)*(n1p+a1)*(np1+b1)
    e_ic11 = np.log2(e_deno/e_nume)
    # calc the variance of IC11
    tmp = ((npp-n11+r-r11)/((n11+r11)*(1+npp+r)) + (npp-np1+a-a1)/((n1p+a1)*(1+npp+a)) + (npp-np1+b-b1)/((np1+b1)*(1+npp+b)))
    head = 1/math.log(2)**2
    v_ic11 = head*tmp
    # return the signal
    lower_CI = e_ic11-2*math.sqrt(v_ic11)
    return lower_CI
