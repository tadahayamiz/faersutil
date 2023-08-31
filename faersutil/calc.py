# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

FARESUtil

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
from tqdm import tqdm
import pickle

from .xml_loader import XMLoader
from .cleanser import Cleanser
from .calculator import Calculator
from .identifier.name_identifier import Chem
from .plot import Plot

SEP = os.sep

class FAERS():
    def __init__(self):
        self.__load = XMLoader()
        self.__clean = Cleanser()
        self.__identify = Chem()
        self.__calc = Calculator()
        self.__plot = Plot()
        self.data = pd.DataFrame()
        self.res = pd.DataFrame()


    def set_data(self,data=None,url=""):
        """ setter for pickle data """
        self.data = data
        if self.data is None:
            if len(url)==0:
                raise ValueError("!! Give data or its url !!")
            else:
                with open(url,"rb") as f:
                    self.data = pickle.load(f)
        self.__clean.set_data(data=self.data)
        self.__calc.set_data(data=self.data)


    def get_data(self):
        """ getter """
        return self.data


    def get_private(self):
        """ getter of private instances """
        print("return (Cleanser, Chem, Calculator)")
        return self.__clean,self.__identify,self.__calc


    ############ data preparation ############
    def data_cleansing(self,url,format="2014Q3-",fileout="",to_pickle=True,to_csv=False):
        """
        a module for converting the large xml files of FAERS data to cleansed one
        take care the date of data (data format), which affects loading algorithm

        1) download zip files from FDA (https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html)
        2) thaw a file into xml
        3) add the path to the below form
        4) run
        
        !! Caution !!
        relatively heavy task

        Parameters
        ----------
        url: str
            a path for xml files from FAERS

        format: str
            indicates the data format

        to_pickle: boolean
            whether result is exported as pickle
            
        to_csv: boolean
            whether result is exported as csv

        """
        if len(url)==0:
            raise ValueError("!! Give url for the directory containing FAERS xml files !!")
        if format=="2014Q3-":
            self.__load.to_2014Q3_()
        elif format=="2012Q4-2014Q2":
            self.__load.to_2012Q4_2014Q2()
        else:
            raise KeyError("!! Indicate correct format: {}!!".format(self.__load.get_format()))
        if (to_pickle==True) or (to_csv==True):
            if len(fileout)==0:
                fileout = url + SEP + "cleansed.pkl"
        p = Path(url)
        filelist = list(map(lambda x: x.as_posix(),list(p.glob("*.xml"))))
        print("xml files (total {}):".format(len(filelist)))
        for f in filelist:
            print("    {}".format(f))
        print("")
        cleansed = []
        n_record = []
        for v in tqdm(filelist):
            self.__clean.data_cleansing(data=self.__load.load_xml(url=v))
            cleansed.append(self.__clean.get_data())
            n_record.append(self.__clean.get_record())
        self.__load = XMLoader() # memory care
        self.data = pd.concat(cleansed)
        del cleansed
        bef = 0
        for v in n_record:
            bef += v[0]
        print("total {} records were obtained".format(self.data.shape[0]))
        print("({} records were dropped by cleansing)".format(bef - self.data.shape[0]))
        if to_pickle:
            self.data.to_pickle(fileout)
        else:
            self.data.to_pickle(fileout)
            # self.data.to_csv(fileout,index=False,encoding="us-ascii")
            """
            something wrong in encoding
            probably due to mixing byte and str type

            """ 

    def narrow_record(self,scope="chemical",fileout=""):
        """
        narrow down records by scope of the analysis

        Parameters
        ----------
        scope: str
            choose one scope from the follwowings:
                - chemical
                    only chemicals listed in PubChem
                - chemical and protein
                    chemical + protein drugs such as antibody drugs
                - all
                    all records detected in FAERS

        """
        self.__clean.set_exception()
        if scope=="all":
            print("!! 'all' was selected: no processing !!")
        elif scope=="chemical":
            self.__clean.exclude_exception(chem_only=True)
            print("narrow down records: chemicals only")            
        elif scope=="chemical and protein":
            self.__clean.exclude_exception(chem_only=False)
            print("narrow down records: chemicals and proteins")            
        else:
            print("!! Wrong scope: choose 'chemical' or 'chemical and protein' !!")        
        n_record = self.__clean.get_record()
        print("BEFORE: {}".format(n_record[0]))
        print("AFTER: {}".format(n_record[1]))
        self.data = self.__clean.get_data()
        self.__calc.set_data(self.data)
        
        # export
        if len(fileout)!=0:
            self.data.to_pickle(fileout)

    # 230727 基本使わない, warning出るが無視 -> ToDo この点を改修
    # def set_identifier(self,data=None,url=""):
    #     """ set identifier for name identification """
    #     if len(url)==0:
    #         url = __file__.replace("faersutil.py",f"identifier{SEP}dwh{SEP}identifier.txt")
    #     self.__identify.set_identifier(data,url)


    def identify_name(self,fileout=""):
        """
        conduct name identification and convert names into IDs
        
        Parameters
        ----------
        fileout: str
            indicates the ouput path

        """
        if len(self.__identify.encoder)==0:
            whole_drug = set(chain.from_iterable(list(self.data["Active Substances"])))
            
            # check the existed data
            self.__identify.set_info()
            self.__identify.prep_ref()
            self.__identify.prep_identifier(list(whole_drug))
            nf = self.__identify.get_not_found()

            # gather information
            if len(nf) > 0:
                temp,nf = self.__identify.get_info(list(nf.keys()))
                self.__identify.update_info(whole_update=True)
            
            # prepare name identification
            new_info = self.__identify.info.get_data().copy()
            self.__identify = Chem() # initialization
            self.__identify.set_info(data=new_info)
            self.__identify.prep_ref()
            self.__identify.prep_identifier(list(whole_drug))

        # name identification
        self.data["Active Substances"] = self.data["Active Substances"].map(lambda x: self.__identify.encode_set(x))
        #self.data["Active Substances"] = self.data["Active Substances"].map(lambda x: self.__identify.decode_set(x))

        # export
        if len(fileout)!=0:
            self.data.to_pickle(fileout)
        self.__calc.set_data(data=self.data) # set data for calculation


    ############ calculation ############
    def check_drug(self,interest=[]):
        """ check whether drugs related to the indicated keywords """
        self.__calc.check_drug(interest)


    def check_rxn(self,interest=[],layer=["SOC","HLGT","HLT"]):
        """
        check whether adverse events related to the indicated keywords
                
        """
        if self.__calc.meddra.empty:
            self.__calc.set_meddra()
        self.__calc.check_rxn(interest,layer)


    def set_rxn(self,interest=[],layer="SOC"):
        """
        set PT (adverse events) based on the indicated SOC, HLGT, or HLT

        Parameters
        ----------
        interest: str
            indicates the adverse events of interest

        """
        if self.__calc.meddra.empty:
            self.__calc.set_meddra()
        self.__calc.set_rxn(interest,layer)


    def calc(self,fileout=""):
        """
        calculate the relationship between drugs and the focused PT
        take a long time

        Parameters
        ----------
        fileout: str
            indicates the output path

        """
        self.__calc.generate_table()
        try:
            self.__calc.calc_stat()
        except Exception as e:
            print("!! an ERROR occurs in calc_stat() !!",e.args)
            print("returns the contingency table")
        self.res = self.__calc.res
        if len(fileout)!=0:
            self.res.to_csv(fileout,sep="\t")
        return self.res


    ############ visualization ############
    def forest_plot(self,data=None,figsize=(6,4),color="darkblue",title="Forest Plot",
                    markersize=15,linewidth=2,fontsize=14,labelsize=14,
                    fileout="",dpi=100,alpha=0.7,log=False,forced=False,xmin=1e-1,xmax=None):
        """
        visualize data with forest plot
        
        Parameters
        ----------
        data: dataframe
            the output of calc()

        forced: boolean
            when the size of data is too large, this method returns an error as default
            forced=True forcibly visualizes data despite size
            
        """
        self.__plot.forest_plot(
            data=data, figsize=figsize, color=color, title=title,
            markersize=markersize, linewidth=linewidth, fontsize=fontsize,
            labelsize=labelsize, fileout=fileout, dpi=dpi, alpha=alpha,
            log=log, forced=forced, xmin=xmin, xmax=xmax
                    )