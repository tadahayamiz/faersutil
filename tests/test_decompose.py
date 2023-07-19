# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

@author: tadahaya
"""
import unittest
import pandas as pd
import os
import sys
import math

from depro.depro import Depro

BASEPATH = os.path.dirname(os.path.abspath(__file__))

class SampleTest(unittest.TestCase):
    SETUP = None # initially loaded data

    # called when test class initialization
    @classmethod
    def setUpClass(cls):
        print('****** setUpClass method is called. ******')
        cls.SETUP = pd.read_csv(BASEPATH.replace("tests","tests\\large_grm.csv"),index_col=0)

    # called when test class end
    @classmethod
    def tearDownClass(cls):
        print('****** setDownClass method is called. ******')

    # called when a test method runs
    def setUp(self):
        self.smpl = Depro() # initialization must be done
        self.smpl.set_data(SampleTest.SETUP)
        self.smpl.preprocess()
        nfactor = self.smpl.estimate(method="parallel",perm=3,scree_plot=True)

    # called when a test method ends
    def tearDown(self):
        pass

    # wrapper for test if necessary
    def wrapper1(self,df):
        if type(df)!=pd.core.frame.DataFrame:
            return False
        elif df.shape[0]==0:
            return False
        else:
            head = df.head(1)
            judge = math.isnan(head.iat[0,0])
            return not judge

    def wrapper2(self,tpl):
        cL = self.wrapper1(tpl[0])
        cF = self.wrapper1(tpl[1])
        return all([cL,cF])


    # main test
    def test_decompose(self):
        # prepare test patterns
        test_patterns = [
            # (initialize,rotate,nfactor,exflag,expected)
            ("principal_component","varimax",None,False,True), # (arg1, arg2, ..., exception flag, expected result)
            ("principal_component","promax",None,False,True), # (arg1, arg2, ..., exception flag, expected result)
            ("principal_component","varimax",5,False,True), # (arg1, arg2, ..., exception flag, expected result)
            ("hoge","varimax",None,True,KeyError), # (arg1, arg2, ..., exception flag, expected result)
            ("principal_component","hoge",None,False,True) # (arg1, arg2, ..., exception flag, expected result)
            ]

        ### loop for sweeping all conditions
        for tinit,trot,tn,flag,exp in test_patterns:
            with self.subTest(initialize=tinit,rotate=trot,nfactor=tn):
                print(tinit,trot,tn,flag,exp)
                if flag:
                    with self.assertRaises(exp):
                        self.smpl.decompose(initialize=tinit,rotate=trot,nfactor=tn)
                else:                   
                    tpl = self.smpl.decompose(initialize=tinit,rotate=trot,nfactor=tn)
                    print(type(tpl[0]),type(tpl[1]),type(tpl[2]))
                    nf = self.smpl.data.nfactor
                    self.assertEqual(self.wrapper2(tpl),exp)