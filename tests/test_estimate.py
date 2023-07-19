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
    SETUP = None

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


    def test_estimate(self):
        # prepare test patterns
        test_patterns = [
            ("parallel",True,False,int), # (arg1, arg2, ..., exception flag, expected result)
            ("smc",True,False,int), # (arg1, arg2, ..., exception flag, expected result)
            ("map",True,False,int), # (arg1, arg2, ..., exception flag, expected result)
            ("accumulation",True,False,int), # (arg1, arg2, ..., exception flag, expected result)
            ("parallel",False,False,int), # (arg1, arg2, ..., exception flag, expected result)
            ("hoge",True,True,KeyError) # (arg1, arg2, ..., exception flag, expected result)
            ]

        ### loop for sweeping all conditions
        for tmethod,tplot,flag,exp in test_patterns:
            with self.subTest(method=tmethod,plot=tplot):
                print(tmethod,tplot,flag,exp)
                if flag:
                    with self.assertRaises(exp):
                        self.smpl.estimate(tmethod,tplot)
                else:
                    self.assertEqual(type(int(self.smpl.estimate(tmethod,tplot))),exp)