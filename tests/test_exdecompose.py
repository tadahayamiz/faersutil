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

    # called when a test method ends
    def tearDown(self):
        pass

    # wrapper for test if necessary
    def wrapper(self,df):
        if type(df)!=pd.core.frame.DataFrame:
            return False
        elif df.shape[0]==0:
            return False
        else:
            head = df.head(1)
            judge = math.isnan(head.iat[0,0])
            return not judge


    # initializer
    def initialize(self):
        nfactor = self.smpl.estimate(method="parallel",perm=3,scree_plot=True)
        L,F,C = self.smpl.decompose(initialize="principal_component",rotate="varimax",nfactor=nfactor)
        df_deg = pd.read_csv(BASEPATH.replace("tests","tests\\gep.csv"),index_col=0)
        deg = self.smpl.prep_deg(df_deg)
        return deg,L


    # main test
    def test_exdecompose(self):
        deg,L = self.initialize()

        # prepare test patterns
        test_patterns = [
            (deg,None,False,True), # (arg1, arg2, ..., exception flag, expected result)
            (deg,L,False,True), # (arg1, arg2, ..., exception flag, expected result)
            (dict(),L,True,ValueError) # (arg1, arg2, ..., exception flag, expected result)
            ]

        ### loop for sweeping all conditions
        for tdeg,tL,flag,exp in test_patterns:
            with self.subTest(deg=tdeg,loading=tL):
                print(tdeg,tL,flag,exp)
                if flag:
                    with self.assertRaises(exp):
                        self.smpl.decompose_ex(deg=tdeg,loading=tL)
                else:                   
                    self.assertEqual(self.wrapper(self.smpl.decompose_ex(deg=tdeg,loading=tL)),exp)