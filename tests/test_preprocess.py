# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

@author: tadahaya
"""
import unittest
import pandas as pd
import numpy as np
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

    # called when a test method ends
    def tearDown(self):
        pass

    # wrapper for test if necessary
    def wrapper(self,X):
        if type(X)!=np.ndarray:
            return False
        elif X.shape[0]==0:
            return False
        else:
            judge = math.isnan(X[0,0])
            return not judge

    # main test
    def test_preprocess(self):
        # prepare test patterns
        test_patterns = [
            (True,0.05,True,False,True), # (arg1, arg2, ..., flag, expected result)
            (False,0.05,True,False,True), # (arg1, arg2, ..., flag, expected result)
            (True,0.05,False,False,True), # (arg1, arg2, ..., flag, expected result)
            (False,0.05,False,False,True), # (arg1, arg2, ..., flag, expected result)
            (True,0.2,True,False,True) # (arg1, arg2, ..., flag, expected result)
            ]

        ### loop for sweeping all conditions
        self.smpl.set_data(SampleTest.SETUP) # take care of how to access class variables
        for tmirror,talpha,tsphere,flag,exp in test_patterns:
            with self.subTest(mirror=tmirror,alpha=talpha,sphere=tsphere):
                print(tmirror,talpha,tsphere,flag,exp)
                if flag:
                    with self.assertRaises(exp): # give the error expected
                        self.smpl.preprocess(mirror=tmirror,alpha=talpha,sphere=tsphere) # process expected to fail
                else:
                    self.smpl.preprocess(mirror=tmirror,alpha=talpha,sphere=tsphere)
                    self.assertEqual(self.wrapper(self.smpl.get_processed()["X"]),exp)