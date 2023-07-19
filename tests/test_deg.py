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
        cls.SETUP = pd.read_csv(BASEPATH.replace("tests","tests\\gep.csv"),index_col=0)

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
    def wrapper(self,dic):
        if type(dic)!=dict:
            return False
        elif len(list(dic.keys()))==0:
            return False
        else:
            return True

    # main test
    def test_deg(self):
        # prepare test patterns
        test_patterns = [
            (3,None,None,True,"_",0,"control",False,True), # (arg1, arg2, ..., exception flag, expected result)
            (3,10,None,True,"_",0,"control",False,True), # (arg1, arg2, ..., exception flag, expected result)
            (3,None,500,True,"_",0,"control",False,True), # (arg1, arg2, ..., exception flag, expected result)
            (3,None,None,False,"_",0,"control",False,True), # (arg1, arg2, ..., exception flag, expected result)
            (3,None,None,True,"-",0,"control",True,ValueError), # (arg1, arg2, ..., exception flag, expected result)
            (3,None,None,True,"_",1,"control",False,True), # (arg1, arg2, ..., exception flag, expected result)
            (3,None,None,True,"_",0,"control",False,True), # (arg1, arg2, ..., exception flag, expected result)
            (3,None,None,True,"_",0,"hoge",False,True), # (arg1, arg2, ..., exception flag, expected result)
            ]

        ### loop for sweeping all conditions
        for tfold,tmin,tmax,traw,tsep,tpos,tcon,flag,exp in test_patterns:
            print(tfold,tmin,tmax,traw,tsep,tpos,tcon,flag,exp)
            with self.subTest(fold=tfold,nmin=tmin,nmax=tmax,raw=traw,sep=tsep,position=tpos,control=tcon):
                if flag:
                    with self.assertRaises(exp):
                        self.smpl.prep_deg(data=SampleTest.SETUP,fold=tfold,nmin=tmin,nmax=tmax,
                                           raw=traw,sep=tsep,position=tpos,control=tcon)
                else:
                    dic = self.smpl.prep_deg(data=SampleTest.SETUP,fold=tfold,nmin=tmin,nmax=tmax,
                                    raw=traw,sep=tsep,position=tpos,control=tcon)
                    self.assertEqual(self.wrapper(dic),exp)
