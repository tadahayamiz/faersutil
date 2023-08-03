# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

Chem o-- Info o-- Editor

@author: tadahaya
"""
import pandas as pd
import numpy as np
import re

__all__ = ["RawEditor", "DeleteEditor", "ReplaceEditor"]

### indicate deletion keys
DEL_CHAR = [",", ";", ":", "-", "\.", "\(", "\)", "\.alpha\.-", "\?"]

DEL_COMP = [
    "acetate", "acetonide", "ammonium", "anhydrous",
    "benzonate", "benzoate", "besylate", "bitartrate", "bromide", 
    "calcium", "carbonate", "chloride",  
    "decanoate", "dequalinium", "diacetate", "diammonium", "dihydrate",
    "dihydrochloride", "dimer", "dipotassium", "dipropionate", "disodium", 
    "fumarate",
    "gadolinium", "gallium",
    "hemifumarate", "hydrobromide", "hydrochloride", "hydroxide", 
    "indium", "iodide", 
    "lactate", "lithium", "lutetium", 
    "malate", "maleate", "magnesium", "methanesulfonate", "methylbromide",
    "mesylate", "monohydrate", "monohydrochloride", 
    "nitrate", 
    "orotate", "oxalate", "oxide", 
    "pamoate", "palmitate", "pentahydrate", "phosphate", "pivalate",
    "pollen", "propionate", "pottasium",
    "rubidium", 
    "saccharate", "salicylate", "samarium", "sodium", "strontium", "succinate", "sulfate",
    "tartrate", "technetium", "tetrasodium", "tosylate", "trifluoroacetate",
    "trihydrate", "trifluoroacetic acid", "trisodium",
    "yttrium", 
    ]

RE = [r"\s?\((.*?)\)\s?"]

### indicate replacement keys
REPLACE = [".", ",", ";", ":", "-"]


# abstract class
class Editor():
    def __init__(self):
        self.__editor = RawEditor()


    def edit(self, word):
        """ edit the word by the indicated editor """
        return self.__editor.edit(word)


    def add(self, lst):
        """ additional elements to be tested """
        self.__editor.add(lst)


    def to_raw(self):
        self.__editor = RawEditor()


    def to_delete(self):
        self.__editor = DeleteEditor()


    def to_replace(self):
        self.__editor = ReplaceEditor()


# concrete class
class RawEditor():
    def __init__(self):
        pass


    def edit(self,word):
        """ return the word without any modification """
        return [word]


# concrete class
class DeleteEditor():
    def __init__(self):
        self.re_list = []
        self.re_list += RE
        ap = self.re_list.append
        for d in DEL_CHAR:
            ap(r"\s?{}\s?(.*)".format(d))
            ap(r"\s?{}\s?".format(d))
        for d in DEL_COMP:
            ap(r"\s?{}\s?(.*)".format(d))
            ap(r"\s?{}\s?".format(d))


    def edit(self, word):
        """ edit the word by deletion """
        w = r'{}'.format(word)
        res = []
        ap = res.append
        for r in self.re_list:
            repat = re.compile(r)
            temp = repat.sub("",w)
            if (temp!=w) and (len(temp) > 0):
                ap(temp)
        return list(set(res))


    def add(self, lst):
        """ add elements """
        ap = self.re_list.append
        for e in lst:
            ap(r"\s?{}\s?(.*)".format(e))
            ap(r"\s?{}\s?".format(e))


# concrete class
class ReplaceEditor():
    def __init__(self):
        self.replace_list = []
        self.replace_list += REPLACE


    def edit(self, word):
        """ edit the word by replacement """
        res = []
        ap = res.append
        for r in self.replace_list:
            temp = word.replace(r," ")
            if temp!=word:
                ap(temp)
        return list(set(res))


    def add(self, lst):
        """ add elements """
        self.replace_list += lst