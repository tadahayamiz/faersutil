# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

Editor of chemicals with regex

@author: tadahaya
"""
import re

class ChemEditor():
    def __init__(self):
        self.pattern = None
        self.strings = []
        self._pre_chars = [
            ".", ",", "?", ":", ";", 
            "(", ")", "[", "]", "{", "}",
            "+", "-", "*", "/", "=", "<", ">", 
            "!", "@", "#", "$", "%", "^", "&", 
            "_", "`",
            ]

        self._pre_compounds = [
            # A
            "acetate", "acetonide", "ammonium", "anhydrous",
            # B
            "benzonate", "benzoate", "besylate", "bitartrate", "bromide", 
            # C
            "calcium", "carbonate", "chloride",  
            # D
            "decanoate", "dequalinium", "diacetate", "diammonium", "dihydrate",
            "dihydrochloride", "dimer", "dipotassium", "dipropionate", "disodium", 
            # E
            # F
            "fumarate",
            # G
            "gadolinium", "gallium",
            # H
            "hemifumarate", "hydrobromide", "hydrochloride", "hydroxide", 
            # I
            "indium", "iodide", 
            # J
            # K
            # L
            "lactate", "lithium", "lutetium", 
            # M
            "malate", "maleate", "magnesium", "methanesulfonate", "methylbromide",
            "mesylate", "monohydrate", "monohydrochloride", 
            # N
            "nitrate", 
            # O
            "orotate", "oxalate", "oxide", 
            # P
            "pamoate", "palmitate", "pentahydrate", "phosphate", "pivalate",
            "pollen", "propionate", "pottasium",
            # Q
            # R
            "rubidium", 
            # S
            "saccharate", "salicylate", "salt", "samarium", "sodium",
            "strontium", "succinate", "sulfate",
            # T
            "tartrate", "technetium", "tetrasodium", "tosylate", "trifluoroacetate",
            "trihydrate", "trifluoroacetic acid", "trisodium",
            # U
            # W
            # X
            # Y
            "yttrium", 
            # Z
            ]
            # 230805 updated

    def set_strings(self, strings:list=[], refresh:bool=True):
        """
        set a list of strings to be compiled
        note these strings are handledsimultaneously
        
        Parameters
        ----------
        strings: list
            a list of strings to be edited

        refresh: bool
            whether to initialize the previous strings or not

        """
        if refresh:
            self.strings = []
        self.strings += strings # to save the previous

    def get_strings(self):
        return self.strings

    def get_special_char_list(self):
        return self._pre_chars

    def get_compound_list(self):
        return self._pre_compounds

    def compile(self, strings:list=[]):
        """
        compile a list of strings
        to handle these strings simultaneously
        
        Parameters
        ----------
        strings: list
            a list of strings to be edited

        """
        if len(strings)==0:
            strings = self.strings
        re_strings = r"{}".format('|'.join(map(re.escape, strings)))
        self.pattern = re.compile(re_strings)


    def delete(self, target:list):
        """
        delete all the matched phrases of the words in the given list
        
        """
        return [self.pattern.sub("", w).strip() for w in target]


    def replace(self, target:list, rep_char:str=" "):
        """
        replace all the matched phrases of the words in the given list
        with the indicated character
        
        """
        return [self.pattern.sub(rep_char, w).strip() for w in target]


    def set_special_chars(self, refresh:bool=False):
        """ set the predefined special characters to be compiled """
        self.set_strings(self._pre_chars, refresh=refresh)


    def set_compounds(self, refresh:bool=False):
        """ set the predefined compounds to be compiled """
        self.set_strings(self._pre_compounds, refresh=refresh)


    def check_ate(self, target:list):
        """ check whether -ate or not """
        pattern = re.compile(r"ate$")
        return [bool(pattern.search(v)) for v in target]
    
    def check_ide(self, target:list):
        """ check whether -ide or not """
        pattern = re.compile(r"ide$")
        return [bool(pattern.search(v)) for v in target]

    def check_ium(self, target:list):
        """
        check whether -ium or not
        note this also catches bacterium
        
        """
        pattern = re.compile(r"ium$")
        return [bool(pattern.search(v)) for v in target]