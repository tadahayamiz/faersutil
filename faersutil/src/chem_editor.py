# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

Editor of chemicals with regex

@author: tadahaya
"""
import re

# deletion
# restriction: leading word with a trailing space, or trailing word with a preceding space
HARD_SET = [
    "hcl", "\, salt",
]


PREFIX_SET = [
    "dl-", "d-", "l-",
]


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
        self._pre_other = [
            "alfa", "capsular", 
        ]

        self._pre_compounds = [
            # A
            "acetate", "aceponate", "acetonide", "adipate", "aluminium",
            "ammonium", "anhydrous", "aspartate", 
            # B
            "benzonate", "benzoate", "benzylate", "besilate", "besylate",
            "bitartrate", "borate", "bromide", "butyrate",
            # C
            "calcium", "camphorsulfonate", "camsylate", "caproate", "carbonate",
            "chloride", "chromium", "citrate", "condensate", "concentrate",
            "conjugate", "cypionate", 
            # D
            "decanoate", "decahydrate", "dequalinium", "diacetate", "diammonium",
            "dihydrate", "dihydrochloride", "dimalonate", "dimer", "dipotassium",
            "diphosphonate", "dipropionate", "disodium", "dodecahydrate", 
            # E
            "edisylate", "embonate", "enantate", "enanthate", "esylate", 
            # F
            "fumarate", "furoate", 
            # G
            "gadolinium", "gallium", "gluceptate", "gluconate", "glycolate",
            "glutamate", "glycinate", 
            # H
            "hemiethanolate", "hemifumarate", "hemihydrate", "heminonahydrate", "heptahydrate",
            "hexahydrate", "hydrate", "hydrobromide", "hydrochloride", "hydroxide", 
            # I
            "indium", "iodide", "isethionate", "isovalerate",
            # J
            # K
            # L
            "lactate", "lactobionate", "lithium", "lutetium", 
            # M
            "malate", "maleate", "magnesium", "magniesium", "methanesulfonate",
            "methylbromide", "mesilate", "mesylate", "monohydrate", "monohydrochloride",
            "monopotassium", "monostearate", 
            # N
            "nitrate", 
            # O
            "octahydrate", "oleate", "orotate", "oxalate", "oxide",
            "oxoglurate",
            # P
            "pamoate", "palmitate", "pentahydrate", "pentetate", "phosphate",
            "phosphide", "pidolate", "pivalate", "pollen", "propionate",
            "potassium", "pottasium",
            # Q
            # R
            "rubidium", "radium", 
            # S
            "saccharate", "salicylate", "salt", "samarium", "silicate",
            "sodium", "stearate", "strontium", "succinate", "sulfate",
            "sulfonate", 
            # T
            "tartrate", "tannnate", "technetium", "teoclate", "terephthalate",
            "tetrabutyrate", "tetrahydrate", "tetrasodium", "thiocyanate", "tosilate",
            "tosylate", "trifluoroacetate", "triflutate", "trihydrate", "trifluoroacetic acid",
            "trisodium",
            # U
            # V
            "valerate", 
            # W
            # X
            # Y
            "yttrium", 
            # Z
            "zinc", "zirconium",
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