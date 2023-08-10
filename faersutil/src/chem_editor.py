# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:25:44 2019

Editor of chemicals with regex

@author: tadahaya
"""
import re
import numpy as np
import pandas as pd

# simple deletion
## characters that can be unexpectedly included
SET_SPECIAL_CHAR = [
    '?', '*', '=', '<', '>', 
    '!', '@', '#', '$', '%', '^', '&', 
    '_', '`',
    ]
## not used because some of them affect meanings of compoudns
# SET_ALL_SPECIAL_CHAR = [
#     '.', ',', '?', ':', ';', 
#     '(', ')', '[', ']', '{', '}',
#     '+', '-', '*', '/', '=', '<', '>', 
#     '!', '@', '#', '$', '%', '^', '&', 
#     '_', '`',
#     ]

# deletion with restriction:
#     a part of leading word
SET_LEADING = [
    ".alpha.-", ".beta.-", ".gamma.-", ".delta.-",
    "d-", "l-", "dl-", 
    "(+)-", "(-)-", "(+/-)-",
]

# deletion with restriction:
#     a part of trailing word
SET_TRAILING = [
    ", d-", ", l-", ", dl-", 
    ", (+)-", ", (-)-", ", (+/-)-",
    ", salt", ", sodium salt",
    ", unspecified",
]

# deletion with restriction:
#     leading word with a trailing space, or trailing word with a preceding space
SET_HARD = [
    "anhydrous",
    "capsular", "concentrate", "condensate", "conjugate", 
    "monomer", "dimer", 
    "hcl", "salt", 
]
SET_SALT = [
    "lithium", "sodium", "disodium", "trisodium", "tetrasodium",
    "magnesium", "magniesium", "aluminium", "aluminum", 
    "potassium", "monopotassium", "dipotassium", "tripotassium",
    "calcium", "zinc", "rubidium", "barium", "gadolinium", "gallium",
    "radium", 
    "ammonium", "diammonium",
]
SET_ANION = [
    # A
    "acetate", "aceponate", "acetonide", "adipate", 
     "aspartate", 
    # B
    "benzonate", "benzoate", "benzylate", "besilate", "besylate",
    "bitartrate", "borate", "bromide", "butyrate",
    # C
    "camphorsulfonate", "camsylate", "caproate", "carbonate",
    "chloride",
    # D
    "decanoate", "decahydrate", "diacetate", 
    "dihydrate", "dihydrochloride", "dimalonate", 
    "diphosphonate", "dipropionate", "dodecahydrate", 
    # E
    "edisylate", "embonate", "enantate", "enanthate", "esylate", 
    # F
    "fumarate", "furoate", 
    # G
    "gluceptate", "gluconate", "glycolate",
    "glutamate", "glycinate", 
    # H
    "hemiethanolate", "hemifumarate", "hemihydrate", "heminonahydrate", "heptahydrate",
    "hexahydrate", "hydrate", "hydrobromide", "hydrochloride", "hydroxide", 
    # I
    "iodide", "isethionate", "isovalerate",
    # J
    # K
    # L
    "lactate", "lactobionate", 
    # M
    "malate", "maleate", "methanesulfonate",
    "methylbromide", "mesilate", "mesylate", "monohydrate", "monohydrochloride",
    "monostearate", 
    # N
    "nitrate", 
    # O
    "octahydrate", "oleate", "orotate", "oxalate", "oxide",
    "oxoglurate",
    # P
    "pamoate", "palmitate", "pentahydrate", "pentetate", "phosphate",
    "phosphide", "pidolate", "pivalate", "propionate",    
    # Q
    # R
    # S
    "saccharate", "salicylate", "samarium", "silicate",
    "stearate", "strontium", "succinate", "sulfate",
    "sulfonate", 
    # T
    "tartrate", "tannnate", "technetium", "teoclate", "terephthalate",
    "tetrabutyrate", "tetrahydrate", "thiocyanate", "tosilate",
    "tosylate", "trifluoroacetate", "triflutate", "trihydrate", "trifluoroacetic acid",
    # U
    # V
    "valerate", 
    # W
    # X
    # Y
    # Z
]

def del_characters(target:list=[]):
    """
    edit the target compounds by deleting special characters
    
    """
    ## comopunds
    dat = ChemEditor()
    dat.set_string(string=SET_SPECIAL_CHAR, method="simple_match")
    ## conversion
    dat.compile()
    return dat.delete(target)


def del_compounds(target:list=[]):
    """
    delete hard coding set, salt set, and anion set
    under locational restriction
    
    """
    ## comopunds
    dat = ChemEditor()
    dat.set_string(string=SET_HARD + SET_SALT + SET_ANION, method="leading_word")
    dat.set_string(string=SET_HARD + SET_SALT + SET_ANION, method="trailing_word")
    ## numbers and characters
    dat.set_regex(regex=[r"[a-zA-Z]", r"[0-9]+"], method="leading_word")
    dat.set_regex(regex=[r"[a-zA-Z]", r"[0-9]+"], method="trailing_word")
    ## conversion
    dat.compile()
    return dat.delete(target)


def replace_compounds(target:list=[]):
    """
    replace hard coding set, salt set, and anion set with a single space
    in the middle positions of whole words

    """
    ## comopunds
    dat = ChemEditor()
    dat.set_string(string=SET_HARD + SET_SALT + SET_ANION, method="middle_word")
    ## numbers and characters
    dat.set_regex(regex=[r"[a-zA-Z]"], method="middle_word")
    dat.set_regex(regex=[r"[0-9]+"], method="middle_word")
    ## conversion
    dat.compile()
    return dat.replace(target)


def del_parts(target:list=[]):
    """
    delete leading and trailing set
    under locational restriction
    
    """
    ## comopunds
    dat = ChemEditor()
    dat.set_string(string=SET_LEADING, method="leading_part")
    dat.set_string(string=SET_TRAILING, method="trailing_part")
    ## conversion
    dat.compile()
    return dat.delete(target)


def main(target:list=[]):
    """
    handle the compounds in the given list
    with a predefined flow
    
    """
    tmp = target.copy()
    result = {"Base":tmp}
    # 1. delete harmless characters
    tmp = del_characters(tmp)
    result["1st_del_harmless_chars"] = tmp.copy()
    # 2. compound deletion with leading/trailing word restriction
    tmp = del_compounds(tmp)
    result["2nd_del_comp_lead_trail_word"] = tmp.copy()
    # 3. compound replacement with middle word restriction
    tmp = replace_compounds(tmp)
    result["3rd_rep_comp_mid_word"] = tmp.copy()
    # 4. compound deletion with leading/trailing part restriction
    tmp = del_parts(tmp)
    result["4th_del_comp_lead_trail_part"] = tmp.copy()
    return pd.DataFrame(result)


class ChemEditor():
    def __init__(self):
        self.pattern = None
        self.conved = ""
        self.regex = {
            "simple_match":[], # no restriction
            "leading_word":[], # restricted to the leading word with a trailing space
            "trailing_word":[], # restricted to the trailing word with a preceding space
            "leading_part":[], # restricted to a part of the leading word
            "trailing_part":[], # restricted to a part of the trailing word
            "middle_word":[], # restricted to the middle words
        }
        self.string = {
            "simple_match":[], # no restriction
            "leading_word":[], # restricted to the leading word with a trailing space
            "trailing_word":[], # restricted to the trailing word with a preceding space
            "leading_part":[], # restricted to a part of the leading word
            "trailing_part":[], # restricted to a part of the trailing word
            "middle_word":[], # restricted to the middle words
        }


    def init_regex(self):
        self.regex = {
            "simple_match":[], "leading_word":[], "trailing_word":[], "leading_part":[], 
            "trailing_part":[], "middle_word":[],
        }


    def init_string(self):
        self.string = {
            "simple_match":[], "leading_word":[], "trailing_word":[], "leading_part":[], 
            "trailing_part":[], "middle_word":[],
        }


    def set_string(
            self, string:list=[], method:str="leading_word", init:bool=False
            ):
        """
        set a list of strings
        note these strings are handled simultaneously based on the indicated method
        
        Parameters
        ----------
        strings: list
            a list of strings to be edited

        method: str
            indicates how to apply matching

        init: bool
            whether to initialize the previous strings or not

        """
        if init:
            self.init_string()
            print("init is true")
        try:
            self.string[method] += string
        except KeyError:
            raise KeyError(
                f"!! Wrong method: choose a method from the following {self.string.keys()} !!"
                )


    def set_regex(
            self, regex:list=[], method:str="leading_word", init:bool=False
            ):
        """
        set a list of regulalized expressions
        note these strings are handled simultaneously based on the indicated method
        
        Parameters
        ----------
        regex: list
            a list of strings to be edited

        method: str
            indicates how to apply matching

        init: bool
            whether to initialize the previous strings or not

        """
        if init:
            self.init_regex()
            print("init is true")
        try:
            self.regex[method] += regex
        except KeyError:
            raise KeyError(
                f"!! Wrong method: choose a method from the following {self.regex.keys()} !!"
                )


    def get_string(self):
        return self.string


    def get_regex(self):
        return self.regex


    def get_predefined_list(self):
        tmp = {
            "special_characters":SET_SPECIAL_CHAR,
            "leading_part":SET_LEADING,
            "tailing_part":SET_TRAILING,
            "hard_coding":SET_HARD,
            "salt":SET_SALT,
            "anion":SET_ANION,
        }
        return tmp


    def compile(self):
        """
        compile a list of strings
        to handle these strings simultaneously
        
        Parameters
        ----------
        strings: list
            a list of strings to be edited

        """
        # prepared converted expression
        conv = r""
        k = "simple_match"
        if len(self.regex[k]) > 0:
            conv += r"{}".format('|'.join(self.regex[k])) + r'|' # no escape
        if len(self.string[k]) > 0:
            conv += r"{}".format('|'.join(map(re.escape, self.string[k]))) + r'|' # char needes escape
        k = "leading_word"
        if len(self.regex[k]) > 0:
            conv += r"^{}".format(' |^'.join(self.regex[k])) + r' |'
        if len(self.string[k]) > 0:
            conv += r"^{}".format(' |^'.join(map(re.escape, self.string[k]))) + r' |'
        k = "trailing_word"
        if len(self.regex[k]) > 0:
            conv += r" {}".format('$| '.join(self.regex[k])) + r'$|'
        if len(self.string[k]) > 0:
            conv += r" {}".format('$| '.join(map(re.escape, self.string[k]))) + r'$|'
        k = "leading_part"
        if len(self.regex[k]) > 0:
            conv += r"^{}".format('|^'.join(self.regex[k])) + r'|'
        if len(self.string[k]) > 0:
            conv += r"^{}".format('|^'.join(map(re.escape, self.string[k]))) + r'|'
        k = "trailing_part"
        if len(self.regex[k]) > 0:
            conv += r"{}".format('$|'.join(self.regex[k])) + r'$|'
        if len(self.string[k]) > 0:
            conv += r"{}".format('$|'.join(map(re.escape, self.string[k]))) + r'$|'
        k = "middle_word"
        if len(self.regex[k]) > 0:
            conv += r" {}".format(' | '.join(self.regex[k])) + r' |'
        if len(self.string[k]) > 0:
            conv += r" {}".format(' | '.join(map(re.escape, self.string[k]))) + r' |'
        # compile
        self.pattern = re.compile(conv)
        self.conved = conv


    def delete(self, target:list):
        """
        delete all the matched phrases of the words in the given list
        
        """
        # strip deletes leading or trailing space
        return [self.pattern.sub("", w).strip() for w in target]


    def replace(self, target:list, rep_char:str=" "):
        """
        replace all the matched phrases of the words in the given list
        with the indicated character
        
        """
        return [self.pattern.sub(rep_char, w).strip() for w in target]


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