#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 14:28:37 2021

FAERSUtil o-- Calculator

@author: tadahaya
"""
import pandas as pd
import numpy as np
import sqlite3
from contextlib import closing
from tqdm.auto import tqdm
import time
from scipy import stats
import statsmodels.stats.multitest as multitest
from tqdm import tqdm
import math
import time

class Calculator:
    """
    calculate FAERS statistics

    ### contingency table ###
                        focusing adverse effect  not focusing AE     sum

    focusing drug                 n11                n12             n1p

    not focusing drug             n21                n22             n2p

    sum                           np1                np2             npp

    """
    def __init__(self, db_path:str):
        self.db_path = db_path
        self.tmake = Cond2Text()
        self.__stats = CommonFxn()
        self.drug_id = set()
        self.rxn_id = set()
        self.case_id = set()
        self.table = None


    def make_table(
        self, drug_list:list=[], category:int=1, 
        rxn_list:list=[], layer:str="soc",
        qualification:int=2,
        ):
        """
        reflect ID preprocessing

        Parameters
        ----------

        """
        # 0. init
        start = time.time()
        ## narrow drugs based on category
        self.prep_drug([{"ctype":"ge", "key":"category", "value":category}])
        ## narrow records based on qualification
        self.prep_case([{"ctype":"ge", "key":"qualification", "value":qualification}])
        ## no pre-selection for reactions
        dic_cond = dict()
        dic_cond["drug_id"] = {"ctype":"in", "key":"drug_id", "value":tuple(self.drug_id)}
        dic_cond["case_id"] = {"ctype":"in", "key":"case_id", "value":tuple(self.case_id)}
        # 1. npp
        print("count npp...", end="")
        conds = [dic_cond["case_id"], dic_cond["drug_id"]]
        npp = self._count_any("drug_rxn_table", "case_id", conds, inner_join="AND")
        print("DONE")
        # 2. prep rxns
        print("count np1 and np2...", end="")
        dic_cond["rxn_id"] = {"ctype":"in", "key":"rxn_id", "value":tuple(self._select_rxn(rxn_list, layer))}
        conds.append(dic_cond["rxn_id"])
        np1 = self._count_any("drug_rxn_table", "case_id", conds, inner_join="AND")
        np2 = npp - np1
        conds = conds[:-1] # restore for n1p and n2p
        print("DONE")
        # 3. count n1p and n2p
        print("count drugs...")
        drugs = self._select_drug(drug_list)
        result = dict() # {drug_id: [n11, n12, n1p, n21, n22 n2p, np1, np2, npp]}
        idx = ["n11", "n12", "n1p", "n21", "n22", "n2p", "np1", "np2", "npp"]
        for d in tqdm(drugs):
            tmp = {"ctype":"equal", "key":"drug_id", "value":d}
            conds.append(tmp)
            n1p = self._count_any("drug_rxn_table", "case_id", conds, inner_join="AND")
            n2p = npp - n1p
            conds.append(dic_cond["rxn_id"]) # drug x rsn
            n11 = self._count_any("drug_rxn_table", "case_id", conds, inner_join="AND")
            n12 = n1p - n11
            n21 = np1 - n11
            n22 = n2p - n21
            result[d] = [n11, n12, n1p, n21, n22, n2p, np1, np2, npp]
            conds = conds[:-2] # restore the base
        result = pd.DataFrame(result, index=idx).T
        # summary
        decoder = dict()
        for d in drugs:
            c0 = {"ctype":"equal", "key":"drug_id", "value":d}
            c1 = {"ctype":"equal", "key":"representative", "value":1}
            tmp = self._select_any("drug_dict", "drug_name", [c0, c1])
            assert len(tmp) == 1 # selected result sould be unique
            decoder[d] = tmp.pop()
        ids = list(result.index)
        ids = [decoder[c] for c in ids]
        result.index = ids
        ## convert drug_id to drug_name
        self.table = result
        print("> completed")
        elapsed = time.time() - start
        h, rem = divmod(elapsed, 3600)
        m, s = divmod(rem, 60)
        print(f"elapsed time: {int(h)} hr {int(m)} min {s:.1f} sec")
        return result


    def calc_stat(self, table=None):
        """
        calculate statistics from the contingency table
        
        """
        if table is None:
            if self.res.empty:
                raise ValueError("!! Give table or run make_table() before this process !!")
            else:
                table = self.table.copy()
        n11, n12, n21, n22 = table["n11"].values, table["n12"].values, table["n21"].values, table["n22"].values
        ror = self.__stats.calc_ror(n11, n12, n21, n22)
        lci, uci = self.__stats.calc_ci(n11, n12, n21, n22)
        chi2, p, q_bh, q_holm = self.__stats.calc_chi2(n11, n12, n21, n22)
        table.loc[:, "lower_CI"] = lci
        table.loc[:, "upper_CI"] = uci
        table.loc[:, "ROR"] = ror
        table.loc[:, "chi-square"] = chi2
        table.loc[:, "adjusted_p_value_BH"] = q_bh
        table.loc[:, "adjusted_p_value_Holm"] = q_holm
        table.loc[:, "p_value"] = p
        table = table.sort_values("lower_CI", ascending=False)
        self.table = table


    def integrate_drug(self, name="", drugs=[]):
        """
        merge counts of the indicated drugs
        
        Parameters
        ----------
        name: str
            indicates the drug name after integration

        drugs: list
            indicates the drugs to be integrated
        
        """
        if len(name)==0:
            name = "integrated_drug"
        merge = [v for v in drugs if v in self.drug]
        if len(merge) > 1:
            rest = [v for v in self.drug if v not in merge]
            temp = self.res.loc[merge,:]
            npp = temp.loc[:, "npp"].values[0]
            np1 = temp.loc[:, "np1"].values[0]
            np2 = temp.loc[:, "np2"].values[0]
            n1p = np.sum(temp["n1p"])
            n2p = npp - n1p
            n11 = np.sum(temp["n11"])
            n12 = np.sum(temp["n12"])
            n21 = np1 - n11
            n22 = np2 - n12
            col = ["n11", "n12", "n1p", "n21", "n22", "n2p", "np1", "np2", "npp"]
            rest = self.res.loc[rest, col]
            df = pd.DataFrame([n11, n12, n1p, n21, n22, n2p, np1, np2, npp],index=col,columns=[name]).T
            rest = pd.concat([rest, df], join="inner", axis=0)
            self.calc_stat(rest)
        else:
            print("CAUTION: no integration because the indicated drugs are not enough")


    def _select_any(
        self, name_table:str, name_field:str, conditions:list, inner_join:str="AND"
        ):
        """ inner method to select records """
        # prepare text
        txt = self.tmake.simple_text(conditions, inner_join)
        # select
        with closing(sqlite3.connect(self.db_path)) as conn:
            cur = conn.cursor()
            order = f"SELECT {name_field} FROM {name_table} WHERE {txt}"
            cur.execute(order)
            content = cur.fetchall()
            content = set(map(lambda x: x[0], content))
            return content


    def _count_any(
        self, name_table:str, name_field:str, conditions:list=[], inner_join:str="AND"
        ):
        """ inner method to count records """
        with closing(sqlite3.connect(self.db_path)) as conn:
            cur = conn.cursor()
            if len(conditions) > 0:
                # prepare text
                txt = self.tmake.simple_text(conditions, inner_join)
                order = f"SELECT COUNT (DISTINCT {name_field}) FROM {name_table} WHERE {txt}"
            else:
                order = f"SELECT COUNT (DISTINCT {name_field}) FROM {name_table}"
            cur.execute(order)
            n = cur.fetchall()[0][0]
            return n


    def _select_drug(self, drug_list:list):
        """ select drug based on the given name list """
        tmp = [d.lower() for d in drug_list]
        condition = {"ctype":"in", "key":"drug_name", "value":tuple(tmp)}
        tmp = self._select_any("drug_dict", "drug_id", [condition])
        if len(tmp) == 0:
            raise ValueError("!! No drugs were selected: check the given list !!")
        if len(self.drug_id) == 0:
            return tmp
        else:
            return self.drug_id & tmp


    def _select_rxn(self, rxn_list:list, layer:str="soc"):
        """
        select reactions based on the given rxn list

        """
        tmp = [d.lower() for d in rxn_list]
        if len(rxn_list) == 1:
            condition = {"ctype":"match", "key":layer, "value":tmp[0]}
        else:
            condition = {"ctype":"in", "key":layer, "value":tuple(tmp)}
        tmp = self._select_any("rxn_table", "rxn_id", [condition])
        assert len(tmp) == len(rxn_list)
        if len(self.rxn_id) == 0:
            return tmp
        else:
            return self.rxn_id & tmp


    def prep_drug(self, conditions:list, inner_join:str="AND", if_exists:str="and"):
        """
        select drug based on the given conditions

        conditions: a list of dict
            dict keys (ctype, key, value)

        inner_join: str
            how to concatenate conditions

        """
        tmp = self._select_any("drug_table", "drug_id", conditions, inner_join)
        if if_exists.lower() == "replace":
            self.drug_id = tmp
        elif if_exists.lower() == "and":
            if len(self.drug_id) == 0:
                self.drug_id = tmp
            else:
                self.drug_id = self.drug_id & tmp
        elif if_exists.lower() == "or":
            self.drug_id = self.drug_id | tmp
        else:
            raise KeyError("!! if_exists should be chozen from {replace, and, or} !!")


    def prep_rxn(self, conditions:list, inner_join:str="AND", if_exists:str="and"):
        """
        select reactions based on the given conditions

        Parameters
        ----------
        conditions: a list of dict
            dict keys (ctype, key, value)

        inner_join: str
            how to concatenate conditions

        """
        tmp = self._select_any("rxn_table", "rxn_id", conditions, inner_join)
        if if_exists.lower() == "replace":
            self.rxn_id = tmp
        elif if_exists.lower() == "and":
            if len(self.rxn_id) == 0:
                self.rxn_id = tmp
            else:
                self.rxn_id = self.rxn_id & tmp
        elif if_exists.lower() == "or":
            self.rxn_id = self.rxn_id | tmp
        else:
            raise KeyError("!! if_exists should be chozen from {replace, and, or} !!")


    def prep_case(self, conditions:list, inner_join:str="AND", if_exists:str="and"):
        """
        select cases based on the given conditions

        Parameters
        ----------
        conditions: a list of dict
            dict keys (ctype, key, value)

        inner_join: str
            how to concatenate conditions

        """
        tmp = self._select_any("case_table", "case_id", conditions, inner_join)
        if if_exists.lower() == "replace":
            self.case_id = tmp
        elif if_exists.lower() == "and":
            if len(self.case_id) == 0:
                self.case_id = tmp
            else:
                self.case_id = self.case_id & tmp
        elif if_exists.lower() == "or":
            self.case_id = self.case_id | tmp
        else:
            raise KeyError("!! if_exists should be chozen from {replace, and, or} !!")


class Cond2Text:
    def __init__(self):
        pass


    def nested_text(
        self, conditions:list, inner_joint:list=[], outer_joint:str="AND"
        ):
        """
        convert the given conditions to a query

        Parameters
        ----------
        conditions: a list of dict
            dict {ctype, key, value}

        inner_joint: a list of str
            indicates each inner condition is concatenated with AND or OR
            default, AND

        outer_joint: str
            indicates nested conditions are concatenated with AND or OR

        """
        assert outer_joint in {"AND", "OR", "and", "or"}
        if len(inner_joint)==0:
            inner_joint = ["AND"] * len(conditions)
        txt = []
        for c, j in zip(conditions, inner_joint):
            t = self.simple_text(c, j)
            t = f"({t})"
            txt.append(t)
        txt = f" {outer_joint.upper()} ".join(txt)
        return txt


    def simple_text(self, conditions:list=[], inner_join:str="AND"):
        """
        convert a set of conditions to a text

        Parameters
        ----------
        conditions: a list of dict
            dict keys (ctype, key, value)

        inner_join: str
            how to concatenate conditions

        """
        assert inner_join in {"AND", "OR", "and", "or"}
        txt = []
        for c in conditions:
            txt.append(self._make_text(c))
        return f" {inner_join.upper()} ".join(txt)


    def _make_text(self, condition:dict):
        """ convert a condition to text """
        # init
        assert {"ctype", "key", "value"} == set(condition.keys())
        assert condition["ctype"] in {
            "in", "match", "greater", "g", "greater_equal", "ge", "lower", "l",
            "lower_equal", "le", "equal", "e",
            }
        ctype = condition["ctype"]
        key = condition["key"]
        val = condition["value"]
        # fxn
        if ctype == "in":
            if type(list(val)[0])==str:
                val = set(map(lambda x: f'{x}', val))
                return f"{key} IN {tuple(val)}"
            else:
                return f"{key} IN {val}"
        if ctype == "match":
            return f"{key} = '{val}'"
        if ctype in {"g", "greater"}:
            return f"{key} > {val}"
        if ctype in {"ge", "greater_equal"}:
            return f"{key} >= {val}"
        if ctype in {"l", "lower"}:
            return f"{key} < {val}"
        if ctype in {"le", "lower_equal"}:
            return f"{key} <= {val}"
        if ctype in {"e", "equal"}:
            return f"{key} = {val}"
        raise KeyError # invalid ctype

class CommonFxn():
    def __init__(self):
        self.ror = np.array([])


    def calc_ror(self, n11, n12, n21, n22):
        """ calculate ROR """
        self.ror = (n11*n22)/(n12*n21)
        return self.ror


    def calc_ci(self, n11, n12, n21, n22):
        """ calculate confidence interval """
        if len(self.ror)==0:
            ror = self.calc_ror(n11,n12,n21,n22)
        temp = 1.96*np.sqrt(1/n11 + 1/n12 + 1/n21 + 1/n22)
        upperCI = np.exp(np.log(self.ror) + temp)
        lowerCI = np.exp(np.log(self.ror) - temp)
        return lowerCI, upperCI


    def calc_chi2(self, n11, n12, n21, n22):
        """ calculate chi-squared values """
        contingency = np.stack([n11,n12,n21,n22],axis=1)
        lst_chi = [stats.chi2_contingency(np.reshape(v,(2,2))) for v in contingency]
        chi = np.array([v[0] for v in lst_chi])
        p = np.array([v[1] for v in lst_chi])
        q_bh = multitest.multipletests(p,alpha=0.05,method="fdr_bh")[1]
        q_holm = multitest.multipletests(p,alpha=0.05,method="holm")[1]
        return chi, p, q_bh, q_holm


    def calc_ic(self, n11, n12, n21, n22):
        """
        return the IC 95% lower confidence interval signal
        if the signal > 0, the signal is statistically significant
        """
        n1p = n11 + n12
        n2p = n21 + n22
        np1 = n11 + n21
        #np2 = n12 + n22
        npp = n1p + n2p
        a1 = 1
        b1 = 1
        a = 2
        b = 2
        r11 = 1
        #r12 = 1
        #r21 = 1
        #r22 = 1
        r = r11*(npp+a)*(npp+b) / ((n1p+a1)*(n1p+b1))
        # calc the IC11 expected value
        e_deno = (n11+r11)*(npp+a)*(npp+b)
        e_nume = (npp+r)*(n1p+a1)*(np1+b1)
        e_ic11 = np.log2(e_deno/e_nume)
        # calc the variance of IC11
        tmp = ((npp-n11+r-r11)/((n11+r11)*(1+npp+r)) + (npp-np1+a-a1)/((n1p+a1)*(1+npp+a)) + (npp-np1+b-b1)/((np1+b1)*(1+npp+b)))
        head = 1/math.log(2)**2
        v_ic11 = head*tmp
        # return the signal
        lower_CI = e_ic11-2*math.sqrt(v_ic11)
        return lower_CI


"""
*** Note ***
len([condition]) was faster than sum(condition) and len(list(filter))

"""