# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

handle sqlite

@author: tadahaya
"""

import pandas as pd
import os
import datetime
import sqlite3
from sqlite3.dbapi2 import OperationalError, ProgrammingError
from contextlib import closing

SEP = os.sep

class DBhandler():
    def __init__(self):
        self.path = ""


    def set_path(self, url:str=""):
        """ set DB path and load it """
        assert len(url) > 0
        self.path = url
        self._logging(table="history", note="init")


    def head(self, name:str="", n:int=5):
        """
        cueing a table
        
        Parameters
        ----------
        n: int
            the number of cueing        
        
        """
        assert len(name) > 0
        with closing(sqlite3.connect(self.path)) as conn:
            cur = conn.cursor()
            try:
                cur.execute(f"SELECT * FROM {name}")
                content = cur.fetchall() # list
                for c in content[:n]:
                    print(c)
                print(f">> total {len(content)} records in {name}")
            except OperationalError:
                raise OperationalError(f"!! {name} table does not exist !!")


    def make_case_table(self, df:pd.DataFrame, if_exists:str=None):
        """
        generate case table from clean_xxx.txt
        - case_id
        - age
        - qualification, etc.

        Parameters
        ----------
        df: pd.DataFrame
            table data that contains below fields:
            - case_id: int, unique ID
            - sex: str (F, M, and N correspond to Female, Male, and NaN)
            - event_date: int (YYYYMMDD or 0)
            - event_country: str (like JP, or NaN as N)
            - patient_age: int (age or 0)
            - qualification: float
            - stored_year: float, like 2014.3
        
        if_exists: str
            indicates the order if the table already exists

        """
        # check df
        col = list(df.columns)
        field = [
            "case_id", "sex", "event_date", "event_country",
            "patient_age", "qualification", "stored_year"
            ]
        if col != field:
            try:
                df = df[field]
            except KeyError:
                raise KeyError(f"!! col of df should be {field} !!")
        # check existence
        if self._check_table("case_table"):
            if if_exists is None:
                raise KeyError(
                    "!! case_table already exists or indicate 'if_exists' (append or replace) !!"
                    )
        else:
            # preparation
            ci = "case_id INTEGER PRIMARY KEY"
            se = "sex TEXT"
            ed = "event_date INTEGER"
            ec = "event_country TEXT"
            pa = "patient_age INTEGER"
            ql = "qualification INT"
            sy = "stored_year REAL"
            constraint = f"{ci}, {se}, {ed}, {ec}, {pa}, {ql}, {sy}"
            # prepare table for indicating primary constraint
            with closing(sqlite3.connect(self.path)) as conn:
                cur = conn.cursor()
                cur.execute(f"CREATE TABLE case_table ({constraint})")
                conn.commit()
            # update if_exists
            if_exists = "append"
        # add record
        with closing(sqlite3.connect(self.path)) as conn:
            df.to_sql("case_table", con=conn, index=False, if_exists=if_exists)


    def make_rxn_table(self, df:pd.DataFrame, if_exists:str=None):
        """
        make reaction table from MedDRA_X_X.txt
        Note this is not used for update the table

        Parameters
        ----------
        df: pd.DataFrame
            table data that contains below fields:
            - rxn_id, int
            - PT, str
            - HLT, str
            - HLGT, str
            - SOC, str

        if_exists: str
            indicates the order if the table already exists

        """
        # check df
        col = list(df.columns)
        field = ["rxn_id", "pt", "hlt", "hlgt", "soc"]
        if col != field:
            try:
                col = [v.lower() for v in df.columns]
                df.columns = col
                df = df[field]
            except KeyError:
                raise KeyError(f"!! col of df should be {field} !!")
        # check existence
        if self._check_table("rxn_table"):
            if if_exists is None:
                raise KeyError(
                    "!! rxn_table already exists or indicate 'if_exists' (append or replace) !!"
                    )
        else:
            # preparation
            ri = "rxn_id INTEGER PRIMARY KEY"
            pt = "pt TEXT"
            hlt = "hlt TEXT"
            hlgt = "hlgt TEXT"
            soc = "soc TEXT"
            constraint = f"{ri}, {pt}, {hlt}, {hlgt}, {soc}"
            # prepare table for indicating primary constraint
            with closing(sqlite3.connect(self.path)) as conn:
                cur = conn.cursor()
                cur.execute(f"CREATE TABLE rxn_table ({constraint})")
                conn.commit()
            # update if_exists
            if_exists = "append"            
        # add record
        with closing(sqlite3.connect(self.path)) as conn:
            df.to_sql("rxn_table", con=conn, index=False, if_exists=if_exists)


    def make_drug_dict(self, df:pd.DataFrame, if_exists:str=None):
        """
        make drug dict table from Drug_dict_updated_xxx.txt

        Parameters
        ----------
        df: pd.DataFrame
            table data that containes below fields:
            - drug_dict_id (unique)
            - drug_id (derived from concept_id of OHDSI)
            - drug_name (derived from concept_name of OHDSI)
            - representative (indicating whether it's representative compound name or not)

        if_exists: str
            indicates the order if the table already exists

        """
        # check df
        col = list(df.columns)
        field = ["drug_dict_id", "drug_name", "drug_id", "representative"]
        if col != field:
            try:
                col = [v.lower() for v in df.columns]
                dic = {"key":"drug_name", "value":"drug_id"}
                col = [dic.get(c, c) for c in col]
                df.columns = col
                df = df[field] # sorting
            except KeyError:
                raise KeyError(f"!! col of df should be {field} !!")
        # check existence
        if self._check_table("drug_dict"):
            if if_exists is None:
                raise KeyError(
                    "!! drug_dict already exists or indicate 'if_exists' (append or replace) !!"
                    )
        else:
            # preparation
            dri = "drug_dict_id INTEGER PRIMARY KEY"
            name = "drug_name TEXT"
            did = "drug_id INTEGER"
            rep = "representative INTEGER"
            constraint = f"{dri}, {name}, {did}, {rep}"
            # prepare table for indicating primary constraint
            with closing(sqlite3.connect(self.path)) as conn:
                cur = conn.cursor()
                cur.execute(f"CREATE TABLE drug_dict ({constraint})")
                conn.commit()
            # update if_exists
            if_exists = "append"  
        # add record
        with closing(sqlite3.connect(self.path)) as conn:
            df.to_sql("drug_dict", con=conn, index=False, if_exists=if_exists)


    def make_drug_table(self, df:pd.DataFrame, if_exists:str=None):
        """
        make drug table from Drug_curated_xxx.txt

        Parameters
        ----------
        df: pd.DataFrame
            table data that containes below fields:
            - concept_name (drug_name, str)
            - concept_id (drug_id, unique int)
            - CID (derived from PubChem, int)
            - CanonicalSMILES (derived from PubChem, str)
            - IUPACName (derived from PubChem, str)
            - MolecularFormula (derived from PubChem, str)
            - MolecularWeight (derived from PubChem, float)
            - TPSA (derived from PubChem, float)
            - XLogP (derived from PubChem, float)
            - category (int)            

        if_exists: str
            indicates the order if the table already exists

        """
        # check df
        col = list(df.columns)
        field = [
            "drug_id", "drug_name", "category", "cid", "smiles", "iupacname",
                "molecular_formula", "mw", "tpsa", "xlogp"
            ]
        if col != field:
            try:
                col = [v.lower() for v in df.columns]
                dic = {
                    "concept_name":"drug_name", "concept_id":"drug_id",
                    "canonicalsmiles":"smiles", "molecularformula":"molecular_formula",
                    "molecularweight":"mw"
                    }
                col = [dic.get(c, c) for c in col]
                df.columns = col
                df = df[field] # sorting
            except KeyError:
                raise KeyError(f"!! col of df should be {field} !!")
        # check existence
        if self._check_table("drug_table"):
            if if_exists is None:
                raise KeyError(
                    "!! drug_table already exists or indicate 'if_exists' (append or replace) !!"
                    )
        else:
            # preparation
            did = "drug_dict_id INTEGER PRIMARY KEY"
            nam = "drug_name TEXT"
            cat = "category INTEGER"
            cid = "cid INTEGER"
            smi = "smiles TEXT"
            iup = "iupacname TEXT"
            mfo = "molecular_formula TEXT"
            mwe = "mw REAL"
            tps = "tpsa REAL"
            log = "xlogp REAL"
            constraint = f"{did}, {nam}, {cat}, {cid}, {smi}, {iup}, {mfo}, {mwe}, {tps}, {log}"
            # prepare table for indicating primary constraint
            with closing(sqlite3.connect(self.path)) as conn:
                cur = conn.cursor()
                cur.execute(f"CREATE TABLE drug_table ({constraint})")
                conn.commit()
            # update if_exists
            if_exists = "append"  
        # add record
        with closing(sqlite3.connect(self.path)) as conn:
            df.to_sql("drug_table", con=conn, index=False, if_exists=if_exists)


    def make_drug_rxn_table(self, df:pd.DataFrame, if_exists:str=None):
        """
        make drug-reaction cross table from drug_rxn_XXX.txt

        Parameters
        ----------
        df: pd.DataFrame
            table data that containes below fields:
            - active_substances
            - reactions
            - case_id
            - stored_year

        if_exists: str
            indicates the order if the table already exists

        """
        # check df
        col = list(df.columns)
        field = ["case_id", "active_substances", "reactions"]
        ## stored_year is unnecessary
        if col != field:
            try:
                col = [v.lower() for v in df.columns]
                df.columns = col
                df = df[field] # sorting
            except KeyError:
                raise KeyError(f"!! col of df should be {field} !!")
        # check existence
        if self._check_table("drug_rxn_table"):
            if if_exists is None:
                raise KeyError(
                    "!! drug_rxn_table already exists or indicate 'if_exists' (append or replace) !!"
                    )
        else:
            # preparation
            cid = "case_id INTEGER PRIMARY KEY AUTOINCREMENT"
            did = "drug_id INTEGER"
            rid = "rxn_id INTEGER"
            constraint = f"{cid}, {did}, {rid}"
            # prepare table for indicating primary constraint
            with closing(sqlite3.connect(self.path)) as conn:
                cur = conn.cursor()
                cur.execute(f"CREATE TABLE drug_rxn_table ({constraint})")
                conn.commit()
        # update if_exists
        if_exists = "append" 
        # add record
        with closing(sqlite3.connect(self.path)) as conn:
            df.to_sql("drug_rxn_table", con=conn, index=False, if_exists=if_exists)


    def make_qulification_table(self, df:pd.DataFrame, if_exists:str=None):
        """
        make drug-reaction cross table from drug_rxn_XXX.txt

        Parameters
        ----------
        df: pd.DataFrame
            table data that containes below fields:
            - active_substances
            - reactions
            - case_id
            - stored_year

        if_exists: str
            indicates the order if the table already exists

        """
        raise NotImplementedError


    def _check_table(self, name:str=""):
        """ whether the indicated table exists or not """
        with closing(sqlite3.connect(self.path)) as conn:
            cur = conn.cursor()
            try:
                cur.execute(f'SELECT * FROM {name}')
                # will fail if the table does not exist
                flag = True
            except OperationalError:
                flag = False
        return flag
    

    def _logging(self, table:str="", note:str=""):
        """
        save history
        
        """
        if self._check_table("history"):
            pass
        else:
            # preparation
            hid = "history_id INTEGER PRIMARY KEY AUTOINCREMENT"
            dat = "date INTEGER"
            nam = "table TEXT"
            des = "description TEXT"
            constraint = f"{hid}, {dat}, {nam}, {des}"
            # prepare table for indicating primary constraint
            with closing(sqlite3.connect(self.path)) as conn:
                cur = conn.cursor()
                cur.execute(f"CREATE TABLE history ({constraint})")
                conn.commit()
        # add record
        with closing(sqlite3.connect(self.path)) as conn:
            now = datetime.datetime.now().strftime('%Y%m%d')
            cur = conn.cursor()
            order = "INSERT INTO history (date, table, description) VALUES (?, ?, ?)"
            tmp = (int(now), table, note)
            cur.execute(order, tmp)
            conn.commit()

# ToDo
# - check_tableにfieldも表示する
# 更新のtableを入れたい