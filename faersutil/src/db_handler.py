# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

handle sqlite

@author: tadahaya
"""

import pandas as pd
import os
import sqlite3
from sqlite3.dbapi2 import OperationalError, ProgrammingError
from contextlib import closing

SEP = os.sep

class DBhandler():
    def __init__(self):
        self.path = ""


    def set_path(self, url:str):
        """ set DB path and load it """
        self.path = url


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
        case table
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
        # check
        if self.check_table("case_table"):
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
            ql = "qualification REAL"
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
        focused = [
            "case_id", "sex", "event_date", "event_country",
            "patient_age", "qualification", "stored_year"
            ]
        with closing(sqlite3.connect(self.path)) as conn:
            df[focused].to_sql("case_table", con=conn, index=False, if_exists=if_exists)


    def make_rxn_table(self, df:pd.DataFrame, if_exists:str=None):
        """
        make reaction table
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
        # check
        if self.check_table("rxn_table"):
            if if_exists is None:
                raise KeyError(
                    "!! rxn_table already exists or indicate 'if_exists' (append or replace) !!"
                    )
        else:
            # preparation
            col = [v.lower() for v in df.columns]
            df.columns = col
            df = df[["rxn_id", "pt", "hlt", "hlgt", "soc"]]
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
        drug dict table

        Parameters
        ----------
        df: pd.DataFrame
            table data that containes below fields:
            - drug_id (derived from concept_id of OHDSI)
            - drug_name (derived from concept_name of OHDSI)
            - representative (indicating whether it's representative compound name or not)

        if_exists: str
            indicates the order if the table already exists

        """
        # check
        if self.check_table("drug_dict"):
            if if_exists is None:
                raise KeyError(
                    "!! drug_dict already exists or indicate 'if_exists' (append or replace) !!"
                    )
        else:
            # preparation
            df.loc[:, "drug_dict_id"] = df.index
            df = df[["drug_dict_id", "key", "value", "representative"]]
            df.columns = ["drug_dict_id", "drug_name", "drug_id", "representative"]
            dri = "drug_dict_id INTEGER PRIMARY KEY AUTOINCREMENT"
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


    def make_drug_rxn_table(self, df:pd.DataFrame, if_exists:str=None):
        """
        drug-reaction cross table

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
        # check
        if self.check_table("drug_rxn_table"):
            if if_exists is None:
                raise KeyError(
                    "!! drug_rxn_table already exists or indicate 'if_exists' (append or replace) !!"
                    )
        else:
            # preparation
            df = df[["case_id", "active_substances", "reactions"]]
            ## stored_year is unnecessary
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


    def check_table(self, name:str=""):
        """ whether the indicated table exists or not """
        with closing(sqlite3.connect(self.path)) as conn:
            cur = conn.cursor()
            try:
                cur.execute(f"SELECT * FROM {name}")
                # will fail if the table does not exist
                flag = True
            except OperationalError:
                flag = False
        return flag
    

# ToDo
# update時はIDがユニークである必要があるため, 入力段階でその辺りコントロールできた方がよい
# make_dbレベルでDBに向けてちょうどよいdfへと加工し, DBHandlerでは与えるだけにする