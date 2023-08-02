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

if os.name == 'nt':
    SEP = "\\"
elif os.name == 'posix':
    SEP = "/"
else:
    raise ValueError("!! Something wrong in OS detection !!")


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

    def make_case_table(self, df:pd.DataFrame):
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
            - stored_year: float, like 2014.3

        """
        # check
        if self.check_table("case_table"):
            raise KeyError("!! case_table already exists !!")
        # preparation
        ci = "case_id INTEGER PRIMARY KEY"
        se = "sex TEXT"
        ed = "event_date INTEGER"
        ec = "event_country TEXT"
        pa = "patient_age INTEGER"
        sy = "stored_year REAL"
        constraint = f"{ci}, {se}, {ed}, {ec}, {pa}, {sy}"
        # prepare table for indicating primary constraint
        with closing(sqlite3.connect(self.path)) as conn:
            cur = conn.cursor()
            cur.execute(f"CREATE TABLE case_table ({constraint})")
            conn.commit()
        # add record
        focused = ["case_id", "sex", "event_date", "event_country", "patient_age", "stored_year"]
        with closing(sqlite3.connect(self.path)) as conn:
            df[focused].to_sql("case_table", con=conn, index=False, if_exists="append")

    def add_case_table(
        self, df:pd.DataFrame, index_label:str="case_id", if_exists:str="append"
        ):
        """ add new records to the case table """
        raise NotImplementedError

    def update_case_table(
        self, df:pd.DataFrame, index_label:str="case_id", if_exists:str="replace"
        ):
        """ update the case table """
        raise NotImplementedError

    def make_rxn_table(self, df:pd.DataFrame):
        """
        make reaction table
        Note this is not used for update the table

        Parameters
        ----------
        df: pd.DataFrame
            table data that contains below fields:
            - PT, str
            - HLT, str
            - HLGT, str
            - SOC, str

        """
        # check
        if self.check_table("rxn_table"):
            raise KeyError("!! rxn_table already exists !!")
        # preparation
        col = [v.lower() for v in df.columns]
        df.columns = col
        df = df[["pt", "hlt", "hlgt", "soc"]]
        ri = "rxn_id INTEGER PRIMARY KEY AUTOINCREMENT"
        pt = "pt TEXT"
        hlt = "hlt TEXT"
        hlgt = "hlgt"
        soc = "soc TEXT"
        constraint = f"{ri}, {pt}, {hlt}, {hlgt}, {soc}"
        # prepare table for indicating primary constraint
        with closing(sqlite3.connect(self.path)) as conn:
            cur = conn.cursor()
            cur.execute(f"CREATE TABLE rxn_table ({constraint})")
            conn.commit()
        # add record
        with closing(sqlite3.connect(self.path)) as conn:
            df.to_sql("rxn_table", con=conn, index=False, if_exists="append")


    def make_drug_table(self, df:pd.DataFrame):
        """
        drug table
        - drug_id
        - drug_name

        Parameters
        ----------

        """
        conn = sqlite3.connect(self.path)

        conn.commit()
        conn.close()


    def drug_rxn_table(self):
        """
        reaction table
        - unique_id
        - record_id
        - drug_id
        - rxn_id

        """
        conn = sqlite3.connect(self.path)

        conn.commit()
        conn.close()

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
# - prepare drug table
# - prepare drug x rxn table