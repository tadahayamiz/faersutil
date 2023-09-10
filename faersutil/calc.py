# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

FARESUtil

@author: tadahaya
"""

import pandas as pd
import numpy as np
import argparse
import os
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multitest
from tqdm import tqdm

from src.calculator import Calculator
from src.plot import Plot

### setup ###
SEP = os.sep

parser = argparse.ArgumentParser(description='convert FAERS data into a sqlite database')
parser.add_argument('--note', type=str, help='convert FAERS data into a sqlite database')
parser.add_argument(
    'workdir',
    type=str,
    help='working directory that contains outputs of preprocess such as cleansed FAERS data'
    )
args = parser.parse_args()

### main ###
def main():
    raise NotImplementedError

### calculatrion ###
def calculate():
    raise NotImplementedError

### visualization ###
def forest_plot(data=None,figsize=(6,4),color="darkblue",title="Forest Plot",
                markersize=15,linewidth=2,fontsize=14,labelsize=14,
                fileout="",dpi=100,alpha=0.7,log=False,forced=False,xmin=1e-1,xmax=None):
    """
    visualize data with forest plot
    
    Parameters
    ----------
    data: dataframe
        the output of calc()

    forced: boolean
        when the size of data is too large, this method returns an error as default
        forced=True forcibly visualizes data despite size
        
    """
    dat = Plot()
    dat.forest_plot(
        data=data, figsize=figsize, color=color, title=title,
        markersize=markersize, linewidth=linewidth, fontsize=fontsize,
        labelsize=labelsize, fileout=fileout, dpi=dpi, alpha=alpha,
        log=log, forced=forced, xmin=xmin, xmax=xmax
        )
    

if __name__ == '__main__':
    main()     