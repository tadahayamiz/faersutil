# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:09:08 2019

FAERSUtil o-- Plot

@author: tadahaya
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multitest

class Plot():
    def __init__(self):
        self.data = pd.DataFrame()


    def set_data(self,data=None,url=""):
        """ setter """
        self.data = data
        if self.data is None:
            if len(url)==0:
                raise ValueError("!! Give data or its url !!")
            else:
                self.data = pd.read_csv(url,index_col=0,sep="\t")


    def forest_plot(
            self, data=None, figsize=(6,4), color="darkblue", title="Forest Plot",
            markersize=15, linewidth=2, fontsize=14, labelsize=14,
            fileout="", dpi=100, alpha=0.7, log=False, forced=False, xmin=1e-1, xmax=None
            ):
        """
        visualize data with forest plot
        
        Parameters
        ----------
        data: dataframe
            dataframe of statistics such as result of "extract_interest"
            
        """
        if data is None:
            if self.data.empty:
                raise ValueError("!! Give data !!")
            else:
                data = self.data
        if data.shape[0] > 100:
            if not forced:
                raise ValueError("!! data is too large: reduce the size or use forced option !!")

        ### data prep
        sample = list(data.index)
        ror = data["ROR"].values
        lower = data["lower_CI"].values
        upper = data["upper_CI"].values
        zipped = zip(sample,ror,lower,upper)
        if xmax is None:
            temp = np.argmax(upper)
            xmax = int((upper[temp] + ror[temp] - 1)*1.1)

        ### visualization
        plt.figure(figsize=figsize)
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        plt.xlim(left=xmin,right=xmax)
        if log:
            plt.xscale("log")
        plt.xlabel('ROR (95% CI)',fontsize=fontsize)
        plt.title(title,fontsize=fontsize)
        plt.tick_params(labelsize=labelsize)
        plt.xticks(list(range(len(sample))), sample)
        for sa, ro, lo, up in zipped:    
            plt.plot([lo, up], [sa, sa], linewidth=linewidth, color=color)
        plt.plot([1.0,1.0], [sample[0], sample[-1]], color='grey', linestyle="dashed")
        plt.plot(
            ror, sample, linestyle='', linewidth=0, color=color,
                marker='o', markersize=markersize, alpha=alpha
                )
        plt.tight_layout()
        if len(fileout) > 0:
            plt.savefig(fileout, dpi=dpi, box_inches='tight')
        plt.show()