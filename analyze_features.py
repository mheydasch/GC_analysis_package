#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 11:30:17 2019

@author: max
"""

import os
import sys
import seaborn as sns
import numpy as np
import pandas as pd
sys.path.append(os.path.realpath(__file__))
import load_data as exp
from load_data import KnockdownFeatures_class as kd


#%%
#add the paths to the experiment folders
path=['/Users/max/Desktop/Office/test/data_test/SiRNA_30/segmented/', '/Users/max/Desktop/Office/test/data_test/SiRNA_31/segmented/']
#add the knockdowns you want to load
Knockdowns=['CTRL', 'DLC1', 'ARHGAP17', 'DOCK10', 'ITSN1']

data=exp.Experiment_data(path, Knockdowns)
data.extract_all()
#%%
def boxplot(feature):
    #makes a boxplot of the values from one feature, grouped by knockdown
    sns.boxplot(x=data.grouped_features[feature]['KD'], y=data.grouped_features[feature]['value'])
#either computes the MAD (robust==True),
#or the standard deviation(robust==False)
robust=False
def MAD_robust(x):
    if robust==True:
        med=np.median(x)
        dif=[np.abs(i-med) for i in x]
        return np.median(dif)
    else:
        return np.std(x)
#either computes the median (robust==True), or the mean (robust==False)
def Mean_robust(x):
    if robust==True:
        return np.median(x)
    else:
        return np.mean(x)
def calc_mean_features():
    '''
    calculates the mean values of each feature grouped by timepoint and by experiment
    excluding the ctrl
    '''
    mean_features=[]
    for f in data.features:
        temp=pd.DataFrame()
        temp=data.grouped_features[f][data.grouped_features[f]['KD']!='CTRL'].groupby(['timepoint', 'experiment'], as_index=False).agg({'value':[MAD_robust, Mean_robust]})    
        temp['feature']=f
        mean_features.append(temp)
    mean_features = pd.concat(mean_features, axis=0, sort=True)
    mean_features.columns = ["_".join(x) for x in mean_features.columns.ravel()]
    return mean_features
def calc_mean_ctrl():
    '''
    calculates the mean values of each feature grouped by timepoint and by experiment
    only for the ctrl
    '''
    mean_ctrl=[]
    for f in data.features:
        print(f)
        temp=pd.DataFrame()
        temp=data.grouped_features[f][data.grouped_features[f]['KD']=='CTRL'].groupby(['timepoint', 'experiment'], as_index=False).agg({'value':[MAD_robust, Mean_robust]})    
        temp['feature']=f
        mean_ctrl.append(temp)
    mean_ctrl = pd.concat(mean_ctrl, axis=0, sort=True)
    mean_ctrl.columns = ["_".join(x) for x in mean_ctrl.columns.ravel()]
    return mean_ctrl    
    
#%%    
