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
import scipy.stats as stats	
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.formula.api import ols
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
def boxplot(feature, value):
    #makes a boxplot of the values from one feature, grouped by knockdown
    ax=sns.catplot(x='experiment', y=value, hue='KD',\
                   data=data.grouped_features[feature], kind='box')
    sig=calc_Bonferroni(feature)
    axes = ax.axes.flatten()
    axes[0].set_xlabel(feature)
    plt.show()   
   # ax=ax.get_figure()
    plt.close()
    return ax
    #axes[1].set_title("External")

    
    #data.grouped_features[feature].boxplot('z_score', by='KD', figsize=(12, 8))
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
        temp=data.grouped_features[f][data.grouped_features[f]['KD']!='CTRL'].groupby(['timepoint', 'experiment', 'KD'], as_index=False).agg({'value':[MAD_robust, Mean_robust]})    
        temp['feature']=f
        mean_features.append(temp)
    mean_features = pd.concat(mean_features, axis=0, sort=True)
    mean_features.columns = ["_".join(x) for x in mean_features.columns.ravel()]
    mean_features=mean_features.reset_index(drop=True)
    return mean_features
def calc_mean_ctrl():
    '''
    calculates the mean values of each feature grouped by timepoint and by experiment
    only for the ctrl
    '''
    mean_ctrl=[]
    for f in data.features:
        temp=pd.DataFrame()
        temp=data.grouped_features[f][data.grouped_features[f]['KD']=='CTRL'].groupby(['timepoint', 'experiment'], as_index=False).agg({'value':[MAD_robust, Mean_robust]})    
        temp['feature']=f
        mean_ctrl.append(temp)
    mean_ctrl = pd.concat(mean_ctrl, axis=0, sort=True)
    mean_ctrl.columns = ["_".join(x) for x in mean_ctrl.columns.ravel()]
    mean_ctrl=mean_ctrl.reset_index(drop=True)
    return mean_ctrl    

def calc_z_score():
    '''
    uses each value 
    '''
    #mean_features=calc_mean_features()
    mean_ctrl=calc_mean_ctrl()
    for f in data.features:    
        data.grouped_features[f]['z_score']=''
        print(f)
        for enum, row in  enumerate(data.grouped_features[f]['value']):
            #populating variables for the identifiers of the df
            k=data.grouped_features[f].iloc[enum]['KD']

            e=data.grouped_features[f].iloc[enum]['experiment']
            t=data.grouped_features[f].iloc[enum]['timepoint']
            #creating boolean mask to subset the correct items from the mean feature df
            #feature_mask=((mean_features['feature_']==f) & (mean_features['KD_']==k) & (mean_features['experiment_']==e) & (mean_features['timepoint_']==t))
            #creating boolean mask for the ctrl, same as for knockdown, but without the knockdown parameter
            ctrl_mask=((mean_ctrl['feature_']==f) & (mean_ctrl['experiment_']==e) & (mean_ctrl['timepoint_']==t))
            #print(f, k, e, t, enum)
            #calculating z_score
            z_score=(row-mean_ctrl['value_Mean_robust'][ctrl_mask].values[0])/(mean_ctrl['value_MAD_robust'][ctrl_mask].values[0])
          
            data.grouped_features[f].loc[[enum],['z_score']]=np.float64(z_score)
        data.grouped_features[f]['z_score']= data.grouped_features[f]['z_score'].astype(float)
         
def calc_ANOVA(f):    
    model_name= ols('value ~ C(KD)', data=data.grouped_features[f]).fit()
    return model_name
#https://pythonfordatascience.org/anova-python/
def calc_Tukey(f):
    mc=MultiComparison(data.grouped_features[f]['value'],data.grouped_features[f]['KD'])
    mc_results=mc.tukeyhsd()
    return mc_results

def calc_Bonferroni(f):
    #creating the groups based on the knockdowns for each experiment. 
    #gives a series with the experiment being the index
    grps=data.grouped_features[f].groupby('experiment')['KD'].unique()
    alpha=0.05/len(grps)
    sig={}
    for enum, experiment in enumerate(grps):
        for group in experiment:
            for g in experiment:
                if g !=group:
                    temp=stats.ttest_ind(data.grouped_features[f]['value'][(data.grouped_features[f]['KD'] == group) & (data.grouped_features[f]['experiment'] == grps.index[enum])], \
                                        data.grouped_features[f]['value'][(data.grouped_features[f]['KD'] == g) & (data.grouped_features[f]['experiment'] == grps.index[enum])])
                    key=grps.index[enum]+group+':'+g
                    sig.update({key:temp})
    return sig

#%%    
# =============================================================================
# def calc_man_ANOVA(f):
#     '''
#     manual anova
#     f: feature to analyze
#     '''
#     data.grouped_features[f].boxplot('value', by='KD', figsize=(12, 8))
#     ctrl=data.grouped_features[f]['value'][data.grouped_features[f]['KD']=='CTRL']
#     grps=pd.unique(data.grouped_features[f]['KD'])
#     k=len(grps) #number of conditions
#     N=len(data.grouped_features[f]['value'])#number of observations
#     n=data.grouped_features[f].groupby('KD').size() #number of observations in each group
#     DFbetween=k-1
#     DFwithin=N-k
#     DFtotal=N-1
#     SSbetween = (sum(data.grouped_features[f].groupby('KD')['value'].sum()**2)/n) \
#     - (data.grouped_features[f]['value'].sum()**2)/N
#     sum_y_squared = sum([value**2 for value in data.grouped_features[f]['value'].values])
#     SSwithin = sum_y_squared - sum(data.grouped_features[f].groupby('KD')['value'].sum()**2)/n
#     SStotal = sum_y_squared - (data.grouped_features[f]['value'].sum()**2)/N
#     MSbetween=SSbetween/DFbetween
#     MSwithin=SSwithin/DFwithin
#     F=MSbetween/MSwithin
#     p=stats.f.sf(F, DFbetween, DFwithin)     	
#     eta_sqrd = SSbetween/SStotal
#     om_sqrd = (SSbetween - (DFbetween * MSwithin))/(SStotal + MSwithin)
#     return F, p, om_sqrd 
# =============================================================================

# =============================================================================
# def calc_z_score():
#     '''
#     uses average value per group
#     '''
#     mean_features=calc_mean_features()
#     mean_ctrl=calc_mean_ctrl()
#    
#     mean_features['z_score']=''
#     for enum, row in  enumerate(mean_features['value_Mean_robust']):
#         #populating variables for the identifiers of the df
#         k=mean_features.iloc[enum]['KD_']
#         e=mean_features.iloc[enum]['experiment_']
#         t=mean_features.iloc[enum]['timepoint_']
#         #creating boolean mask to subset the correct items from the mean feature df
#         #feature_mask=((mean_features['KD_']==k) & (mean_features['experiment_']==e) & (mean_features['timepoint_']==t))
#         #creating boolean mask for the ctrl, same as for knockdown, but withouxt the knockdown parameter
#         ctrl_mask=((mean_ctrl['experiment_']==e) & (mean_ctrl['timepoint_']==t))
#         print(k, e, t, enum, row)
#         #calculating z_score
#         z_score=(row-mean_ctrl['value_Mean_robust'][ctrl_mask].values[0])/(mean_ctrl['value_MAD_robust'][ctrl_mask].values[0])
#       
#         mean_features.loc[[enum],['z_score']]=z_score 
#     return mean_features
#         
# =============================================================================