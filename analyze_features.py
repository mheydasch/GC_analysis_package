# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 11:30:17 2019

@author: max
"""

#!/usr/bin/env python3
import os
import sys
#import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as stats	
#import matplotlib.pyplot as plt
import argparse
import plotly 
#import plotly.plotly as py
import plotly.graph_objs as go
#init_notebook_mode(connected=True)
#import statsmodels.api as sm
#from xattr import xattr
import time
#import subprocess
from plotly import __version__
#from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
print(__version__) # requires version >= 1.9.0
#from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.formula.api import ols
sys.path.append(os.path.realpath(__file__))
import load_data as exp
#from load_data import KnockdownFeatures_class as kd



#%%
#add the paths to the experiment folders
# =============================================================================
# path=['/Users/max/Desktop/Office/test/data_test/SiRNA_31/segmented/']
# #add the knockdowns you want to load
# knockdowns=['CTRL', 'ARHGAP17', 'DOCK10', 'ITSN1']
# =============================================================================

def parseArguments():
  # Define the parser and read arguments
  parser = argparse.ArgumentParser(description='a function including various statistical tools to be applied to the data objects.')
  parser.add_argument('-d','--dir', nargs='+', help='add the directories with spaces between them', required=True)
  parser.add_argument('-k','--kd', nargs='+', help='add the directories with spaces between them', required=True)
  parser.add_argument('-t','--TSNE', help='set True for TSNE output', required=False)

  args = parser.parse_args()
  return(args)
  
#%%
# =============================================================================
# def boxplot(feature, value):
#     #makes a boxplot of the values from one feature, grouped by knockdown
#     ax=sns.catplot(x='experiment', y=value, hue='KD',\
#                    data=data.grouped_features[feature], kind='box')
#     sig, alpha=calc_Bonferroni(feature)
#     plot_median=data.grouped_features[feature].groupby(['experiment', 'KD'])[value].median()
#     nobs=[sig[x][1] for x in sig]
#     axes = ax.axes.flatten()
#     axes[0].set_xlabel(feature)
#     pos=range(len(nobs))
#     sns.FacetGrid.set_xticklabels(ax, nobs)
#     plt.show()   
#    # ax=ax.get_figure()
#     plt.close()
#     return ax
# =============================================================================
    #axes[1].set_title("External")
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' + directory) 

def pyplot(feature, value):
    to_tag=False
    
    #gets the keys to the groups for indexing
    x_data=list(data.grouped_features[feature].groupby(['experiment', 'KD']).groups.keys())
    #gets a group object the groups are referring to
    y_index=data.grouped_features[feature].groupby(['experiment', 'KD'])[value]
    #y_data=data.grouped_features[feature].iloc[list(y_index.groups[x_data[0]])]['value']
    #y_data=data.grouped_features[feature].groupby(['experiment', 'KD']).groups[x_data[1]]
    traces=[]
    sig, alpha=calc_Bonferroni(feature)
    #https://stackoverflow.com/questions/26536899/how-do-you-add-labels-to-a-plotly-boxplot-in-python
    for enum, xd in enumerate(x_data):              
        traces.append(go.Box(
        #list(y_index.groups[xd]) applies the index of one group to the grouped dataframe to obtain
        # a list of indices for that group. This list of indeces is used to index the dataframe, and obtain
        #the value column of it.        
        y=data.grouped_features[feature].iloc[list(y_index.groups[xd])][value],
        name=str(xd),
        #adds the points for each value next to the box
        boxpoints='all',
        #boxpoint='all',
        jitter=0.5,
        whiskerwidth=0.2,
        marker=dict(
            size=2,
        ),
        line=dict(width=1),
        ))

        layout = go.Layout(              
        title=feature,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            dtick=5,
            gridcolor='rgb(255, 255, 255)',
            gridwidth=1,
            zerolinecolor='rgb(255, 255, 255)',
            zerolinewidth=2,
            ),
        margin=dict(
            l=40,
            r=30,
            b=80,
            t=100,
        ),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)',
        showlegend=True
        )
 
    fig = go.Figure(data=traces, layout=layout)
    #counts the number of observations for each group
    count=data.grouped_features[feature].groupby(['experiment', 'KD'])['value'].count()
    for enum, xd in enumerate(x_data):
        sig_index=xd[0]+xd[1]
        #gets the number of observations for the current box
        n=str(count[xd])

        
        #getting the title for each column the following way:
        #getting the sig_index by concatonating the two strings of xd
        #and using this as the key for the bonferrony corrected t_test
        #to obtain the second value, which is the p value
        try:
            p=round(sig[sig_index][1], 4)
            #adds a star if the p value is significant
            if p < alpha:
                p=str(p)
                p=p+'*'
                #marks the plot as being significant
                to_tag=True
            p=str(p)
        #exception, if no p value exists (i.e. for control)
        except:
            p=''   
        
        fig['layout']['annotations']+=tuple([dict(
                    #positions on x axis based on current box
                    x=enum,
                    #positions text based on y axis based on the median of current box
                    y=data.grouped_features[feature].iloc[list(y_index.groups[xd])]['value'].median(),
                    yref='y',                
                    xref='x',
                    text='p: {}<br>n: {}'.format(p, n),
                    showarrow=True,
                    #determines the length of the arrow for the annotation text
                    arrowhead=0,
                    ax=0,
                    ay=-10
                    )])
    if to_tag==True:
        #saves the plot in a different folder, if one or more groups show significance
        sig_folder=os.path.join(path[0], 'significant')
        createFolder(sig_folder)
        file='{}/{}.html'.format(sig_folder,feature)
    else:
        file='{}{}.html'.format(path[0],feature)
    plotly.offline.plot(fig, filename = file, auto_open=False)
        
    return fig

def loop_graph(function, value):
    '''
    creates a graph for each feature
    '''
    for f in data.features:
        function(f, value)
        time.sleep(1)
        
        
#%%    
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
            for g in experiment:
                if g != 'CTRL':
                    temp=stats.ttest_ind(data.grouped_features[f]['value'][(data.grouped_features[f]['KD'] == 'CTRL') & (data.grouped_features[f]['experiment'] == grps.index[enum])], \
                                        data.grouped_features[f]['value'][(data.grouped_features[f]['KD'] == g) & (data.grouped_features[f]['experiment'] == grps.index[enum])])
                    key=grps.index[enum]+g
                    sig.update({key:temp})
    return sig, alpha
#%%
if __name__ == '__main__':
    args=parseArguments()
    path=args.dir  
    knockdowns=args.kd
    TSNE=args.TSNE
    data=exp.Experiment_data(path, knockdowns)
    data.extract_all()
    loop_graph(pyplot, 'value')
    print(args)
    if TSNE=='True':
        data.pca_feature_data()
        data.pca_attribute_data()
        data.save_df(data.wide_feature, path[0], 'wide_feature')
        data.save_df(data.wide_feature, path[0], 'wide_time')
        data.save_df(data.wide_attribute, path[0], 'wide_attribute')
#%%
#pyplot(feature, 'value')
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