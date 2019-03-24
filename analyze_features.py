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

'''
This script is taking a list of folders and a list of knockdowns as the input.
These will be used as input to the load_data.py module.
This script is containing several statistical functions to analyze the data object
created by load_data.py.
It can be run from the command line and is giving the options to run a bonferroni corrected
t-test to compare the knockdowns with the respective control and print figures for all features
with significance annotations. Figures can either be based on raw data, or on the z_score
In addition it gives the option to print .csv files of median feature values to be fed 
into the PCA analysis app.
Dependencies:
    KnockdownFeatures_class.py
    load_data.py
'''

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
  parser.add_argument('-k','--kd', nargs='+', help='add the knockdown folder names with spaces between them', required=True)
  parser.add_argument('-t','--TSNE', help='set True for TSNE output. leave empty to skip this output', required=False)
  parser.add_argument('-f','--figures', help='set True for figure printing of raw data, z_score for figure printing of z_scores. leave empty to skip this output', required=False)


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
def median(a, l, r): 
    n = r - l + 1
    n = (n + 1) // 2 - 1
    return n + l 

# Function to calculate IQR 
def IQR(a, n): 

    a.sort() 

    # Index of median of entire data 
    mid_index = median(a, 0, n) 

    # Median of first half 
    Q1 = a[median(a, 0, mid_index)] 

    # Median of second half 
    Q3 = a[median(a, mid_index + 1, n)] 

    # IQR calculation 
    return Q3

def pyplot(feature, value):
    to_tag=False
    
    #gets the keys to the groups for indexing
    x_data=list(data.grouped_features[feature].groupby(['experiment', 'KD']).groups.keys())
    #gets a group object the groups are referring to
    y_index=data.grouped_features[feature].groupby(['experiment', 'KD'])[value]
    #y_data=data.grouped_features[feature].iloc[list(y_index.groups[x_data[0]])]['value']
    #y_data=data.grouped_features[feature].groupby(['experiment', 'KD']).groups[x_data[1]]
    traces=[]
   # Q3=[]
    sig, alpha=calc_Bonferroni(feature)
    #https://stackoverflow.com/questions/26536899/how-do-you-add-labels-to-a-plotly-boxplot-in-python
    for enum, xd in enumerate(x_data):     
        #Q3.append(IQR(list(data.grouped_features[feature].iloc[list(y_index.groups[xd])][value]), len(data.grouped_features[feature].iloc[list(y_index.groups[xd])][value])))         
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
        autosize=True,
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            dtick=5,
            gridcolor='rgb(255, 255, 255)',
            gridwidth=1,
            zerolinecolor='rgb(255, 255, 255)',
            zerolinewidth=2,
            automargin=True,
            ),
           
# =============================================================================
#         margin=dict(
#             l=40,
#             r=30,
#             b=80,
#             t=3*max(Q3),
#         ),
# =============================================================================
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
def MAD_robust(x, robust=False):
    if robust==True:
        med=np.median(x)
        dif=[np.abs(i-med) for i in x]
        return np.median(dif)
    else:
        return np.std(x)
#either computes the median (robust==True), or the mean (robust==False)
def Mean_robust(x, robust=False):
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
def calc_mean_ctrl(all_features=False):
    '''
    calculates the mean values of each feature grouped by timepoint and by experiment
    only for the ctrl
    '''
    #all features==False used fors tandard z_score. Only calculates the mean and standard deviation for the control
    if all_features==False:
        mean_ctrl=[]
        for f in data.features:
            temp=pd.DataFrame()
            temp=data.grouped_features[f][data.grouped_features[f]['KD']=='CTRL'].groupby(['timepoint', 'experiment'], as_index=False).agg({'value':[MAD_robust, Mean_robust]})    
            temp['feature']=f
            mean_ctrl.append(temp)
        mean_ctrl = pd.concat(mean_ctrl, axis=0, sort=True)
        mean_ctrl.columns = ["_".join(x) for x in mean_ctrl.columns.ravel()]
        mean_ctrl=mean_ctrl.reset_index(drop=True)
    #if all features==True used for internal z_score, computes the mean and standard deviation
    #for all knockdowns
    if all_features==True:
        mean_ctrl=[]
        for k in data.knockdowns:
            for f in data.features:
                temp=pd.DataFrame()
                temp=data.grouped_features[f][data.grouped_features[f]['KD']==k].groupby(['timepoint', 'experiment'], as_index=False).agg({'value':[MAD_robust, Mean_robust]})    
                temp['feature']=f
                temp['knockdown']=k
                mean_ctrl.append(temp)
        mean_ctrl = pd.concat(mean_ctrl, axis=0, sort=True)
        mean_ctrl.columns = ["_".join(x) for x in mean_ctrl.columns.ravel()]
        mean_ctrl=mean_ctrl.reset_index(drop=True)
    return mean_ctrl    

def calc_z_score(internal=False):
    '''
    if internal==False calculates standard z_score using calc_mean_ctrl(all_features=False)
    calculates the z_score as the distance of the values from the mean control.
    
    if internal==True calculates the z_score using calc_mean_ctrl(all_features=True)
    calculates the z_score as the distance of values of each knockdown from it's own mean.
    The sum of these values can be used as a measure for heterogeneity.
    '''
    #mean_features=calc_mean_features()
    if internal==False:
        mean_ctrl=calc_mean_ctrl()
        for f in data.features:    
            data.grouped_features[f]['z_score']=''
            print(f)
            for enum, row in  enumerate(data.grouped_features[f]['value']):
                #populating variables for the identifiers of the df
                #k=data.grouped_features[f].iloc[enum]['KD']
    
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
    if internal==True:
        mean_ctrl=calc_mean_ctrl(all_features=True)
        
        for f in data.features: 
            print(f)
            data.grouped_features[f]['z_score_int']=''
            for k in data.knockdowns:
                print(k)   
                #creating a boolean mask for the knockdowns
                knockdown_mask=data.grouped_features[f]['KD']==k              
                #iterating through the values of current featurer only for the current knockdown
                for enum, row in  enumerate(data.grouped_features[f][knockdown_mask]['value']):
                    #populating variables for the identifiers of the df
                    e=data.grouped_features[f][knockdown_mask].iloc[enum]['experiment']
                    t=data.grouped_features[f][knockdown_mask].iloc[enum]['timepoint']
                    #creating a boolean mask to subset the correct items by feature, experiment, timepoint and knockdown
                    ctrl_mask=((mean_ctrl['feature_']==f) & (mean_ctrl['experiment_']==e) & (mean_ctrl['timepoint_']==t) & 
                               (mean_ctrl['knockdown_']==k))
                    #is getting the index in the original dataframe, as we here work on a boolean subset which changes indexing
                    original_index=data.grouped_features[f][knockdown_mask].iloc[enum][:].name
                    #calculates the z_score as subtracting from the current value, the corresponding mean and dividing this by the corresponding
                    #standard deviation
                    z_score=(row-mean_ctrl['value_Mean_robust'][ctrl_mask].values[0])/(mean_ctrl['value_MAD_robust'][ctrl_mask].values[0])
                    #adds the z_score value to the column 'z_score_int' as a float64
                    data.grouped_features[f].loc[[original_index],['z_score_int']]=np.float64(z_score)
        data.grouped_features[f]['z_score_int']= data.grouped_features[f]['z_score_int'].astype(float)


        
def calc_hetero(norm=True, internal=True):
    '''
    calculates the hetereogeneity.
    Proper use should be with norm=True and internal=True.
    In that case the heterogeneity is expressed as the sum of absolute values
    of z_scores that are calculated of the distance from the mean of the same knockdown.
    if internal=False. Z_score is calculated as the distance from the mean of the control.
    This leads to the output being a measure of the distance from the control, rather
    than of heterogeneity.
    If norm=False the standarddeviation of the unnormalized values will be used instead,
    leading to bias if larger values differ.
    '''
    if norm==False:
        #calculates the total standard deviation for each knock down per feature and the sum of it.
        standev=[]
        for f in data.feature_list:
            temp=pd.DataFrame(data.grouped_features[f].groupby('KD')['value'].std())
            temp=temp.rename(columns={'value':f})
            standev.append(temp)
        standev=pd.concat(standev, axis=1, sort=True)
        standev['total']=standev.sum(axis=1)    
        return standev
    if norm==True:
        if internal==True:
            score='z_score_int'
            if score not in data.grouped_features[data.features[1]]:
                calc_z_score(internal=True)
        if internal==False:
             #calculates the z_score using the calc_z_score() function.
             score='z_score'
             if score not in data.grouped_features[data.features[1]]:
                 calc_z_score()
                    
        #calculates the total z_score (absolute values) for each knock down per feature and the sum of it.
        z_sum=[]
        #loops through the features
        for f in data.feature_list:
            #calculates the sum of the absolute values, of the z_score for each knockdown and counts the
            #amount of values
            temp=pd.DataFrame(data.grouped_features[f].groupby('KD').agg({score:[(lambda c: c.abs().sum()), 'count']}))
            #drops outer index value created by .agg
            temp.columns=temp.columns.droplevel(0)
            #divides the cumulative z_score by the number of z_scores
            temp=pd.DataFrame(temp['<lambda>']/temp['count'], columns=[f])
            #appends temp to a list
            z_sum.append(temp)
        #concatonates the list z_sum into a data frame    
        z_sum=pd.concat(z_sum, axis=1, sort=True)
        #sums up the average z_scores of the features and divides them by the 
        #number of features
        z_sum['total']=z_sum.sum(axis=1)/(z_sum.shape[1]-1)  
        return z_sum
    
def calc_ANOVA(f):
    '''
    f= name of feature
    '''    
    model_name= ols('value ~ C(KD)', data=data.grouped_features[f]).fit()
    return model_name
#https://pythonfordatascience.org/anova-python/
def calc_Tukey(f):
    '''
    f= name of feature
    '''
    mc=MultiComparison(data.grouped_features[f]['value'],data.grouped_features[f]['KD'])
    mc_results=mc.tukeyhsd()
    return mc_results

def calc_Bonferroni(f):
    '''
    f= name of feature
    '''
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
    figures=args.figures
    data=exp.Experiment_data(path, knockdowns)
    print(args)
    data.extract_all()
    if figures=='True':
        loop_graph(pyplot, 'value')
        print('figures are saved at {}'.format(path[0]))
    if figures=='z_score':
        calc_z_score()
        loop_graph(pyplot, 'z_score')
        print('figures are saved at {}'.format(path[0]))
    
    
    if TSNE=='True':
        data.pca_feature_data()
        data.pca_attribute_data()
        data.save_df(data.wide_feature, path[0], 'wide_feature')
        data.save_df(data.wide_feature, path[0], 'wide_time')
        data.save_df(data.wide_attribute, path[0], 'wide_attribute')
        print('csv files for TSNE are saved at{}'.format(path[0]))
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