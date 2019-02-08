#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 17:08:39 2019

@author: max
"""
#%%
import os
import argparse
import sys
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go

sys.path.append(os.path.realpath(__file__))
import KnockdownFeatures_class

def parseArguments():
  # Define the parser and read arguments
  parser = argparse.ArgumentParser(description='Get tags from all files in a directory.')
  parser.add_argument('-d', '--dir', type=str, help='The folder where your experiments are located', required=True)  
  parser.add_argument('-f', '--folders', nargs='+', help='list of folders you want to take into account', required=True)
  args = parser.parse_args()
  return(args)

#add the paths to the experiment folders
path=['/Users/max/Desktop/Office/test/data_test/SiRNA_30/segmented/']
#add the knockdowns you want to load
Knockdowns=['CTRL', 'DLC1', 'ARHGAP17']
class Experiment_data:
    def __init__(self, path, knockdowns):
        self.knockdowns=knockdowns
        self.path=path
        #the upcoming functions will become elements of the class
def load_groups():
    experiment={}
    #for each of the specified knockdowns 
    #create an object from the KnockdownFeatures class at the first path instance
    for i in Knockdowns:
        temp=KnockdownFeatures_class.KnockdownFeatures(path[0], i)
        #for the current object call the objects load_all function to load the features
        temp.load_all()
        #adds the object to a dictionary with the objects experiment identifier and the
        #current knockdown as the key
        experiment.update({temp.experiment_identifier+'_'+i:temp})
    return experiment
#%%

def feature_extraction(feature):
    '''
    feature: input for the feature to extract
    '''
    l=[]
    #calls the load groups function to get all the groups of the experiment
    experiment=load_groups()
    #gets the features based on the features of the first object in the dict
    #features=next(iter(experiment.values())).features
    #for each object in the dict
    for i in experiment:
        #creates a list with each element being a dataframe for the same feature
        #for a different group
        temp=experiment[i].all_features[feature]
        l.append(temp)
    #concatonates the list to a dataframe    
    cross_group_feature = pd.concat(l, axis=0, sort=True)
    return cross_group_feature

def stat_test(feature):
    cross_group_feature=feature_extraction(features[1])
    py.iplot(cross_group_feature['value'][cross_group_feature['KD']=='CTRL'])
    