#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 17:08:39 2019

@author: max
"""
#%%
import os
import sys
import pandas as pd
import seaborn as sns

sys.path.append(os.path.realpath(__file__))
import KnockdownFeatures_class

#add the paths to the experiment folders
#path=['/Users/max/Desktop/Office/test/data_test/SiRNA_30/segmented/', '/Users/max/Desktop/Office/test/data_test/SiRNA_31/segmented/']
#add the knockdowns you want to load
#Knockdowns=['CTRL', 'DLC1', 'ARHGAP17']
class Experiment_data:
    '''
    Initialize with the pathname and a list of all knockdowns.
    Creates an object holding all objects of an experiment with their given features.
    Create a dictionary of Knockdowns and features by launching load_groups()
    Create a dataframe of Knockdowns and one feature by launching feature_extraction(feature)
    Create a dictionary mapping each feature to such a dataframe with extract_all()
    '''
    def __init__(self, path, knockdowns):
        self.knockdowns=knockdowns
        self.path=path
        #the upcoming functions will become elements of the class
    
    def info(self):
        print(self.knockdowns, self.path, self.features)
    
    def load_groups(self):
        experiment={}
        #for each of the specified knockdowns 
        #create an object from the KnockdownFeatures class at the first path instance
        for p in self.path:
            print('loading experiment', p)
            for i in self.knockdowns:
                if os.path.isdir(os.path.join(p, i)):
                    print('loading group: ', i)
                    temp=KnockdownFeatures_class.KnockdownFeatures(p, i)
                    #for the current object call the objects load_all function to load the features
                    temp.load_all()
                    #adds the object to a dictionary with the objects experiment identifier and the
                    #current knockdown as the key
                    experiment.update({temp.experiment_identifier+'_'+i:temp})
                    self.features=next(iter(experiment.values())).features
        return experiment

    def feature_extraction(self, feature):
        '''
        feature: input for the feature to extract
        '''
        l=[]
        if self.experiment is None:
            self.experiment=self.load_groups()
        #for each object in the dict
        for i in self.experiment:
            #print('extracting feature: ', feature, 'for group: ', i)
            #creates a list with each element being a dataframe for the same feature
            #for a different group
            temp=self.experiment[i].all_features[feature]
            l.append(temp)
        #concatonates the list to a dataframe    
        cross_group_feature = pd.concat(l, axis=0, sort=True)
        cross_group_feature=cross_group_feature.reset_index(drop=True)
        return cross_group_feature
    
    def extract_all(self):
        '''
        extract all features for the given experiment
        '''
        #calls the load groups function to get all the groups of the experiment
        self.experiment=self.load_groups()
        self.grouped_features={}
        for feature in self.features:       
            #print('extracting feature: ', feature)
            self.grouped_features.update({feature:self.feature_extraction(feature)})
     
        
        
#%%
