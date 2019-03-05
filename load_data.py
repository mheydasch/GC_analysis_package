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
        '''
        loads the objects for each individual feature
        '''
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
    def pca_feature_data(self):
        '''
        creates a wide format dataframe with the feature data for each cell
        to use for PCA analysis
        '''
        temp=[]
        #loops through features
        for enum, f in enumerate(self.features):
            #computes the median value of the current feature for each group and appends the 
            #resulting dataframe consisting of variavle and median value
            #to the list
            temp.append(self.grouped_features[f].groupby('variable').agg({'value':'median'}))
            #renames the column of the dataframe to the current feature
            temp[enum].rename(columns = {'value':'{}'.format(f)}, inplace = True)
        #concatonates the list to one data frame adding it as an attribute to the object
        self.wide_feature=pd.concat(temp, axis=1, sort=True)
        self.wide_feature=self.wide_feature.fillna(0)

    def pca_attribute_data(self):
        '''
        creates a wide format dataframe with the attributes, experiment and knockdown,
        for each cell.
        to use for PCA analysis
        '''
        kd={}
        exp={}
        exp_kd={}
        #loops through the features
        for f in self.features:
            #prints the current feature to show progress
            print('collecting attributes of feature {}'.format(f))
            #loops through the variables
            for enum, i in enumerate(self.grouped_features[f]['variable']):
                #if the current variable is not already in the dictionary
                if f not in kd:
                    #updates the dictionary with the variable and the knockdown
                    kd.update({i:self.grouped_features[f].loc[enum]['KD']})
                    #updates the dictionary with the variable and the experiment
                    exp.update({i:self.grouped_features[f].loc[enum]['experiment']})
                    comb=str(self.grouped_features[f].loc[enum]['experiment']+self.grouped_features[f].loc[enum]['KD'])
                    exp_kd.update({i:comb})
        #makes dataframes from the two dictionaries            
        temp1=pd.DataFrame.from_dict(kd, orient='index', columns=['knockdown'])
        temp2=pd.DataFrame.from_dict(exp, orient='index', columns=['experiment'])
        temp3=pd.DataFrame.from_dict(exp_kd, orient='index', columns=['exp_kd'])
        #joins the two dataframes and adds them as an attribute to the object
        self.wide_attribute=temp1.join(temp2)
        self.wide_attribute=self.wide_attribute.join(temp3)
        
        
        
        
        
        
        
    def save_df(self, df, path, name):
        '''
        saves a dataframe to a csv
        df= DataFrame
        path= full path where to save
        name= name of the csv file.
        '''
        df.to_csv('{}{}.csv'.format(path, name), index_label='ID')
#%%
