#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 11:19:25 2019

@author: max
"""
import os
import re
import functools
import pandas as pd
import numpy as np

class KnockdownFeatures:
    '''
    class that holds all the features of one knockdown
    Usage: Initialize with the path to the experiment and a list of the folder names
    you want to search for features
    either use find_csv() and then load_features(feature) to load an
    individual feature, or use load_all() to load all features.
    
    '''
    def __init__(self, kd_path, KD):
        self.kd_path=kd_path
        self.KD=KD
        self.KD_pattern=re.compile('({}){}[0-9]+[^{}]'.format(self.KD, os.sep, os.sep))
        self.Experiment_pattern=re.compile('SiRNA_[0-9]+')
        self.experiment_identifier=re.search(self.Experiment_pattern, self.kd_path).group()
    def info(self):
        '''
        prints info about this instance of the class
        '''
        print('experiment folder: \n', self.kd_path)
        print('experiment name: \n', self.experiment_identifier)
        print('Knockdown: \n', self.KD)
        


    def find_csv(self):
        '''
        returns a list with all csv files in the KD folder
        '''
        #pattern: must end with '.csv'
        csv_find=re.compile('\.csv$')
        #finds the directory with that specific name
        find_dir='GCAFeatureExtraction'
        self.i_dirs=[]
        #Knockdown pattern, must match inserted variable self.__KD followed by '/'
        #and one or more digits
        for root, dirs, files in os.walk(self.kd_path):
            #looks in each folder from given path
            #if folder is matching the KD pattern, the find_dir pattern and
            #if there are csv files
            if re.search(self.KD_pattern, root) and find_dir in root and len([x for x in files if re.search(csv_find, x)])!=0:
                #finds the csv files in the folder and adds them to a list
                csv_files=[x for x in files if re.search(csv_find, x)]
                for csv in csv_files:
                    #each csv file is added to the path and appended to a list
                    self.i_dirs.append(os.path.join(root, csv))            
        return self.i_dirs
   
    def get_features(self):
        '''
        creates a list of all the feature names
        calls find_csv to get i_dirs variable
        for interactive feature selection, change here!
        '''
        #pattern: must start with a '/' followed by characters that are not '/'
        #and end with '.csv'
        self.features=[]
        csv_pattern=re.compile('[{}][^{}]+\.csv$'.format(os.sep, os.sep))
        self.i_dirs=self.find_csv()
        for file in self.i_dirs:
            #making tuples for what to replace with what
            repls = ('.csv', ''), (os.sep, '')     
            filename=re.search(csv_pattern, file).group()
            #applies repls to replace name parts
            filename=functools.reduce(lambda a, kv:a.replace(*kv), repls, filename)
            if filename not in self.features:
                self.features.append(filename)
                
        return self.features
  

  
    def load_feature(self, feature):
        '''
        loads all csvs of a single feature
        needs to be called by load_all
        '''
        GC_list=[]
        time_pattern=re.compile('n[0-9]+')
        for file in self.i_dirs:
            if feature in file:
                identifier=re.search(self.KD_pattern, file).group() 
                try:
                    temp=pd.read_csv(file, header=None)
                except:
                    break
                rows, columns=temp.shape
                num_identifier=[]
                #creates a list of identfiers to match the columns
                for i in range(0, columns):                
                    num_identifier.append(self.experiment_identifier+'/'+identifier+'n'+str(i))  
                #renames columns with identifier    
                temp.columns=num_identifier
                #adding each loaded file to a list
                GC_list.append(temp)
        #concatonating the list to a dataframe    
        full_feature = pd.concat(GC_list, axis=1, sort=True)
        #creating an index for melting
        rows, columns=full_feature.shape
        full_feature['meltid']= range(0, rows)
        #melting  to long format
        long_feature=pd.melt(full_feature, id_vars='meltid')
        #dropping NAs
        long_feature=long_feature.dropna()
        long_feature=long_feature.reset_index(drop=True)
        #dropping all values that cannot be converted to a float
        for n, i in enumerate(long_feature['value']):
            try:
                long_feature.loc[n, 'value']=float(i)
            except ValueError:
                long_feature=long_feature.drop(n)
        #reset index again        
        long_feature=long_feature.reset_index(drop=True)
        long_feature['experiment']=self.experiment_identifier
        long_feature['KD']=self.KD
        long_feature['item']='placeholder'
        long_feature['timepoint']='placeholder'
        long_feature['value']=long_feature['value'].astype('float')
        for n, var in enumerate(long_feature['variable']):
            long_feature.loc[[n],['item']]=re.search(self.KD_pattern, var).group()
            long_feature.loc[[n],['timepoint']]=re.search(time_pattern, var).group()
        return long_feature 
            
    def load_all(self):
        '''
        loops over load_feature for each feature
        calls get_features to get the features and i_dirs
        '''
        self.get_features()
        self.all_features={}
        for i in self.features:
                #for each element in feature the load_feature function is called
                #and its output is added to a dictionary with the feature as a key
                feature=self.load_feature(i)
                self.all_features.update({i:feature})                 
        return self.all_features
    




























