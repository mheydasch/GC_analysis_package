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
    '''
    def __init__(self, kd_path, KD):
        self.__kd_path=kd_path
        self.__KD=KD
        self.__KD_pattern=re.compile('({})/[0-9]+'.format(self.__KD))
    def info(self):
        '''
        prints info about this instance of the class
        '''
        print('experiment folder: \n', self.__kd_path)
        print('Knockdown: \n', self.__KD)
    def find_csv(self):
        '''
        returns a list with all csv files in the KD folder
        '''
        #pattern: must end with '.csv'
        csv_find=re.compile('\.csv$')
        #finds the directory with that specific name
        find_dir='GCAFeatureExtraction'
        i_dirs=[]
        #Knockdown pattern, must match inserted variable self.__KD followed by '/'
        #and one or more digits
        for root, dirs, files in os.walk(self.__kd_path):
            #looks in each folder from given path
            #if folder is matching the KD pattern, the find_dir pattern and
            #if there are csv files
            if re.search(self.__KD_pattern, root) and find_dir in root and len([x for x in files if re.search(csv_find, x)])!=0:
                #finds the csv files in the folder and adds them to a list
                csv_files=[x for x in files if re.search(csv_find, x)]
                for csv in csv_files:
                    #each csv file is added to the path and appended to a list
                    i_dirs.append(os.path.join(root, csv))            
        return i_dirs
   
    def get_features(self):
        '''
        creates a list of all the feature names
        calls find_csv to get i_dirs variable
        for interactive feature selection, change here!
        '''
        #pattern: must start with a '/' followed by characters that are not '/'
        #and end with '.csv'
        features=[]
        csv_pattern=re.compile('[/][^/]+\.csv$')
        i_dirs=self.find_csv()
        for file in i_dirs:
            #making tuples for what to replace with what
            repls = ('.csv', ''), ('/', '')      
            filename=re.search(csv_pattern, file).group()
            #applies repls to replace name parts
            filename=functools.reduce(lambda a, kv:a.replace(*kv), repls, filename)
            if filename not in features:
                features.append(filename)
        return features, i_dirs
  

  
    def load_feature(self, feature, i_dirs):
        '''
        loads all csvs of a single feature
        needs to be called by load_all
        '''
        GC_list=[]
        for file in i_dirs:
            if feature in file:
                identifier=re.search(self.__KD_pattern, file).group() 
                temp=pd.read_csv(file, header=None)
                rows, columns=temp.shape
                num_identifier=[]
                #creates a list of identfiers to match the columns
                for i in range(0, columns):                
                    num_identifier.append(identifier+'_'+str(i))  
                #renames columns with identifier    
                temp.columns=num_identifier
                #adding each loaded file to a list
                GC_list.append(temp)
                print(GC_list)
        #concatonating the list to a dataframe    
        full_feature = pd.concat(GC_list, axis=1, sort=True)
        #creating an index for melting
        rows, columns=full_feature.shape
        full_feature['meltid']= range(0, rows)
        #melting  to long format
        long_feature=pd.melt(full_feature, id_vars='meltid')
        #dropping NAs
        long_feature=long_feature.dropna()
        return long_feature, GC_list
    
    def load_all(self):
        '''
        loops over load_feature for each feature
        calls get_features to get the features and i_dirs
        '''
        features, i_dirs = self.get_features()
        self.all_features={}
        for i in features:
                #for each element in feature the load_feature function is called
                #and its output is added to a dictionary with the feature as a key
                feature, GC_list=self.load_feature(i, i_dirs)
                self.all_features.update({i:feature})                 
        return self.all_features
    




























