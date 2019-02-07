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
    






























#%%    
def file_load(path, file_name, parameter, nandrop=False):
"""
file_load(path, file_name, parameter, nandrop=False):
    path= path where the files are located
    file_name= name of the csv file to be opened
    parameter= name of the measured parameter, this will rename the value column
    nandrop= False by default, will replaces NaNs by 0s and interpolates 0 values
            if True, all nan values will be dropped from the data
    
"""
    for filepath in glob.glob(path + '*{}'.format(identifier)):
        file = os.path.join(filepath, file_name)	
        filename = vars()['filepath'].split('/')[-1]    
        print(filename)
        try: 
            temp=pd.read_csv(file)
            rows, columns = temp.shape
            timer=0 
            temp=temp.dropna(axis=0, how='all', subset=[x for x in list(temp) if x.startswith('cell ')]) 
            for x in range(len(temp)):
                temp.at[x, 'time']=timer
                timer+=15 #Adjust for respective time interval  
            well=re.search(wellpattern, file).group().strip('_')
            temp['well']=well
            fov = re.search(fovpattern, file).group().strip('_' + identifier)
            temp['fov']=fov
            temp=temp.drop(['frame'], axis=1)
            length_l.append(temp) #appends current dataframe to list length_l
        except: 
            print('{} could not be found'.format(file)) 
            next 
    allwells = pd.concat(length_l, sort=True)
    allwells= pd.melt(allwells, id_vars=['well', 'fov', 'time'])
    allwells['unique']=allwells['well']+allwells['fov']+allwells['variable']
    allwells.rename(index=str, columns={'value' : parameter}, inplace=True)
    if nandrop == False:
        allwells=allwells.fillna(0) #filling all NA values
        for cell in allwells.groupby('unique'):
            #for each parameter, for each row, if the value is 0 it will be replaced by the average of it's two
            #adjacent values    
            for i in range(0, len(cell[1])):
                if i < (len(cell[1])-1):
                    i_plus=i+1
                else:
                    i_plus=i
                if i >0:
                    i_minus=i-1
                else:
                    i_minus=i
                    #if the current value of the measurement is 0, but the two adjacent ones are not
                if cell[1][parameter][cell[1][parameter].index[i]]==0 and cell[1][parameter][cell[1][parameter].index[i_minus]]!=0 and cell[1][parameter][cell[1][parameter].index[i_plus]]!=0:
                    new_value=(cell[1][parameter][cell[1][parameter].index[i_minus]]+cell[1][parameter][cell[1][parameter].index[i_plus]])/2
                    #the current value in the dataframe (not the grouped object) will be replaced
                    #by the mean of the two adjacent values
                    allwells.loc[[cell[1][parameter].index[i]],[parameter]]=new_value
    if nandrop == True:
        allwells=allwells.dropna() #removing all lines where no value could be obtained, this could probably already be implemented at the "allwells" level
        allwells=allwells[allwells[parameter]!=0] #activate if you want to ignore zero values        
    return(allwells)