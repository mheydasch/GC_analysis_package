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
sys.path.append(os.path.realpath(__file__))
import KnockdownFeatures_class
def parseArguments():
  # Define the parser and read arguments
  parser = argparse.ArgumentParser(description='Get tags from all files in a directory.')
  parser.add_argument('-d', '--dir', type=str, help='The folder where your experiments are located', required=True)  
  parser.add_argument('-f', '--folders', nargs='+', help='list of folders you want to take into account', required=True)
  args = parser.parse_args()
  return(args)
