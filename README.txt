

Goal of this package is to create a python pipeline to analyze the data output from the GrowthConeAnalyzer MATLAB software: 
"Bagonis, Maria M., et al. "Automated profiling of growth cone heterogeneity defines relations between morphology and motility." J Cell Biol 218.1 (2019): 350-379."


Contents of this package:
analyze_features.py
load_data.py
KnockdownFeatures_class.py

analyze_features.py:
launch as python3 analyze_features.py -h 
from command line for help.


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

load_data.py:
This module is calling KnockdownFeatures_class.py to create an object holding all knockdowns
of an experiment with their respective features.
It can be either run by the analyze_features.py script from the command line, which allows
some automated data analysis options or from an editor to explore and analyze data more
interactively. 
It takes a list of paths and knockdown names as an input.
Folder structure must be the following:
    Experiment1 |-Knockdown1|-Cell1|-GCAFeatureExtraction|-...|-feature1.csv
                                                              |-feature2.csv
                |-Cell2|-GCAFeatureExtraction|-...
                |-Knockdown2|-...

    Experiment2 |-Knockdown2|-...
                |-Knockdown3|-...
GCAFeatureExtraction is hard coded, every feature needs to be contained as a .csv
file somewhere in a folder exactly fitting this name. The path to this file
is allowed to have an undefined number of subfolders within GCAFeatureExtraction.
Not each experiment has to have all the knockdowns. If a knockdown is not found
for an experiment it is simply skipped.
The object, when initialized will only contain the information about path and knockdowns
and contains function to extract and organize the given data into one object.

self.extract_all() will use the input to extract the data from the different features 
by using the KnockdownFeatures_class.py module on each individual knockdown.
It will create the following data structures:
    self.path:       A list of the experiment folders to extract data from (first input argument)
    self.knockdowns: A list of the knockdowns to extract (second input argument)
    self.experiment: A dictionary mapping the experiment+knockdown name to the 
                     respective object
                     created by the KnockdownFeatures_class.py module
    self.features: A list of strings of all the features that could be extracted
    grouped_features: A dictionary mapping the feature names to long format 
                      DataFrames containing the following information combined for
                      all input data:
                     
                      KD: Name of the respective knockdown (second input argument)
                      experiment: name of the experiment (first input argument)
                      item: name of the knockdown+the cell number within that group
                      meltid: row number of the original csv file, can be for example
                              the ID of an individual filopodium
                      timepoint: column number of the original csv file.
                                 usually the timepoint
                      value: cell value of the original csv file. This is the measured
                             value for the feature
                      variable: experiment+item+timepoint. This is the unique ID
                      for the cell and timepoint.
   
Dependencies:
    KnockdownFeatures_class.py

KnockdownFeatures_class.py:
    creates a class that holds all the features of one knockdown
    Usage: Initialize with the path to the experiment and a list of the folder names
    you want to search for features
    either use find_csv() and then load_features(feature) to load an
    individual feature, or use load_all() to load all features.