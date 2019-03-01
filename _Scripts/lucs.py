# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 19:39:43 2018

@author: Tvon-
"""

import numpy as np
import pandas as pd

def read_luc(filename):
    raw_data = pd.read_csv(filename,header=None)
    #read out the three tables corresponding to factors, FLuc data and RLuc data.
    #Flatten the two-dimensonal data into linear vectors.
    Factors = raw_data.iloc[2:10,1:13].values.flatten()
    FLuc = pd.Series(raw_data.iloc[13:21,1:13].values.flatten())
    RLuc = pd.Series(raw_data.iloc[24:32,1:13].values.flatten())
    #check if the Factors table contains more than one factor level
    #if so, split into two factor columns at ':'
    if ':' in Factors[0]:
        luc_df = pd.DataFrame([entries.split(':') for entries in Factors])
        luc_df.columns = ['Factor 1','Factor 2']
        luc_df['FLuc'] = FLuc
        luc_df['RLuc'] = RLuc
    else:
        Factors = pd.Series(Factors)
        luc_df = pd.concat([Factors,FLuc,RLuc], axis=1)
        luc_df.columns = ['Factor','FLuc','RLuc']

    #replace the FLuc and RLuc values with the F/R ratio
    luc_df['FLuc'] = pd.to_numeric(luc_df['FLuc'])
    luc_df['RLuc'] = pd.to_numeric(luc_df['RLuc'])
    luc_df['FR_ratio'] = luc_df.FLuc / luc_df.RLuc
    del luc_df['RLuc'], luc_df['FLuc']

    #remove rows where factors indicate failed runs
    remove_factors = ['cont','high','other']
    for word in remove_factors:
        luc_df = luc_df[luc_df['Factor 1'] != word]

    return luc_df

#==============================================================================

def describe_luc(filename):
    luc_df = read_luc(filename)

    print('File ' + filename + ' has ' + str(luc_df.shape[1]-2) + ' factor levels.')
    for factor in range(luc_df.shape[1]-2):
        levels = luc_df.iloc[:,factor].value_counts().index
        print('Factor ' + str(factor + 1) + ' has ' + str(len(levels)) + ' values: ' + str(list(levels)))

#==============================================================================

def identify_controls(df, column_index):

    #compare against a list of routinely used standard controls
    standard_controls = ['C','RLucC','control','Control']

    #extract the factor levels from the dataframe
    if type(column_index) == int:
        control_candidates = list(df.iloc[:,column_index].value_counts().index)
    else:
        control_candidates = list(df.loc[:,column_index].value_counts().index)

    #determine which factor levels are in the standard_controls list
    comparison_index = [level in standard_controls for level in control_candidates]

    #extract levels corresponding to standard control list entries from control_candidates
    identified_controls = []
    for counter in range(len(control_candidates)):
        if comparison_index[counter]:
                identified_controls.append(control_candidates[counter])

    #if there is exactly one level corresponding to a standard_control,
    #use this as the returned control
    if len(identified_controls)  == 1:
        return identified_controls[0]
    #otherwise ask for the control level
    else:
        user_input = input('No unambiguous control level could be identified for column ' + str(column_index) + '\nPlease enter the control level for this column.')
        #check that user_input exists among the factor levels
        if user_input in control_candidates:
            return user_input
        else:
            raise Exception('Entered value does not exist in this column.')

#==============================================================================

def luc_double_normalise(df, df_controls):

    #ensure that two controls have been defined
    if len(df_controls) != 2: raise Exception('luc_double_normalise requires exactly two controls to be defined for normalisation')

    #normalise all values to factor 1 control
    #slice the data so only data for one plasmid are retained
    for plasmid in ['RLucI','RLucII','RLucC']:
        #extract all entries for plasmid
        this_slice = df[df['Factor 2']==plasmid]
        this_control_set = this_slice[this_slice['Factor 1']==df_controls[0]]
        norm_factor = np.mean(this_control_set['FR_ratio'])
        this_slice['norm1'] = [value/norm_factor for value in list(this_slice['FR_ratio'])]
        #create df_norm1 from this_slice if it does not exist, or else append
        #this_slice to existing df_norm1
        if 'df_norm1' in locals():
            df_norm1 = pd.concat([df_norm1, this_slice])
        else:
            df_norm1 = this_slice

    #now normalise all values to factor 2 control
    for condition in list(df_norm1['Factor 1'].value_counts().index):
        this_slice = df_norm1[df_norm1['Factor 1']==condition]
        this_control_set = this_slice[this_slice['Factor 2']==df_controls[1]]
        norm_factor = np.mean(this_control_set['norm1'])
        this_slice['norm2'] = [value/norm_factor for value in list(this_slice['norm1'])]
        #create df_norm2 from this_slice if it does not exist, or else append
        #this_slice to existing df_norm2
        if 'df_norm2' in locals():
            df_norm2 = pd.concat([df_norm2, this_slice])
        else:
            df_norm2 = this_slice

    return df_norm2

#==============================================================================

def prep_luc(*filenames, controls):

    #ensure that controls is a list
    controls = list(controls)
    #for each of the filenames passed to prep_luc:
    for name in filenames:

        #read in the file and split factors if required via read_luc()
        new_read = read_luc(name)

        #double normalise the FR ratios in new_read
        new_norm = luc_double_normalise(new_read, controls)

        #if luc_df already exists, append new_read to the existing variable
        if 'luc_df' in locals():
            try:
                luc_df = pd.concat([luc_df, new_norm])
            except:
                print('Could not concatenate multiple arrays.')
        #otherwise copy new_read into read_luc
        else:
            luc_df = new_norm

        #remove RLucC data which are only for normalisation
        luc_df = luc_df[luc_df['Factor 2'] != controls[1]]

    return luc_df
