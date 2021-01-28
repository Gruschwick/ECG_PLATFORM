# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:30:50 2019

@author: x

AIDMED QRS detector

"""

import pandas as pd
import wfdb
import glob
import os
from os.path import join
from biosppy.signals import ecg
from Helpers.Physionet import LoadFile, PhysionetConstants
from Helpers.Metrices import *
from subprocess import call
from collections import namedtuple

base_dir_raw = join(os.getcwd(),"Datasets/circ1/data_raw")
ret_qrs_scores_list = []
ret_scores_by_class_list = []

ChannelDescription = namedtuple('ChannelDescription',['blocks','annotations'])


def processECGFile(data):
    exec_path = join(os.getcwd(), "Binaries/App/demo" )
    ecg_path = join(os.getcwd(), "Temp/AppTemp/ecg_file.csv" )
    beat_path = join(os.getcwd(), "Temp/AppTemp/beat_file.txt" )
    block_path = join(os.getcwd(), "Temp/AppTemp/block_file.txt" )
    
    #for each channel
    retDict = {}

    #right now only channel [0]
    for channel in [0]:
        
        pd.DataFrame({'ecg': data.data[:,channel]}).to_csv(ecg_path)
           
        prev_wd = os.getcwd()
        
        try:
            os.chdir(os.path.dirname(exec_path))
            execution_result = call([exec_path, str(data.getFreq()), ecg_path, beat_path, block_path])
                
        finally:
            os.chdir(prev_wd)
        
        if  execution_result == 0:
            #read output files
            
            blocksAnn = pd.read_csv(block_path, sep = ";", 
                                names = ["start_T","end_T","start","end","bpm","type"])
            beatsAnnRaw = pd.read_csv(beat_path, sep = "\t",
                               names = ["raw"])
        
            beatsAnn = beatsAnnRaw.raw.str.strip().str.split(" ", expand = True)
            beatsAnn.columns = ["idx","annotation"]                
            beatsAnn = beatsAnn[beatsAnn['annotation'].isnull() == False]
            beatsAnn['idx'] = pd.to_numeric(beatsAnn['idx'], errors = 'coerce')
            retDict[channel] = ChannelDescription(blocksAnn, beatsAnn)            
        else:
            retDict[channel] = None
        
        #except ValueError:
               # print('There was an error when reading file')
              #  print('File skipped')
    return retDict

#for file in PhysionetConstants.DS2:

for file in ["158"]:
    
    try:
        print("current file: {0}".format(file))
        data = LoadFile(join(base_dir_raw, file))
    #rpeaks = ecg.hamilton_segmenter(data.data[:,0], data.freq)[0]
    
        print("freq: {0}".format(data.freq))
    
        ret = processECGFile(data)
    
        annotations_df = ret[list(ret.keys())[0]].annotations
    
        rpeaks = annotations_df.loc[annotations_df.annotation == "R","idx"].values    
    
        ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
        anntype_selected = np.array(data.annotations.anntype)[ann_selector]
        annsamp_selected = data.annotations.annsamp[ann_selector]
    
        ref_peaks = sample_to_time(annsamp_selected, data.freq)
        pred_peaks = sample_to_time(rpeaks, data.freq)
        peaks_matching = match_peaks(ref_peaks, pred_peaks)
            
        print("ref# {0}, pred# {1}, match# {2}".format(len(ref_peaks), len(pred_peaks), len(peaks_matching)))
        
        ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
        print("mean error %.4f, error std %.4f, TPR %.4f, VPP %.4f"%ret_qrs_scores)
    
        ret_qrs_scores_list += [(file,) + ret_qrs_scores]
    
        ret_scores_by_class = qrs_detection_by_class(anntype_selected, peaks_matching)
        print(ret_scores_by_class)
    
        ret_scores_by_class_list += [(file,) + ret_scores_by_class]
    
    except ValueError:
                print('There was an error when reading file')
                print('File skipped')

pd.DataFrame.from_records(ret_qrs_scores_list, columns = ["file", "error", "sd", "TPR", "VPP"])
