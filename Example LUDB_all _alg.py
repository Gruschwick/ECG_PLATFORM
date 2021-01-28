#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 02:48:56 2020

@author: mateusz
"""
import platform
import os
import wfdb
from wfdb import processing
import glob
import numpy as np
import pandas as pd
from os.path import join
from biosppy.signals import ecg
from Helpers.Physionet import LoadFile, PhysionetConstants
from Helpers.Metrices import *
from subprocess import call
from collections import namedtuple
import traceback
import csv
import h5py
import Helpers.ECGprocessing as fcg
import matplotlib.pyplot as plt


ChannelDescription = namedtuple('ChannelDescription',['blocks','annotations'])

# segmenters = [ecg.hamilton_segmenter, ecg.engzee_segmenter]
segmenters = [ecg.hamilton_segmenter]

databases = ["ludb2"]*200
patients = [str(item) for item in list(range(1,201))]

def processECGFile(data):
    exec_path = join(os.getcwd(), "Binaries/App/demo.exe" )
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
            
    return retDict

def check_segmentation(dbcode, patient_nr):

    ret = []
    
    if(platform.system() == "Linux"):
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
    else:
        base_dir_raw = os.path.join(os.getcwd(),"Datasets\\"+dbcode+"\\data_raw")
    
    if (dbcode=='ludb2'):
        data = LoadFile(join(base_dir_raw, patient_nr), "atr_ii")
    else:
        data = LoadFile(join(base_dir_raw, patient_nr))

# ecg segmenters

    # for segmenter in segmenters:
    #     rpeaks = segmenter(data.data[:,0], data.freq)[0]

    #     ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
    #     anntype_selected = np.array(data.annotations.anntype)[ann_selector]
    #     annsamp_selected = data.annotations.annsamp[ann_selector]

    #     ref_peaks = sample_to_time(annsamp_selected, data.freq)
    #     pred_peaks = sample_to_time(rpeaks, data.freq)
        
    #     try:
    #         peaks_matching = match_peaks(ref_peaks, pred_peaks)
            
    #         ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
    #         ret.append([dbcode] + [patient_nr] + [segmenter.__name__] + [len(ref_peaks)] + 
    #                   [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
        
    #     except Exception:
    #         traceback.print_exc()
    #         print("---")
    #         ret.append([dbcode] + [patient_nr] + [segmenter.__name__] + [-1]*7)
        
# # WTDelin        
    # rpeaks = fcg.delineateMultiLeadECG(data.data[:,0][:,np.newaxis],data.freq)[0][:,7]
    # try: 
    #     ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
    #     anntype_selected = np.array(data.annotations.anntype)[ann_selector]
    #     annsamp_selected = data.annotations.annsamp[ann_selector]
        
    #     ref_peaks = sample_to_time(annsamp_selected, data.freq)
    #     pred_peaks = sample_to_time(rpeaks, data.freq)
        
    #     peaks_matching = match_peaks(ref_peaks, pred_peaks)
            
    #     ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
    #     ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [len(ref_peaks)] + 
    #                     [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
    # except Exception:
    #     traceback.print_exc()
    #     print("---")
    #     ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [-1]*7)
    
    # return ret

#ECGPUWAVE
    try:
        data2 = LoadFile(join(base_dir_raw, patient_nr), "qrs")
        ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
        ann_selector2 = [x in PhysionetConstants.T_beats for x in data2.annotations.anntype]
        anntype_selected = np.array(data.annotations.anntype)[ann_selector]
        annsamp_selected = data.annotations.annsamp[ann_selector]
        annsamp_selected2 = data2.annotations.annsamp[ann_selector2]
        fs = data.freq
        
        ref_peaks = sample_to_time(annsamp_selected, data.freq)
        pred_peaks = sample_to_time(annsamp_selected2, data.freq)
        
        peaks_matching = match_peaks(ref_peaks, pred_peaks)
            
        ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
        # ret.append([dbcode] + [patient_nr] + [fs] + ["WTDelin"] +[anntype_selected] + [ppeaks] + [len(ref_peaks)] + 
        #                 [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
        
        ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [len(ref_peaks)] + 
              [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
    except Exception:
        traceback.print_exc()
        print("---")
        ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [-1]*7)
    
    return ret

    # rpeaks = processing.XQRS(data.data[:,0],data.freq)[0]
    # rpeaks.detect()
    # rpeaks = processing.xqrs_detect(sig=data.data[:,0], fs=data.freq)
    # # rpeaks.detect()
    # # XQRS = processing.XQRS(sig, fs)
    # # rpeaks.detect()
    # try: 
    #     ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
    #     anntype_selected = np.array(data.annotations.anntype)[ann_selector]
    #     annsamp_selected = data.annotations.annsamp[ann_selector]
        
    #     ref_peaks = sample_to_time(annsamp_selected, data.freq)
    #     pred_peaks = sample_to_time(rpeaks, data.freq)
        
    #     peaks_matching = match_peaks(ref_peaks, pred_peaks)
            
    #     ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
    #     ret.append([dbcode] + [patient_nr] + ["XQRS"] + [len(ref_peaks)] + 
    #                    [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
    # except Exception:
    #     traceback.print_exc()
    #     print("---")
    #     ret.append([dbcode] + [patient_nr] + ["XQRS"] + [-1]*7)
    
    # return ret

#print(from_3dbs("cinc1", "101"))
#print(from_3dbs("mitdb", "100"))
#print(from_3dbs("cinc2", "1003"))

resultsdf = pd.DataFrame(columns=['dbcode', 'patient_nr', 'segmenter', 
                                  'ref_peaks', 'pred_peaks', 'peaks_matching',
                                  'mu', 'sigma', 'TPR', 'PPV'])

#for j in range(len(databases)):
#for j in range(330,333):
#for j in [10, 58, 74, 109, 122, 142, 154, 190, 196, 197, 250]:
#for j in ([190, 196, 197, 249, 250] + list(range(251,248+82))):
#for j in range(248, 248+82):
#for j in range(248+82+34, 248+82+200):
#for j in range(248+82+200):
for j in range(len(databases)):
    print(j, databases[j], patients[j])
    a=check_segmentation(databases[j], patients[j])
    for i in range(len(a)):
        resultsdf.loc[len(resultsdf)] = a[i]

#print(databases[58], patients[58])
