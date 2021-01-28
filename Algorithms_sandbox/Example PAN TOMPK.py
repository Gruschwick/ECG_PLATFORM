#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 01:59:56 2020

@author: mateusz
"""
import platform
import os
import wfdb
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
import scipy.signal as signal
import pathlib


ChannelDescription = namedtuple('ChannelDescription',['blocks','annotations'])

segmenters = [ecg.hamilton_segmenter, ecg.engzee_segmenter]

databases = ["cinc1"]*100 +["cinc2"]*100+["mitdb"]*48+["qtdb"]*82+["ludb"]*200
patients = [str(i) for i in list(range(100,200))] + [
        '1003', '1009', '1016', '1019', '1020', '1022', '1023', '1028', '1032',
        '1033', '1036', '1043', '1069', '1071', '1073', '1077', '1169', '1195',
        '1242', '1284', '1354', '1376', '1388', '1447', '1456', '1485', '1503',
        '1522', '1565', '1584', '1683', '1686', '1715', '1742', '1774', '1804',
        '1807', '1821', '1858', '1866', '1900', '1906', '1954', '1993', '1998',
        '2041', '2063', '2132', '2164', '2174', '2201', '2203', '2209', '2247',
        '2277', '2279', '2283', '2296', '2327', '2370', '2384', '2397', '2469',
        '2527', '2552', '2556', '2602', '2639', '2664', '2714', '2728', '2732',
        '2733', '2798', '2800', '2812', '2839', '2850', '2879', '2885', '2886',
        '2907', '2923', '2970', '3188', '3266', '41024', '41025', '41081', 
        '41164', '41173', '41180', '41566', '41778', '41951', '42228', '42511',
        '42878', '42961', '43247' 
] + ['100', '101', '102', '103', '104', '105', '106', '107', '108', '109',
     '111', '112', '113', '114', '115', '116', '117', '118', '119', '121',
     '122', '123', '124', '200', '201', '202', '203', '205', '207', '208',
     '209', '210', '212', '213', '214', '215', '217', '219', '220', '221',
     '222', '223', '228', '230', '231', '232', '233', '234' 
] + ['sel14046', 'sel14157', 'sel14172', 'sel15814', 'sel16265', 'sel16272', 
     'sel16273', 'sel16420', 'sel16483', 'sel16539', 'sel16773', 'sel16786', 
     'sel16795', 'sel17152', 'sel17453', 'sele0104', 'sele0106', 'sele0107', 
     'sele0110', 'sele0111', 'sele0112', 'sele0114', 'sele0116', 'sele0121', 
     'sele0122', 'sele0124', 'sele0126', 'sele0129', 'sele0133', 'sele0136', 
     'sele0166', 'sele0170', 'sele0203', 'sele0210', 'sele0211', 'sele0303', 
     'sele0405', 'sele0406', 'sele0409', 'sele0411', 'sele0509', 'sele0603', 
     'sele0604', 'sele0606', 'sele0607', 'sele0609', 'sele0612', 'sele0704', 
     'sel100', 'sel102', 'sel103', 'sel104', 'sel114', 'sel116', 'sel117', 
     'sel123', 'sel213', 'sel221', 'sel223', 'sel230', 'sel231', 'sel232', 
     'sel233', 'sel301', 'sel302', 'sel306', 'sel307', 'sel308', 'sel310', 
     'sel803', 'sel808', 'sel811', 'sel820', 'sel821', 'sel840', 'sel847', 
     'sel853', 'sel871', 'sel872', 'sel873', 'sel883', 'sel891'
] + [str(item) for item in list(range(1,201))]

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


def MWA_from_name(function_name):
    if function_name == "cumulative":
        return MWA_cumulative
    elif function_name == "convolve":
        return MWA_convolve
    elif function_name == "original":
        return MWA_original
    else: 
        raise RuntimeError('invalid moving average function!')

#Fast implementation of moving window average with numpy's cumsum function 
def MWA_cumulative(input_array, window_size):
    
    ret = np.cumsum(input_array, dtype=float)
    ret[window_size:] = ret[window_size:] - ret[:-window_size]
    
    for i in range(1,window_size):
        ret[i-1] = ret[i-1] / i
    ret[window_size - 1:]  = ret[window_size - 1:] / window_size
    
    return ret

#Original Function 
def MWA_original(input_array, window_size):

    mwa = np.zeros(len(input_array))
    mwa[0] = input_array[0]
    
    for i in range(2,len(input_array)+1):
        if i < window_size:
            section = input_array[0:i]
        else:
            section = input_array[i-window_size:i]        
        
        mwa[i-1] = np.mean(section)

    return mwa

#Fast moving window average implemented with 1D convolution 
def MWA_convolve(input_array, window_size):
    
    ret = np.pad(input_array, (window_size-1,0), 'constant', constant_values=(0,0))
    ret = np.convolve(ret,np.ones(window_size),'valid')
    
    for i in range(1,window_size):
        ret[i-1] = ret[i-1] / i
    ret[window_size-1:] = ret[window_size-1:] / window_size
    
    return ret


def normalise(input_array):

    output_array = (input_array-np.min(input_array))/(np.max(input_array)-np.min(input_array))

    return output_array


def panPeakDetect(detection, fs):    

    min_distance = int(0.25*fs)

    signal_peaks = [0]
    noise_peaks = []

    SPKI = 0.0
    NPKI = 0.0

    threshold_I1 = 0.0
    threshold_I2 = 0.0

    RR_missed = 0
    index = 0
    indexes = []

    missed_peaks = []
    peaks = []

    for i in range(len(detection)):

        if i>0 and i<len(detection)-1:
            if detection[i-1]<detection[i] and detection[i+1]<detection[i]:
                peak = i
                peaks.append(i)

                if detection[peak]>threshold_I1 and (peak-signal_peaks[-1])>0.3*fs:
                        
                    signal_peaks.append(peak)
                    indexes.append(index)
                    SPKI = 0.125*detection[signal_peaks[-1]] + 0.875*SPKI
                    if RR_missed!=0:
                        if signal_peaks[-1]-signal_peaks[-2]>RR_missed:
                            missed_section_peaks = peaks[indexes[-2]+1:indexes[-1]]
                            missed_section_peaks2 = []
                            for missed_peak in missed_section_peaks:
                                if missed_peak-signal_peaks[-2]>min_distance and signal_peaks[-1]-missed_peak>min_distance and detection[missed_peak]>threshold_I2:
                                    missed_section_peaks2.append(missed_peak)

                            if len(missed_section_peaks2)>0:           
                                missed_peak = missed_section_peaks2[np.argmax(detection[missed_section_peaks2])]
                                missed_peaks.append(missed_peak)
                                signal_peaks.append(signal_peaks[-1])
                                signal_peaks[-2] = missed_peak   

                else:
                    noise_peaks.append(peak)
                    NPKI = 0.125*detection[noise_peaks[-1]] + 0.875*NPKI

                threshold_I1 = NPKI + 0.25*(SPKI-NPKI)
                threshold_I2 = 0.5*threshold_I1

                if len(signal_peaks)>8:
                    RR = np.diff(signal_peaks[-9:])
                    RR_ave = int(np.mean(RR))
                    RR_missed = int(1.66*RR_ave)

                index = index+1      
    
    signal_peaks.pop(0)

    return signal_peaks

def pan_tompkins_detector(unfiltered_ecg, freq, MWA_name='cumulative'):
        """
        Jiapu Pan and Willis J. Tompkins.
        A Real-Time QRS Detection Algorithm. 
        In: IEEE Transactions on Biomedical Engineering 
        BME-32.3 (1985), pp. 230–236.
        """
        fs = freq
        f1 = 5/fs
        f2 = 15/fs

        b, a = signal.butter(1, [f1*2, f2*2], btype='bandpass')

        filtered_ecg = signal.lfilter(b, a, unfiltered_ecg)        

        diff = np.diff(filtered_ecg) 

        squared = diff*diff

        N = int(0.12*fs)
        mwa = MWA_from_name(MWA_name)(squared, N)
        mwa[:int(0.2*fs)] = 0

        mwa_peaks = panPeakDetect(mwa, fs)
    
        return mwa_peaks

def check_segmentation(dbcode, patient_nr):

    ret = []
    
    if(platform.system() == "Linux"):
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
    else:
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
    
    if (dbcode=='ludb'):
        data = LoadFile(join(base_dir_raw, patient_nr), "atr_ii")
    else:
        data = LoadFile(join(base_dir_raw, patient_nr))

#ecg segmenters
#
#    for segmenter in segmenters:
#        rpeaks = segmenter(data.data[:,0], data.freq)[0]
#
#        ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
#        anntype_selected = np.array(data.annotations.anntype)[ann_selector]
#        annsamp_selected = data.annotations.annsamp[ann_selector]
#
#        ref_peaks = sample_to_time(annsamp_selected, data.freq)
#        pred_peaks = sample_to_time(rpeaks, data.freq)
#        
#        try:
#            peaks_matching = match_peaks(ref_peaks, pred_peaks)
#            
#            ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
#            ret.append([dbcode] + [patient_nr] + [segmenter.__name__] + [len(ref_peaks)] + 
#                       [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
#        
#        except Exception:
#            traceback.print_exc()
#            print("---")
#            ret.append([dbcode] + [patient_nr] + [segmenter.__name__] + [-1]*7)
#        
        
#WTDelin        
#    rpeaks = fcg.delineateMultiLeadECG(data.data[:,0][:,np.newaxis],data.freq)[0][:,4]
#    try: 
#        ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
#        anntype_selected = np.array(data.annotations.anntype)[ann_selector]
#        annsamp_selected = data.annotations.annsamp[ann_selector]
#        
#        ref_peaks = sample_to_time(annsamp_selected, data.freq)
#        pred_peaks = sample_to_time(rpeaks, data.freq)
        
#        peaks_matching = match_peaks(ref_peaks, pred_peaks)
            
#        ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
#        ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [len(ref_peaks)] + 
#                       [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
#    except Exception:
#        traceback.print_exc()
#        print("---")
#        ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [-1]*7)
    
#    return ret

    rpeaks = pan_tompkins_detector(data.data[:,0],data.freq)[0]
    # print (rpeaks)
    try: 
        ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
        anntype_selected = np.array(data.annotations.anntype)[ann_selector]
        annsamp_selected = data.annotations.annsamp[ann_selector]
        
        ref_peaks = sample_to_time(annsamp_selected, data.freq)
        pred_peaks = sample_to_time(rpeaks, data.freq)
        
        peaks_matching = match_peaks(ref_peaks, pred_peaks)
            
        ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
        ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [len(ref_peaks)] + 
                        [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
    except Exception:
        traceback.print_exc()
        print("---")
        ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [-1]*7)
    
    return ret

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
