import platform
import os
import wfdb
import glob
import numpy as np
from os.path import join
from biosppy.signals import ecg
from Helpers.Physionet import LoadFile, PhysionetConstants
from Helpers.Metrices import *

def example_calc_qrs_detection_quality():
    
    if(platform.system() == "Linux"):
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/mitdb/data_raw")
    else:
        base_dir_raw = os.path.join(os.getcwd(),"Datasets\\mitdb\\data_raw")
    
    data = LoadFile(join(base_dir_raw, "100"))

    rpeaks = ecg.hamilton_segmenter(data.data[:,0], data.freq)[0]

    ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
    anntype_selected = np.array(data.annotations.anntype)[ann_selector]
    annsamp_selected = data.annotations.annsamp[ann_selector]

    ref_peaks = sample_to_time(annsamp_selected, data.freq)
    pred_peaks = sample_to_time(rpeaks, data.freq)
    peaks_matching = match_peaks(ref_peaks, pred_peaks)
        
    print("ref# {0}, pred# {1}, match# {2}".format(len(ref_peaks), len(pred_peaks), len(peaks_matching)))
    
    ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)    
    print("mean error %.4f, error std %.4f, TPR %.4f, VPP %.4f"%ret_qrs_scores)
    
    ret_scores_by_class = qrs_detection_by_class(anntype_selected, peaks_matching)
    print(ret_scores_by_class)
    
    data.plot(1000,3000)

example_calc_qrs_detection_quality()

