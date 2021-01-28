import platform
import os
import wfdb
import glob
import numpy as np
from os.path import join
from biosppy.signals import ecg
from Helpers.Physionet import LoadFile, PhysionetConstants
from Helpers.Metrices import *

segmenters = [ecg.hamilton_segmenter, ecg.engzee_segmenter]

def from_3dbs(dbcode, patient_nr):

    ret = []
    
    if(platform.system() == "Linux"):
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
    else:
        base_dir_raw = os.path.join(os.getcwd(),"Datasets\\"+dbcode+"\\data_raw")
    
    data = LoadFile(join(base_dir_raw, patient_nr))

    for segmenter in segmenters:
        rpeaks = segmenter(data.data[:,0], data.freq)[0]

        ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
        anntype_selected = np.array(data.annotations.anntype)[ann_selector]
        annsamp_selected = data.annotations.annsamp[ann_selector]

        ref_peaks = sample_to_time(annsamp_selected, data.freq)
        pred_peaks = sample_to_time(rpeaks, data.freq)
        peaks_matching = match_peaks(ref_peaks, pred_peaks)
        
        ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
        ret.append([dbcode] + [patient_nr] + [segmenter.__name__] + [len(ref_peaks)] + 
                   [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))

    return ret

print(from_3dbs("cinc1", "101"))
print(from_3dbs("mitdb", "100"))
print(from_3dbs("cinc2", "1022"))