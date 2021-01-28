import matplotlib.pyplot as plt
import h5py
import wfdb
import numpy as np
import csv
import platform
import os
import glob
import numpy as np
import pandas as pd
from os.path import join
from biosppy.signals import ecg
from Helpers.Physionet import LoadFile, PhysionetConstants
from Helpers.Metrices import *
import Helpers.ECGprocessing as fcg
import traceback

dbcode= "mitdb"; patient_nr = '101'

ret = []

if(platform.system() == "Linux"):
    base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
else:
    base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")


data = LoadFile(join(base_dir_raw, patient_nr))
sig = data.data
freq = data.getFreq()
ECGdelin = fcg.delineateMultiLeadECG(sig[:,0][:,np.newaxis],freq)

#rpeaks = [y[:,5] for y in ECGdelin]
rpeaks = ECGdelin[0][:,4]
       
ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
anntype_selected = np.array(data.annotations.anntype)[ann_selector]
annsamp_selected = data.annotations.annsamp[ann_selector]

ref_peaks = sample_to_time(annsamp_selected, data.freq)
pred_peaks = sample_to_time(rpeaks, data.freq)

peaks_matching = match_peaks(ref_peaks, pred_peaks)
    
ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
ret.append([dbcode] + [patient_nr] + ["WTDelin"] + [len(ref_peaks)] + 
               [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))

plt.plot(np.arange(1, 3000, 1)/data.getFreq(), data.data[1:3000,0])
plt.show()
plt.clf()
tmp3=500
tmp4=1399
plt.scatter(pred_peaks[tmp3:tmp4], ref_peaks[tmp3:tmp4])

print ("CD")
