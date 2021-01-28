# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:30:50 2019

@author: x

hamilton QRS detector

"""

import pandas as pd
import wfdb
import glob
import os
from os.path import join
from biosppy.signals import ecg
from Helpers.Physionet import LoadFile, PhysionetConstants
from Helpers.Metrices import *

base_dir_raw = join(os.getcwd(),"Datasets\\mitdb\\data_raw")



ret_qrs_scores_list = []
ret_scores_by_class_list = []

for file in PhysionetConstants.DS2:
    print("current file: {0}".format(file))
    data = LoadFile(join(base_dir_raw, file))
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
    
    ret_qrs_scores_list += [(file,) + ret_qrs_scores]
    
    ret_scores_by_class = qrs_detection_by_class(anntype_selected, peaks_matching)
    print(ret_scores_by_class)
    
    ret_scores_by_class_list += [(file,) + ret_scores_by_class]

pd.DataFrame.from_records(ret_qrs_scores_list, columns = ["file", "error", "sd", "TPR", "VPP"])
#   file     error        sd       TPR       VPP
#0   100 -0.001585  0.001549  0.999560  1.000000
#1   103 -0.002009  0.001255  0.998560  1.000000
#2   105 -0.000035  0.018473  0.995334  0.986893
#3   111  0.001944  0.012216  0.998588  1.000000
#4   113 -0.001646  0.001365  0.999443  1.000000
#5   117  0.010147  0.008639  1.000000  1.000000
#6   121 -0.001988  0.004858  0.998926  0.999463
#7   123 -0.003229  0.001054  0.998024  1.000000
#8   200 -0.001904  0.006938  0.998462  0.998846
#9   202 -0.002255  0.002189  0.989700  1.000000
#10  210 -0.000320  0.004668  0.969434  0.998834
#11  212 -0.002161  0.001212  0.999636  1.000000
#12  213 -0.000219  0.006216  0.993540  1.000000
#13  214 -0.004218  0.010419  0.996905  1.000000
#14  219 -0.000972  0.002396  0.998143  1.000000
#15  221 -0.002054  0.001808  0.902761  1.000000
#16  222 -0.001060  0.003928  0.997986  0.999597
#17  228 -0.000280  0.015100  0.942523  0.995370
#18  231  0.000370  0.001102  0.999363  1.000000
#19  232 -0.003393  0.002233  1.000000  0.998878
#20  233 -0.002293  0.010560  0.995453  1.000000
#21  234 -0.000806  0.001261  0.998547  1.000000

all_events = sum([x[2] for x in ret_scores_by_class_list], Counter())
hit_events = sum([x[3] for x in ret_scores_by_class_list], Counter())

[(key, hit_events[key]/all_events[key]) for key in all_events.keys()]
#[('N', 0.9957743387114477),
# ('A', 0.9994239631336406),
# ('V', 0.8913043478260869),
# ('Q', 0.5714285714285714),
# ('L', 0.9983034415899176),
# ('a', 0.7),
# ('F', 0.9716494845360825),
# ('E', 1.0),
# ('R', 0.9994246260069045),
# ('j', 1.0),
# ('J', 1.0)]
