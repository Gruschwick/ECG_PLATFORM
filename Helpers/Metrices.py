# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:56:51 2019

@author: x
"""

import numpy as np
from collections import Counter

class MetricesConstants(object):
    #qrs_cutoff_distance = 0.2
    qrs_cutoff_distance = 0.120 #https://www.sciencedirect.com/science/article/abs/pii/S1746809417300216

def sample_to_time(samples, freq):
    return samples/freq

def match_peaks( ref_peaks, pred_peaks, cutoff_distance = None):
    '''
    calc best matching between ref_peaks and pred_peaks with cutoff (error time distance no longer than cutoff_distance)
    [(ref_peaks[r], pred_peaks[c]) for r, c in zip(row_ind, col_ind)
    '''
    from scipy.optimize import linear_sum_assignment
    assert np.all(ref_peaks >= 0), "positive time"
    assert np.all(pred_peaks >= 0), "positive time"
    
    if cutoff_distance is None:
        cutoff_distance = MetricesConstants.qrs_cutoff_distance
    
    max_ref_peaks = np.max(ref_peaks)
    len_ref_peaks = len(ref_peaks)
    max_pred_peaks = np.max(pred_peaks)
    len_pred_peaks = len(pred_peaks)
    
    max_len = max(len_ref_peaks, len_pred_peaks)
    max_peaks = max(max_ref_peaks, max_pred_peaks)
    max_distance = max_peaks*10000    
    
    ref_peaks = np.pad(ref_peaks, ((0,max_len - len_ref_peaks),), 'constant', constant_values=(0, max_distance))        
    pred_peaks = np.pad(pred_peaks, ((0,max_len - len_pred_peaks),), 'constant', constant_values=(0, max_distance))        

    distance_matrix = np.abs(ref_peaks[:,np.newaxis] - pred_peaks[np.newaxis,:])
    
    distance_matrix[distance_matrix > cutoff_distance] = max_distance
    
    row_ind, col_ind= linear_sum_assignment(distance_matrix)
    
    matching_filtered = [(r,c) for r, c in zip(row_ind, col_ind) if distance_matrix[r,c] <= cutoff_distance]
    
    #ref_peaks[r], pred_peaks[c]
    return matching_filtered

def qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching):
    deltas = [(ref_peaks[r] - pred_peaks[c]) for r, c in peaks_matching]
    tpr = len(peaks_matching)/len(ref_peaks)
    ppv = len(peaks_matching)/len(pred_peaks)
    
    return np.mean(deltas), np.std(deltas), tpr, ppv

def qrs_detection_by_class(ref_peaks_class, peaks_matching):
    ref_counts = Counter(ref_peaks_class)
    detected_counts = Counter(ref_peaks_class[r] for r, c in peaks_matching)
    
    return {(k, detected_counts.get(k,0)/ref_counts[k]) for k in ref_counts.keys()}, ref_counts, detected_counts