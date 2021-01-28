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

ChannelDescription = namedtuple('ChannelDescription',['blocks','annotations'])

segmenters = [ecg.hamilton_segmenter, ecg.engzee_segmenter]

databases = ["cinc1"]*100 +["cinc2"]*100+["mitdb"]*48
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
]

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
    
    data = LoadFile(join(base_dir_raw, patient_nr))

#biosppy ecg segmenters

    for segmenter in segmenters:
        rpeaks = segmenter(data.data[:,0], data.freq)[0]

        #physionet specifics
        ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
        anntype_selected = np.array(data.annotations.anntype)[ann_selector]
        annsamp_selected = data.annotations.annsamp[ann_selector]
        #end of physionet specifics

        ref_peaks = sample_to_time(annsamp_selected, data.freq)
        pred_peaks = sample_to_time(rpeaks, data.freq)
        peaks_matching = match_peaks(ref_peaks, pred_peaks)
        
        ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
        
        ret.append([dbcode] + [patient_nr] + [segmenter.__name__] + [len(ref_peaks)] + 
                   [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
    

#print(from_3dbs("cinc1", "101"))
#print(from_3dbs("mitdb", "100"))
#print(from_3dbs("cinc2", "1003"))

resultsdf = pd.DataFrame(columns=['dbcode', 'patient_nr', 'segmenter', 
                                  'ref_peaks', 'pred_peaks', 'peaks_matching',
                                  'mu', 'sigma', 'TPR', 'PPV'])
#for j in range(len(databases)):
for j in range(1,3):
    print(j, databases[j], patients[j])
    a=check_segmentation(databases[j], patients[j])
    for i in range(len(a)):
        resultsdf.loc[len(resultsdf)] = a[i]

#print(databases[58], patients[58])
