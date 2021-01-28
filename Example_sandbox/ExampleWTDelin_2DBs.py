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
import Helpers.ECGprocessing as ecg
import traceback

#ChannelDescription = namedtuple('ChannelDescription',['blocks','annotations'])

tests = ['BR', 'BC1', 'BC2', 'BI1', 'BI2', 'BI3', 'BI4', 'BI5', 'PC1', 'PC2', 'PR1', 'PR2']
leadNames = ['V1','V2','V3','V4','V5','V6','DI','DII','DIII','aVR','aVL','aVF']
h5attr = [['Annotation names'],['Pon,P,Pend, QRSon,R,QRSend,Ton,Tp,Tend']]

databases = ["mitdb"]*48+["qtdb"]*82
patients = ['100', '101', '102', '103', '104', '105', '106', '107', '108', '109',
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
]

#dbcode= "mitdb"; patient_nr = '100'
def check_segmentation(dbcode, patient_nr):

    ret = []
    
    if(platform.system() == "Linux"):
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
    else:
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
    
  
    data = LoadFile(join(base_dir_raw, patient_nr))
    
    
    try:
        
        s = data.data
        freq = data.getFreq()
        #s, att = wfdb.rdsamp(rec.zfill(4),pb_dir=dbase) # Read from Physionet
                
        # Read the annotations
       # a1,a2,_ = annot[patient][measurement+1].split(';')
      #  a1 = int(int(a1)*att['fs']) # Start of balloon inflation
      #  a2 = int(int(a2)*att['fs']) # End of balloon inflation
                
        # Only keep the ischemic part of the signal
      #  s = s[a1:a1+a2,:]
        # Calculate augmented limb leads and append them to the signals
      #  aVR, aVL, aVF = ecg.augmentedLimbs(s[:,-3], s[:,-2])
     #   s = np.concatenate((s, aVR, aVL, aVF), axis=1) 
                
        # Delineate all the ECG leads using the WT and fusion techniques
        ECGdelin = ecg.delineateMultiLeadECG(s,freq)
      #  ECGdelin = ecg.delineateMultiLeadECG(patient_nr,)
                
        # Create subgroup for the patient and save all the leads
       # grp = outFile.create_group(tests[measurement] + '/' + str(patient_nr).zfill(3))
      #  for idx, ECG in enumerate(ECGdelin):
           # dsetName =  leadNames[idx]
         #   dset = grp.create_dataset(dsetName,ECG.shape,ECG.dtype)
          #  dset[...] = ECG
            #dset.attrs.create(h5attr[0],h5attr[1])                   
                    
            #delinDB[measurement] += [ECGdelin]
    except ValueError:
        print('There was an error when reading file', patient_nr+',', 'file skipped')
          #delinDB[measurement] += ['Error reading record ' + rec]
                
#        else: # If the record doesn't exist, add an empty list
#            delinDB[measurement] += [[]]
                
def segmenter(sig, freq):
    
    #sig = s[:,0]
    #sig[:,np.newaxis].shape
    ECGdelin = ecg.delineateMultiLeadECG(sig[:,np.newaxis],freq)
    return [x[:,5] for x in ECGdelin]
                                      
#aidmed mywave
    ret2 = segmenter(sig, freq)
    try: 
        annotations_df = ret2[list(ret2.keys())[0]].annotations
        rpeaks = annotations_df.loc[annotations_df.annotation == "R","idx"].values  
        #rpeaks = segmenter(data.data[:,5], data.freq)[0]

        ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
        anntype_selected = np.array(data.annotations.anntype)[ann_selector]
        annsamp_selected = data.annotations.annsamp[ann_selector]
        ref_peaks = sample_to_time(annsamp_selected, data.freq)
        pred_peaks = sample_to_time(rpeaks, data.freq)
        peaks_matching = match_peaks(ref_peaks, pred_peaks)
        ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks, peaks_matching)
        ret.append([dbcode] + [patient_nr] + ["WTdelineator"] + [len(ref_peaks)] + 
              [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores))
    except Exception:
        traceback.print_exc()
        print("---")
        ret.append([dbcode] + [patient_nr] + ["WTdelineator"] + [-1]*7)
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
for j in range(48+82):
    print(j, databases[j], patients[j])
    a=check_segmentation(databases[j], patients[j])
    for i in range(len(a)):
        resultsdf.loc[len(resultsdf)] = a[i]

#print(databases[58], patients[58])        
