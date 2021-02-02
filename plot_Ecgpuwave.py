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
from Helpers import ecgpu_detector
from subprocess import call
from collections import namedtuple
import traceback
import csv
import h5py
import Helpers.ECGprocessing as fcg
import matplotlib.pyplot as plt
from IPython.display import display
import shutil
import posixpath

#ChannelDescription = namedtuple('ChannelDescription',['blocks','annotations'])

#segmenters = [ecg.hamilton_segmenter, ecg.engzee_segmenter]

# databases = ["cinc1"]*100 +["cinc2"]*100+["mitdb"]*48++["ludb"]*200+
# patients = [str(i) for i in list(range(100,200))] + [
#         '1003', '1009', '1016', '1019', '1020', '1022', '1023', '1028', '1032',
#         '1033', '1036', '1043', '1069', '1071', '1073', '1077', '1169', '1195',
#         '1242', '1284', '1354', '1376', '1388', '1447', '1456', '1485', '1503',
#         '1522', '1565', '1584', '1683', '1686', '1715', '1742', '1774', '1804',
#         '1807', '1821', '1858', '1866', '1900', '1906', '1954', '1993', '1998',
#         '2041', '2063', '2132', '2164', '2174', '2201', '2203', '2209', '2247',
#         '2277', '2279', '2283', '2296', '2327', '2370', '2384', '2397', '2469',
#         '2527', '2552', '2556', '2602', '2639', '2664', '2714', '2728', '2732',
#         '2733', '2798', '2800', '2812', '2839', '2850', '2879', '2885', '2886',
#         '2907', '2923', '2970', '3188', '3266', '41024', '41025', '41081', 
#         '41164', '41173', '41180', '41566', '41778', '41951', '42228', '42511',
#         '42878', '42961', '43247' 
# ] + ['100', '101', '102', '103', '104', '105', '106', '107', '108', '109',
#      '111', '112', '113', '114', '115', '116', '117', '118', '119', '121',
#      '122', '123', '124', '200', '201', '202', '203', '205', '207', '208',
#      '209', '210', '212', '213', '214', '215', '217', '219', '220', '221',
#      '222', '223', '228', '230', '231', '232', '233', '234' 
# ] + ['sel14046', 'sel14157', 'sel14172', 'sel15814', 'sel16265', 'sel16272', 
#       'sel16273', 'sel16420', 'sel16483', 'sel16539', 'sel16773', 'sel16786', 
#       'sel16795', 'sel17152', 'sel17453', 'sele0104', 'sele0106', 'sele0107', 
#       'sele0110', 'sele0111', 'sele0112', 'sele0114', 'sele0116', 'sele0121', 
#       'sele0122', 'sele0124', 'sele0126', 'sele0129', 'sele0133', 'sele0136', 
#       'sele0166', 'sele0170', 'sele0203', 'sele0210', 'sele0211', 'sele0303', 
#       'sele0405', 'sele0406', 'sele0409', 'sele0411', 'sele0509', 'sele0603', 
#       'sele0604', 'sele0606', 'sele0607', 'sele0609', 'sele0612', 'sele0704', 
#       'sel100', 'sel102', 'sel103', 'sel104', 'sel114', 'sel116', 'sel117', 
#       'sel123', 'sel213', 'sel221', 'sel223', 'sel230', 'sel231', 'sel232', 
#       'sel233', 'sel301', 'sel302', 'sel306', 'sel307', 'sel308', 'sel310', 
#       'sel803', 'sel808', 'sel811', 'sel820', 'sel821', 'sel840', 'sel847', 
#       'sel853', 'sel871', 'sel872', 'sel873', 'sel883', 'sel891'
# ] + 

# patients = [str(item) for item in list(range(1,201))]
     
# #databases = ["qtdb2"]*105
# databases = ["mitdbp"]*12
# patients = ['100', '101', '103', '106','117', '119', 
#       '122', '207', '214',
#       '222', '223', '231'
#       ]

#databases = ["qtdb"]*82
# patients = ['sel14046', 'sel14157', 'sel14172', 'sel15814', 'sel16265', 'sel16272', 
#       'sel16273', 'sel16420', 'sel16483', 'sel16539', 'sel16773', 'sel16786', 
#       'sel16795', 'sel17152', 'sel17453', 'sele0104', 'sele0106', 'sele0107', 
#       'sele0110', 'sele0111', 'sele0112', 'sele0114', 'sele0116', 'sele0121', 
#       'sele0122', 'sele0124', 'sele0126', 'sele0129', 'sele0133', 'sele0136', 
#       'sele0166', 'sele0170', 'sele0203', 'sele0210', 'sele0211', 'sele0303', 
#       'sele0405', 'sele0406', 'sele0409', 'sele0411', 'sele0509', 'sele0603', 
#       'sele0604', 'sele0606', 'sele0607', 'sele0609', 'sele0612', 'sele0704', 
#       'sel100', 'sel102', 'sel103', 'sel104', 'sel114', 'sel116', 'sel117', 
#       'sel123', 'sel213', 'sel221', 'sel223', 'sel230', 'sel231', 'sel232', 
#       'sel233', 'sel301', 'sel302', 'sel306', 'sel307', 'sel308', 'sel310', 
#       'sel803', 'sel808', 'sel811', 'sel820', 'sel821', 'sel840', 'sel847', 
#       'sel853', 'sel871', 'sel872', 'sel873', 'sel883', 'sel891'
# ]

#databases = ["qtdb2"]*105
# #patients = ['sel100','sel102','sel103','sel104','sel114','sel116','sel117','sel123',
#             'sel14046','sel14157','sel14172','sel15814','sel16265','sel16272','sel16273',
#             'sel16420','sel16483','sel16539','sel16773','sel16786','sel16795','sel17152',
#             'sel17453','sel213','sel221','sel223','sel230','sel231','sel232','sel233',
#             'sel301','sel302','sel306','sel307','sel308','sel31','sel310','sel32',
#             'sel33','sel34','sel35','sel36','sel37','sel38','sel39','sel40','sel41',
#             'sel42','sel43','sel44','sel45','sel46','sel47','sel48','sel49','sel50','sel51',
#             'sel52','sel803','sel808','sel811','sel820','sel821','sel840','sel847','sel853',
#             'sel871','sel872','sel873','sel883','sel891','sele0104','sele0106','sele0107',
#             'sele0110','sele0111','sele0112','sele0114','sele0116','sele0121','sele0122',
#             'sele0124','sele0126','sele0129','sele0133','sele0136','sele0166','sele0170',
#             'sele0203','sele0210','sele0211','sele0303','sele0405','sele0406','sele0409',
#             'sele0411','sele0509','sele0603','sele0604','sele0606','sele0607','sele0609',
#             'sele0612','sele0704']


# databases = ["ludb"] * 200
# patients = [str(item) for item in list(range(1,201))]

# databases = ["but-pdb"] * 50

# patients = [str(i) for i in list(range(1,51))]

if(platform.system() == "Linux"):
    base_dir_raw = os.path.join(os.getcwd(),"Datasets/mitdbp/data_raw/214")
else:
    base_dir_raw = os.path.join(os.getcwd(),"Datasets\\"+dbcode+"\\data_raw")
    
    if (dbcode=='mitdbp'):
        data = LoadFile(join(base_dir_raw, patient_nr), "pwave")
    else:
        data = LoadFile(join(base_dir_raw, patient_nr))
        
        
#############################################################
# record = wfdb.rdrecord(base_dir_raw, sampto = 500)
# # #record = wfdb.rdsamp(base_dir_raw, sampto = 1000)
# annotation = wfdb.rdann(base_dir_raw, 'pwave', sampto = 500)
# annotation2 = wfdb.rdann(base_dir_raw, 'qrs', sampto = 500)

# wfdb.plot_wfdb(record=record, annotation=annotation,
#                   title='Załamek P MIT-BIH rekord 214 ',
#                   time_units='seconds')
# wfdb.plot_wfdb(record=record, annotation=annotation2,
#                   title='Ecgpuwave MIT-BIH rekord 214 ',
#                   time_units='seconds')

# # # wfdb.show_ann_labels()
#############################################################
sig, fields = wfdb.rdsamp(base_dir_raw, channels=[0], sampto=5000)
ann_ref = wfdb.rdann(base_dir_raw,'atr', sampto=5000)


xqrs = ecgpu_detector.XQRS(sig=sig[:,0], fs=fields['fs'])
xqrs.detect()
# Alternatively, use the gateway function to get the QRS indices directly
# qrs_inds = ecgpu_detector.xqrs_detect(sig=sig[:,0], fs=fields['fs'])

# Compare detected QRS complexes to reference annotation.
# Note, first sample in 100.atr is not a QRS.
comparitor = processing.compare_annotations(ref_sample=ann_ref.sample[1:],
                                            test_sample=xqrs.qrs_inds,
                                            window_width=int(0.1 * fields['fs']),
                                            signal=sig[:,0])

# Print and plot the results
comparitor.print_summary()
#comparitor.plot(title='Ecgpuwave detected QRS vs reference annotations')
comparitor.plot(title='Ecgpuwave CINC2 nagranie 2201')


############################################################
# def peaks_hr(sig, peak_inds, fs, title, figsize=(10, 5), saveto=None):
#     "Plot a signal with its peaks and heart rate"
#     # Calculate heart rate
#     hrs = processing.hr.compute_hr(sig_len=sig.shape[0], qrs_inds=peak_inds, fs=fs)
    
#     N = sig.shape[0]
    
#     fig, ax_left = plt.subplots(figsize=figsize)
#     ax_right = ax_left.twinx()
    
#     ax_left.plot(sig, color='#3979f0', label='Signal')
#     ax_left.plot(peak_inds, sig[peak_inds], 'rx', marker='x', 
#                   color='#8b0000', label='Peak', markersize=12)
#     ax_right.plot(np.arange(N), hrs, label='Heart rate', color='m', linewidth=2)

#     ax_left.set_title(title)

#     ax_left.set_xlabel('Czas (ms)')
#     ax_left.set_ylabel('ECG (mV)', color='#3979f0')
#     ax_right.set_ylabel('HR (bpm)', color='m')
#     # Make the y-axis label, ticks and tick labels match the line color.
#     ax_left.tick_params('y', colors='#3979f0')
#     ax_right.tick_params('y', colors='m')
#     if saveto is not None:
#         plt.savefig(saveto, dpi=100)
#     plt.show()

# # Load the WFDB record and the physical samples
# record = wfdb.rdrecord(base_dir_raw, sampfrom=0, sampto=5000, channels=[0])

# # Use the XQRS algorithm to detect QRS locations in the first channel
# qrs_inds = ecgpu_detector.xqrs_detect(sig=record.p_signal[:,0], fs=record.fs)
# p_ref = wfdb.rdann(base_dir_raw,'qrs', sampto=5000)

# # Plot results
# # peaks_hr(sig=record.p_signal, peak_inds=qrs_inds, fs=record.fs,
# #          title="Ecgpuwave peak detection on record 100")
# peaks_hr(sig=record.p_signal, peak_inds=qrs_inds, fs=record.fs,
#           title="Ecgpuwave CINC2 nagranie 2728")