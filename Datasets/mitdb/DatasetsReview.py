# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:16:48 2019

@author: x
"""

import os
import wfdb
import glob

###############################################################################
#
#   import mit-bih
#

def plot_wfdb_rec(base_dir_raw, file, sampfrom = 0, sampto = 2000):
    example_file = os.path.join(base_dir_raw, file)
    record = wfdb.rdsamp(example_file, sampfrom=sampfrom, sampto=sampto)    # record = (signals, fields)
    annotations = wfdb.rdann(example_file, "atr", sampfrom=sampfrom, sampto=sampto)

    wfdb.plotrec(record=record, annotation=annotations,
             figsize=(20, 5), timeunits="seconds", title='MIT-BIH Record 100')    

# Conventional training set
DS1 = ["101", "106", "108", "109", "112", "114", "115", "116", "118", "119", "122",
       "124", "201", "203", "205", "207", "208", "209", "215", "220", "223", "230"]

# Conventional testing set
DS2 = ["100", "103", "105", "111", "113", "117", "121", "123", "200", "202", "210",
       "212", "213", "214", "219", "221", "222", "228", "231", "232", "233", "234"]

wfdb.show_ann_labels()

base_dir_raw = os.path.join(os.getcwd(),"Datasets\\mitdb\\data_raw")

all_recordings = [x.split("\\")[-1].split(".")[0] for x in glob.glob(os.path.join(base_dir_raw,"*.dat"))]

len(set(all_recordings))
len(DS1)
assert len(set(DS1) & set(DS2)) == 0, "DS do not overlab"

for rec in DS1:
    annotations = wfdb.rdann(os.path.join(base_dir_raw, rec), "atr")
    annotations.get_contained_labels()
    print(annotations.contained_labels)
    print()
    
    
plot_wfdb_rec(base_dir_raw, "100")
plot_wfdb_rec(base_dir_raw, "119")

###############################################################################




