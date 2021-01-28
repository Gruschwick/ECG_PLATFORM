# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 16:22:54 2019

@author: x
"""

import platform
import os
from os.path import join
import shutil
import glob
from Helpers.Physionet import LoadFile as PhysionetLoadFile 

###############################################################################
#
#   import mit-bih
#
if(platform.system() == "Linux"):
    base_dir_temp = join(os.getcwd(),"Datasets/Temp/mitdb")
    base_dir_raw = join(os.getcwd(),"Datasets/mitdb/data_raw")

    if not os.path.isdir(base_dir_temp):
        os.mkdir(base_dir_temp)

    get_cmd = "wget64 -r -l1 --no-parent {0} -P {1}".format(
        "https://www.physionet.org/physiobank/database/mitdb/",
        base_dir_temp
        )

    if os.system(get_cmd) == 0:
       print("wget success")
    else:
        print("!!! wget faild !!!")

    if not os.path.isdir(join(os.getcwd(),"Datasets/mitdb")):
        os.mkdir(join(os.getcwd(),"Datasets/mitdb"))
    
    if not os.path.isdir(join(os.getcwd(),"Datasets/mitdb/data")):
        os.mkdir(join(os.getcwd(),"Datasets/mitdb/data"))
    
    if not os.path.isdir(os.path.join(os.getcwd(),"Datasets/mitdb/data_raw")):
        os.mkdir(join(os.getcwd(),"Datasets/mitdb/data_raw"))
    
        for filename in sum( [glob.glob(os.path.join(base_dir_temp, "www.physionet.org/physiobank/database/mitdb/*.{0}".format(ext))) for ext in ['atr','dat','hea']],[]):
            shutil.copy(filename, base_dir_raw)
    
#for systems different than Linux
else:
    base_dir_temp = join(os.getcwd(),"Datasets\\Temp\\mitdb")
    base_dir_raw = join(os.getcwd(),"Datasets\\mitdb\\data_raw")

    if not os.path.isdir(base_dir_temp):
        os.mkdir(base_dir_temp)

    get_cmd = "wget64 -r -l1 --no-parent {0} -P {1}".format(
        "https://www.physionet.org/physiobank/database/mitdb/",
        base_dir_temp
        )

    if os.system(get_cmd) == 0:
        print("wget success")
    else:
        print("!!! wget faild !!!")

    if not os.path.isdir(join(os.getcwd(),"Datasets\\mitdb")):
        os.mkdir(join(os.getcwd(),"Datasets\\mitdb"))
    
    if not os.path.isdir(join(os.getcwd(),"Datasets\\mitdb\\data")):
        os.mkdir(join(os.getcwd(),"Datasets\\mitdb\\data"))
    
    if not os.path.isdir(os.path.join(os.getcwd(),"Datasets\\mitdb\\data_raw")):
        os.mkdir(join(os.getcwd(),"Datasets\\mitdb\\data_raw"))
    
        for filename in sum( [glob.glob(os.path.join(base_dir_temp, "www.physionet.org\physiobank\database\mitdb\*.{0}".format(ext))) for ext in ['atr','dat','hea']],[]):
            shutil.copy(filename, base_dir_raw)
    
    
    
