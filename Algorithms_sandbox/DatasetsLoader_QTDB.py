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

#C:\Users\blaze\.spyder-py3\ECG_ALG_EXAMPLE\Datasets

if(platform.system() == "Linux"):
    base_dir_temp = join(os.getcwd(),"Datasets/Temp/qtdb")
    base_dir_raw = join(os.getcwd(),"Datasets/qtdb/data_raw")

    if not os.path.isdir(base_dir_temp):
        os.mkdir(base_dir_temp)

    get_cmd = "wget64 -r -l1 --no-parent {0} -P {1}".format(
        "https://www.physionet.org/physiobank/database/qtdb/",
        base_dir_temp
        )

    if os.system(get_cmd) == 0:
       print("wget success")
    else:
        print("!!! wget faild !!!")

    if not os.path.isdir(join(os.getcwd(),"Datasets/qtdb")):
        os.mkdir(join(os.getcwd(),"Datasets/qtdb"))
    
    if not os.path.isdir(join(os.getcwd(),"Datasets/qtdb/data")):
        os.mkdir(join(os.getcwd(),"Datasets/qtdb/data"))
    
    if not os.path.isdir(os.path.join(os.getcwd(),"Datasets/qtdb/data_raw")):
        os.mkdir(join(os.getcwd(),"Datasets/qtdb/data_raw"))
    
        for filename in sum( [glob.glob(os.path.join(base_dir_temp, "www.physionet.org/physiobank/database/qtdb/*.{0}".format(ext))) for ext in ['atr','dat','hea']],[]):
            shutil.copy(filename, base_dir_raw)
    
#for systems different than Linux
else:
    base_dir_temp = join(os.getcwd(),"Datasets\\Temp\\qtdb")
    base_dir_raw = join(os.getcwd(),"Datasets\\qtdb\\data_raw")

    if not os.path.isdir(base_dir_temp):
        os.mkdir(base_dir_temp)

    get_cmd = "wget64 -r -l1 --no-parent {0} -P {1}".format(
        "https://www.physionet.org/physiobank/database/qtdb/",
        base_dir_temp
        )

    if os.system(get_cmd) == 0:
        print("wget success")
    else:
        print("!!! wget faild !!!")

    if not os.path.isdir(join(os.getcwd(),"Datasets\\qtdb")):
        os.mkdir(join(os.getcwd(),"Datasets\\qtdb"))
    
    if not os.path.isdir(join(os.getcwd(),"Datasets\\qtdb\\data")):
        os.mkdir(join(os.getcwd(),"Datasets\\qtdb\\data"))
    
    if not os.path.isdir(os.path.join(os.getcwd(),"Datasets\\qtdb\\data_raw")):
        os.mkdir(join(os.getcwd(),"Datasets\\qtdb\\data_raw"))
    
        for filename in sum( [glob.glob(os.path.join(base_dir_temp, "www.physionet.org\physiobank\database\qtdb\*.{0}".format(ext))) for ext in ['atr','dat','hea']],[]):
            shutil.copy(filename, base_dir_raw)
    
    
    
