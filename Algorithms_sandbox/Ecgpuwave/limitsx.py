import platform
import os
import wfdb
import glob
from scipy import signal
import math
import numpy as np
import pandas as pd
from os.path import join
from biosppy.signals import ecg
from Helpers.Physionet2 import LoadFile, PhysionetConstants
from Helpers.Metrices import *
from subprocess import call
from collections import namedtuple
import traceback
import csv
import h5py
import Helpers.ECGprocessing as fcg
import matplotlib.pyplot as plt
from wfdb import processing
from .matlab import main


databases = ["qtdb2"]*105

patients = ['sel100','sel102','sel103','sel104','sel114','sel116','sel117','sel123',
            'sel14046','sel14157','sel14172','sel15814','sel16265','sel16272','sel16273',
            'sel16420','sel16483','sel16539','sel16773','sel16786','sel16795','sel17152',
            'sel17453','sel213','sel221','sel223','sel230','sel231','sel232','sel233',
            'sel301','sel302','sel306','sel307','sel308','sel31','sel310','sel32',
            'sel33','sel34','sel35','sel36','sel37','sel38','sel39','sel40','sel41',
            'sel42','sel43','sel44','sel45','sel46','sel47','sel48','sel49','sel50','sel51',
            'sel52','sel803','sel808','sel811','sel820','sel821','sel840','sel847','sel853',
            'sel871','sel872','sel873','sel883','sel891','sele0104','sele0106','sele0107',
            'sele0110','sele0111','sele0112','sele0114','sele0116','sele0121','sele0122',
            'sele0124','sele0126','sele0129','sele0133','sele0136','sele0166','sele0170',
            'sele0203','sele0210','sele0211','sele0303','sele0405','sele0406','sele0409',
            'sele0411','sele0509','sele0603','sele0604','sele0606','sele0607','sele0609',
            'sele0612','sele0704']

def fread(fid, nelements, dtype):

    if dtype is np.str:
        dt = np.uint8  
    else:
        dt = dtype

    data_array = np.fromfile(fid, dt, nelements)
    data_array.shape = (nelements, 1)

    return data_array


def limits (dbcode, patient_nr):
    
    ret = []
    
    if(platform.system() == "Linux"):
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
    else:
        base_dir_raw = os.path.join(os.getcwd(),"Datasets/"+dbcode+"/data_raw")
    
    if (dbcode=='ludb'):
        data = LoadFile(join(base_dir_raw, patient_nr), "atr_ii")
    else:
        data = LoadFile(join(base_dir_raw, patient_nr))
    
    res = 0
    nl =1
    nbo_flag = 0
    Kq = 1.5
    Kr = 5
    Ks = 3
    Krr = 5
    Kpb = 1.35
    Kpe = 2
    Ktb = 2
    Kte = 3.5
    pco = 8
    
    rmax=0
    index=0
    ilastqrs=0
    fin=0
    ns=0
    prewindt=0
    a1=1
    qt1=1
    q1=1
    r1=1
    s1=1
    qrs1=1
    nt1=1
    QTpos=[]
    QTval=[] 
    QTCpos=[]
    QTCval=[]
    QWpos=[]
    QWval=[]
    RWpos=[]
    RWval=[]
    SWpos=[]
    SWval=[]
    QRSpos=[]
    QRSval=[]
    anntPonset=[]
    anntP=[]
    anntPoffset=[]
    anntQRSonset=[]
    anntQ=[]
    anntR=[]
    anntfiducial=[]
    anntS=[]
    anntR2=[]
    anntQRSoffset=[]
    anntTonset=[]
    anntT=[]
    anntT2=[]
    anntToffset=[]
    annAmpP=[]
    annAmpQ=[]
    annAmpR=[]
    annAmpS=[]
    annAmpR2=[]
    annAmpT=[]
    annAmpT2=[]
    banottime=[]
    banotanntyp=[]
    banotsubtyp=[]
    banotnum=[]
    banotchan=[]
    
    ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
    anntype_selected = np.array(data.annotations.anntype)[ann_selector]
    annsamp_selected = data.annotations.annsamp[ann_selector]
    annsubtype = np.array(data.annotations.subtype)[ann_selector]
    annnum = data.annotations.num[ann_selector]
    annum = data.annotations.num[ann_selector]
    annaux = data.annotations.aux[ann_selector]
    annchan = data.annotations.chan[ann_selector]
    anottime = []
    
        #rpeaks = processing.gqrs_detect(sig=data.data[:,0], fs=data.freq)
    #ann_time = rpeaks.Annotation.
    # ann_time = processing.
        
    
    #rpeaks = processing.xqrs_detect(sig=data.data[:,0], fs=data.freq)
    
    #nqrs = np.size(rpeaks)
    
    sig=data.data[:,0]
    fs = data.getFreq
    # nsamp = np.size(data.data[:,0])
    # nsamp = np.size(sig)
    nsamp = sig.shape[0]
    dim = nsamp/ (fs*2)
    anottime = np.zeros(dim, 1)
    pos = 0
    i = 0
    
    with open(join(base_dir_raw, patient_nr), 'rb') as fid:
        data_array = np.fromfile(fid, str, 2)
        data_array.shape = (2, 1)
        
        for i in fid:
            A = data_array[0] + data_array[1] * 256
            II = A % 1024
            A = math.floor(A/1024)
        
            if A == 59:
                skip = np.fromfile(fid,str,4)
                skip.shape = (4,1)
                I = skip[2]+256*(skip[3]+256*(skip[0]+256*skip[1]))
                pos=pos+I
                I = 0
                
            pos = pos + II
            anottime[i] = pos
            
            i += i
            
    if i > dim:
        anottime[i-1: dim] = []
      
    
    if anottime[0]<0.5*fs:
        anottime[0] = []
        
        
    no_beats = np.size(anottime)

    
    nlat = 100
    iqrs = np.zeros(1, nlat)
    ilat = 1
    nlat = min(nlat, np.size(anottime))
    
    while ilat < no_beats-3:
        flat = ilat +nlat-3
        if nlat+3>=np.size(anottime):
            flat = np.size(anottime)
            nlat = nlat-3
        
        ta = anottime[ilat-1:flat]
        antyp=anntype_selected[ilat-1:flat]
        
        Xpa, Xbp, D, F, Der = main2.lynfilt(sig, fs, ns)
        
        n = 0
        nqrs = 0
        i = 0
        iqrs = []
        for i in np.size(ta):
            if nbo_flag == 0 or (nbo+flag ==1 and antyp[i] == "N"):
                ibe = max(1, ta[i-1] - round(0.2*fs))
                ien = ta[i-1] + round (0.17*fs)
                if ien<np.size(Xpa):
                    
                    ymax=max(Xpa[ibe-1:ien])
                    imax=max(Xpa[ibe-1:ien])
                    imax=ibe+imax-1
                    ymin=min(Xpa[ibe-1:ien])
                    imin=min(Xpa[ibe-1:ien])
                    imin=ibe+imin-1 
 
                if abs(ymin)>1.3*abs(ymax):
                    iqrs[n-1]=imin
                else:
                    iqrs[n-1]=imax
                
                atyp[n] = antyp[i]
                nqrs=nqrs+1
                
                n=n+1;
     

    index=1;
   
    if iqrs[0]>ilastqrs+0.4*Fs:
        iprimerqrs=1
    elif iqrs[1]>ilastqrs+0.4*Fs:
        iprimerqrs=2
    elif iqrs[2]>ilastqrs+0.4*Fs:
        iprimerqrs=3
    elif iqrs[4]<=ilastqrs+0.4*Fs:
        iprimerqrs=4 


#%Detection of limits and peaks.
# [POS,AMP,ANNOT,POS_ANNOT,NUMFIELD,SUBTYPEFIELD,CHANFIELD,POS_QT,VAL_QT,VAL_QTC,AMP_Q,POS_Q,AMP_R,POS_R,AMP_S,POS_S,VAL_QRS,POS_QRS,prewindt]=proces...
# (fidan,X,Xpa,Xpb,D,F,Der,ti,tf,iprimerqrs,nqrs,iqrs,atyp,ns,Fs,nl,res,prewindt,Kq,Kr,Ks,Krr,Kpb,Kpe,Ktb,Kte,pco); 

    Pwav, QRS, Twav = proces(sig, Xpa, Xpb, D,F,Der,iprimerqrs, iqrs, antyp, ns, fs, nl, res,prewindt,Kq,Kr,Ks,Krr,Kpb,Kpe,Ktb,Kte,pco)
    
    
for j in range(len(databases)):
    print(j, databases[j], patients[j])
    a=limits(databases[j], patients[j])

# %Recording of results.
# if (~isempty(POS_ANNOT))
# a2=a1+length(POS_ANNOT)-1;
# banot.time(a1:a2,1)=POS_ANNOT;
# banot.anntyp(a1:a2,1)=setstr(ANNOT); 
# banot.subtyp(a1:a2,1)=num2str(SUBTYPEFIELD);
# banot.chan(a1:a2,1)=num2str(CHANFIELD);
# banot.num(a1:a2,1)=num2str(NUMFIELD);     
# a1=a2+1;
# if ~isempty(POS_QT)
# qt2=qt1+length(POS_QT)-1;
# QT.pos(qt1:qt2,1)=POS_QT;
# QT.val(qt1:qt2,1)=VAL_QT;
# QTC.pos(qt1:qt2,1)=POS_QT;
# QTC.val(qt1:qt2,1)=VAL_QTC;
# qt1=qt2+1;
# end
# if ~isempty(POS_Q)
# q2=q1+length(POS_Q)-1;
# QW.pos(q1:q2,1)=POS_Q;
# QW.val(q1:q2,1)=AMP_Q;
# q1=q2+1;
# end
# if ~isempty(POS_R)
# r2=r1+length(POS_R)-1;
# RW.pos(r1:r2,1)=POS_R;
# RW.val(r1:r2,1)=AMP_R;
# r1=r2+1;
# end
# if ~isempty(POS_S)
# s2=s1+length(POS_S)-1;
# SW.pos(s1:s2,1)=POS_S;
# SW.val(s1:s2,1)=AMP_S;
# s1=s2+1;
# end
# if ~isempty(POS_QRS)
# qrs2=qrs1+length(POS_QRS)-1;
# QRS.pos(qrs1:qrs2,1)=POS_QRS;
# QRS.val(qrs1:qrs2,1)=VAL_QRS;
# qrs1=qrs2+1;
# end

# nt2=nt1+length(POS.Ponset)-1;
# annt.Ponset(nt1:nt2)=POS.Ponset; 
# annt.P(nt1:nt2)=POS.P;
# annt.Poffset(nt1:nt2)=POS.Poffset;
# annt.QRSonset(nt1:nt2)=POS.QRSonset;
# annt.Q(nt1:nt2)=POS.Q;
# annt.R(nt1:nt2)=POS.R;
# annt.fiducial(nt1:nt2)=POS.fiducial;
# annt.S(nt1:nt2)=POS.S;
# annt.R2(nt1:nt2)=POS.R2;
# annt.QRSoffset(nt1:nt2)=POS.QRSoffset;
# annt.Tonset(nt1:nt2)=POS.Tonset;
# annt.T(nt1:nt2)=POS.T;
# annt.T2(nt1:nt2)=POS.T2;
# annt.Toffset(nt1:nt2)=POS.Toffset;
# annAmp.P(nt1:nt2)=AMP.P;
# annAmp.Q(nt1:nt2)=AMP.Q;
# annAmp.R(nt1:nt2)=AMP.R;
# annAmp.S(nt1:nt2)=AMP.S;
# annAmp.R2(nt1:nt2)=AMP.R2;
# annAmp.T(nt1:nt2)=AMP.T;
# annAmp.T2(nt1:nt2)=AMP.T2;
# nt1=nt2+1;
# end
            
            
            
    
