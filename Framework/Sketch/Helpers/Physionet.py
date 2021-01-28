# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:07:16 2019

@author: x
"""

import wfdb
import numpy as np
from .Plot import PlotWithAnnotations
from .DataSeriesBase import DataSeriesBase

def LoadFile(file, ext = 'atr'):
    record = wfdb.rdsamp(file)
    annotations = wfdb.rdann(file,ext)
    return PhysionetData(file, record, annotations)

class PhysionetConstants(object):
    
    
    #beats classification
    N_beats = ["N", "L", "R", "e", "j"]
    
    SVEB_beats = ["A", "a", "J", "S"]
    #supraventricular ectopic
    
    VEB_beats = ["V", "E"]
    #ventricular ectopic
    
    F_beats = ["F"]
    #fusion beat
    
    Q_beats = ["/", "f", "Q"]
    #unknown beat
    
    all_beats = N_beats + SVEB_beats  + VEB_beats + F_beats + Q_beats
        
    # Conventional training set
    DS1 = ["101", "106", "108", "109", "112", "114", "115", "116", "118", "119", "122",
       "124", "201", "203", "205", "207", "208", "209", "215", "220", "223", "230"]

    # Conventional testing set
    DS2 = ["100", "103", "105", "111", "113", "117", "121", "123", "200", "202", "210",
       "212", "213", "214", "219", "221", "222", "228", "231", "232", "233", "234"]
    
    
    
    

class PhysionetAnnotationDecorator(object):
    def __init__(self, ann):
        #print(dir(ann))
        
        self.annsamp = ann.annsamp if 'annsamp' in dir(ann) else ann.sample 
        self.aux = ann.aux if 'aux' in dir(ann) else ann.aux_note
        self.anntype = ann.anntype if 'anntype' in dir(ann) else ann.symbol
        self.chan = ann.chan 
        self.ann_org = ann

class PhysionetData(DataSeriesBase):
    def __init__(self, file, record, annotations):
        super(PhysionetData,self).__init__(file)
        self.record = record
        self.annotations = PhysionetAnnotationDecorator(annotations)
        self.freq = record[1]['fs']
        self.data = record[0]
        
    def getFreq(self):
        return self.freq
        
    def getData(self):
        return self.data
        
    def plot(self, start = None, to = None, addRPeaks = False):
        annIdx = self.annotations.annsamp #if 'annsamp' in dir(self.annotations) else self.annotations.sample 
        annLabel = self.annotations.aux #if 'aux' in dir(self.annotations) else self.annotations.aux_note
        annType = self.annotations.anntype #if 'anntype' in dir(self.annotations.anntype) else self.annotations.symbol
        
        annList = list(map(lambda x: x[0]+"\n"+x[1], zip(annType,annLabel)))
        
        annColor = np.where(np.asarray(annType)=='N','r','b')
                
        for i in range(self.data.shape[1]):
            
            title = self.file + ' ' + self.record[1]['sig_name'][i]
                        
            s = self.data[:,i]
            t = np.arange(len(s))
                    
            PlotWithAnnotations(s,t, annIdx, annList, annColor, title, rpeaks = (self.features[i].getRPeaks() if addRPeaks else None), start=start, to=to)