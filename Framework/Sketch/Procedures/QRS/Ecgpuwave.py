#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 03:59:01 2020

@author: mateusz
"""
from ..Abstract.ProcedureBase import Base
from biosppy.signals import ecg
from .Resume.Results import *
import wfdb
from wfdb import processing
from .Ecgpuwave import ecgpu_detector

class Ecgpuwave(Base):

	def requiredDataFields(self):
		#to be checked
		return ['signal', 'freq'] + resDataFields()

	def AttributesError(self):
		return "REQUIRED ATTRIBUTES: signals  freq  QRSann  dbcode  patient_nr  signal_index "


	def resultsFor(self, data):
		return results(ecgpu_detector.xqrs_detect(sig=data.signal, fs=data.freq), data,"Ecgpuwave") 
