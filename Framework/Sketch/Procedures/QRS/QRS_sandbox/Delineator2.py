import numpy as np
from ...Helpers.ECGprocessing import delineateMultiLeadECG as fcg
from ..Abstract.ProcedureBase import Base
from biosppy.signals import ecg
from .Resume.Results import *

class Delineator(Base):

	def requiredDataFields(self):
		#to be checked
		return ['signal', 'freq'] + resDataFields()

	def AttributesError(self):
		return "REQUIRED ATTRIBUTES: signals  freq  QRSann  dbcode  patient_nr  signal_index "


	def resultsFor(self, data):
        	ECGdelin = fcg.delineateMultiLeadECG(data.signal[:,np.newaxis],data.freq)
        return results(ECGdelin[0][:,4], data, "Delineator")
        #return results(ECGdelin[0][:,4], data, "Delineator") 
