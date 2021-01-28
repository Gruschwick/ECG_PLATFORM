import numpy as np
#from ...Helpers.ECGprocessing import delineateMultiLeadECG
#from QRS.ECGprocessing import delineateMultiLeadECG
#from .ECGprocessing import delineateMultiLeadECG
from .ECGprocessing import delineateMultiLeadECG
from ..Abstract.ProcedureBase import Base
from biosppy.signals import ecg
from .Resume.Results import *

class Delineator(Base):
        
        def requiredDataFields(self):
                
                return ['signal', 'freq'] + resDataFields()
        
        def AttributesError(self):
                
                return "REQUIRED ATTRIBUTES: signals freq QRSann dbcode patient_nr signal_index "
        
        #ECGdelin = delineateMultiLeadECG(data.data[:,0][:,np.newaxis],data.freq)
        
        def resultsFor(self, data):
                
                ECGdelin = delineateMultiLeadECG(data.signal[:,np.newaxis],data.freq)
                # ECGdelin = delineateMultiLeadECG(data.signal[:,0],data.freq)
                
                return results(ECGdelin[0][:,4], data, "Delineator")

