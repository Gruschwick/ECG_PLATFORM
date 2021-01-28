from ..Abstract.ProcedureBase import Base
from biosppy.signals import ecg
from .Resume.Results import *

class EngzeeSegmenter(Base):

	def requiredDataFields(self):
		#to be checked
		return ['signal', 'freq'] + resDataFields()

	def AttributesError(self):
		return "REQUIRED ATTRIBUTES: signals  freq  QRSann  dbcode  patient_nr  signal_index "


	def resultsFor(self, data):
		return results(ecg.engzee_segmenter(data.signal, data.freq)[0], data,"EngzeeSgmenter") 

