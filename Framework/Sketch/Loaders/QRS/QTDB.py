from ..Abstract.LoaderBase import Base
from ...Helpers.Physionet import LoadFile, PhysionetConstants
from ...DataStruct.BasicData import BasicData
import re
import numpy as np

class Loader(Base):
	def DataFrom(self,path):
		data = LoadFile(path)
		
		ann_selector = [x in PhysionetConstants.all_beats for x in data.annotations.anntype]
		anntype_selected = np.array(data.annotations.anntype)[ann_selector]
		annsamp_selected = data.annotations.annsamp[ann_selector]

		ret = BasicData()
		ret.signal = data.data[:,0]
		ret.freq = data.freq
		ret.QRSann = annsamp_selected
		ret.dbcode = "QTDB"
		ret.patient_nr = re.findall("([a-zA-Z0-9]+)$",path)
		ret.signal_index = 0
		
		return [ret]
	
	def DataFields(self):
		#signal - array containing raw ECG signal
		#freq - sampling frequency
		#QRSann - array of indices on which QRS peaks were anotated
		#signal_index - nr of signal slice, here is only one, hence 0
		#patient_nr - name of file from which data were loaded
		#dbcode - database name from which data are loaded
		return ["freq", "QRSann", "signal", "dbcode", "patient_nr", "signal_index"]