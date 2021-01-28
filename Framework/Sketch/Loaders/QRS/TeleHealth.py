from ..Abstract.LoaderBase import Base
from ...DataStruct.BasicData import BasicData

import re
import numpy as np

def MaskBoundariesFrom(mask_change):
	indices = []
	low = mask_change[0]

	for i in range(1,len(mask_change)):
		if (mask_change[i-1]+1 != mask_change[i]):
			indices.append({'low':low, 'high':mask_change[i-1]})
			low = mask_change[i]
		
		if (i == len(mask_change)-1):
			indices.append( {'low':low, 'high':mask_change[i]} )
			
	return indices		



class Loader(Base):

	def __init__(self, frequecy = 500):
		self.frequecy = frequecy
		super().__init__()
		

	def DataFrom(self,path):
		data = np.loadtxt(path + ".dat", delimiter=',')
		ret = []

		#countinuos indices from data where mask changes
		mask_change = np.where(data[0:,2]==0)[0]
		indices_pairs = MaskBoundariesFrom(mask_change)

		#getting rid of bad data (with annotated artifacts)
		i = 0
		for index in indices_pairs:
			tmp_data = data[index['low']:index['high'],0:]
			qrs_indices = np.nonzero(tmp_data[0:,1])[0]

			tmp = BasicData()
			tmp.signal = tmp_data[0:,0]
			tmp.freq = self.frequecy
			tmp.QRSann = qrs_indices
			tmp.dbcode = "TeleHealth"
			tmp.patient_nr = re.findall("([a-zA-Z0-9_]+)$",path)
			tmp.signal_index = i			

			ret.append( tmp )
			i += 1

		return ret

	def DataFields(self):
		return ["signal", "freq", "QRSann", "dbcode", "patient_nr", "signal_index"]
	
