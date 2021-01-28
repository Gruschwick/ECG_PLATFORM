import numpy as np

class BasicData(object):
	def __init__(self,raw_signal, QRS_ann, frequency, dbcode, patient_nr):
		self.signals = raw_signal
		self.peaks_annotations = QRS_ann
		self.freq = frequency
		self.dbcode = dbcode
		self.patient = patient_nr     

'''create dictionaries of {'low':low_index_val, 'high':high_index_val}
	and put them into array
input data is an array of indices, where values goes like:
	[1,2,(...),199,200,   400,401,(...),560,561 ]
	[<n, n+1, n+2, ...>,  <m, m+1, m+2,... >	]	
ouput is an array of dictionaries, where:
	[{'low': low_index_val1, 'high':high_index_val1},
	 {'low': low_index_val2, 'high':high_index_val2},
	 ...
	 {'low': low_index_valn, 'high':high_index_valn}]
	
	for exaple above:
	[{'low':1,   'high':200}, 
	 {'low':400, 'high':561}]
'''
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


def LoadFile(filename,frequency, dbcode):
	#for CSV format
	data = np.loadtxt(filename, delimiter=',')
	ret = []

	#countinuos indices from data where mask changes
	mask_change = np.where(data[0:,2]==0)[0]
	indices_pairs = MaskBoundariesFrom(mask_change)

	#getting rid of bad data (with annotated artifacts)
	for index in indices_pairs:
		tmp_data = data[index['low']:index['high'],0:]
		qrs_indices = np.nonzero(tmp_data[0:,1])[0]
		ret.append( BasicData(tmp_data[0:, 0], qrs_indices, frequency, dbcode, filename) )
	
	return ret




