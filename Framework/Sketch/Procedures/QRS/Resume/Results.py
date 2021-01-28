from ....Helpers.Metrices import *

#return standarized results, takes peaks found by segmenters and BasicData object
def results(rpeaks, data, segmenter):

	ref_peaks = sample_to_time(data.QRSann, data.freq)
	pred_peaks = sample_to_time(rpeaks, data.freq)
	#wyobrębnić statystykę do osobej klasy/modułu
	try:
		peaks_matching = match_peaks(ref_peaks, pred_peaks)
		ret_qrs_scores = qrs_detection_scores( ref_peaks, pred_peaks , peaks_matching)
		return [data.dbcode] + [data.patient_nr[0]+ "_" + str(data.signal_index)] + [segmenter] + [len(ref_peaks)] + [len(pred_peaks)] + [len(peaks_matching)] + list(ret_qrs_scores)
        
	except Exception as e:
		print(segmenter +" : "+ e)
		return [data.dbcode] + [data.patient_nr[0]+ "_" +str(data.signal_index)] + [segmenter] + [-1]*7 


#fields required in BasicData object for results(...) function above
def resDataFields():
	return ['QRSann', "dbcode", "patient_nr", "signal_index"]
