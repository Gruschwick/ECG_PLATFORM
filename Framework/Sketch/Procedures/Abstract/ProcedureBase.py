"""
Abstract base class for Procedures, public methods implement interface for other classes
"""

class Base(object):

	#return list of atributes required by Run() (and Function indirectly)
	def requiredDataFields(self):
		pass

	#called by Comparator, if requiredDataFields differ from ones given by used Loader
	def AttributesError(self):
		pass

	#returns results of whatever is to be done on data using self.func (see constructor)
	def resultsFor(self,data):
		pass

	#Patients is a list of files, from which data will be loaded into BasicData structure
	#assumption - there may be different sets of patients for each segmenter
	def __init__(self, Patients):	
		self.patients = Patients


	#just checks if all fields required by Function are loaded by loader to BasicData
	def _canRunOn(self, DataFields):
		for field in self.requiredDataFields():
			if not field in DataFields:
				return False;
		return True;	
		

