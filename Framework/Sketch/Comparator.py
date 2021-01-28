from .DataStruct.Database import Database
import os

class Comparator(object):

	def __init__(self):
		#lists of:
		self.Results = []
		self.Errors = []
		self.Databases = []

	def AddDatabase(self,loader, path):
		self.Databases.append( Database(path, loader) )

	#saves results to file (CSV for ease?)
	#def Save(values_to_be_saved, filename_under_which_they_are_to_be_saved):
		#no idea for now on how to do this
		

	def Run(self):
		results = []
		for database in self.Databases:
			for Procedure in database.Procedures:
				if(Procedure._canRunOn(database.DataFields())):
					for patient in Procedure.patients:
						for data in database.load(os.path.join(database.path, patient)):
							try:
								results.append(Procedure.resultsFor(data))
							except Exception as e:
								print(e)
								self.Errors.append(e)
				else:
					self.Errors.append(Procedure.AttributesError())				
		self.Results = results


"""
Usage:

from Loaders.TeleHealth.teleECG import LoadFile as teleECG
from Loaders.Physionet.MITBIH import LoadFile as MITBIH
from Loaders.Physionet.CINC import LoadFile as CINC

from Procedures.QRS.HamiltonSegmenter import HamiltonSegmenter
from Procedures.QRS.MyWave import MyWave

from Framework import Comparator


Comparator cmp

cmp.AddDatabase(teleECG, "Datasets/telehealth/data_raw")
cmp.Database[-1].addProcedure(MyWave( MyWave_patients ))
cmp.Database[-1].addProcedure(HamiltonSegmenter( HamSeg_patients ) )

cmp.AddDatabase(MITBIH, "Datasets/mitdb/data_raw")
cmp.Database[-1].addProcedure(MyWave( MyWave_patients ))
cmp.Database[-1].addProcedure(HamiltonSegmenter( HamSeg_patients ) )

cmp.AddDatabase(CINC, "Datasets/cinc1/data_raw")
cmp.Database[-1].addProcedure(MyWave( MyWave_patients ))
cmp.Database[-1].addProcedure(HamiltonSegmenter( HamSeg_patients ) )

cmp.Run()
"""
