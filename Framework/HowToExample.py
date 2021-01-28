from Sketch.Loaders.QRS.MITBIH import Loader as MITBIH
from Sketch.Loaders.QRS.TeleHealth import Loader as TeleHealth
from Sketch.Loaders.QRS.CINC1 import Loader as CINC1
from Sketch.Loaders.QRS.CINC2 import Loader as CINC2
from Sketch.Loaders.QRS.LuDB import Loader as LuDB
from Sketch.Loaders.QRS.QTDB import Loader as QTDB


from Sketch.Procedures.QRS.HamiltonSegmenter import HamiltonSegmenter
from Sketch.Procedures.QRS.EngzeeSegmenter import EngzeeSegmenter
#from Sketch.Procedures.QRS.Delineator import Delineator

from Sketch.Comparator import Comparator

def Run():
	on_patients1 = ["100","101","102"]
	on_patients2 = ["1_1","2_2","3_3"]
	on_patients3 = ["100","101","102"]
	on_patients4 = ["1003","1009","1016"]
  
	CMP = Comparator()
    
	CMP.AddDatabase(MITBIH(), "../Datasets/mitdb/data_raw")
	CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_patients1) )
	CMP.Databases[-1].addProcedure( EngzeeSegmenter(on_patients1) )
	CMP.Databases[-1].addProcedure( Delineator(on_patients1) )

	CMP.AddDatabase(TeleHealth(frequecy=500), "../Datasets/telehealth/data_raw")
	CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_patients2) )
	CMP.Databases[-1].addProcedure( Delineator(on_patients2) )
    
	CMP.AddDatabase(CINC1(), "../Datasets/cinc1/data_raw")
	CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_patients3) )
	CMP.Databases[-1].addProcedure( EngzeeSegmenter(on_patients3) )
	CMP.Databases[-1].addProcedure( Delineator(on_patients3) )
    
	CMP.AddDatabase(CINC2(), "../Datasets/cinc2/data_raw")
	CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_patients4) )
	CMP.Databases[-1].addProcedure( EngzeeSegmenter(on_patients4) )
	CMP.Databases[-1].addProcedure( Delineator(on_patients4) )
    

	CMP.Run()

	for res in CMP.Results:
		print(res)
    


Run()
