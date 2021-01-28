from Sketch.Loaders.QRS.MITBIH4 import Loader as MITBIH
from Sketch.Loaders.QRS.TeleHealth import Loader as TeleHealth
# from Sketch.Loaders.QRS.CINC1 import Loader as CINC1
# from Sketch.Loaders.QRS.CINC2 import Loader as CINC2
#from Sketch.Loaders.QRS.LuDB import Loader as LuDB
from Sketch.Loaders.QRS.QTDB import Loader as QTDB


from Sketch.Procedures.QRS.HamiltonSegmenter import HamiltonSegmenter
from Sketch.Procedures.QRS.EngzeeSegmenter import EngzeeSegmenter


from Sketch.Comparator import Comparator

#def test_TeleHealth(path):
	#T = TeleHealth()
##	print(T.DataFrom(path)[0].patient_nr)

#test_MITBIH("../Datasets/mitdb/data_raw/101")

def test_HamiltonSegmenter():
	H = HamiltonSegmenter([])
	print(H.requiredDataFields())
	print(MITBIH.DataFields())
	print(H._canRunOn(MITBIH.DataFields()))

#test_HamiltonSegmenter()

def test_Comparator():
    on_HSpatients = ["100","101","102"]
    # on_HSpatients2 = ["101","102","103"]
    on_HSpatients3 = ["1003","1009","1016"]
    on_HSpatients4 = ["sel30","sel31","sel32"]
    on_HSpatients2 = [str(i) for i in list(range(1,251))]
    
    #patients = [str(i) for i in list(range(100,200))]
    #on_HSpatients3 = ["100","101","102"]
    #on_HSpatients4 = ["1003","1009","1016"]
    #on_HSpatients5 = ["sel14046","sel14157", "sel14172"]
    CMP = Comparator()
    
    # CMP.AddDatabase(MITBIH(), "../Datasets/mitdb/data_raw")
    # CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_HSpatients) )
    # CMP.Databases[-1].addProcedure( EngzeeSegmenter(on_HSpatients) )
    
    # CMP.AddDatabase(CINC1(), "../Datasets/cinc1/data_raw")
    # CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_HSpatients2) )
    # CMP.Databases[-1].addProcedure( EngzeeSegmenter(on_HSpatients2) )

    # CMP.AddDatabase(CINC2(), "../Datasets/cinc2/data_raw")
    # CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_HSpatients3) )
    # CMP.Databases[-1].addProcedure( EngzeeSegmenter(on_HSpatients3) )

    #CMP.AddDatabase(QTDB(), "../Datasets/qtdb/data_raw")
    #CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_HSpatients4) )    
    #CMP.Databases[-1].addProcedure( EngzeeSegmenter(on_HSpatients5) )
    
    CMP.AddDatabase(TeleHealth(frequecy=500), "../Datasets/telehealth/data_raw")
    CMP.Databases[-1].addProcedure( HamiltonSegmenter(on_HSpatients2) )
    # CMP.Databases[-1].addProcedure( EngzeeSegmenter(on_HSpatients2) )

    
    CMP.Run()

    for res in CMP.Results:
        print(res)
        
        
    print(CMP.Errors)	


test_Comparator()


