from Sketch.Loaders.QRS.MITBIH4 import Loader as MITBIH
from Sketch.Loaders.QRS.TeleHealth import Loader as TeleHealth

from Sketch.Procedures.QRS.Delineator import Delineator

from Sketch.Comparator import Comparator

def test_MITBIH(path):
	L = MITBIH()
	print(L.DataFields)
	print(L.DataFrom(path)[0].patient_nr)
    
def test_Delineator():
	D = Delineator([])
	print(D.requiredDataFields())
	print(MITBIH.DataFields())
	print(D._canRunOn(MITBIH.DataFields()))

def test_Comparator():
    #on_Dpatients = ["100","101","102"]
    on_DSpatients2 = [str(i) for i in list(range(1,251))]
    
    CMP = Comparator()

    #CMP.AddDatabase(MITBIH(), "../Datasets/mitdb/data_raw")
    #CMP.Databases[-1].addProcedure( Delineator(on_Dpatients) )

    CMP.AddDatabase(TeleHealth(frequecy=500), "../Datasets/telehealth/data_raw")
    CMP.Databases[-1].addProcedure( Delineator(on_DSpatients2) )
    
    CMP.Run()

    for res in CMP.Results:
        print(res)

    print(CMP.Errors)	


test_Comparator()    