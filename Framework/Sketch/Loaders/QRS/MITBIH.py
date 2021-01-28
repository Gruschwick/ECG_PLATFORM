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
                ret.signal = data.data
                ret.freq = data.freq
                ret.QRSann = annsamp_selected
                ret.dbcode = "MIT_BIH"
                ret.patient_nr = re.findall("([a-zA-Z0-9]+)$",path)
                ret.signal_index = 0
                
                return [ret]
        
        def DataFields(self):
                
                return ["freq", "QRSann", "signal", "dbcode", "patient_nr", "signal_index"]