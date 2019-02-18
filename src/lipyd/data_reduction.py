import pyopenms as pms


class MSProcessing(object):
    '''
    '''
    def __init__(self,
                in_file_name=None,
                out_file_name=None):
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.in_map = pms.MSExperiment()
        self.out_map = pms.MSExperiment()
        pms.MzMLFile().load(self.in_file_name, self.in_map)

    def peakPicking(self):
        pp = pms.PeakPickerHiRes()
        pp.pickExperiment(self.in_map, self.out_map)
        self.out_map.updateRanges()
        pms.MzMLFile().store(self.out_file_name, self.out_map)

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    
        


if __name__ == "__main__":
    app = MSProcessing(in_file_name = "/home/igor/Documents/Black_scripts/Raw__data/STARD10_invivo_raw/mzml/150310_Popeye_MLH_AC_STARD10_A10_pos.mzML",
                         out_file_name = "STARD10_A10_pos_picked_HiRes2.mzML")

    app.peakPicking()
