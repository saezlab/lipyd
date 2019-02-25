import pyopenms as pms
#from pyopenms import *

class Data_Processing(object):
    '''
    Main class for ms data processing using pyopenms library.

    Parameters
    ----------
    
    '''
    def __init__(self,
                raw_data_file = None,
                peaked_file = None,
                feature_file = None,
                seeds = None,
                fh = None,
                options = None,
                picked_input_map = None,
                ff = None,
                features = None,
                name = None,
                ff_params = None,
                featureXML_file_name=None,
                raw_map = None,
                aligned_map = None,
                input_fXML_1 = None,
                input_fXML_2 = None,
                ):
        
        #peak picking
        self.raw_data_file = raw_data_file
        self.peaked_file = peaked_file
        self.raw_map = pms.MSExperiment()
        self.picked_out_map = pms.MSExperiment()
        self.pp_params = pms.PeakPickerHiRes().getParameters() 


        #feature finder
        self.seeds = pms.FeatureMap()
        self.fh = pms.MzMLFile()
        self.options = pms.PeakFileOptions()
        self.picked_input_map = pms.MSExperiment()
        self.ff = pms.FeatureFinder()
        self.features = pms.FeatureMap()
        self.name = pms.FeatureFinderAlgorithmPicked().getProductName()
        self.ff_params = pms.FeatureFinder().getParameters(self.name)
        self.featureXML_file_name = featureXML_file_name

        #consensus map
        self.input_fXML_1 = input_fXML_1
        self.input_fXML_2 = input_fXML_2

        #map alignment
        self.alignment_map = pms.MSExperiment()
        self.aligned_map = aligned_map


    def peak_picking(self):

        pms.MzMLFile().load(self.raw_data_file, self.raw_map)
        pp = pms.PeakPickerHiRes()
        pp.setParameters(self.pp_params)
        pp.pickExperiment(self.raw_map, self.picked_out_map)
        self.picked_out_map.updateRanges()
        pms.MzMLFile().store(self.peaked_file, self.picked_out_map)


    def feature_detection(self):

        self.options.setMSLevels([1,1])
        self.fh.setOptions(self.options)
        self.fh.load(self.peaked_file, self.picked_input_map)
        self.picked_input_map.updateRanges()
        self.ff.setLogType(pms.LogType.CMD)
        
        # Run the feature finder
        self.ff.run(self.name, self.picked_input_map, self.features, self.ff_params, self.seeds)
        self.features.setUniqueIds()
        self.fh = pms.FeatureXMLFile()
        self.fh.store(self.featureXML_file_name, self.features)
     

    def consensus_mapping(self):

        self.features = pms.ConsensusFeature()   


    def map_alignment(self):

        self.fh.load(self.input_fXML_1, self.input_fXML_2)
        ma = pms.MapAlignmentAlgorithmIdentification()
        ma.align(self.input_fXML_1, self.input_fXML_2)
        self.fh.store(self.aligned_map)


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


if __name__ == "__main__":
    
    a = Data_Processing(raw_data_file = "/home/igor/Documents/Black_scripts/Data_Processing/Filter_preproc/STARD10_A10_pos_Gauss_removed_base_line.mzML",
                         peaked_file = "STARD10_A10_pos_picked_HiRes_RBL.mzML",
                        featureXML_file_name = "STARD10_A10_pos_HiRes_RBL.featureXML ",)
                        #input_fXML_1 = "STARD10_A09_pos.featureXML",
                        #input_fXML_2 = "STARD10_A10_pos.featureXML",
                        #aligned_map = "test_align_map.featureXML")
                        
    a.peak_picking()
    a.feature_detection()
    #a.map_alignment()