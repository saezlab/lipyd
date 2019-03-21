#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://denes.omnipathdb.org/
#

from future.utils import iteritems

import os
import imp

import pyopenms as oms




def utf8_to_bin(s):
    return s.encode("utf-8")


class BaseEntity(object):
    """
    Base class for setup and check parameter for all derived classes
    """
    def __init__(self,
                entity,
                **kwargs
                ):
        
        if not hasattr(entity, "getDefaults"):
            raise RuntimeError( "BaseEntity: the entity has no getDefaults attr" )
        
        self.entity = entity

    def set_parameters(self, **kwargs):
        """
        Check type of value: if it`s str - convert to binary string
        """
        param = self.entity.getDefaults()   #get default param;
        for k, v in kwargs.items():
            if isinstance(v, str):              
                param.setValue(k, utf8_to_bin(v))   
            else:
                param.setValue(k, v)               
        
        self.entity.setParameters(param)


class MapAlignment(object):

    def __init__(self,
                input_file1 = None,  #
                input_file2 = None,  #
                output_file1 = None, #
                output_file2 = None,  #
                ):
        
        if not input_file1 or not input_file2 or \
            not output_file1 or not output_file2:
                raise RuntimeError( "You don`t specify all necessary files" )
        
        self.input_file1 = input_file1  
        self.input_file2 = input_file2
        self.output_file1 = output_file1
        self.output_file2 = output_file2

        self.reference = reference
        self.toAlign = toAlign
        self.algorithm = algorithm


    def load_maps(self):
        self.reference = oms.FeatureMap()
        self.toAlign = oms.FeatureMap()
        oms.FeatureXMLFile().load(self.input_file1, self.reference)
        oms.FeatureXMLFile().load(self.input_file2, self.toAlign)
        
   
    def set_algorithm(self): 
        self.algorithm = BaseEntity(oms.MapAlignmentAlgorithmPoseClustering())

    def apply_align(self, **kwargs):
        self.load_maps()
        #self.set_algorithm() # recreate self.algorithm - may be overload;
        #self.set_parameters(**kwargs) #set default and set custom param;
        
        self.algorithm.set_parameters(**kwargs)
        self.algorithm.entity.setReference(self.reference)
        #create object for the computed transformation
        transformation = oms.TransformationDescription()  
        #align
        self.algorithm.entity.align(self.toAlign, transformation)

    def store_to_file(self):
        #store results
        oms.FeatureXMLFile().store(self.aligned_file1, self.reference)
        oms.FeatureXMLFile().store(self.aligned_file2, self.toAlign)

"""
Since I think the good way to overload all parameters of every process,
we can do that way. It`s removes repetitive methods.
"""


class MtdEntity(BaseEntity):
    """ oms.MassTraceDetection() """
    def __init__(self, **kwargs):
        super(MtdEntity, self).__init__(oms.MassTraceDetection(),
                                                **kwargs)

class EpdEntity(BaseEntity):
    """ oms.ElutionPeakDetection() """
    def __init__(self, **kwargs):
        super(EpdEntity, self).__init__(oms.ElutionPeakDetection(),
                                                **kwargs)

class FfmEntity(BaseEntity):
    """ oms.FeatureFindingMetabo() """
    def __init__(self, **kwargs):
        super(FfmEntity, self).__init__(oms.FeatureFindingMetabo(),
                                                **kwargs)

class FeatureFindingMetabo(object):
    
    def __init__(self,
                in_file_name = None,
                out_file_name = None,
                ):
        
                                                
        if not in_file_name:
            raise RuntimeError( "You don`t specify all necessary files" )

        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.output_mt = []
        self.splitted_mt = []
        self.filtered_mt = []
        self.chromatograms = [[]]

        self.mtd = MtdEntity()
        self.epd = EpdEntity()
        self.ffm = FfmEntity()
                

    def load_map(self):
        self.input_map = oms.PeakMap()
        self.fm = oms.FeatureMap()
        oms.MzMLFile().load(self.in_file_name, self.input_map)


    def apply_ffm(self):

        self.mtd.entity.run(self.input_map, self.output_mt)
        self.epd.entity.detectPeaks(self.output_mt, self.splitted_mt)
        self.ffm.entity.run(self.splitted_mt, self.fm, self.chromatograms)

    def store(self):
        oms.FeatureXMLFile().store(self.out_file_name, self.fm)
    
    def main(self):

        self.load_map()  
        self.mtd.set_parameters(**param)
        self.epd.set_parameters(**param)
        self.ffm.set_parameters(**param)
        self.apply_ffm()
        self.store()

class Convert2mgf():

    def __init__(self,
                mzml_file = None,
                mgf_file = None,

                ):
        
        self.mzml_file = mzml_file
        self.mgf_file = mgf_file

    def convert(self):

        file = pyopenms.MzMLFile()
        msdata = pyopenms.MSExperiment()
        file.load(self.mzml_file, msdata)

        outfile = open(self.mgf_file, "w")

        # Create header
        outfile.write("COM=Testfile\n")
        outfile.write("ITOL=1\n")
        outfile.write("ITOLU=Da\n")
        outfile.write("CLE=Trypsin\n")
        outfile.write("CHARGE=1,2,3\n")

        # Iterate through all spectra, skip all MS1 spectra and then write mgf format
        nr_ms2_spectra = 0
        for spectrum in msdata:
            if spectrum.getMSLevel() == 2:
                continue
            nr_ms2_spectra += 1
            outfile.write("\nBEGIN IONS\n")
            outfile.write("TITLE=%s\n" % spectrum.getNativeID())
            outfile.write("RTINSECONDS=%s\n" % spectrum.getRT())
            try:
                outfile.write("PEPMASS=%s\n" % spectrum.getPrecursors()[0].getMZ())
                ch = spectrum.getPrecursors()[0].getCharge()
                if ch > 0:
                    outfile.write("CHARGE=%s\n" % ch)
            except IndexError:
                outfile.write("PEPMASS=unknown\n")
            for peak in spectrum:
                outfile.write("%s %s\n" % (peak.getMZ(), peak.getIntensity() ))
            outfile.write("END IONS\n")

        if nr_ms2_spectra == 0:
            print("Did not find any MS2 spectra in your input, thus the output file is empty!")

        outfile.close()
    
if __name__ == "__main__":


    in_file_name = "/home/igor/Documents/lipyd/src/lipyd_ms_preproc/150310_Popeye_MLH_AC_STARD10_A09_pos/150310_Popeye_MLH_AC_STARD10_A09_pos__peaks.mzML"
    out_file_name = "150310_Popeye_MLH_AC_STARD10_A09_pos.featureXML"
    ffm = FeatureFindingMetabo(in_file_name = in_file_name,
                                out_file_name = out_file_name)
    
    param = {"noise_threshold_int":10.0,
            "chrom_peak_snr":3.0,
            "chrom_fwhm":5.0,
            "mz_scoring_13C":      "true",
            "report_convex_hulls":  "true",
            "remove_single_traces": "true"
            }

    ffm.main()
