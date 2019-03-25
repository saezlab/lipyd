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
import re

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
        self.set_parameters(**kwargs)

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
                src = ".+\.mzML$",               #"/path/to/src/.+\.mzML"
                dst = None,                     #/path/to/dst
                suffix_dst_files = "_feature",          #for example : "_feature"
                ext_dst_files = "featureXML",#the string may begin with a dot
                **kwargs
                ):
        
                                                
        if not src:
            raise RuntimeError( "You don`t specify all necessary files" )

        self.src = src
        self.dst = dst
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs
        self.init_entity(**self.kw)
        self.path_parsing()

    
    def init_entity(self, **kwargs):
        self.mtd = MtdEntity(**kwargs)
        self.epd = EpdEntity(**kwargs)
        self.ffm = FfmEntity(**kwargs)
        self.output_mt = []
        self.splitted_mt = []
        self.filtered_mt = []
        self.chromatograms = [[]]       


    def convert_src_to_dst_file_name(self, src):
        file_name = os.path.splitext( os.path.basename(src) )[0] # get file name only;
        file_name += self.suffix_dst_files # add suffix;
        if self.ext_dst_files:
            #add dot in front of ext if dot is not present:
            self.ext_dst_files = ("."+self.ext_dst_files) \
                if self.ext_dst_files[0] != "." else self.ext_dst_files
        else:
            self.ext_dst_files = "" #empty ext;
        dst_file_name = file_name + self.ext_dst_files
        
        return dst_file_name    #dst file name without dst dir name;

    def path_parsing(self):
        self.src_full_name_list = []                #list of src full file name;
        src_dir = os.path.dirname(self.src)         #get dir name from src file name;
        src_dir = src_dir if src_dir else os.getcwd() #if src dir name is empty get current dir name;
        self.dst = src_dir if not self.dst else self.dst #if dst dir name is empty get src dir name or dst dir name;
        pattern = os.path.basename(self.src)        #get file name as pattern;
        pattern = pattern if pattern != "" else ".*" #if pattern is empty set any char pattern;
        #print("src_dir=", src_dir)
        #print("pattern=", pattern)
        
        for file_name in os.listdir(src_dir): #for all file name in src dir:
            full_file_name = os.path.join(src_dir, file_name) #to build full file name;
            if os.path.isfile( full_file_name ): #only for files, except dir name;
                #print("file_name=", file_name)
                match = None    # result re match
                try:            # try to compile patetrn:
                    rec = re.compile(pattern)
                    match = rec.match(file_name)    # apply pattern to file name
                except re.error as e:
                    raise RuntimeError("Match pattern error.") # raise exept if pattern error compile;
                
                if match: # if result of re match is ok
                    self.src_full_name_list.append(full_file_name) #take file name anf put into result list
                    #print("full_file_name=", full_file_name) #debug only;

    def main(self): #the 3 step in one;
        #after path_parsing method we have self.src_full_name_list
        for f in self.src_full_name_list:

            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)

            print("source file=", f)
            self.input_map = None
            self.input_map = oms.PeakMap() # the 1st step: load map;
            self.fm = None
            self.fm = oms.FeatureMap()
            oms.MzMLFile().load(f, self.input_map)
            # the 2nd step: apply_ffm;
            self.mtd.entity.run(self.input_map, self.output_mt)
            self.epd.entity.detectPeaks(self.output_mt, self.splitted_mt)
            self.ffm.entity.run(self.splitted_mt, self.fm, self.chromatograms)
            # the 3d step: is store result into file;
            dst_full_file_name = os.path.join(self.dst,\
                self.convert_src_to_dst_file_name(f) )
            #print("dst=",dst_full_file_name)
            oms.FeatureXMLFile().store(dst_full_file_name, self.fm)


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


    param = {"noise_threshold_int":10.0,
            "chrom_peak_snr":3.0,
            "chrom_fwhm":5.0,
            "mz_scoring_13C":"true",
            "report_convex_hulls":"true",
            "remove_single_traces":"true"
            }
 
    ffm = FeatureFindingMetabo(
                src = "/home/igor/Documents/Scripts/Data/Picked_data/.+\.mzML$",
                dst = "/home/igor/Documents/Scripts/Data/Picked_data/",
                suffix_dst_files = "_feature",          #for example : "_feature"
                ext_dst_files = "featureXML",           #the string may begin with a dot
                **param)
    
    ffm.main()
