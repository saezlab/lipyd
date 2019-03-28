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


class PPEntity(BaseEntity):
    #oms.PeakPickingHiRes()
    def __init__(self, **kwargs):
        super(PPEntity, self).__init__(oms.PeakPickerHiRes(),
                                                **kwargs)

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

class MAEntity(BaseEntity):
    """oms.MapAlignmentAlgorithmPoseClustering"""
    def __init__(self, **kwargs):
        super(MAEntity, self).__init__(oms.MapAlignmentAlgorithmPoseClustering(),
                                                **kwargs)

class PeakPicking(object):

    def __init__(self,
                src = ".+\.mzML$",               #"/path/to/src/.+\.mzML"
                dst = None,                     #/path/to/dst
                suffix_dst_files = "_centroided",          #for example : "_feature"
                ext_dst_files = "mzML",#the string may begin with a dot
                centroid_out_map = None,
                input_map = None,
                **kwargs
                ):

        if not src:
            raise RuntimeError( "You don`t specify all necessary files" )

        self.src = src
        self.dst = dst
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs
        self.input_map = input_map
        self.centroid_out_map = centroid_out_map

        self.init_entity(**self.kw)
        self.path_parsing()

    def init_entity(self, **kwargs):

        self.pp = PPEntity(**kwargs)

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
                    self.src_full_name_list.append(full_file_name)


    def main(self): #the 3 step in one;
        #after path_parsing method we have self.src_full_name_list
        for f in self.src_full_name_list:

            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)

            print("source file=", f)
            
            self.input_map = oms.MSExperiment() # the 1st step: load map;

            oms.MzMLFile().load(f, self.input_map)

            self.centroid_out_map = oms.MSExperiment()

            # the 2nd step: apply_ffm;
            self.pp.entity.pickExperiment(self.input_map, self.centroid_out_map)
            
            self.centroid_out_map.updateRanges()

            # the 3d step: is store result into file;
            dst_full_file_name = os.path.join(self.dst,\
                self.convert_src_to_dst_file_name(f) )
            #print("dst=",dst_full_file_name)
            oms.MzMLFile().store(dst_full_file_name, self.centroid_out_map)


class FeatureFindingMetabo(object):
    
    def __init__(self,
                src = ".+\.mzML$",               #"/path/to/src/.+\.mzML"
                dst = None,                     #/path/to/dst
                suffix_dst_files = "_feature",          #for example : "_feature"
                ext_dst_files = "featureXML",#the string may begin with a dot
                input_map = None,
                fm = None,
                **kwargs
                ):
        
                                                
        if not src:
            raise RuntimeError( "You don`t specify all necessary files" )

        self.src = src
        self.dst = dst
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs
        self.input_map = input_map
        self.fm = fm
        
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
                    #print("full_file_name=", full_file_name) 


    def main(self): #the 3 step in one;
        #after path_parsing method we have self.src_full_name_list
        for f in self.src_full_name_list:

            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)

            #print("source file=", f)
            
            self.input_map = oms.PeakMap() # the 1st step: load map;
            
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
                src = ".+\.featureXML$",               #"/path/to/src/.+\.mzML"
                dst = None,                     #/path/to/dst
                suffix_dst_files = "_aligned",          
                ext_dst_files = "featureXML",
                reference_file = None,
                reference_map = None,
                toAlign_map = None,
                **kwargs
                ):
        
        if not (src and reference_file):
           raise RuntimeError( "You don`t specify all necessary files" )
        
        self.src = src
        self.dst = dst
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs
        
        self.init_entity(**self.kw)
        self.path_parsing()

        self.reference_file = reference_file
        self.reference_map = reference_map
        self.toAlign_map = toAlign_map


    def init_entity(self, **kwargs):

        self.ma = MAEntity(**kwargs)

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

            self.reference_map = oms.FeatureMap()
            self.toAlign_map = oms.FeatureMap()
            
            oms.FeatureXMLFile().load(self.reference_file, self.reference_map)
            oms.FeatureXMLFile().load(f, self.toAlign_map)
            
            #Set reference_map file
            self.ma.entity.setReference(self.reference_map)
            
            #3rd step create object for the computed transformation
            transformation = oms.TransformationDescription()

            # the 4rd step:
            self.ma.entity.align(self.toAlign_map, transformation)
            # the 5th step: is store result into file;
            dst_full_file_name = os.path.join(self.dst,\
                self.convert_src_to_dst_file_name(f) )
            
            #print("dst=",dst_full_file_name)
            oms.FeatureXMLFile().store(dst_full_file_name, self.toAlign_map)

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

    
    param = {"signal_to_noise": 1.0 }
    
    a = PeakPicking(
                src = "/home/igor/Documents/Scripts/Data/Raw_data_STARD10/.+\.mzML$",
                dst = "/home/igor/Documents/Scripts/Data/Picked_STARD10/",
                #reference_file = "/home/igor/Documents/Scripts/Data/feature_STARD10/A09_pos_picked_feature.featureXML",
                suffix_dst_files = "_picked",          #for example : "_feature"
                ext_dst_files = "mzML",           #the string may begin with a dot
                **param)
    a.main()
