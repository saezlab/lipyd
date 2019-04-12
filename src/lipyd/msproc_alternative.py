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

import lipyd.log as log
import lipyd.common as common
import lipyd.settings as settings


def convert_src_to_dst_file_name(src, dst, suffix_dst_files, ext_dst_files):
    """
    Global function to transform source directory to destination directory

    """
    file_name = os.path.splitext( os.path.basename(src) )[0] # get file name only;
    file_name += suffix_dst_files # add suffix;
    if ext_dst_files:
        #add dot in front of ext if dot is not present:
        ext_dst_files = ("."+ext_dst_files) \
            if ext_dst_files[0] != "." else ext_dst_files
    else:
        ext_dst_files = "" #empty ext;
    dst_file_name = file_name + ext_dst_files
    
    return dst_file_name    #dst file name without dst dir name;


def get_list_full_names(src):
    """
    Global function for getting full names from list
    """
    src_full_name_list = []                #list of src full file name;
    src_dir = os.path.dirname(src)         #get dir name from src file name;
    src_dir = src_dir if src_dir else os.getcwd() #if src dir name is empty get current dir name;
    pattern = os.path.basename(src)        #get file name as pattern;
    pattern = pattern if pattern != "" else ".*" #if pattern is empty set any char pattern;
    
    for file_name in os.listdir(src_dir): #for all file name in src dir:
        full_file_name = os.path.join(src_dir, file_name) #to build full file name;
        if os.path.isfile( full_file_name ): #only for files, except dir name;
            match = None    # result re match

            try:            # try to compile patern:
                rec = re.compile(pattern)
                match = rec.match(file_name)    # apply pattern to file name
            except re.error as e:
                raise RuntimeError("Match pattern error.") # raise exept if pattern error compile;
            
            if match:       # if result of re match is ok
                src_full_name_list.append(full_file_name)

    return src_full_name_list


def utf8_to_bin(s):
    return s.encode("utf-8")


class BaseEntity(object):
    """
    Base class for setup and check parameter for all derived classes
    
    Parameters
    ----------
    entity : :obj:
        Variable takes parameters from class which has method setParameters (pyyopenms)

    Attributes
    ----------
    entity : obj
        Variable takes parameters from class which has method setParameters (pyyopenms)

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
    """oms.PeakPickingHiRes()"""
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
    """
    Class for peak picking implementation.
    
    This class implements a fast peak-picking algorithm best suited for high resolution MS data
    (FT-ICR-MS, Orbitrap). In high resolution data, the signals of ions with similar
    mass-to-charge ratios (m/z) exhibit little or no overlapping and therefore allow 
    for a clear separation. Furthermore, ion signals tend to show well-defined peak shapes
    with narrow peak width.

    This peak-picking algorithm detects ion signals in profile data and reconstructs
    the corresponding peak shape by cubic spline interpolation. Signal detection depends
    on the signal-to-noise ratio which is adjustable by the user (see parameter signal_to_noise).
    A picked peak's m/z and intensity value is given by the maximum of the underlying peak spline.

    So far, this peak picker was mainly tested on high resolution data. 
    With appropriate preprocessing steps (e.g. noise reduction and baseline subtraction),
    it might be also applied to low resolution data.

    Parameters
    ----------
    src : str
        Source directory consists source file(s)
    dst :  str, optional
        Destination directory
    suffix_dst_files : str, optional
        Additional part of result file name
    ext_dst_files: str, optional
        Extension of resulting files
    logger:   
        System variable for tracking

    Attributes
    ----------
    src : str
        Source directory consists source file(s)  
    dst :  str, optional
        Destination directory
    suffix_dst_files : str, optional
        Additional part of result file name
    ext_dst_files: str, optional
        Extension of resulting files
    kw : obj
        Additional arguments


    """

    def __init__(self,
                src = ".+\.mzML$",               
                dst = None,                     
                suffix_dst_files = "_centroided",        
                ext_dst_files = "mzML",
                logger = None,
                **kwargs
            ):

        self.log = (
            logger or
            log.new_logger(
                'lipyd.msproc',
                logdir = 'lipyd_log',
            )
        )

        if not src:
            raise RuntimeError( "You don`t specify all necessary files" )

        self.src = src
        self.dst = dst
        
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs

        self.init_entity(**self.kw)
  

    def init_entity(self, **kwargs):

        self.pp = PPEntity(**kwargs)


    def main(self):
        #after path_parsing method we have self.src_full_name_list
        for f in get_list_full_names(self.src):   #call 'global' function;

            self.log.msg(
            'Performing peak picking on experiment `%s`.' % self.src
        )
        
            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)

            print("source file=", f)
            
            input_map = oms.MSExperiment() # the 1st step: load map;

            oms.MzMLFile().load(f, input_map)

            centroid_out_map = oms.MSExperiment()

            # the 2nd step: apply_ffm;
            self.pp.entity.pickExperiment(input_map, centroid_out_map)
            
            centroid_out_map.updateRanges()

            # the 3d step: is store result into file:
            #convert_src_to_dst_file_name(src, dst, suffix_dst_files, ext_dst_files)
            dst_full_file_name = os.path.join(self.dst,\
                convert_src_to_dst_file_name(f,
                                            self.dst,
                                            self.suffix_dst_files,
                                            self.ext_dst_files
                                            ) )   #call 'global' function;
            #print("dst=",dst_full_file_name)
            oms.MzMLFile().store(dst_full_file_name, centroid_out_map)

            self.log.msg(
            'Peak picking finished. Centroid data has been '
            'written to `%s`.' % self.dst
        )

class FeatureFindingMetabo(object):
    """
    Class for feature detection implementation.
    
    Method for the assembly of mass traces belonging to the same isotope pattern, i.e.,
    that are compatible in retention times, mass-to-charge ratios, and isotope abundances.

    Parameters
    ----------
    src : str
        Source directory consists source file(s)
    dst :  str, optional
        Destination directory
    suffix_dst_files : str, optional
        Additional part of result file name
    ext_dst_files: str, optional
        Extension of resulting files
    logger:   
        System variable for tracking

    Attributes
    ----------
    src : str
        Source directory consists source file(s)  
    dst :  str, optional
        Destination directory
    suffix_dst_files : str, optional
        Additional part of result file name
    ext_dst_files: str, optional
        Extension of resulting files
    kw : obj
        Additional arguments

    """
    def __init__(self,
                src = ".+\.mzML$",               
                dst = None,                     
                suffix_dst_files = "_feature",  
                ext_dst_files = "featureXML",
                logger = None,
                **kwargs
                ):
            
        self.log = (
            logger or
            log.new_logger(
                'lipyd.msproc',
                logdir = 'lipyd_log',
            )
        )
                                                
        if not src:
            raise RuntimeError( "You don`t specify all necessary files" )

        self.src = src
        self.dst = dst
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs
        
        self.init_entity(**self.kw)
        

    def init_entity(self, **kwargs):
        self.mtd = MtdEntity(**kwargs)
        self.epd = EpdEntity(**kwargs)
        self.ffm = FfmEntity(**kwargs)
        self.output_mt = []
        self.splitted_mt = []
        self.filtered_mt = []
        self.chromatograms = [[]]       

    def main(self): #the 3 step in one;
        #after path_parsing method we have self.src_full_name_list
        for f in get_list_full_names(self.src):

            self.log.msg(
            'Performing feature finding metabo on experiment `%s`.' % self.src
            )
        
            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)
            
            self.input_map = oms.PeakMap() # the 1st step: load map;
            fm = oms.FeatureMap()
            oms.MzMLFile().load(f, self.input_map)
            # the 2nd step: apply_ffm;
            self.mtd.entity.run(self.input_map, self.output_mt)
            self.epd.entity.detectPeaks(self.output_mt, self.splitted_mt)
            self.ffm.entity.run(self.splitted_mt, fm, self.chromatograms)
            # the 3d step: is store result into file;
            dst_full_file_name = os.path.join(self.dst,\
                self.convert_src_to_dst_file_name(f,
                                            self.dst,
                                            self.suffix_dst_files,
                                            self.ext_dst_files) )
            #print("dst=",dst_full_file_name)
            oms.FeatureXMLFile().store(dst_full_file_name, fm)

            self.log.msg(
            'Feature finding finished. Centroided data has been '
            'written to `%s`.' % self.dst
            )

class MapAlignment(object):
    """
    Class for map alignment process.

    A map alignment algorithm based on pose clustering.

    Pose clustering analyzes pair distances to find the most probable
    transformation of retention times.
    The algorithm chooses the x most intensity peaks/features per map. 
    This is modeled via the parameter 'max_num_peaks_considered',
    which in turn influences the runtime and stability of the results. 
    Bigger values prolong computation, smaller values might lead to no or
    unstable trafos. Set to -1 to use all features (might take very long for large maps).

    Parameters
    ----------
    src : str
        Source directory consists source file(s)
    dst :  str, optional
        Destination directory
    suffix_dst_files : str, optional
        Additional part of result file name
    ext_dst_files: str, optional
        Extension of resulting files
    logger:   
        System variable for tracking
    reference_file: obj
        The file by which other files will be aligned


    Attributes
    ----------
    src : str
        Source directory consists source file(s)  
    dst :  str, optional
        Destination directory
    suffix_dst_files : str, optional
        Additional part of result file name
    ext_dst_files: str, optional
        Extension of resulting files
    kw : obj
        Additional arguments
    reference_file: obj
        The file by which other files will be aligned
    """

    def __init__(self,
                src = ".+\.featureXML$",               #"/path/to/src/.+\.mzML"
                dst = None,                     #/path/to/dst
                suffix_dst_files = "_aligned",          
                ext_dst_files = "featureXML",
                reference_file = None,
                logger = None,
                **kwargs
                ):
        
        self.log = (
            logger or
            log.new_logger(
                'lipyd.msproc',
                logdir = 'lipyd_log',
            )
        )

        if not (src and reference_file):
           raise RuntimeError( "You don`t specify all necessary files" )
        
        self.src = src
        self.dst = dst
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs
        
        self.init_entity(**self.kw)

    def init_entity(self, **kwargs):

        self.ma = MAEntity(**kwargs)


    def main(self): #the 3 step in one;
        #after path_parsing method we have self.src_full_name_list
        
        for f in get_list_full_names(self.src):

            self.log.msg(
            'Performing feature finding metabo on experiment `%s`.' % self.src
            )
            
            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)

            reference_map = oms.FeatureMap()
            toAlign_map = oms.FeatureMap()
            
            oms.FeatureXMLFile().load(self.reference_file, reference_map)
            oms.FeatureXMLFile().load(f, toAlign_map)
            
            #Set reference_map file
            self.ma.entity.setReference(reference_map)
            
            #3rd step create object for the computed transformation
            transformation = oms.TransformationDescription()

            # the 4rd step:
            self.ma.entity.align(toAlign_map, transformation)
            # the 5th step: is store result into file;
            dst_full_file_name = os.path.join(self.dst,\
                self.convert_src_to_dst_file_name(f) )
            
            #print("dst=",dst_full_file_name)
            oms.FeatureXMLFile().store(dst_full_file_name, toAlign_map)

            self.log.msg(
            'Map alignmnet finished. Identification data has been '
            'written to `%s`.' % self.dst
            )

class Convert2mgf():
    """
    Class for convertation mzml data to MGF format (MS2 data)
    
    Parameters
    ----------
    mzml_file : obj
        File with mzml extension
    mgf_file : obj
        File with mgf extension
    
    Attributes
    ----------
    mzml_file : obj
        File with mzml extension
    mgf_file : obj
        File with mgf extension
    """
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
            if spectrum.getMSLevel() == 1:
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

       
    param = {}
    
    a = PeakPicking(
                src = "/home/igor/Documents/Scripts/Data/Raw_data_STARD10/.+\.mzML$",
                dst = "/home/igor/Documents/Scripts/Data/Picked_STARD10/",
                #reference_file = "/home/igor/Documents/Scripts/Data/feature_STARD10/A09_pos_picked_feature.featureXML",
                suffix_dst_files = "_picked",          #for example : "_feature"
                ext_dst_files = "mzML",           #the string may begin with a dot
                **param)
    a.main()
    """
    get_list_full_names(src = "/home/igor/Documents/Scripts/Data/Test/.+\.featureXML$")
        #dst = "/home/igor/Documents/Scripts/Data/Test/Glob")
    
    
    convert_src_to_dst_file_name(
        src = "/home/igor/Documents/Scripts/Data/Test/.+\.featureXML$",
        dst = "/home/igor/Documents/Scripts/Data/Test/Glob",
        suffix_dst_files = "_glob",
        ext_dst_files = "featureXML"
        )
    """
