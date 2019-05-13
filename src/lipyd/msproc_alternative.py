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
import numpy as np

import pyopenms as oms

import lipyd.session as session
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
    
    src_full_name_list = []                # list of src full file name;
    src_dir = os.path.dirname(src)         # get dir name from src file name;
    src_dir = src_dir if src_dir else os.getcwd() # if src dir name is empty
                                                  # get current dir name;
    pattern = os.path.basename(src)        # get file name as pattern;
    pattern = pattern if pattern != "" else ".*" # if pattern is empty
                                                 # set any char pattern;
    
    for file_name in os.listdir(src_dir): # for all file name in src dir:
        full_file_name = os.path.join(src_dir, file_name) # to build full
                                                          # file name;
        if os.path.isfile( full_file_name ): # only for files,
                                             # except dir name;
            match = None    # results re match

            try:            # try to compile patern:
                rec = re.compile(pattern)
                match = rec.match(file_name)    # apply pattern to file name
            except re.error as e:
                raise RuntimeError("Match pattern error.") # raise exeption
                                                           # if regex does
                                                           # not match
            
            if match:       # if result of re match is ok
                src_full_name_list.append(full_file_name)

    return src_full_name_list


class BaseEntity(session.Logger):
    """
    Base class for setup and check parameters for all derived classes.
    
    Parameters
    ----------
    entity : object
        Pyopenms object having setParameters method.

    Attributes
    ----------
    entity : object
        Pyopenms object having setParameters method.
    """
    
    
    def __init__(self,
                entity,
                **kwargs
                ):
        
        session.Logger.__init__(self, name = 'msproc')
        
        if not hasattr(entity, 'getDefaults'):
            
            raise RuntimeError(
                'BaseEntity: the entity has no getDefaults attr'
            )
        
        self.entity = entity
        #TODO: get default parameters from settings
        self.set_parameters(**kwargs)
        #TODO: check if no parameter is out of acceptable range
    
    
    def set_parameters(self, **kwargs):
        """
        Check type of value: if it`s str - convert to binary string
        """
        
        param = self.entity.getDefaults() # get default param;
        
        kwargs = common.dict_ensure_bytes(kwargs)
        
        for name, value in iteritems(kwargs):
            
            param.setValue(name, value)
            
            self._log(
                'Parameter `{}` of `{}` set to `{}`.'.format(
                    common.ensure_unicode(name),
                    self.entity.__class__.__name__,
                    common.ensure_unicode(value),
                ),
                1,
            )
        
        self.entity.setParameters(param)


class PPEntity(BaseEntity):
    """
    Wrapper around ``oms.PeakPickingHiRes``.
    """
    
    
    def __init__(self, **kwargs):
        
        super(PPEntity, self).__init__(
            oms.PeakPickerHiRes(),
            **kwargs,
        )


class MtdEntity(BaseEntity):
    """
    Wrapper around ``oms.MassTraceDetection``.
    """
    
    
    def __init__(self, **kwargs):
        
        super(MtdEntity, self).__init__(
            oms.MassTraceDetection(),
            **kwargs,
        )


class EpdEntity(BaseEntity):
    """
    Wrapper around ``oms.ElutionPeakDetection``.
    """
    
    
    def __init__(self, **kwargs):
        
        super(EpdEntity, self).__init__(
            oms.ElutionPeakDetection(),
            **kwargs,
        )


class FfmEntity(BaseEntity):
    """
    Wrapper around ``oms.FeatureFindingMetabo()``.
    """
    
    
    def __init__(self, **kwargs):
        
        super(FfmEntity, self).__init__(
            oms.FeatureFindingMetabo(),
            **kwargs,
        )


class MAEntity(BaseEntity):
    """
    Wrapper around ``oms.MapAlignmentAlgorithmPoseClustering``.
    """
    
    
    def __init__(self, **kwargs):
        
        super(MAEntity, self).__init__(
            oms.MapAlignmentAlgorithmPoseClustering(),
            **kwargs,
        )

class PeakPicking(session.Logger):
    """
    Class for peak picking implementation.
    
    This class implements a fast peak-picking algorithm best suited for high
    resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the
    signals of ions with similar mass-to-charge ratios (m/z) exhibit little
    or no overlapping and therefore allow for a clear separation. Furthermore,
    ion signals tend to show well-defined peak shapes with narrow peak width.

    This peak-picking algorithm detects ion signals in profile data and
    reconstructs the corresponding peak shape by cubic spline interpolation.
    Signal detection depends on the signal-to-noise ratio which is adjustable
    by the user (see parameter signal_to_noise). A picked peak's m/z and
    intensity value is given by the maximum of the underlying peak spline.

    So far, this peak picker was mainly tested on high resolution data. With
    appropriate preprocessing steps (e.g. noise reduction and baseline
    subtraction), it might be also applied to low resolution data.

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
                **kwargs
            ):

        session.Logger.__init__(self, name = 'peak_picking')

        if not src:
            raise RuntimeError( "You don`t specify all necessary files" )

        self.src = src
        self.dst = dst
        #self.log = log

        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs

        self.init_entity(**self.kw)
  

    def init_entity(self, **kwargs):

        self.pp = PPEntity(**kwargs)


    def main(self):
        #after path_parsing method we have self.src_full_name_list
        print("Peak Picking implementation")

        for f in get_list_full_names(self.src):

            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)

            print("source file:", f)
            
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
            
            print("Picked data stored into:", dst_full_file_name)


class FeatureFindingMetabo(session.Logger):
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
    
    
    def __init__(
            self,
            src = ".+\.mzML$",
            dst = None,
            suffix_dst_files = "_feature",
            ext_dst_files = "featureXML",
            **kwargs,
        ):
        
        session.Logger.__init__(self, name = 'feature_finder_metabo')
        
        if not src:
            raise RuntimeError( "You don`t specify all necessary files" )

        self.src = src
        self.dst = dst
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.kw = kwargs
        # set default value to self.param_*
        # if user will not call set_param_*
        self.param_mtd = self.kw
        self.param_epd = self.kw
        self.param_ffm = self.kw
        
        self.init_entity(**self.kw)
    

    def init_entity(self, **kwargs):
        self.mtd = MtdEntity(**kwargs)
        self.epd = EpdEntity(**kwargs)
        self.ffm = FfmEntity(**kwargs)
        self.output_mt = []
        self.splitted_mt = []
        self.filtered_mt = []
        self.chromatograms = [[]]   

    def set_param_mtd(self, **kwargs):
        self.param_mtd = kwargs

    def set_param_epd(self, **kwargs):
        self.param_epd = kwargs
    
    def set_param_ffm(self, **kwargs):
        self.param_ffm = kwargs   
    
    
    def main(self):
        #after path_parsing method we have self.src_full_name_list
        print("FeatureFindingMetabo implementation")
        
        for f in get_list_full_names(self.src):

            print("Source file:", f)
            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)
            
            input_map = oms.PeakMap() # the 1st step: load map;
            fm = oms.FeatureMap()
            oms.MzMLFile().load(f, input_map)

            self.mtd.set_parameters( **self.param_mtd )
            self.epd.set_parameters( **self.param_epd )
            self.ffm.set_parameters( **self.param_ffm )

            # the 2nd step: apply_ffm;
            self.mtd.entity.run(input_map, self.output_mt)
            self.epd.entity.detectPeaks(self.output_mt, self.splitted_mt)
            self.ffm.entity.run(self.splitted_mt, fm, self.filtered_mt)
            # the 3d step: is store result into file;
            dst_full_file_name = os.path.join(self.dst,\
                convert_src_to_dst_file_name(f,
                                            self.dst,
                                            self.suffix_dst_files,
                                            self.ext_dst_files) )
           
            oms.FeatureXMLFile().store(dst_full_file_name, fm)
            
            print("Centroided data stored into:", dst_full_file_name)


class MapAlignment(session.Logger):
    """
    Class for map alignment process.

    A map alignment algorithm based on pose clustering.

    Pose clustering analyzes pair distances to find the most probable
    transformation of retention times.
    The algorithm chooses the x most intensity peaks/features per map.
    This is modeled via the parameter 'max_num_peaks_considered',
    which in turn influences the runtime and stability of the results.
    Bigger values prolong computation, smaller values might lead to no or
    unstable trafos. Set to -1 to use all features (might take very long for
    large maps).

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
                dst_full_file_name = None,
                reference_map = None,
                toAlign_map = None,
                **kwargs
                ):
        
        session.Logger.__init__(self, name = 'map_alignment')

        if not (src and reference_file):
            
            raise RuntimeError( "You don`t specify all necessary files" )
        
        self.src = src
        self.dst = dst
        self.suffix_dst_files = suffix_dst_files
        self.ext_dst_files = ext_dst_files
        self.dst_full_file_name = dst_full_file_name
        self.kw = kwargs

        self.reference_map = reference_map
        self.toAlign_map = toAlign_map
        self.reference_file = reference_file
        self.init_entity(**self.kw)
        

    def init_entity(self, **kwargs):

        self.ma = MAEntity(**kwargs)


    def main(self): 
        #after path_parsing method we have self.src_full_name_list
        
        for f in get_list_full_names(self.src):
            print("Map Alignment implementation")
            print("Source file:", f)
            # to prepare(init) empty list and entity;
            self.init_entity(**self.kw)

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
            self.dst_full_file_name = os.path.join(self.dst,\
                convert_src_to_dst_file_name(f,
                                            self.dst,
                                            self.suffix_dst_files,
                                            self.ext_dst_files) )
            
            #print("dst=",dst_full_file_name)
            oms.FeatureXMLFile().store(self.dst_full_file_name, self.toAlign_map)
            oms.FeatureXMLFile().store(self.dst_full_file_name, self.reference_map)

            print("Aligned data stored into:", self.dst_full_file_name)
        

    

class Convert2mgf(session.Logger):
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
    
    
    def __init__(
            self,
            mzml_file = None,
            mgf_file = None,
        ):
        
        session.Logger.__init__(self, name = 'mgf_export')
        
        self.mzml_file = mzml_file
        self.mgf_file = mgf_file
    
    
    def convert(self):
        """
        Generates MGF format MS2 spectra and writes them into the output file.
        """
        
        file = oms.MzMLFile()
        msdata = oms.MSExperiment()
        file.load(self.mzml_file, msdata)
        
        outfile = open(self.mgf_file, "w")
        
        # Create header
        outfile.write("COM=Testfile\n")
        outfile.write("ITOL=1\n")
        outfile.write("ITOLU=Da\n")
        outfile.write("CLE=Trypsin\n")
        outfile.write("CHARGE=1,2,3\n")
        
        # Iterate through all spectra,
        # skip all MS1 spectra and then write mgf format
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
            
            self._log(
                'Could not find any MS2 spectra in the input, '
                'thus the output MGF file is empty!',
                -1,
            )
        
        outfile.close()


class Preprocessing(session.Logger):
    """
    
    Constructor class for all preprocessing stages implementation
    
    Usage:

    Firstly, you need to define variable as preproc = Preprocessing()
    Secondly, you have to specify desirable parameters according to pyopenms
    library. All parameters you can see on https://abibuilder.informatik.uni-tuebingen.de/
    /archive/openms/Documentation/release/2.4.0/html/index.html or TOPP application
    in vocabulary view, e.g param_pp = {"signal_to_noise": 1.0}
    Note: Boolean value should be in quoted, e.g. param_epd = {"enabled": "true"}
    
    param_pp - variable for Peak Picking parameters
    param_ff_com - variable for common Feature Finding Metabo parameters
    param_mtd - variable for Mass Trace Detection (1st step of FFM) parameters
    param_epd - variable for Elution Peak Detection (2nd step) parameters
    param_ffm - variable for Feature Finder (3rd step) parameters
    param_ma - variable for Map Alignment process
    Thirdly you have to call next methods consistently:
    preproc.peak_picking(src = "/your source directory/",
                          dst = "/destination directory",
                          suffix_dst_files = "",
                          ext_dst_files = "mzML" or "featureXML")
    prerpoc.feature_finding_metabo (the same arguments)
    preproc.map_alignment (the same srguments)
    preproc.run()
    
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
    
    """

    def __init__(
            self,
            src = None,
            dst = None,
            mzs = None,
            intensities = None,
            rts = None
        ):

        session.Logger.__init__(self, name = 'Preprocessing')

        self.pp = None
        self.ff = None
        self.ma = None

        self.src = src
        self.dst = dst

        self.mzs = mzs
        self.intensities = intensities
        self.rts = rts
        
    def peak_picking(self,
                src = None,
                dst = None,
                suffix_dst_files = "_picked",
                ext_dst_files = "mzML",
                **param):
        
        self.pp = PeakPicking(
                src = src,
                dst = dst,
                suffix_dst_files = "_picked",          #for example : "_feature"
                ext_dst_files = "mzML",           #the string may begin with a dot
                **param)

    def feature_finding_metabo(self,
                src = None,
                dst = None,
                suffix_dst_files = "_feature",
                ext_dst_files = "featureXML",
                **param):
        
        self.ff = FeatureFindingMetabo(
                src = src,
                dst = dst,
                suffix_dst_files = suffix_dst_files, #for example : "_feature"
                ext_dst_files = ext_dst_files,       #the string may begin with a dot
                **param)

    def map_alignment(self,
                src = None,
                dst = None,
                suffix_dst_files = "",
                ext_dst_files = "featureXML",
                reference_file = None,
                **param):
        
        self.ma = MapAlignment(
                src = src,
                dst = dst,
                reference_file = reference_file,
                suffix_dst_files = suffix_dst_files, #for example : "_feature"
                ext_dst_files = ext_dst_files,       #the string may begin with a dot
                **param)


    def setup_pp_params(self, **param):
        self.pp.pp.set_parameters(**param)


    def setup_ff_mtd_params(self, **param):
        self.ff.set_param_mtd(**param)

    def setup_ff_epd_params(self, **param):
        self.ff.set_param_epd(**param)

    def setup_ff_ffm_params(self, **param):
        self.ff.set_param_ffm(**param)

    def setup_ma_params(self, **param ):
        self.ma.ma.set_parameters(**param)

    def run(self):
        
        self.pp.main()
        self.ff.main()
        self.ma.main()


    def get_sampleset(self, src):
        """
        Methods for extracting mzs, rts, intensities from all files
        as 2 dimensional arrays

        """

        if not(src or self.src):
            raise RuntimeError("you have to point src pattern.")
        
        src_pattern = src if src else self.src
        
        self.mzs = []
        self.intensities = []
        self.rts = []
        
        for f in get_list_full_names(src_pattern):
            
            xml_file = oms.FeatureXMLFile()
            fmap = oms.FeatureMap()
            xml_file.load(f, fmap)
            
            print("get_xml_data f=", f)
            
            rts_tmp = []
            mzs_tmp = []
            intensities_tmp = []
            
            for n in fmap:
                _rt = n.getRT()
                _mz = n.getMZ()
                _intensities = n.getIntensity()
                rts_tmp.append(_rt)
                mzs_tmp.append(_mz)
                intensities_tmp.append(_intensities)
            
            self.mzs.append(mzs_tmp)
            self.intensities.append(intensities_tmp)
            self.rts.append(rts_tmp)

        #print(rts. mzs. intensities)

        return {
            'mzs': self.mzs,
            'rt_means': self.rts,
            'intensities': self.intensities
        }  

    def get_sampleset_2(self, src):  
        """Another methods for extracting data"""  
        
        if not (src or self.src):
            raise RuntimeError("you have to point src pattern.")
        
        #override path pattern name to src;
        src_pattern = src if src else self.src

        self.mzs = []
        self.intensities = []
        self.rts = []
        
        for f in get_list_full_names(src_pattern):
            
            #print("get_data.f = ", f)
            mzml_file = oms.MzMLFile()
            exp = oms.MSExperiment()
            mzml_file.load(f, exp)  # f is current file name from list of names;
            
            rt_for_file = []        # rt list for this file;
            mzs_for_rt = []         # mz list of list for this rt;
            intensities_for_rt = []        # ints list of list for this rt;
            
            for spec in exp:        # spec is MSSpectrum type object;
                
                mzs_tmp = []
                intensities_tmp = []
                
                _rt = spec.getRT()
                rt_for_file.append(_rt)
                
                for p in spec:
                    _mz = p.getMZ()
                    _ints = p.getIntensity()
                    mzs_tmp.append(_mz)
                    intensities_tmp.append(_ints)
                    #print("get_data: _rt={}, _mz={}, _ints={} ".format(_rt, _mz, _ints) )

                mzs_for_rt.append(mzs_tmp)
                ints_for_rt.append(intensities_tmp)
            
            self.rts.append(rt_for_file) # add rt list for current file to global rt list;
            self.mzs.append(mzs_for_rt)  # same for mz;
            self.ints.append(ints_for_rt)# same for ints;
        
        return  {
            'mzs': self.mzs,
            'rt_means': self.rts,
            'intensities': self.intensities
        } 



if __name__ == "__main__":
    
    
    a = Preprocessing() 

    param_pp = {"signal_to_noise": 1.0}

    param_ff_com = {"noise_threshold_int": 9.0,
                "chrom_peak_snr": 2.0,
                "chrom_fwhm": 4.0}
    
    param_mtd = {"mass_error_ppm" : 19.0,
                "reestimate_mt_sd": "true",
                "quant_method": "area"}
    
    
    param_epd = {"enabled": "true"}
    
    
    param_ffm = {"mz_scoring_13C": "true",
                "report_convex_hulls": "true",
                "remove_single_traces": "true"}
    
    #We need to apply and update common parameters for each 3 processes in Feature Finding
    param_mtd.update(param_ff_com)
    param_epd.update(param_ff_com)
    param_ffm.update(param_ff_com)

    param_ma = {"max_num_peaks_considered": 990}
    
    
    a.peak_picking(src = "/home/igor/Documents/Black_scripts/Src_data/.+\.mzML$",
                                dst = "/home/igor/Documents/Black_scripts/Src_data/picked/",
                                suffix_dst_files = "_picked",
                                ext_dst_files = "mzML",
                                )
    

    a.feature_finding_metabo(
                                src = "/home/igor/Documents/Black_scripts/Src_data/picked/.+\.mzML$",
                                dst = "/home/igor/Documents/Black_scripts/Src_data/featured/",
                                suffix_dst_files = "_feature",
>>>>>>> 3096b0cbf61c140290233c85b3345a7e90e93f2b
                                ext_dst_files = "featureXML",
                                )
    
    
    a.map_alignment(src = "/home/igor/Documents/Black_scripts/Src_data/featured/.+\.featureXML$",
                                dst = "/home/igor/Documents/Black_scripts/Src_data/aligned",
                                suffix_dst_files = "_aligned",
                                ext_dst_files = "featureXML",
                                reference_file = "/home/igor/Documents/Black_scripts/Src_data/featured/150310_Popeye_MLH_AC_STARD10_A09_pos_picked_feature.featureXML",
                                ) 
    

    a.setup_pp_params(**param_pp)
    
    a.setup_ff_mtd_params(**param_mtd)
    a.setup_ff_epd_params(**param_epd)
    a.setup_ff_ffm_params(**param_ffm)
    a.setup_ma_params()
    
  
    
    a.run()
    
    p = "/home/igor/Documents/Black_scripts/Src_data/aligned/.+\.featureXML$"
    #a.get_data( p )
    a.get_sampleset(p)
