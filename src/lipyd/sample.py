#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems
from past.builtins import xrange, range


import os
import re
import imp
import warnings
import itertools
import operator
import numpy as np


import lipyd.common as common
basestring = common.basestring
from lipyd import reader
import lipyd.reader.peaks
import lipyd.moldb as moldb
import lipyd.ms2 as ms2
import lipyd.mgf as mgf
import lipyd.settings as settings
import lipyd.progress as progress


remgf  = re.compile(r'(\w+)_(pos|neg)_([A-Z])([0-9]{1,2})\.mgf')
remgf2 = re.compile(r'(\w+)_([A-Z])([0-9]{1,2})_(pos|neg)\.mgf')


class SampleReader(object):
    
    reader_classes = {
        'peaks': reader.peaks.PeaksReader,
    }
    
    def __init__(self, input_type, ionmode = None, **kwargs):
        """
        Reads data from files and creates ``Sample``, ``SampleSet`` and
        ``FeatureAttributes`` objects.
        
        :param str input_type:
            Format of the input files. Possible values: ``peaks``, ``mzml``.
            Reader for ``mzml`` not yet available.
        :param str ionmode:
            Ion mode of the experiment: ``pos`` or ``neg``.
            By deafult is ``None`` because the reading is possible without
            the reader being aware of the ion mode. However it is most
            convenient to provide it here.
        :param **kwargs:
            Arguments for the reader. Depends on the input format, please
            refer to classes in ``lipyd.reader`` modules.
        """
        
        if input_type not in self.reader_classes:
            
            raise RuntimeError('Unknown input type: %s' % input_type)
        
        self.ionmode = ionmode
        self.reader_class = self.reader_classes[input_type]
        self.reader_args  = kwargs
        
        self.read()
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def read(self):
        
        self.reader = self.reader_class(**self.reader_args)
    
    def get_attributes(self):
        """
        Returns ``lipyd.sample.FeatureAttributes`` object.
        This object contains variables describing series of features
        across all samples. E.g. Quality, significance, mean RT,
        centroid m/z, etc
        """
        
        attrs = self.reader.get_attributes()
        attrs['ionmode'] = self.ionmode
        
        return FeatureAttributes(**attrs)
    
    def get_samples(self, bind = True, sample_args = None):
        """
        Yields ``lipyd.sample.Sample`` objects for each sample read.
        
        To extract all data from ``PeaksReader`` the ``get_sampleset`` method
        is more convenient.
        
        :param bool bind:
            Bind samples to each other. This way they will sort together, i.e.
            if any of them is sorted all the others follow the same order.
        """
        
        feature_attrs = self.get_attributes()
        
        sorter = (
            (feature_attrs.sorter or FeatureIdx(len(self.reader.mzs)))
                if bind else
            None
        )
        
        for _sample_args in self.reader.get_samples():
            
            if not bind:
                
                feature_attrs = self.get_attributes()
            
            _sample_args['sorter']        = feature_attrs.sorter
            _sample_args['feature_attrs'] = feature_attrs
            _sample_args['ionmode']       = self.ionmode
            
            if sample_args:
                
                _sample_args.update(sample_args)
            
            yield Sample(**_sample_args)
    
    def get_sampleset(self, sampleset_args = None):
        """
        Returns a ``SampleSet`` and a ``FeatureAttributes`` object.
        """
        
        _sampleset_args = self.reader.get_sampleset()
        feature_attrs  = self.get_attributes()
        _sampleset_args['feature_attrs'] = feature_attrs
        _sampleset_args['sorter']        = feature_attrs.sorter
        _sampleset_args['ionmode']       = self.ionmode or self.reader.ionmode
        
        if sampleset_args:
            
            _sampleset_args.update(sampleset_args)
        
        return SampleSet(**_sampleset_args)


class FeatureBase(object):
    
    def __init__(self, sorter = None, **kwargs):
        """
        Serves as a base class for various classes handling arrays of
        features. Some of its predecessors represent one sample others
        more than one, or maybe data or annotations but all represent
        a series of features detected in one LC MS/MS run or more than
        one runs aligned with each other.
        """
        
        self.var = self.var if hasattr(self, 'var') else set()
        self.sorted_by = None
        
        for attr, data in iteritems(kwargs):
            
            self._add_var(data, attr)
            
        self.sorter = (
            FeatureIdx(len(self))
            if (
                sorter is None and
                self.__class__.__name__ != 'FeatureIdx' and
                len(self) is not None
            ) else
            sorter
        )
        
        if self.sorter is not None:
            
            self.sorter.register(self)
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def _add_var(self, data, attr):
        """
        Registers a variable (array of data). If an array added this way
        it will be always sorted the same way as all the other arrays in
        order to keep the data of features together.
        This method should not be called by users it is to be called
        automatically at creation of the object.
        """
        
        setattr(self, attr, data)
        
        if data is None:
            
            self.missing.add(attr)
            
        else:
            
            self.var.add(attr)
    
    def sort_all(
            self,
            by = None,
            desc = False,
            resort = False,
            propagate = True,
            return_isort = False,
            indices = (),
        ):
        """
        Sorts all data arrays according to an index array or values in one of
        the data arrays.
        
        :param np.ndarray,str by:
            Either an array of indices or the attribute name of any data
            array of the object.
        :param bool propagate:
            Whether to propagate the sorting to bound objects. This is to
            avoid loops of sorting.
        :param bool return_isort:
            Return the argsort vector. A vector of indices which can be used
            to apply the same sorting on other arrays.
        """
        
        isort = None
        
        if isinstance(by, basestring):
            
            if (
                (by == self.sorted_by and not resort) or
                not hasattr(self, by) or
                getattr(self, by) is None
            ):
                
                return None
            
            byarray = getattr(self, by)
            
            if len(byarray.shape) > 1:
                
                if len(indices) < len(byarray.shape) - 1:
                    
                    warnings.warn(
                        'You requested sort by array of %u dimensions but '
                        'selected only %u indices.\n'
                        'Selecting 0 for all other axes by default.' % (
                            len(byarray.shape),
                            len(indices),
                        )
                    )
                    
                    indices = tuple(
                        itertools.chain(
                            indices,
                            ((0,) * (len(byarray.shape) - len(indices)))
                        )
                    )
                
                if len(indices) > len(byarray.shape) - 1:
                    
                    warnings.warn(
                        'You requested sort by array of %u dimensions but'
                        'selected %u indices.\n'
                        'Dropping indices from the end.' % (
                            len(byarray.shape),
                            len(indices),
                        )
                    )
                    
                    indices = indices[
                        :(len(byarray.shape) -  len(indices) - 1)
                    ]
                
                for i in indices:
                    
                    byarray = byarray[:,i]
            
            isort = byarray.argsort()
            
            if desc:
                isort = isort[::-1]
            
            self.sorted_by = by
            
        if isinstance(by, np.ndarray):
            
            if len(by) != len(self):
                
                raise RuntimeError(
                    'Can not sort object of length %u by '
                    'index array of length %u.' (len(self), len(by))
                )
            
            isort = by
            
            self.sorted_by = None
        
        if isort is not None:
            
            for var in self.var:
                
                setattr(self, var, getattr(self, var)[isort])
            
            if propagate and self.sorter is not None:
                
                self.sorter.sort_all(by = isort, origin = id(self))
    
    def __len__(self):
        
        return (
            getattr(self, next(self.var.__iter__())).shape[0]
            if self.var else
            0
        )
    
    def filter(self, idx, negative = False, propagate = True):
        """
        Filters features at indices in ``idx``.
        
        :param list,np.array idx:
            A list or array of inidices along axis 0.
        :param bool negative:
            Drop the features at indices instead of keeping them and removing
            all others.
        :param bool propagate:
            Perform filtering on all object bound to the same sorter.
            This is to avoid propagation in loops.
        """
        
        selection = self._idx_to_selection(idx, negative = negative)
        
        self._filter(selection, propagate = propagate)
    
    @staticmethod
    def _idx_to_selection(idx, negative = False):
        
        if isinstance(idx, list):
            
            idx = np.array(idx)
        
        selection = np.in1d(np.arange(len(self)), idx)
        
        if negative:
            
            selection = np.logical_not(selection)
        
        return selection
    
    def _filter(self, selection, propagate = True):
        
        for var in self.var:
            
            setattr(self, var, getattr(self, var)[selection])
        
        if propagate:
            
            self.sorter._filter(selection = selection, origin = id(self))
    
    def threshold_filter(self, var, threshold, op = operator.gt):
        """
        Filters the features by any numeric attribute by simply applying
        a threshold. If the attribute does not exist silently does nothing.
        
        :param str var:
            Name of the attribute.
        :param float threshold:
            Value of the threshold.
        :param operator op:
            The operator to use. By default ``operator.gt`` (greater than),
            which means the features having attribute value greatre than
            the threshold will be kept and all others removed.
        """
        
        if var not in self.var:
            
            return
        
        selection = op(getattr(self, var), threshold)
        
        self._filter(selection)


class FeatureAttributes(FeatureBase):
    
    def __init__(
            self,
            ionmode = None,
            attrs = None,
            **kwargs
        ):
        """
        Features are m/z values accompanied by other attributes.
        Features might have been detected in one sample (one LC MS/MS run),
        or across multiple runs (``SampleSet``).
        ``FeatureAttributes`` handles the attributes for either a ``Sample``
        or a ``SampleSet``.
        """
        
        FeatureBase.__init__(self, **kwargs)
        
        self.attrs = attrs or {}
        self.ionmode = ionmode
    
    def __len__(self):
        
        var = next(iter(self.var))
        
        if var:
            
            return len(getattr(self, var))
    
    def charges(self):
        """
        Returns a set of ion charges observed in the sample(set).
        """
        
        return sorted(set(self.charge))


class Sample(FeatureBase):
    
    ms2_collection_methods = {
        'mgf':  'collect_mgf',
        'mzml': 'collect_mzml',
    }
    
    def __init__(
            self,
            mzs,
            ionmode = None,
            intensities = None,
            rts = None,
            sample_id = None,
            attrs = None,
            feature_attrs = None,
            sorter = None,
            ms2_format = 'mgf',
            ms2_param = None,
            silent = False,
            **kwargs,
        ):
        """
        Represents one LC MS/MS run.
        Has at least a vector of m/z's, optionally vector of intensites,
        retention times and other metadata.
        
        :param str,tuple,callable sample_id:
            An identifier for the sample or a callable method which takes
            the sample attributes as an argument and returns such an ID.
            If ``None`` the default method will be applied which attempts
            to use the fraction ID from the attributes. If fails it uses
            a random string as a unique identifier for the sample.
        :param str ms2_format:
            Format of the MS2 data. At the moment ``mgf`` is accepted.
        :param dict ms2_param:
            - ``mgf_files``:
            A list of MGF file paths.
            - ``mgf_match_method`` (callable):
            A method which finds the MGF files. It accepts 2 arguments:
            the directory path with the MGF files and the attributes used
            to match the file names.
            - ``mgfdir`` (path):
            Path to the directory containing the MGF files.
            - ``mgf_charge``:
            If a value provided here only the scans with this precursor
            charge will be processed from the MGF files. If ``None``, charge
            won't be considered.'
        """
        
        self.var     = set()
        self.missing = set()
        
        if mzs is None:
            
            raise RuntimeError('Sample object must have at least m/z values.')
        
        self._add_var(mzs, 'mzs')
        
        FeatureBase.__init__(
            self,
            intensities = intensities,
            rts = rts,
            sorter = sorter,
            **kwargs,
        )
        
        self.silent     = silent
        self.attrs      = attrs or {}
        self.ionmode    = ionmode
        self.feattrs    = feature_attrs
        self.ms2_format = ms2_format
        self.ms2_param  = ms2_param or {}
        self._sample_id = sample_id
        self._set_sample_id()
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def _add_var(self, data, attr):
        
        # mzs is the first variable so it is added without
        # checking its length. All others must be the same
        # length as mzs.
        if attr != 'mzs' and data is not None and len(data) != len(self):
            
            raise RuntimeError(
                'Sample object: length of `%s` (%u) is not the same '
                'as number of features in the sample (%u).' % (
                    attr, len(data), len(self),
                )
            )
        
        FeatureBase._add_var(self, data, attr)
    
    def sort_all(
            self,
            by = 'mzs',
            desc = False,
            resort = False,
            propagate = True,
            return_isort = False,
        ):
        """
        Sorts all data arrays according to an index array or values in one of
        the data arrays.
        
        :param np.ndarray,str by:
            Either an array of indices or the attribute name of any data
            array of the object.
        :param bool desc:
            Sort descending. Applies only if ``by`` is the attribute name
            of a data array.
        :param bool resort:
            Sort even if the object is already sorted according to the
            value stored in ``sorted_by``.
        :param bool propagate:
            Whether to propagate the sorting to bound objects. This is to
            avoid loops of sorting.
        :param bool return_isort:
            Return the argsort vector. A vector of indices which can be used
            to apply the same sorting on other arrays.
        """
        
        isort = FeatureBase.sort_all(
            self,
            by = by,
            desc = desc,
            resort = resort,
            propagate = propagate,
            return_isort = return_isort,
        )
        
        if return_isort:
            
            return isort
    
    def _set_sample_id(self):
        
        self.sample_id = self._get_sample_id(self._sample_id, self.attrs)
    
    @classmethod
    def _get_sample_id(cls, sample_id = None, attrs = None):
        
        if sample_id is None:
            
            return cls._default_sample_id_method(attrs)
            
        elif hasattr(sample_id, '__call__'):
            
            return sample_id(attrs)
            
        else:
            
            return sample_id
    
    @staticmethod
    def _default_sample_id_method(attrs):
        
        if (
            attrs is not None and
            'label' in attrs and
            'fraction' in attrs['label']
        ):
            
            return attrs['label']['fraction']
        
        return common.random_string()
    
    def __len__(self):
        """
        Returns number of MS1 ions (m/z's) detected in the sample.
        """
        
        return len(self.mzs)
    
    #
    # Methods for processing
    #
    
    def process(
            self,
            basic_filters = True,
            database_lookup = True,
            ms2_analysis = True,
            basic_filters_args = None,
            database_lookup_args = None,
            ms2_analysis_args = None,
        ):
        """
        This method implements the workflow of processing the samples.
        At the moment it does trivial filtering (see ``basic_filters``),
        performs database lookups (see ``database_lookup``) and analyses
        the MS2 data (see ``ms2_analysis``).
        """
        
        if basic_filters:
            self.basic_filters(**(basic_filters_args or {}))
        
        if database_lookup:
            self.database_lookup(**(database_lookup_args or {}))
        
        if ms2_analysis:
            self.ms2_analysis(**(ms2_analysis_args or {}))
    
    def quality_filter(self, threshold = .2):
        """
        If the features has an attribute named `quality` it removes the
        ones with quality below the threshold.
        
        """
        
        self.feattrs.threshold_filter('quality', threshold = threshold)
    
    def rt_filter(self, threshold = 1.):
        """
        Removes the features with mean retention time below the threshold.
        Works only if the features has an attribute named `rt_means`.
        """
        
        self.feattrs.threshold_filter('rt_means', threshold = threshold)
    
    def charge_filter(self, charge = 1):
        """
        Keeps only the features with the preferred charge.
        Does not care about sign.
        """
        
        self.feattrs.threshold_filter(
            'charge',
            threshold = abs(charge),
            op = lambda vc, c: np.abs(vc) == c,
        )
    
    def intensity_filter(self, threshold = 10000.):
        """
        Removes the features with total intensities across samples lower than
        the threshold.
        
        Works by the ``total_intensities`` attribute.
        """
        
        self.feattrs.threshold_filter(
            'total_intensities',
            threshold = threshold,
        )
    
    def basic_filters(
            self,
            quality_min = .2,
            charge = 1,
            rt_min = 1.,
            intensity_min = 10000.
        ):
        """
        Performs 4 trivial filtering steps which discard a large number of
        features at the beginning of the analysis. It removes the features
        with quality lower than 0.2, non single charge, retention time lower
        than 1 minute or total intensity lower than 10,000.
        Different threshold values can be provided to this method.
        """
        
        self.quality_filter(threshold = quality_min)
        self.charge_filter(charge = charge)
        self.intensity_filter(threshold = intensity_min)
        self.rt_filter(threshold = rt_min)
    
    def database_lookup(
            self,
            database_args = None,
            reinit_db = False,
            adduct_constraints = True,
            charge = None,
            tolerance = None,
        ):
        """
        Performs database lookups of all m/z's.
        
        Creates an array variable named ``records`` in the
        ``FeatureAttributes`` object (``feattrs``) with the records retrieved
        from the database.
        """
        
        if not hasattr(moldb, 'db') or reinit_db:
            
            moldb.init_db(**(database_args or {}))
        
        # we add an array variable with the recodrs resulted
        # in the lookup
        self.feattrs._add_var(
            moldb.adduct_lookup_many(
                self.mzs,
                ionmode = self.ionmode,
                adduct_constraints = adduct_constraints,
                charge = charge,
                tolerance = tolerance,
            ),
            'records',
        )
    
    def set_ms2_sources(self, attrs = None):
        """
        Collects the MS2 scans belonging to this sample.
        
        The collected resources will be in the ``ms2_source`` attribute.
        """
        
        self.ms2_source = []
        
        ms2_source = self.collect_ms2(attrs = attrs)
        
        if self.sample_id in ms2_source:
            
            self.ms2_source = ms2_source[self.sample_id]
    
    def collect_ms2(self, sample_id = None, attrs = None):
        """
        Collects the MS2 scans belonging to this sample.
        
        The collected resources returned as ``dict`` which can be passed
        to ``ms2.MS2Feature``.
        """
        
        attrs = attrs or self.attrs
        sample_id = sample_id or self._get_sample_id(attrs)
        
        if self.ms2_format is None:
            
            raise RuntimeError(
                '`ms2_format` is necessary to run MS2 analysis.'
            )
        
        if self.ms2_format in self.ms2_collection_methods:
            
            method = self.ms2_collection_methods[self.ms2_format]
            
            return getattr(self, method)(sample_id = sample_id, attrs = attrs)
    
    def collect_mgf(self, sample_id = None, attrs = None):
        """
        Collects MGF files containing the MS2 spectra.
        
        :param dict attrs:
            A ``dict`` with attributes necessary to find the appropriate
            MGF files.
        """
        
        # if list of filenames provided we simply use those
        if 'mgf_files' in self.ms2_param:
            
            mgf_files = self.ms2_param['mgf_files']
            
            # if directory provided create paths
            if 'mgfdir' in self.ms2_param:
                
                mgf_files = [
                    os.path.join(self.ms2_param['mgfdir'], fname)
                    for fname in mgf_files
                ]
            
        # otherwise we need to collect the files
        else:
            
            # check if directory with MGF files provided
            mgfdir = (
                self.ms2_param['mgfdir']
                    if 'mgfdir' in self.ms2_param else
                # fall back to this default
                'mgf'
            )
            
            # if dir does not exist we can do nothing
            if not os.path.isdir(mgfdir):
                
                raise FileNotFoundError(
                    'Please provide a directory with MGF files.'
                )
            
            # ok, now we can find the MGF files by matching against
            # the sample attributes
            attrs = attrs or self.attrs
            sample_id = sample_id or self.sample_id
            
            # see if a matching method provided
            mgf_match_method = (
                self.ms2_param['mgf_match_method']
                    if 'mgf_match_method' in self.ms2_param else
                # otherwise we use our default
                self._default_mgf_match_method
            )
            
            mgf_files = []
            
            for fname in os.listdir(mgfdir):
                
                mgf_path     = os.path.join(mgfdir, fname)
                mgf_matches  = mgf_match_method(mgf_path, attrs)
                
                if mgf_matches:
                    
                    mgf_files.append(mgf_path)
        
        # check if mgf_charge provided in params
        mgf_charge = (
            self.ms2_param['mgf_charge']
                if 'mgf_charge' in self.ms2_param else
            None
        )
        
        # dict with one key and list of files probably only one element:
        # do like this to be able to directly pass to ms2.MS2Feature
        return  {
            sample_id: [
                mgf.MgfReader(fname, charge = mgf_charge)
                for fname in mgf_files
            ]
        }
    
    @staticmethod
    def _default_mgf_match_method(path, attrs):
        """
        The default method for matching names of MGF files against attributes.
        """
        
        match  = remgf.search(path)
        match2 = remgf2.search(path)
        
        if match:
            
            main, ionmode, fracrow, fraccol = match.groups()
        
        if match2:
            
            main, fracrow, fraccol, ionmode = match2.groups()
        
        return (
            (match or match2) and
            attrs['label']['main'].lower() in main.lower() and
            (fracrow, int(fraccol)) == attrs['label']['fraction'] and
            ionmode == attrs['label']['ionmode']
        )
    
    def collect_mzml(self, attrs = None):
        
        raise NotImplementedError
    
    def ms2_analysis(self, resources = None):
        """
        Runs MS2 identification methods on all features.
        """
        
        if resources is None:
            
            if isinstance(self.ms2_source, list):
                
                resources = {self.sample_id: self.ms2_source}
                
            elif isinstance(self.ms2_source, dict):
                
                resources = self.ms2_source
                
            else:
                
                raise RuntimeError(
                    'No resources provided for MS2 identification.'
                )
        
        ms2_identities = []
        
        if not self.silent:
            
            prg = progress.Progress(len(self), 'Analysing MS2 spectra', 1)
        
        for i in xrange(len(self)):
            
            if not self.silent:
                
                prg.step()
            
            # MS2 identifications:
            ms2_fe = ms2.MS2Feature(
                mz = self.mzs[i],
                ionmode = self.ionmode,
                resources = resources,
                rt = tuple(self.feattrs.rt_ranges[i]),
                ms1_records = self.feattrs.records[i],
            )
            
            ms2_fe.main()
            
            ms2_identities.append(ms2_fe.identities)
        
        if not self.silent:
            
            prg.terminate()
        
        ms2_identities = np.array(ms2_identities)
        
        self.feattrs._add_var(ms2_identities, 'ms2_identities')
    
    def ms2_identify(self):
        
        self.set_ms2_sources()
        self.ms2_analysis()
    
    #
    # Methods for exporting the results
    #
    
    def get_database_records(self, i, database = None, adduct = None):
        """
        Yields database records for one feature, optionally only records
        of certain adduct or database.
        """
        
        if adduct:
            
            records = self.feattrs.records[i][adduct][1]
            
        else:
            
            records =  itertools.chain(
                *(add[1] for add in self.feattrs.records[i].values())
            )
        
        for rec in records:
            
            if not database or database == rec.lab.db:
                
                yield rec
    
    def table(
            self,
            variables = None,
            headers = None,
        ):
        """
        Returns results as a header and a table as list of lists.
        """
        
        hdr = [
            'm/z',
            'average area',
            'RT range',
            'quality',
            'significance',
            'database lookups',
            'MS2 top',
            'MS2 all',
        ]
        
        if variables:
            
            for i, var in enumerate(variables):
                
                hdr.append(var if not headers else headers[i])
        
        yield hdr
        
        for i in xrange(len(self)):
            
            ms2_id = self.feattrs.ms2_identities[i]
            ms2_max_score = (
                max(
                    (
                        ms2i.score_pct
                        for ms2iii in ms2_id
                        for ms2ii in ms2iii.values()
                        for ms2i in ms2ii
                    ),
                    default = 0
                )
            )
            ms2_best = ';'.join(
                ms2i.full_str()
                for ms2iii in ms2_id
                for ms2ii in ms2iii.values()
                for ms2i in ms2ii
                if ms2i.score_pct == ms2_max_score
            )
            ms2_all = ';'.join(
                ms2i.full_str()
                for ms2iii in ms2_id
                for ms2ii in ms2iii.values()
                for ms2i in ms2ii
                if ms2i.score_pct > 0.0
            )
            
            line = [
                '%08f' % self.mzs[i],
                '%u' % self.feattrs.total_intensities[i],
                '%.02f - %.02f' % tuple(self.feattrs.rt_ranges[i]),
                '%.02f' % self.feattrs.quality[i],
                '%.02f' % self.feattrs.significance[i],
                moldb.records_string(
                    records = self.feattrs.records[i],
                    show_ppm = True,
                    show_adduct = True,
                    show_db = True,
                ),
                ms2_best,
                ms2_all,
            ]
            
            if variables:
                
                for var in variables:
                    
                    line.append(str(getattr(self.feattrs, var)[i]))
            
            yield line
    
    def export_table(self, fname, **kwargs):
        
        table = self.table(**kwargs)
        
        hdr = next(table)
        
        with open(fname, 'w') as fp:
            
            _ = fp.write('\t'.join(hdr))
            
            for row in table:
                
                _ = fp.write('\n')
                _ = fp.write('\t'.join(row))


class FeatureIdx(FeatureBase):
    
    def __init__(self, length):
        """
        Helps the sorting of features across multiple samples
        with keeping track of feature IDs.
        """
        
        FeatureBase.__init__(self)
        
        self.var = {'_original', '_current'}
        
        self._original = np.arange(length)
        self._current  = np.arange(length)
        
        self.clients = {}
    
    def _sort(self, argsort):
        """
        Sorts the two index arrays by an index array.
        Do not call this because it does sync with clients.
        """
        
        if len(argsort) != len(self):
            
            raise RuntimeError(
                'FeatureIdx: can not sort by index array of different length.'
            )
        
        self._current  = self._current[argsort]
        self._original = self._current.argsort()
    
    def current(self, o):
        """
        Tells the current index for the original index ``o``.
        
        :param int o:
            An original index.
        """
        
        return self._original[o]
    
    def original(self, c):
        """
        Tells the original index for the current index ``c``.
        
        :param int c:
            An index in the current ordering.
        """
        
        return self._current[c]
    
    def __len__(self):
        
        return len(self._original)
    
    def acurrent(self, ao):
        """
        For a vector of original indices ``ao`` returns a vector of
        corresponding current indices.
        
        :param int ao:
            Vector of original indices.
        """
        
        return self._original[ao]
    
    def aoriginal(self, co):
        """
        For a vector of current indices ``co`` returns a vector of
        corresponding original indices.
        
        :param int co:
            Vector of current indices.
        """
        
        return self._current[ac]
    
    def convert(self, other):
        """
        Converts from the ordering of an other ``FeatureIdx`` instance
        to the ordering of this one.
        """
        
        raise NotImplementedError
    
    def register(self, sortable):
        """
        Binds a sortable object of the same length to this one hence all
        sorting operations will be applied also to this object.
        Sortables clients assumed to share the same original indices.
        It means should have the same ordering at time of registration.
        """
        
        if len(sortable) != len(self):
            
            raise RuntimeError(
                'FeatureIdx: Objects of unequal length can not be co-sorted.'
            )
        
        self.clients[id(sortable)] = sortable
    
    def unregister(self, id_sortable):
        """
        Removes an object from the list of clients and also makes it forget
        about this object to be its sorter. Simply it makes the this object
        and the client independent.
        
        :param int,sortable id_sortable:
            Either the ``id`` of a client or the object itself.
        """
        
        if not isinstance(id_sortable, int):
            
            id_sortable = id(id_sortable)
        
        if id_sortable in self.clients:
            
            self.clients[id_sortable].sorter = None
            del self.clients[id_sortable]
    
    def sort_all(self,
            by = None,
            origin = None,
        ):
        """
        Applies custom sort to the object if ``by`` is an index array.
        If ``by`` is ``None``, it sorts by the original indices, which
        means the restoration of the original order.
        It calls ``sort_all`` by the same argsort on all registered clients
        hence propagating the sorting. But omits the client if the sort
        request originates from it in order to avoid infinite sorting loops.
        Clients don't propagate further the sorting request also to avoid
        loops.
        In summary, this method is called either from outside and then
        it applies a sorting to all object of the system (if called without
        ``by`` restores the original order); or called by the ``sort_all``
        method of one of the clients propagating the sorting to this object
        and to all other clients.
        """
        
        if by is None:
            
            by = self._current.argsort()
        
        self._sort(by)
        
        for id_, client in iteritems(self.clients):
            
            if id_ != origin:
                
                client.sort_all(by = by, propagate = False)
    
    def _filter(self, selection, origin = None):
        """
        Applies filtering to all registered clients.
        """
        
        FeatureBase._filter(self, selection, propagate = False)
        
        for id_, client in iteritems(self.clients):
            
            if id_ != origin:
                
                client._filter(selection, propagate = False)


class SampleSet(Sample):
    
    def __init__(
            self,
            mzs,
            ionmode = None,
            intensities = None,
            rts = None,
            attrs = None,
            feature_attrs = None,
            sorter = None,
            ms2_format = 'mgf',
            ms2_param = None,
            sample_ids = None,
        ):
        """
        This class represents an ordered set of samples.
        It means all arrays have one more dimension than ``Sample``, though
        the arrays can be multidimensional, i.e. in one data array more
        than one value can belong to a feature in a sample.
        It inherits from ``Sample``, see details in docs of ``Sample``.
        
        In ``ms2_param`` in addition to the parameters used in ``Sample``
        you can provide ``ms2_use_samples`` and ``ms2_use_samples_method``,
        a data and a function, respectively which will be used to select
        only certain samples for MS2 analysis. The method should accept
        sample attributes as first argument and the ``ms2_use_samples``
        or ``None`` will be passed as second. It should return ``bool``.
        By deafult MS2 scans from all samples are used.
        
        :param list,callable sample_ids:
            Either a list of sample identifiers with the same length as
            number of samples or a method which generates sample identifiers
            from sample attributes.
        """
        
        centr_mzs = (
            feature_attrs.centr_mzs.copy()
                if (
                    feature_attrs is not None and
                    hasattr(feature_attrs, 'centr_mzs')
                ) else
            mzs
                if len(mzs.shape) == 1 else
            mzs.mean(axis = tuple(range(1, len(mzs.shape))))
        )
        
        # this is the default
        self.numof_samples = 1
        
        if attrs:
            
            # if we have attributes number of samples should be obvious
            self.numof_samples = len(attrs)
            
        else:
            
            # otherwise try to guess from shape of the data arrays
            for var in (mzs, intensities, rts):
                
                if hasattr(var, 'shape') and len(var.shape) > 1:
                    
                    self.numof_samples = var.shape[1]
                    break
        
        # by default None, which means the default method will be
        # called from Sample to get IDs either from the attributes
        # or random strings in worst case
        sample_id = None
        self.sample_ids = [None] * self.numof_samples
        
        if hasattr(sample_ids, '__call__'):
            
            # if it is callable we register it as a sample ID
            # provider method
            sample_id = sample_ids
            
        elif sample_ids and len(sample_ids) == self.numof_samples:
            
            # if it is a list with same length as number of samples
            # we just use it to assign an ID to each sample
            self.sample_ids = sample_ids
        
        Sample.__init__(
            self,
            mzs = centr_mzs,
            ionmode = ionmode,
            intensities = intensities,
            mzs_by_sample = mzs,
            rts = rts,
            attrs = attrs,
            feature_attrs = feature_attrs,
            sorter = sorter,
            ms2_format = ms2_format,
            ms2_param = ms2_param,
            sample_id = sample_id,
        )
    
    @classmethod
    def combine_samples(cls, attrs, samples):
        """
        Initializes the object by combining a series of ``Sample`` objects.
        
        All samples must be of the same length and have the same ordering.
        It means the corresponding elements must belong to the same feature.
        The ``FeatureAttributes`` and the ``sorter`` (``FeatureIdx``) will
        be used from the first sample.
        """
        
        if len(set(len(s) for s in samples)) > 1:
            
            # TODO: call aligner for different length samples
            # TODO: check if at least sorters show the same ordering
            raise RuntimeError(
                'SampleSet: it\'s possible combine only equal length samples '
                'into a set. Also samples assumed to have the same ordering.'
            )
        
        mzs = np.hstack([s.mzs for s in samples])
        
        var = {}
        
        for var_name in sample[0].var:
            
            var[var_name] = np.hstack([getattr(s, var_name) for s in samples])
        
        sample_attrs = [s.attrs for s in samples]
        
        # TODO: combine FeatureAttributes
        feature_attrs = samples[0].feattrs
        sorter = samples[0].sorter
        
        return cls.__init__(
            mzs = mzs,
            attrs = sample_attrs,
            feature_attrs = feature_attrs,
            sorter = sorter,
            **var,
        )
    
    def get_sample_attrs(self, i):
        """
        Returns the sample attributes (dict of metadata) of the ``i``th
        sample in the set.
        """
        
        return self.attrs[i]
    
    def get_sample(self, i):
        """
        Returns the ``i``th sample as a ``Sample`` object.
        """
        
        var = {}
        
        for var_name in self.var:
            
            var[var_name] = getattr(self, var_name)[:,i]
        
        sample_attrs = self.get_sample_attrs(i)
        
        return Sample(
            mzs = self.mzs_by_sample[:,i],
            attrs = sample_attrs,
            feature_attrs = self.feattrs,
            sorter = self.sorter,
            **var,
        )
    
    def get_sample_id(self, i):
        
        if hasattr(self.sample_id, '__call__'):
            
            sample_id = self.sample_id
            
        else:
            
            sample_id = self.sample_ids[i]
        
        return self._get_sample_id(sample_id, self.attrs[i])
    
    def collect_ms2(self, attrs = None):
        """
        Collects the MS2 resources for each sample a similar way as the same
        method does in ``Sample``.
        """
        
        ms2_source = {}
        
        attrs = attrs or self.attrs
        
        for i, sample_attrs in enumerate(attrs):
            
            if 'ms2_use_samples_method' in self.ms2_param:
                
                select_samples = self.ms2_param['ms2_use_samples_method']
                
                ms2_use_samples = (
                    self.ms2_param['ms2_use_samples']
                        if 'ms2_use_samples' in self.ms2_param else
                    None
                )
                
                if not select_samples(sample_attrs, ms2_use_samples):
                    
                    continue
            
            sample_id = self.get_sample_id(i)
            # adding MS2 sources for this sample to the dict
            ms2_source.update(Sample.collect_ms2(
                self,
                sample_id = sample_id,
                attrs = sample_attrs
            ))
        
        return ms2_source
    
    def set_ms2_sources(self, attrs = None):
        """
        Collects the MS2 scans for all samples.
        
        The collected resources will be in the ``ms2_source`` attribute.
        """
        
        self.ms2_source = self.collect_ms2(attrs = attrs)
    
    def peak_size_filter(self, backg, prot, threshold = 2):
        
        def ps(a0, a1):
            
            if np.all(np.isnan(a1)) or np.nanmax(a1) == 0:
                
                return 0.0
            
            if (
                (np.all(np.isnan(a0)) or np.nanmax(a0) == 0) and
                np.nanmax(a1) > 0
            ):
                
                return np.inf
            
            return np.nanmin(a1) / np.nanmax(a0)
        
        backg = set((f[0], int(f[1:])) for f in backg)
        prot  = set((f[0], int(f[1:])) for f in prot)
        
        ibackg = np.array([
            i
            for i, attr in enumerate(self.attrs)
            if attr['label']['fraction'] in backg
        ])
        iprot = np.array([
            i
            for i, attr in enumerate(self.attrs)
            if attr['label']['fraction'] in prot
        ])
        
        peaksize = np.array([
            ps(self.intensities[i,ibackg], self.intensities[i,iprot])
            for i in xrange(len(self))
        ])
        
        self.feattrs._add_var(peaksize, 'peaksize')
        
        self.feattrs.threshold_filter('peaksize', threshold = threshold)
