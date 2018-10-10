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


from lipyd.common import basestring
from lipyd import reader
import lipyd.reader.peaks
import lipyd.moldb as moldb
import lipyd.ms2 as ms2


remgf = re.compile(r'([^\W_]+)_(pos|neg)_([A-Z])([0-9]{1,2})\.mgf')


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
            
            raise ValueError('Unknown input type: %s' % input_type)
        
        self.ionmode = ionmode
        self.reader_class = reader_classes[input_type]
        self.reader_args  = kwargs
        
        self.read()
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def read(self):
        
        self.reader = self.readerclass(**self.reader_args)
    
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
    
    def get_samples(self, bind = True):
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
        
        for sample_args in self.reader.get_samples():
            
            if not bind:
                
                feature_attrs = self.get_attributes()
            
            sample_args['sorter']        = feature_attrs.sorter
            sample_args['feature_attrs'] = feature_attrs
            sample_args['ionmode']       = self.ionmode
            
            yield Sample(**sample_args)
    
    def get_sampleset(self):
        """
        Returns a ``SampleSet`` and a ``FeatureAttributes`` object.
        """
        
        sampleset_args = self.reader.get_sampleset()
        sampleset_args['feature_attrs'] = self.get_attributes()
        sampleset_args['sorter']        = feature_attrs.sorter
        sampleset_args['ionmode']       = ionmode
        
        return SampleSet(**sampleset_args)


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
                
                raise ValueError(
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
            
            self.sorter.filter(idx = idx, origin = id(self))
    
    def threshold_filter(self, attr, threshold, op = operator.gt):
        """
        Filters the features by any numeric attribute by simply applying
        a threshold. If the attribute does not exist silently does nothing.
        
        :param str attr:
            Name of the attribute.
        :param float threshold:
            Value of the threshold.
        :param operator op:
            The operator to use. By default ``operator.gt`` (greater than),
            which means the features having attribute value greatre than
            the threshold will be kept and all others removed.
        """
        
        if attr not in self.var:
            
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
            **kwargs,
        ):
        """
        Represents one LC MS/MS run.
        Has at least a vector of m/z's, optionally vector of intensites,
        retention times and other metadata.
        
        :param str ms2_format:
            Format of the MS2 data. At the moment ``mgf`` is accepted.
        :param dict ms2_param:
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
            
            raise ValueError('Sample object must have at least m/z values.')
        
        self._add_var(mzs, 'mzs')
        
        FeatureBase.__init__(
            self,
            intensities = intensities,
            rts = rts,
            sorter = sorter,
            **kwargs,
        )
        
        self.attrs      = attrs or {}
        self.ionmode    = ionmode
        self.feattrs    = feature_attrs
        self.ms2_format = ms2_format
        self.ms2_param  = ms2_param or {}
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def _add_var(self, data, attr):
        
        if attr != 'mzs' and data is not None and len(data) != len(self):
            
            raise ValueError(
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
                ionmode = ionmode,
                adduct_constraints = adduct_constraints,
                charge = charge,
                tolerance = tolerance,
            ),
            'records',
        )
    
    def collect_ms2(self, attrs = None):
        """
        Collects the MS2 scans belonging to this sample.
        
        The collected resources will be in the ``ms2_source`` attribute.
        """
        
        attrs = attrs or self.attrs
        
        if self.ms2_format in self.ms2_collection_methods:
            
            self.ms2_collection_methods[self.ms2_format](attrs = attrs)
    
    def collect_mgf(self, attrs = None):
        """
        Collects MGF files containing the MS2 spectra.
        
        :param dict attrs:
            A ``dict`` with attributes necessary to find the appropriate
            MGF files.
        """
        
        mgfdir = (
            self.ms2_param['mgfdir']
                if 'mgfdir' in self.ms2_param else
            'mgf'
        )
        
        if not os.isdir(mgfdir):
            
            raise FileNotFoundError(
                'Please provide a directory with MGF files.'
            )
        
        attrs = attrs or self.attrs
        
        mgf_match_method = (
            self.ms2_param['mgf_match_method']
                if 'mgf_match_method' in self.ms2_param else
            self._default_mgf_match_method
        )
        
        mgf_files = mgf_match_method(mgfdir, attrs)
        
        mgf_charge = (
            self.ms2_param['mgf_charge']
                if 'mgf_charge' in self.ms2_param else
            None
        )
        
        self.ms2_source = [
            mgf.MgfReader(fname, charge = mgf_charge)
            for fname in mgf_files
        ]
    
    @staticmethod
    def _default_mgf_match_method(mgfdir, attrs):
        """
        The default method for matching names of MGF files against attributes.
        """
        
        result = set()
        
        for fname in os.listdir(mgfdir):
            
            remgf.search(fname)
            
            if not match:
                
                continue
            
            main, ionmode, fracrow, fraccol = match.groups()
            
            if (
                main == attrs['label']['main'] and
                (fracrow, int(fraccol)) == attrs['label']['fraction'] and
                ionmode == attrs['label']['ionmode']
            ):
                
                result.add(os.path.join(mgfdir, fname))
        
        return result
    
    def ms2_analysis(self):
        """
        Runs MS2 identification methods on all features.
        """
        
        for i in xrange(len(self)):
            
            # MS2 identifications:
            ms2_fe = ms2.MS2Feature(
                self.mzs[i],
                self.ionmode,
                self.mgf_readers[self.ionmode],
                rt = rt,
                ms1_records = ms1_records,
            )
            ms2_fe.main()
            ms2_id = ms2_fe.identity_summary(
                sample_ids = True,
                scan_ids = True,
            )
            ms2_max_score = max(i[1] for i in ms2_id) if ms2_id else 0
            ms2_best = self.identities_str(
                i for i in ms2_id if i[1] > 0 and i[1] == ms2_max_score
            )
            ms2_all = self.identities_str(ms2_id)
            
            data['ms2_new_top'].append(ms2_best)
            data['ms2_new_all'].append(ms2_all)


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
            
            raise ValueError(
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
            
            raise ValueError(
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
            
            if id_ !+ origin:
                
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
        """
        
        centr_mzs = (
            feature_attrs.centr_mzs
                if (
                    feature_attrs is not None and
                    hasattr(feature_attrs, 'centr_mzs')
                ) else
            mzs
                if len(mzs.shape) == 1 else
            mzs.mean(axis = tuple(range(1, len(mzs.shape))))
        )
        
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
            raise ValueError(
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
    
    def collect_ms2(self, attrs = None):
        """
        Collects the MS2 resources for each sample a similar way as the same
        method does in ``Sample``.
        """
        
        ms2_source = []
        
        attrs = attrs or self.attrs
        
        for sample_attrs in attrs:
            
            if 'ms2_use_samples_method' in self.ms2_param:
                
                select_samples = self.ms2_param['ms2_use_samples_method']
                
                ms2_use_samples = (
                    self.ms2_param['ms2_use_samples']
                        if 'ms2_use_samples' in self.ms2_param else
                    None
                )
                
                if not select_samples(sample_attrs, ms2_use_samples):
                    
                    continue
            
            Sample.collect_ms2(self, attrs = sample_attrs)
            
            ms2_source.append(self.ms2_source)
        
        self.ms2_source = ms2_source
    
    def ms2_analysis(self):
        """
        Runs MS2 identification methods on all features.
        """
        
        for i_fe in xrange(len(self)):
            
            
