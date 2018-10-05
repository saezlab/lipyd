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


import imp
import warnings
import itertools
import numpy as np


from lipyd.common import basestring


class FeatureBase(object):
    
    def __init__(self, sorter = None, **kwargs):
        """
        Serves as a base class for various classes handling arrays of
        features. Some of its predecessors represent one sample others
        more than one, or maybe data or annotations but all represent
        a series of features detected in one LC MS/MS run or more than
        one runs aligned with each other.
        """
        
        self.sorted_by = None
        
        for attr, data in iteritems(kwargs):
            
            self._add_var(data, attr)
            
        self.sorter = (
            FeatureIdx(len(self))
            if (
                sorter is None and
                self.__class__.__name__ != 'FeatureIdx'
            ) else
            sorter
        )
        
        if self.sorter is not None:
            
            self.sorter.register(self)
    
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
                        'You requested sort by array of %u dimensions but'
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


class FeatureAttributes(FeatureBase):
    
    def __init__(
            self,
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
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def __len__(self):
        
        return self.data.shape[0]
    
    def charges(self):
        """
        Returns a set of ion charges observed in the sample(set).
        """
        
        return sorted(set(self.charge))


class Sample(FeatureBase):
    
    def __init__(
            self,
            mzs,
            intensities = None,
            rts = None,
            attrs = None,
            feature_attrs = None,
            sorter = None,
        ):
        """
        Represents one LC MS/MS run.
        Has at least a vector of m/z's, optionally vector of intensites,
        retention times and other metadata.
        """
        
        self.var     = set()
        self.missing = set()
        
        if mzs is None:
            
            raise ValueError('Sample object must have at least m/z values.')
        
        self._add_var(mzs, 'mzs')
        
        FeatureBase.__init__(
            self, intensities = intensities, rts = rts, sorter = sorter
        )
        
        self.attrs   = attrs or {}
        self.feattrs = feature_attrs
    
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
        
        pass
    
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
