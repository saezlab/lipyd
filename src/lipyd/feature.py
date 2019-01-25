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
from past.builtins import xrange, range

import numpy as np

import lipyd.sampleattrs as sampleattrs


class FeatureAnalyzer(object):
    """ """
    
    def __init__(self, name, samples, method, **variables):
        """
        Serves for analysis of features using data in the feature vs. sample
        data arrays in ``SampleSet`` objects and feature variables in
        their ``FeatureAttrs`` object and various background variables
        about samples in ``SampleData`` objects.
        
        ``SampleData`` variables can be provided as keyword arguments.
        ``SampleSet`` object must be provided as well as the analysis method.
        
        name : str
            A name for the feature variable created by this analysis object.
        samples : sample.SampleSet
            A ``lipyd.sample.SampleSet`` object.
        method : callable
            A method which accepts keyword arguments and explicitely
            accepts keyword arguments necessary for its operation.
            All feature variables, sampleset data and background variables
            will be provided as keyword arguments. The method called
            for each sample one by one. At the end the resulted array
            will be registered as a new variable in the ``FeatureAttrs``
            object of the ``SampleSet`` under the attribute name ``name``.
        **variables :
            Custom ``SampleData`` or derived objects. Will be provided to
            ``method`` by their argument name.
        """
        
        self.name = name
        self.samples = samples
        self.method = method
        self.variables = variables
        
        result = self.run()
        self.samples.feattrs._add_var(result, self.name)
    
    def run(self):
        """
        Applies the method for each feature and returns an array of the
        results.

        Parameters
        ----------

        Returns
        -------

        """
        
        result = []
        
        for i in xrange(len(self.samples)):
            
            featurevars = {}
            
            for var in self.samples.feattrs.var:
                
                featurevars[var] = getattr(self.samples.feattrs, var)[i]
            
            for var in self.samples.var:
                
                # safe to do this as first axis is always the features
                featurevars[var] = getattr(self.samples, var)[i]
            
            result.append(self.method(**featurevars, **self.variables))
        
        return np.array(result)


class ProfileAnalyzer(FeatureAnalyzer):
    
    def __init__(
            self,
            samples,
            protein = None,
            condition = None,
            profile_filter_args = (),
        ):
        
        self.protein = protein
        self.condition = condition
        
        feature.FeatureAnalyzer.__init__(
            self,
            name = 'profile',
            samples = samples,
            method = self.profile_method,
            _samples = samples,
            profile_filter_args = profile_filter_args,
            protein = protein,
            condition = condition,
        )
    
    @staticmethod
    def profile_method(
            intens_norm,
            _samples,
            profile_filter_args = (),
            protein = None,
            condition = None,
            **kwargs,
        ):
        
        # intensities in fraction categories
        
        # intensities of fractions in peak
        high_protein = [
            intens_norm[_samples.attrs.get_sample_index(sample_id)]
            for sample_id in profile_filter_args['peak'][0]
        ]
        
        # intensities of fractions with small or tiny protein
        some_protein = [
            intens_norm[_samples.attrs.get_sample_index(sample_id)]
            for sample_id in
            itertools.chain(
                profile_filter_args['small'][0],
                profile_filter_args['tiny'][0],
            )
        ]
        # intensities of fractions with no protein
        no_protein = [
            intens_norm[_samples.attrs.get_sample_index(sample_id)]
            for sample_id in
            profile_filter_args['none'][0]
        ]
        
        # exclude those completely missing from the peak
        
        if np.all(np.isnan(high_protein)):
            
            return False
        
        # apply single fraction operators
        
        for (
            label,
            (sample_ids, threshold, op)
        ) in iteritems(profile_filter_args):
            
            for sample_id in sample_ids:
                
                i = _samples.attrs.get_sample_index(sample_id)
                
                if not (
                    op(intens_norm[i], threshold) or (
                        op == operator.le and
                        np.isnan(intens_norm[i])
                    )
                ):
                    
                    return False
        
        # additional constraints
        
        if not nanle(
            np.nanmax(no_protein) * .6,
            (
                np.nanmin(some_protein)
                    if not np.any(np.isnan(some_protein)) else
                np.nan
            )
        ):
            
            return False
        
        # exclude hollow shape profiles
        
        peak_area = [
            _samples.attrs.get_sample_index(sample_id)
            for sample_id in
            itertools.chain(
                profile_filter_args['peak'][0],
                profile_filter_args['small'][0]
            )
        ]
        
        for i in peak_area[1:-1]:
            
            if (
                intens_norm[i] * 2 < intens_norm[i - 1] and
                intens_norm[i] * 2 < intens_norm[i + 1]
            ):
                
                return False
        
        return True


class PeakSize(FeatureAnalyzer):
    """ """
    
    def __init__(
            self,
            samples,
            protein_samples,
            threshold = 2.0,
            name = 'peaksize',
        ):
        
        self.threshold = threshold
        
        FeatureAnalyzer.__init__(
            self,
            name = name,
            samples = samples,
            method = self.peak_size,
            protein_samples = protein_samples,
        )
    
    def peak_size(self, intensities, protein_samples, **kwargs):
        """

        Parameters
        ----------
        intensities :
            
        protein_samples :
            
        **kwargs :
            

        Returns
        -------

        """
        
        protein   = intensities[protein_samples]
        noprotein = intensities[np.logical_not(protein_samples)]
        
        if np.any(np.isnan(protein)):
            
            return False
        
        if np.all(np.isnan(noprotein)):
            
            return True
        
        return np.min(protein) > np.nanmax(noprotein) * self.threshold


class Slope(FeatureAnalyzer):
    """ """
    
    def __init__(
            self,
            samples,
            protein_profile,
            protein_samples,
            highest_samples,
            name = 'slope',
        ):
        
        self.threshold = threshold
        
        FeatureAnalyzer.__init__(
            self,
            name = name,
            samples = samples,
            method = self.slope,
            protein_samples = protein_samples,
        )
    
    def slope(
            self,
            intensities,
            protein_profile,
            protein_samples,
            highest_samples,
        ):
        """

        Parameters
        ----------
        intensities :
            
        protein_profile :
            
        protein_samples :
            
        highest_samples :
            

        Returns
        -------

        """
        
        pass
