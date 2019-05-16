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

import itertools

import numpy as np


class Variable(object):
    
    
    def __init__(
            self,
            y,
            x = None,
            label = None,
            color = None,
            alpha = None,
            shape = None,
            typ   = 'line',
            
        ):
        """
        Represents a variable to be plotted along the features.
        
        y : numpy.array
            Values of the variable. If different length than the features
            ``x`` must be provided.
        x : numpy.array
            Coordinates along the x-axis.
        label : str
            Name of the variable.
        color : str,numpy.array
            Color of the variable, either a single value or an array of
            values for each data point.
        alpha : float,numpy.array
            Alpha values for the variable, either a single value or an array
            of values for each data point.
        shape : str
            Shape in case of plotting line with dots.
        typ : str
            Type of the plot: `line`, `bar` or `range`. The latter highlights
            a range provided in ``x``.
        """
        
        for name, val in iteritems(locals()):
            
            setattr(self, name, val)
        
        self.scale()
    
    
    def scale(self):
        
        self.y = self.y / np.nanmax(self.y)


class ProfileData(object):
    
    
    def __init__(
            self,
            features,
            features_x = None,
            variables = None,
            labels = None,
            colors = None,
            alphas = None,
            linewidths = None,
            shapes = None,
            title = None,
            ylab = r'Intensity [%]',
            xlab = 'Samples',
        ):
        """
        Represents data for making a plot of feature profiles across
        samples (e.g. fractions).
        
        Parameters
        ----------
        features : numpy.array
            An array with columns representing samples and rows features,
            values are intensitites.
        features_x : numpy.array
            If ``x`` provided for any of the variables it must be provided
            also here in order to correctly align the features with the
            variables.
        variables : lipyd.profile.Variable
            Additional variables to plot along the features.
        labels : list
            A list of strings for tick labels. If not provided labels will
            be a sequence of integers starting from 1.
        colors : str,numpy.array
            Either a single color code or a color for each feature.
        alphas : float,numpy.array
            Either a single alpha or a value for each feature.
        shapes : str,numpy.array
            Either a single shape type or a shape for each feature.
        title : str
            Main title for the plot.
        """
        
        for name, val in iteritems(locals()):
            
            setattr(self, name, val)
        
        self.scale()
        self.set_xaxis()
    
    
    def scale(self):
        """
        Scales features to the 0-1 range.
        """
        
        if len(self.features):
            
            self.features = self.features / np.nanmax(self.features)
            self.features = np.nan_to_num(self.features)
    
    
    def set_xaxis(self):
        """
        Sets the scale and tick locations of the x-axis.
        """
        
        self._variables_have_x()
        self._set_xrange()
        self._set_xlabels()
    
    
    def _variables_have_x(self):
        
        self._custom_x = (
            self.variables is not None and
            any(v.x is not None for v in self.variables)
        )
    
    
    def _set_xrange(self):
        
        if self._custom_x:
            
            self.xmax = max(
                max(
                    map(
                        np.nanmax,
                        (v.x for x in self.variables if v.x is not None)
                    )
                ),
                np.nanmax(self.features_x)
            )
            
            self.xmin = min(
                min(
                    map(
                        np.nanmin,
                        (v.x for x in self.variables if v.x is not None)
                    )
                ),
                np.nanmin(self.features_x)
            )
            
        else:
            
            self.xmax = len(self)
            self.xmin = 0
        
        step = (self.xmax - self.xmin) / len(self)
        common_x = np.arange(self.xmin + step, self.xmax + step, step)
        
        for var in self.variables:
            
            if var.x is None:
                
                var.x = common_x
        
        if self.features_x is None:
            
            self.features_x = common_x
    
    
    def _set_xlabels(self):
        
        if self.labels is None:
            
            self.labels = [str(x + 1) for x in xrange(len(self))]
    
    
    def __len__(self):
        
        return self.features.shape[1]
