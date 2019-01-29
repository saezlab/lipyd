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

from past.builtins import xrange, range

import imp
import re
import collections
import mimetypes

import numpy as np
# for peak detection fits
import scipy as sp
import scipy.optimize
import scipy.spatial

import lipyd.settings as settings
import lipyd.common as common
import lipyd.reader.xls as xls
import lipyd.lookup as lookup


refrac = re.compile(r'([A-Z])([0-9]{1,2})')


Fraction = collections.namedtuple(
    'Fraction',
    ['row', 'col', 'start', 'end', 'mean']
)
Fraction.__new__.__defaults__ = (None,)


class SECReader(object):
    """ """
    
    def __init__(self, path):
        
        self.path = path
        self.read()
    
    def reload(self, children = False):
        """

        Parameters
        ----------
        children :
             (Default value = False)

        Returns
        -------

        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def read(self):
        """ """
        
        self.guess_format()
        
        if self.format == 'asc':
            
            self.read_asc()
            
        elif self.format == 'xls':
            
            self.read_xls()
    
    def read_asc(self):
        """Reads SEC UV absorbance profile from asc file output produced by
        the Unicorn software from GE Healthcare.

        Parameters
        ----------

        Returns
        -------

        """
        
        start      = None
        end        = None
        frac       = None
        volume     = []
        absorbance = []
        fractions  = []
        
        with open(self.path, 'r') as fp:
            
            for l in fp:
                
                l = l.strip().split('\t')
                
                if len(l) < 2 or '.' not in l[0] or '.' not in l[1]:
                    
                    continue
                
                vol = common.to_float(l[0])
                ab_ = common.to_float(l[1])
                
                if isinstance(vol, float) and isinstance(ab_, float):
                    
                    volume.append(vol)
                    absorbance.append(ab_)
                
                if len(l) > 3:
                    
                    start = end
                    end   = common.to_float(l[2])
                    
                    if start and end and frac:
                        
                        fractions.append(
                            Fraction(frac[0], int(frac[1]), start, end)
                        )
                    
                    m = refrac.search(l[3])
                    
                    frac = m.groups() if m else None
        
        self.volume = np.array(volume)
        self.absorbance = np.array(absorbance)
        self.fractions = fractions
    
    def read_xls(self):
        """Reads SEC UV absorbance profile from MS Excel XLS file output
        produced by ???.

        Parameters
        ----------

        Returns
        -------

        """
        
        volume     = []
        absorbance = []
        
        tab = xls.read_xls(self.path)
        
        for l in tab:
            
            if len(l) < 2 or '.' not in l[0] or '.' not in l[1]:
                
                continue
            
            vol = common.to_float(l[0])
            ab_ = common.to_float(l[1])
            
            if isinstance(vol, float) and isinstance(ab_, float):
                
                volume.append(vol)
                absorbance.append(ab_)
        
        self.volume = np.array(volume)
        self.absorbance = np.array(absorbance)
    
    def auto_fractions(
            self,
            start_volume = .6,
            size = .15,
            start_row = 'A',
            start_col = 5,
            length = 9,
        ):
        """Autogenerates fraction volume boundaries according to the parameters
        provided.

        Parameters
        ----------
        start_volume :
             (Default value = .6)
        size :
             (Default value = .15)
        start_row :
             (Default value = 'A')
        start_col :
             (Default value = 5)
        length :
             (Default value = 9)

        Returns
        -------

        """
        
        fractions = []
        
        for i in xrange(length):
            
            start = start_volume + size * i
            end   = start + size
            well  = (ord(start_row) - 65) * 12 + start_col + i - 1
            row   = chr(well // 12 + 65)
            col   = well % 12 + 1
            
            fractions.append(Fraction(row, col, start, end))
        
        return fractions
    
    def guess_format(self):
        """ """
        
        mime = mimetypes.guess_type(self.path)[0]
        
        self.format = (
            'xls' if 'excel' in mime or 'openxml' in mime else 'asc'
        )
    
    def _fractions_dict(self):
        
        if hasattr(self, 'fractions'):
            
            self.fractions_by_well = dict(
                ((fr.row, fr.col), fr)
                for fr in self.fractions
            )
    
    def get_fraction(self, row, col):
        
        return self._get_fraction(self.fraction_by_well(row, col))
    
    def _get_fraction(self, frac):
        """
        Returns absorbances measured within a fraction.

        Parameters
        ----------
        frac :
            

        Returns
        -------

        """
        
        return (
            self.absorbance[
                np.logical_and(
                    self.volume  < frac.end,
                    self.volume >= frac.start
                )
            ]
        )
    
    def fraction_by_well(self, row, col):
        """
        Returns fraction data by plate row letter and column number.
        """
        
        return self.fractions_by_well[(row, col)]
    
    def fraction_mean(self, row, col):
        """
        Returns the mean absorbance from a fraction.
        """
        
        return self._fraction_mean(self.fraction_by_well(row, col))
    
    def _fraction_mean(self, frac):
        """
        Returns the mean absorbance from a fraction.

        Parameters
        ----------
        frac :
            

        Returns
        -------

        """
        
        return self._get_fraction(frac).mean()
    
    def background_correction_by_other_chromatograms(
            self,
            others,
            start = None,
            end = None,
        ):
        """
        Subtracts corresponding values of other chromatograms.
        
        Parameters
        ----------
        other : list,SECReader
            One or more instances of ``SECReader``.
        start : float
            Start at this volume. By default the beginning of the
            chromatogram.
        end : float
            Do until this volume. By default the end of the chromatogram.
        """
        
        def get_closest(vol):
            
            return lookup.find(self.volume, vol, np.inf)
        
        istart = 0 if start is None else get_closest(start)
        iend = 0 if end is None else get_closest(end)
        
        others = (
            (others,)
                if not isinstance(others, (tuple, list, np.ndarray)) else
            others
        )
        
        for i in xrange(istart, iend + 1):
            
            vol = self.volume[i]
            i_other = (
                lookup.find(other.volume, vol, np.inf)
                for other in others
            )
            background = np.median([
                others[j].absorbance[io]
                for j, io in enumerate(i_other)
            ])
            self.absorbance[i] = self.absorbance[i] - background
    
    
    def normalize(self):
        """
        Simply subtracts the minimum and divides by the maximum across the
        chromatogram.
        """
        
        self.absorbance = self.absorbance - self.absorbance.min()
        self.absorbance = self.absorbance / self.absorbance.max()
    
    
    def profile(self, **kwargs):
        """Iterates fractions with their mean absorbance values.

        Parameters
        ----------
        **kwargs :
            

        Returns
        -------

        """
        
        fractions = (
            self.fractions
                if hasattr(self, 'fractions') else
            self.auto_fractions(**kwargs)
        )
        
        for frac in fractions:
            
            yield Fraction(*frac[:-1], self._fraction_mean(frac))
    
    def baseline_correction(self):
        """
        Performs a baseline correction by selecting points in the lowest
        decile and subtracting.
        """
        
    
    def baseline_correction_rubberband(self):
        """
        Performs a rubberband baseline correction following
        https://dsp.stackexchange.com/a/30723.
        """
        
        self.absorbance_baselinecorr = self.baseline_correction_rubberband(
            self.volume,
            self.absorbance,
        )
    
    @classmethod
    def rubberband_correction(cls, x, y):
        """
        Performs a rubberband baseline correction following
        https://dsp.stackexchange.com/a/30723.
        """
        
        rb = cls.rubberband(x, y)
        # Rotate convex hull vertices until they start from the lowest one
        rb = np.roll(rb, -rb.argmin())
        # Leave only the ascending part
        rb = rb[:rb.argmax()]
        # Create baseline using linear interpolation between vertices
        rb = np.interp(x, x[rb], y[rb])
        # The baseline is corrected like this
        return y - rb
    
    @staticmethod
    def rubberband(x, y):
        """
        Returns a convex hull of curve defined by ``x`` and ``y``.
        From https://dsp.stackexchange.com/a/30723.
        """
        
        return sp.spatial.ConvexHull(np.array(zip(x, y))).vertices
    
    def peak_detections(self):
        """
        Detects peaks by fitting and deconvoluting asymmetric decaying
        Gaussian functions.
        
        Following http://kitchingroup.cheme.cmu.edu/blog/2013/01/29/
        Curve-fitting-to-get-overlapping-peak-areas/
        """
        
        pass
