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

from past.builtins import xrange

import os
import sys
import imp
import re
import csv
import mimetypes
import warnings
import numpy as np

import lipyd.sample
import lipyd.reader.xls
import lipyd.reader.common as common


relabel = re.compile(r'(.*)_([A-Z])([0-9]{1,2})_(neg|pos)')


class PeaksReader(object):
    
    rehdr = re.compile(r'(.*)(m/z|RT mean|Normalized Area)')
    rertr = re.compile(r'([\d\.]+) - ([\d\.]+)')
    hdr_7 = [
        'Peptide',
        'Quality',
        'Significance (-10lgP)',
        'm/z',
        'RT range',
        'z',
        'Avg. Area',
    ]
    ignore = {
        'Sample Profile (Ratio)',
        'Control Normalized Area',
        'Protein peaks Normalized Area',
        'Group Profile (Ratio)',
        'RT mean',
        'Accession',
        'PTM',
    }
    
    def __init__(
            self,
            fname,
            ionmode = None,
            format = None,
            label_processor = None,
            sample_sorter = None,
        ):
        """
        Reads data from an output file of the PEAKS software.
        
        PEAKS is a proprietary software for preprocessing of LC MS/MS data.
        Among many other things it is able to align peaks across multiple
        samples. Its output is a ``csv`` or ``xlsx`` file.
        See more at http://www.bioinfor.com/
        
        Args
        ----
        :param callable label_processor:
            A method to process headers of coulumns 7+.
        """
        
        self.set_file(fname)
        self.ionmode = common.guess_ionmode(ionmode, self.fname)
        self.guess_format(format)
        self.label_processor = label_processor or self.default_label_processor
        self.sample_sorter = sample_sorter or self.default_sample_sorter
        
        self.read()
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def __repr__(self):
        
        return '<%s.%s fname=\'%s\'>' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.fname,
        )
    
    def read(self):
        
        # feature data
        quality           = []
        significance      = []
        z                 = []
        total_intensities = []
        centr_mzs         = []
        rt_ranges         = []
        total_rt_means    = []
        
        # sample data
        mzs               = []
        rt_means          = []
        intensities       = []
        
        lines = self.iterlines()
        self.hdr_raw = next(lines)
        self.process_header()
        
        for i, line in enumerate(lines):
            
            rtrange = self.rertr.search(line[4])
            
            if not rtrange:
                
                warnings.warn(
                    'Could not parse RT range: %s\n'
                    'File `%s`, line %u' % (line[4], self.fname, i)
                )
                rtrange = (np.nan, np.nan)
            
            quality.append(common.to_float(line[1]))
            significance.append(common.to_float(line[2]))
            centr_mzs.append(common.to_float(line[3]))
            rt_ranges.append((
                common.to_float(rtrange[0]),
                common.to_float(rtrange[1]),
            ))
            z.append(common.to_int(line[5]))
            total_intensities.append(common.to_float(line[6]))
            
            if self.rt_mean_idx:
                total_rt_means.append(common.to_float(line[self.rt_mean_idx]))
            
            mzs.append([
                common.to_float(line[sample['m/z']])
                for sample in self.samples
            ])
            intensities.append([
                common.to_float(line[sample['Normalized Area']])
                for sample in self.samples
            ])
            rt_means.append([
                common.to_float(line[sample['RT mean']])
                for sample in self.samples
            ])
        
        # feature data
        self.quality           = np.array(quality)
        self.significance      = np.array(significance)
        self.z                 = np.array(z)
        self.total_intensities = np.array(total_intensities)
        self.centr_mzs         = np.array(centr_mzs)
        self.rt_ranges         = np.array(rt_ranges)
        self.total_rt_means    = np.array(total_rt_means)
        
        # sample data
        self.mzs               = np.array(mzs)
        self.rt_means          = np.array(rt_means)
        self.intensities       = np.array(intensities)
    
    def process_header(self):
        
        self.samples = []
        
        for i in xrange(7):
            
            if self.hdr_raw[i] != self.hdr_7[i]:
                
                warnings.warn(
                    'Unexpected column header while reading '
                    'PEAKS output file:\n'
                    'column %u expected to be `%s` but `%s` found\n'
                    'in file `%s` ' % (
                        i + 1,
                        self.hdr_7[i],
                        self.hdr_raw[i],
                        self.fname,
                    )
                )
        
        for i in xrange(7, len(self.hdr_raw), 3):
            
            if self.hdr_raw[i] in self.ignore or i + 1 > len(self.hdr_raw):
                
                continue
            
            sample = {}
            
            match = self.rehdr.search(self.hdr_raw[i])
            
            if match:
                
                label, field = match.groups()
                sample['label_raw'] = label
                sample['label'] = self.label_processor(label)
            
            else:
                
                self._column_label_warning(i)
            
            # repeating column triplets
            for j in xrange(3):
                
                col_idx = i + j
                
                match = self.rehdr.search(self.hdr_raw[col_idx])
                
                if match:
                    
                    label, field  = match.groups()
                    sample[field] = col_idx
                    
                else:
                    
                    self._column_label_warning(col_idx)
            
            self.samples.append(sample)
        
        self.samples = self.sample_sorter(self.samples)
        
        try:
            self.rt_mean_idx = self.hdr_raw.index('RT mean')
        except ValueError:
            self.rt_mean_idx = None
    
    def iterlines(self):
        
        if self.format == 'xls':
            
            for line in lipyd.reader.xls.read_xls(self.fname):
                
                yield line
            
        elif self.format == 'csv':
            
            with open(self.fname, 'r') as fp:
                
                dialect = csv.Sniffer().sniff(fp.read(20000))
                fp.seek(0)
                
                for line in csv.reader(fp, dialect):
                    
                    yield line
    
    def guess_format(self, format):
        
        if hasattr(format, 'lower'):
            
            self.format = format.lower()
            self.format = 'xls' if self.format == 'xlsx' else self.format
            
        else:
            
            mime = mimetypes.guess_type(self.fname)
            self.format = (
                'xls' if 'excel' in mime or 'openxml' in mime else 'csv'
            )
    
    def set_file(self, fname):
        
        if os.path.exists(fname):
            
            self.fname = fname
            
        else:
            
            raise FileNotFoundError(fname)
    
    def _column_label_warning(self, col_idx):
        
        warnings.warn(
            'Could not recognize column label `%s`.\nIn PEAKS '
            'output file coulumn triplets expected to end '
            '`m/z`, `RT mean` and `Normalized Area`.\n'
            'In file `%s`' % (
                self.hdr_raw[col_idx],
                self.fname,
            )
        )
    
    @staticmethod
    def default_label_processor(label):
        
        match = relabel.search(label)
        
        if match:
            
            main, row, col, ionmode = match.groups()
            col = int(col)
            
            return {
                'main': main,
                'fraction': (row, col),
                'ionmode': ionmode,
            }
    
    @staticmethod
    def default_sample_sorter(samples):
        
        return sorted(samples, key = lambda s: s['label']['fraction'])
    
    def get_attributes(self):
        """
        Returns ``lipyd.sample.FeatureAttributes`` object.
        This object contains variables describing series of features
        across all samples. Quality, significance, mean RT, centroid m/z, etc
        """
        
        return lipyd.sample.FeatureAttributes(
            quality = self.quality,
            significance = self.significance,
            charge = self.z,
            total_intensities = self.total_intensities,
            centr_mzs = self.centr_mzs,
            rt_ranges = self.rt_ranges,
            rt_means = self.total_rt_means,
        )
    
    def get_samples(self, bind = True):
        """
        Yields ``lipyd.sample.Sample`` objects for each sample read from
        the PEAKS output file.
        
        To extract all data from ``PeaksReader`` the ``get_sampleset`` method
        is more convenient.
        
        :param bool bind:
            Bind samples to each other. This way they will sort together, i.e.
            if any of them is sorted all the others follow the same order.
        """
        
        feature_attrs = self.get_attributes()
        sorter = attrs.sorter if bind else None
        
        for i, sample_attrs in enumerate(self.samples):
            
            if not bind:
                
                feature_attrs = self.get_attributes()
            
            yield lipyd.sample.Sample(
                mzs = self.mzs[:,i],
                intensities = self.intensities[:,i],
                rts = self.rt_means[:,i],
                attrs = sample_attrs,
                feature_attrs = feature_attrs,
                sorter = sorter,
            )
    
    def get_sampleset(self):
        """
        Returns a ``SampleSet`` and a ``FeatureAttributes`` object.
        """
        
        return lipyd.sample.SampleSet(self.get_samples()), None
