#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This code is not for public use
#  For permission please contact the authors
#
#  Copyright (c) 2019 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Website: http://saezlab.github.io/lipyd
#

import imp
import csv
import itertools
import copy

from lxml import etree
import numpy as np
import pyopenms as oms

from lipyd import lipid
from lipyd import settings
from lipyd import lookup
from lipyd import sample
from lipyd import sampleattrs


class PeakPickingEvaluation(object):
    
    def __init__(
            self,
            feature_xml_fname,
            peaks_fname,
            examples_fname,
            outfile,
            progenesis_fname = None,
            sample_id = ('A', 10),
        ):
        
        self.peaks_fname = peaks_fname
        self.feature_xml_fname = feature_xml_fname
        self.progenesis_fname = progenesis_fname
        self.examples_fname = examples_fname
        self.outfile = outfile
        self.sample_id = sample_id
    
    
    def main(self):
        
        self.pc_masses()
        self.read_examples()
        self.read_peaks()
        self.lookup_peaks()
        self.read_progenesis()
        self.lookup_progenesis()
        self.collect_convex_hulls()
        self.export()
    
    
    def read_examples(self):
        
        
        with open(self.examples_fname, 'r') as fp:
            
            self.examples = fp.read().split('\n')
        
        self.examples = [ex.strip() for ex in self.examples if ex.strip()]
    
    
    def read_peaks(self):
        
        reader = sample.SampleReader(
            input_type = 'peaks',
            fname = self.peaks_fname,
        )
        
        self.samples = reader.get_sampleset(
            sampleset_args = {
                'sample_id_proc': sampleattrs.plate_sample_id_processor(),
            }
        )
        
        self.sample_selected = self.samples.mzs_by_sample[
            :,self.samples.attrs.sample_id_to_index[self.sample_id]
        ]
        
        idx = self.sample_selected.argsort()
        
        self.samples.sort_all(by = idx)
        
        self.sample_selected = self.samples.mzs_by_sample[
            :,self.samples.attrs.sample_id_to_index[self.sample_id]
        ]
    
    
    def read_progenesis(self):
        
        # later RTs and isotopes could be read,
        # no only iterating over all m/z's
        # just to see the peaks picked by Progenesis
        
        if not self.progenesis_fname:
            
            return
        
        parser = etree.iterparse(
            self.progenesis_fname,
            events = ('start', 'end'),
        )
        root = next(parser)
        used_elements = []
        mzs = []
        
        for ev, elem in parser:
            
            if ev == 'end' and elem.tag == 'mz':
                
                mzs.append(float(elem.text))
            
            used_elements.append(elem)
            
            # removing used elements to keep memory low
            if len(used_elements) > 1000:
                
                for _ in range(500):
                    
                    e = used_elements.pop(0)
                    e.clear()
        
        self.progenesis_mzs = np.array(mzs)
        self.progenesis_mzs.sort()
    
    
    def lookup_progenesis(self):
        
        self.progenesis_peaks = {}
        
        if not hasattr(self, 'progenesis_mzs'):
            
            return
        
        for ex in self.examples:
            
            mz_theoretical = self.pc_adduct_masses[(ex, 0)]
            
            i = lookup.find(
                self.progenesis_mzs, # all masses in the sample
                mz_theoretical, # mass to search for
                t = 10, # tolerance in ppm
            )
            
            if i:
                
                self.progenesis_peaks[ex] = self.progenesis_mzs[i]
    
    
    def lookup_peaks(self):
        
        self.peaks_peaks = {}
        
        for ex in self.examples:
            
            mz_theoretical = self.pc_adduct_masses[(ex, 0)]
            
            i = lookup.find(
                self.sample_selected, # all masses in the sample
                mz_theoretical, # mass to search for
                t = 10, # tolerance in ppm
            )
            
            if i is not None:
                
                self.peaks_peaks[ex] = (
                    mz_theoretical,
                    self.sample_selected[i],
                )
                
            else:
                
                print('Example not found: %s %.04f' % (ex, mz_theoretical))
    
    
    def pc_masses(self):
        
        
        self.pc_adduct_masses = dict(
            (
                (pc.name, isotope),
                pc.add_h()
            )
            for isotope in range(5)
            for pc in lipid.PC(
                fa_args = {'c': (15, 20), 'u': (0, 4)},
                isotope = isotope, sum_only = True
            )
        )
    
    
    def collect_convex_hulls(self):
        
        self.convex_hulls = []
        
        # opening featureXML
        xml_file = oms.FeatureXMLFile()
        self.fmap = oms.FeatureMap()
        xml_file.load(self.feature_xml_fname, self.fmap)
        feature_mzs = []
        
        for i, fe in enumerate(self.fmap):
            
            feature_mzs.append([i, fe.getMZ()])
        
        feature_mzs = np.array(feature_mzs)
        feature_mzs = feature_mzs[feature_mzs[:,1].argsort(),:]
        
        # looking up the example features in the featureXML
        self.examples_oms_features = dict(
            (
                feature_mzs[lookup.find(feature_mzs[:,1], mz_m, t = 10), 0],
                ex,
            )
            for ex, (mz_t, mz_m) in self.peaks_peaks.items()
        )
        
        # collecting convex hulls
        for ife, fe in enumerate(self.fmap):
            
            if ife in self.examples_oms_features:
                
                hull_list = fe.getConvexHulls()
                
                self.extend_hulls(hull_list, ife, 0)
                
                subord_feature = fe.getSubordinates()
                
                if subord_feature:
                    
                    for subfe in subord_feature:
                        
                        hull_list = subfe.getConvexHulls()
                        
                        self.extend_hulls(hull_list, ife, 1)
        
        # columns: rt, mz, feature index, hull index, is sub-feature
        self.convex_hulls = np.vstack(self.convex_hulls)
        self.oms_feature_mzs = feature_mzs[feature_mzs[:,0].argsort(),:]
    
    
    def extend_hulls(self, hull_list, ife, sub):
        
        for ihull, hull in enumerate(hull_list):
            
            hull_points = hull.getHullPoints() # hull_points is numpy.ndarray
            hull_points = copy.copy(hull_points)
            hull_points = np.hstack((
                hull_points,
                np.full((hull_points.shape[0], 1), ife),
                np.full((hull_points.shape[0], 1), ihull),
                np.full((hull_points.shape[0], 1), sub),
            ))
            self.convex_hulls.append(hull_points)
    
    
    def export(self, fname = None):
        
        num_str = lambda num: '%.012f' % num
        
        fname = fname or self.outfile
        
        with open(fname, 'w') as fp:
            
            _ = fp.write('\t'.join([
                'lipid_species',
                'mz_theoretical',
                'mz_feature_peaks',
                'mz_feature_progensis',
                'mz_feature_oms',
                'mz_scan_oms',
                'rt_scan_oms',
                'isotope',
                'feature_id_oms',
                'convex_hull_id_oms',
                'is_sub_feature',
            ]))
            
            _ = fp.write('\n')
            
            for rt_scan, mz_scan, ife, ihull, is_sub in self.convex_hulls:
                
                ex = self.examples_oms_features[ife]
                
                mz_theoretical, mz_feature_peaks = (
                    self.peaks_peaks[ex]
                        if ex in self.peaks_peaks else
                    (
                        self.pc_adduct_masses[(ex, int(ihull))],
                        np.nan,
                    )
                )
                
                mz_feature_progenesis = (
                    self.progenesis_peaks[ex]
                        if (
                            hasattr(self, 'progenesis_peaks') and
                            ex in self.progenesis_peaks
                        ) else
                    np.nan
                )
                
                mz_feature_oms = self.oms_feature_mzs[int(ife)][1]
                
                _ = fp.write('\t'.join([
                    ex,
                    num_str(mz_theoretical),
                    num_str(mz_feature_peaks),
                    num_str(mz_feature_progenesis),
                    num_str(mz_feature_oms),
                    num_str(mz_scan),
                    num_str(rt_scan),
                    '%u' % int(ihull),
                    '%u' % int(ife),
                    str(is_sub > 0),
                ]))
                
                _ = fp.write('\n')
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
