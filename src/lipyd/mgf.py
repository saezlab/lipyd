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

import os
import re
import imp
import numpy as np

import lipyd.lookup as lookup
import lipyd.session as session
import lipyd.settings as settings


class MgfReader(object):
    
    stRrtinseconds = 'RTINSECONDS'
    stRtitle = 'TITLE'
    stRbe = 'BE'
    stRen = 'EN'
    stRch = 'CH'
    stRpepmass = 'PEPMASS'
    stRempty = ''
    stRcharge = 'CHARGE'
    reln = re.compile(r'^([A-Z]+).*=([\d\.]+)[\s]?([\d\.]*)["]?$')
    
    def __init__(
            self,
            fname,
            label = None,
            charge = 1,
            rt_tolerance = 1.0,
            drift = 1.0
        ):
        """
        Provides methods for looking up MS2 scans from an MGF file.
        """
        
        self.fname  = fname
        self.label  = label
        self.charge = charge
        self.rt_tolerance = rt_tolerance
        self.drift  = drift
        self.index()
        self.log = session.get_log()
        self.ms2_within_range = settings.get('ms2_within_range')
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def index(self):
        """
        Indexing offsets in one MS2 MGF file.
        
        Columns:
            -- pepmass
            -- intensity
            -- retention time
            -- scan num
            -- offset in file
            -- fraction num
        """
        
        features = []
        offset = 0
        cap_next = False
        
        with open(self.fname, 'rb', 8192) as fp:
            
            for l in fp:
                
                l = l.decode('ascii')
                
                if (
                    not l[0].isdigit() and
                    not l[:2] == self.stRbe and
                    not l[:2] == self.stRen
                ):
                    
                    if not l[:2] == self.stRch:
                        
                        try:
                            
                            m = self.reln.match(l).groups()
                            
                        except:
                            
                            sys.stdout.write(
                                'Line in MGF file `%s`'
                                'could not be processed: '
                                '\n\t`%s`\n' % (fl, l)
                            )
                            continue
                        
                        if m[0] == self.stRtitle:
                            
                            scan = float(m[1])
                        
                        if m[0] == self.stRrtinseconds:
                            
                            rtime = float(m[1]) / 60.0
                        
                        if m[0] == self.stRpepmass:
                            
                            pepmass = float(m[1])
                            intensity = (
                                0.0
                                    if m[2] == self.stRempty else
                                float(m[2])
                            )
                            
                            if self.charge is None:
                                
                                cap_next = True
                        
                    else:
                        _charge = int(l[7]) if len(l) >= 8 else None
                        if self.charge is None or _charge == self.charge:
                            cap_next = True
                
                elif cap_next:
                    
                    features.append([
                        pepmass, # precursor ion mass
                        intensity, # intensity
                        rtime, # retention time
                        scan, # scan ID
                        offset, # byte offset in file
                        self.label # fraction ID
                    ])
                    # reset all values
                    scan = None
                    rtime = None
                    intensity = None
                    pepmass = None
                    _charge = None
                    cap_next = False
                
                offset += len(l)
        
        # sorted by precursor mass
        self.mgfindex = np.array(
            sorted(features, key = lambda x: x[0]),
            dtype = np.object
        )
    
    def lookup(self, mz, rt = None):
        """
        Looks up an MS1 m/z and returns the indices of MS2 scans in the
        MGF file.
        
        Returns 2 numpy arrays of the same length: first one with the indices,
        second one with the RT differences.
        """
        
        if self.log.verbosity > 4:
            
            self.log.msg(
                'Recalibrated m/z: %.08f; drift = %.08f; '
                'measured m/z: %.08f' % (mz, self.drift, mz / self.drift)
            )
        
        rt = rt or np.nan
        mz_uncorr = mz / self.drift
        
        idx    = np.array(lookup.findall(self.mgfindex[:,0], mz_uncorr))
        rtdiff = np.array([self.mgfindex[i,2] - rt for i in idx])
        
        if self.log.verbosity > 4:
            
            self.log.msg(
                'Looking up MS1 m/z %.08f. '
                'MS2 scans with matching precursor mass: %u' % (
                    mz_uncorr,
                    len(idx),
                )
            )
        
        if self.ms2_within_range:
            
            if np.isnan(rt):
                
                self.log.msg(
                    'No MS1 RT provided, could not check RT '
                    'difference of MS2 scans.'
                )
                
            else:
                
                idx = idx[
                    np.logical_or(
                        np.isnan(rtdiff),
                        np.abs(rtdiff) < self.rt_tolerance
                    )
                ]
                
                if self.log.verbosity > 4:
                    
                    self.log.msg(
                    'RT range: %.03f--%.03f; '
                    'Matching MS2 scans within this range: %u' % (
                        rtlower,
                        rtupper,
                        len(idx),
                    )
                )
            
        elif self.log.verbosity > 4:
            
            self.log.msg('Not checking RT.')
        
        return idx, rtdiff
    
    def get_scan(self, i):
        """
        Reads MS2 fragment peaks from one scan.
        
        Returns m/z's and intensities in 2 columns array.
        """
        
        scan = []
        
        self.get_file()
        # jumping to offset
        self.fp.seek(int(self.mgfindex[i, 4]), 0)
        
        # zero means no clue about charge
        charge = 0
        
        for l in self.fp:
            
            if l[:6] == self.stRcharge:
                # one chance to obtain the charge
                charge = int(l.strip()[-2])
                continue
            if not l[0].isdigit():
                # end of scan
                break
            else:
                # reading fragment masses
                mi = l.strip().split()
                scan.append([
                    float(mi[0]), # mass
                    float(mi[1]) if len(mi) > 1 else np.nan # intensity
                ])
        
        if self.log.verbosity > 4:
            
            self.log.msg(
                'Read scan #%u from file `%s`;'
                '%u peaks retrieved.' % (
                    self.mgfindex[i, 3],
                    self.fname,
                    len(scan),
                )
            )
        
        return np.array(scan)
    
    def get_scans(self, mz, rt = None):
        """
        Looks up all scans for one precursor mass and RT, yields 2 column
        arrays of MS2 m/z's and intensities and delta RTs.
        
        Calls `get_scan` for all indices returned by `lookup`.
        """
        
        idx, rtdiff = self.lookup(mz, rt = rt)
        
        for i, r in zip(idx, rtdiff):
            
            yield self.get_scan(i), r
    
    def get_file(self):
        """
        Returns the file pointer, opens the file if necessary.
        """
        
        if not hasattr(self, 'fp') or self.fp.closed:
            
            self.fp = open(self.fname, 'r')
        
    def __del__(self):
        
        if hasattr(self, 'fp') and not self.fp.closed:
            
            self.fp.close()
