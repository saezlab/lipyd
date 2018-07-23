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


class MgfReader(object):
    
    stRrtinseconds = 'RTINSECONDS'
    stRtitle = 'TITLE'
    stRbe = 'BE'
    stRen = 'EN'
    stRch = 'CH'
    stRpepmass = 'PEPMASS'
    stRempty = ''
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
        self.mgfindex = np.array(sorted(features, key = lambda x: x[0]))
    
    def lookup(self, mz, rt = None):
        """
        Looks up an MS1 m/z and returns MS2 scans.
        """
        
        mz_uncorr = mz / self.drift
        
        idx = lookup.findall(self.mgfindex[:,0], mz)
