#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `emese` python module
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

class MgfReader(object):
    
    def __init__(self, fname, name_callback = lambda x: {}):
        
        self.fname = fname
        self.process_name = name_callback
        
        for k, v in self.process_name(os.path.split(self.fname)[-1]):
            
            setattr(self, k, v)
    
    def index(self, fl, fr, charge = 1):
        """
        Looking up offsets in one MS2 mgf file.
        
        Columns:
            -- pepmass
            -- intensity
            -- retention time
            -- scan num
            -- offset in file
            -- fraction num
        """
        stRrtinseconds = 'RTINSECONDS'
        stRtitle = 'TITLE'
        stRbe = 'BE'
        stRen = 'EN'
        stRch = 'CH'
        stRpepmass = 'PEPMASS'
        stRempty = ''
        reln = re.compile(r'^([A-Z]+).*=([\d\.]+)[\s]?([\d\.]*)["]?$')
        features = []
        offset = 0
        cap_next = False
        
        with open(fl, 'rb', 8192) as fp:
            
            for l in fp:
                
                l = l.decode('ascii')
                
                if not l[0].isdigit() and \
                    not l[:2] == stRbe and not l[:2] == stRen:
                    
                    if not l[:2] == stRch:
                        
                        try:
                            m = reln.match(l).groups()
                        except:
                            print(fl, l)
                            continue
                        
                        if m[0] == stRtitle:
                            scan = float(m[1])
                        
                        if m[0] == stRrtinseconds:
                            rtime = float(m[1]) / 60.0
                        
                        if m[0] == stRpepmass:
                            pepmass = float(m[1])
                            intensity = 0.0 if m[2] == stRempty else float(m[2])
                            if charge is None:
                                cap_next = True
                        
                    else:
                        _charge = int(l[7]) if len(l) >= 8 else None
                        if charge is None or _charge == charge:
                            cap_next = True
                
                elif cap_next:
                    
                    features.append([pepmass, intensity,
                        rtime, scan, offset, fr])
                    scan = None
                    rtime = None
                    intensity = None
                    pepmass = None
                    _charge = None
                    cap_next = False
                
                offset += len(l)
        
        return features
