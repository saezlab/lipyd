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

class MgfReader(object):
    
    def __init__(self, fname, charge = 1, name_callback = lambda x: {}):
        
        self.fname = fname
        self.process_name = name_callback
        self.charge_to_read = charge
        
        for k, v in self.process_name(os.path.split(self.fname)[-1]):
            
            setattr(self, k, v)
        
        self.var = ['mz', 'intensity', 'rtime', 'scan', 'charge', 'offset']
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def read(self):
        """
        Reads MS1 data from MGF file.
        Indexes offsets of MS2 scans with purpose of later reading.
        """
        
        stRrtinseconds = 'RTINSECONDS'
        stRtitle = 'TITLE'
        stRbe = 'BE'
        stRen = 'EN'
        stRch = 'CH'
        stRpepmass = 'PEPMASS'
        stRempty = ''
        reln = re.compile(r'^([A-Z]+).*=([\d\.]+)[\s]?([\d\.]*)["]?$')
        renondigit = re.compile(r'[^\d\.-]+')
        for var in self.var:
            
            setattr(self, var, [])

        offset = 0
        cap_next = False
        
        with open(self.fname, 'rb', 8192) as fp:
            
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
                            scan = int(m[1])
                        
                        if m[0] == stRrtinseconds:
                            rtime = float(m[1]) / 60.0
                        
                        if m[0] == stRpepmass:
                            pepmass = float(m[1])
                            intensity = 0.0 if m[2] == stRempty else float(m[2])
                        
                    else:
                        charge = int(l[7:-2]) if len(l) >= 8 else None
                        if (self.charge_to_read is None or
                            charge == self.charge_to_read):
                            
                            cap_next = True
                
                elif cap_next:
                    
                    self.mz.append(pepmass)
                    self.intensity.append(intensity)
                    self.rtime.append(rtime)
                    self.scan.append(scan)
                    self.charge.append(
                        charge if charge is not None else np.nan
                    )
                    self.offset.append(offset)
                    scan = None
                    rtime = None
                    intensity = None
                    pepmass = None
                    charge = None
                    cap_next = False
                
                offset += len(l)
        
        for var in self.var:
            
            setattr(self, var, np.array(getattr(self, var)))
