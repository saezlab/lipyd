#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2015-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
    # filtering features accross samples and fractions
    # to find those relevant ones 
    # belonging to LTP bound lipids
#

import numpy as np
import scipy as sp
import pandas as pd
import time
import os
import sys
import copy
import cPickle as pickle

from future.utils import iteritems

from ltp import ltp2

antof = 'highest_antonella.tab'
outf = 'expected_features.tab'
tlr2 = 0.02
tlr1 = 0.01
tlr05 = 0.005
only_closest = False
hdr = ['protein', 'mode', 'mzo', 'mz', 'i', 'fr', 'fe']

l = ltp2.Screening()

hdata = {}

with open(antof, 'r') as f:
    
    for r in f:
        
        r = r.strip().split('\t')
        
        if r[0] not in hdata:
            hdata[r[0]] = {}
        
        if r[1] not in hdata[r[0]]:
            hdata[r[0]][r[1]] = []
        
        hdata[r[0]][r[1]].append(float(r[2]))

l.get_filenames()
l.read_data()

out = []

for protein, ha in iteritems(hdata):
    
    for mode, h in iteritems(ha):
        
        if protein not in l.data:
            
            continue
        
        tbl = l.data[protein][mode]
        raw = tbl['raw']
        ann = tbl['ann']
        
        for mz in h:
            
            ui = ann[:,2].searchsorted(mz)
            
            if only_closest:
                
                tlr = tlr2
                
                i = ui if ann[ui,2] - mz < mz - ann[ui - 1,2] < tlr else ui - 1 if mz - ann[ui - 1,2] < tlr else None
                
                if i:
                    
                    for ifr, fr in enumerate(l.pfracs[protein]):
                        
                        out.append([protein, mode, mz, ann[ui + u,2], raw[ui + u,ifr], fr, 'F%u' % i])
                
                continue
            
            tlr = tlr05
            
            if ann[ui,2] - mz > tlr and mz - ann[ui-1,2] > tlr:
                
                tlr = tlr1
                
                if ann[ui,2] - mz > tlr and mz - ann[ui-1,2] > tlr:
                    
                    tlr = tlr2
            
            u = 0
            
            while True:
                
                if ui + u >= raw.shape[0] or ann[ui + u,2] - mz > tlr:
                    
                    break
                
                for ifr, fr in enumerate(l.pfracs[protein]):
                    
                    out.append([protein, mode, mz, ann[ui + u,2], raw[ui + u,ifr], fr, 'F%u' % (ui + u)])
                
                u += 1
            
            d = 1
            
            while True:
                
                if ui - d < 0 or mz - ann[ui - d,2] > tlr:
                    
                    break
                
                for ifr, fr in enumerate(l.pfracs[protein]):
                    
                    out.append([protein, mode, mz, ann[ui - d,2], raw[ui + u,ifr], fr, 'F%u' % (ui - d)])
                
                d += 1

with open(outf, 'w') as f:
    
    f.write('%s\n' % '\t'.join(hdr))
    
    f.write(
        '\n'.join(
            map(
                lambda r:
                    '\t'.join(
                        map(
                            lambda i:
                                '%.08f' % i if type(i) is float else str(i),
                            r
                        )
                    ),
                out
            )
        )
    )
