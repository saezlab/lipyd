#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `emese` python module
#
#  Copyright (c) 2015-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import sys
import re
import numpy as np
import itertools
import imp

from emese.common import *

class Feature(object):
    """
    Provides additional, more sophisticated methods
    for identification of a single feature.
    
    In the original concept all methods for identification
    based on MS1 and MS2 took place in class Screening(),
    as those could simply iterate through the arrays.
    
    Later more complex methods became necessary, so
    I created this class to group them.
    """
    
    def __init__(self, main, protein, mode, oi, log = True):
        """
        @main : ltp.Screening() instance
            One Screening() instance with MS1 and MS2 processing already done.
        
        @protein : str
            Protein name
        
        @mode : str
            MS mode (`pos` or `neg`)
        
        @oi : int
            Original index of one feature.
        
        @log : bool
            Whether output verbose messages to logfile.
        """
        self.main = main
        self.log = log
        self.protein = protein
        self.mode = mode
        self.oi = oi
        self.ifracs = self.main.fraction_indices(self.protein)
        self.fracsi = dict(map(lambda fr: (fr[1][0], fr[0]),
                               iteritems(self.ifracs)))
        self.tbl = self.main.valids[self.protein][self.mode]
        self.ms2 = self.tbl['ms2'][self.oi]
        self.i = self.main.oi2i(self.protein, self.mode, self.oi)
        self.fa = {}
        self.scans_fractions = map(
            lambda tpl: tuple(map(int, tpl)),
            uniqList(
                map(
                    tuple,
                    # scan ID, fraction ID
                    self.ms2[:,[12,14]]
                )
            )
        )
        self.classes = ['PA', 'PC', 'PE', 'PG', 'PS']
        self.classes2 = ['PA', 'PC', 'PE', 'PG', 'PS', 'PI', 'SM', 'BMP',
                         'Cer', 'Cer1P', 'HexCer', 'DAG', 'TAG', 'FA', 'VA']
        self.identities = set([])
        self.identities2 = {}
        # get carbon counts from MS1
        self.ms1fa = self.tbl['ms1fa'][oi]
        # sorting by fractions/scans
        self.scans = dict(
            map(
                lambda sc_fr:
                    (
                        # scan ID, fraction ID: key
                        (sc_fr[0], sc_fr[1]),
                        # MS2 array slice: value
                        self.ms2[
                            np.where(
                                np.logical_and(
                                    self.ms2[:,12] == sc_fr[0],
                                    self.ms2[:,14] == sc_fr[1]
                                )
                            )
                        ]
                    ),
                self.scans_fractions
            )
        )
        # sorting by intensity desc
        self.scans = dict(
            map(
                lambda i:
                    (
                        i[0],
                        i[1][i[1][:,2].argsort()[::-1],:]
                    ),
                iteritems(self.scans)
            )
        )
        self.deltart = dict(
            map(
                lambda i:
                    (
                        i[0],
                        self.tbl['rtm'][self.i] - i[1][0,11]
                    ),
                iteritems(self.scans)
            )
        )
        self._scans = dict(
            map(
                lambda i:
                    (
                        i[0],
                        # i[0]: (scan ID, fraction ID)
                        # i[1]: MS2 array slice
                        MS2Scan(i[1], i[0], self)
                    ),
                iteritems(self.scans)
            )
        )
        self.maxins = dict(
            map(
                lambda i:
                    (
                        i[0],
                        i[1][0,2]
                    ),
                iteritems(self.scans)
            )
        )
        self.medins = dict(
            map(
                lambda i:
                    (
                        i[0],
                        np.median(i[1][:,2])
                    ),
                iteritems(self.scans)
            )
        )
        self.sort_scans()
        self.select_best_scan()
        self.msg('\n::: Analysing feature: %s :: %s :: index = %u ::'\
                ' m/z = %.03f :: number of MS2 scans: %u\n' % \
            (self.protein, self.mode, self.oi, self.tbl['mz'][self.i],
                len(self._scans))
        )
        self.msg('\n::: Database lookup resulted '\
            'the following species: %s\n' % self.print_db_species())
        self.msg('\n::: Intensities:\n%s%s\n' % \
            (' ' * 24, '          '.join(['A09', 'A10', 'A11', 'A12', 'B01'])))
        self.msg('%s%s' % (' ' * 16, '=' * 63))
        self.msg('\n    - absolute:  %s' % '   '.join(
            map(lambda x: '%10.01f' % x, self.tbl['fe'][self.i,:]))
        )
        self.msg('\n    - relative: %s\n' % \
            '  '.join(
                map(
                    lambda xx:
                        '%10.02f%%' % (xx * 100.0),
                    map(
                        lambda x:
                            x / np.nanmax(self.tbl['fe'][self.i,:]),
                        self.tbl['fe'][self.i,:]
                    )
                )
            )
        )
        self.msg('\n::: MS2 scans available (%u):\n\n' % len(self.scans))
        
        for sc in self._scans.values():
            sc.print_scan()
    
    def sort_scans(self):
        """
        Groups the scans in 3 groups: highest consists of those from the
        fractions with the highest protein level (there might be more than
        one the highest, because the fraction offset limits); the secondary
        contains scans from other protein containing fractions; while the
        other contains the scans from non protein containing fractions.
        Within the groups the scans are sorted from lowest to highest
        deltaRT.
        """
        self.highest = []
        self.secondary = []
        self.other = []
        with_protein = self.main.protein_containing_fractions(self.protein)
        for scan_num, fr in self.scans.keys():
            fr_name = 'a%u' % fr if fr != 13 and fr != 1 else 'b1'
            if fr_name in with_protein:
                if fr_name == self.main.fracs_orderL[self.protein][0][0] or \
                    fr_name == self.main.fracs_orderU[self.protein][0][0]:
                    self.highest.append((scan_num, fr))
                else:
                    self.secondary.append((scan_num, fr))
            else:
                self.other.append((scan_num, fr))
        self.highest = sorted(self.highest, key = lambda sc: abs(self._scans[sc].deltart))
        self.secondary = sorted(self.secondary, key = lambda sc: abs(self._scans[sc].deltart))
        self.other = sorted(self.other, key = lambda sc: abs(self._scans[sc].deltart))
    
    def select_best_scan(self):
        self.best_scan = \
            self.highest[0] if len(self.highest) else \
            self.secondary[0] if len(self.secondary) else \
            self.other[0] if len(self.other) else \
            None
    
    def print_db_species(self):
        return ', '.join(
            map(
                lambda hg:
                    '%s' % (
                        hg \
                            if hg not in self.tbl['ms1fa'][self.oi] \
                            or not len(self.tbl['ms1fa'][self.oi][hg]) \
                            else \
                        ', '.join(
                            map(
                                lambda fa:
                                    '%s(%s)' % (hg, fa),
                                self.tbl['ms1fa'][self.oi][hg]
                            )
                        )
                    ),
                self.tbl['ms1hg'][self.oi]
            )
        ) \
        if len(self.tbl['ms1hg'][self.oi]) \
        else 'none'
    
    def reload(self, children = False):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
        
        if children:
            
            for sc in self._scans.values():
                
                sc.reload()
    
    def __str__(self):
        return ', '.join(
            map(
                lambda hgfas:
                    ', '.join(
                        map(
                            lambda fa:
                                '%s(%s)' % (hgfas[0], fa),
                            hgfas[1]
                        )
                    ),
                iteritems(self.fa)
            )
        )
    
    def get_header_div(self):
        return '\t\t<div class="ms2hdr">\n\t\t'\
            'MS2 scans of feature %.04f'\
            '\n\t\t\t|<span class="scansbutton morescans1"'\
                ' title="Show/hide scans from fractions with '\
                'highest protein concentration">scans+</span>\n'\
            '\n\t\t\t|<span class="scansbutton morescans2"'\
                ' title="Show/hide scans from other protein '\
                'containing fractions">scans++</span>\n'\
            '\n\t\t\t|<span class="scansbutton morescans3"'\
                ' title="Show/hide scans from non protein '\
                'containing fractions">scans+++</span>\n'\
            '\n\t\t\t|<span class="scansbutton morefrags"'\
                ' title="Show/hide fragments after 5 highest'\
                '">frags+</span>\n'\
            '\n\t\t\t|<span class="scansbutton remove"'\
                ' title="Remove scans of this feature'\
                '">remove</span>\n'\
            '\t\t</div>\n' % \
            self.tbl['mz'][self.i]
    
    def html_table(self):
        container = '\t<div id="%s" class="ms2tblcontainer">\n%s%s\n\t</div>'
        header = self.get_header_div()
        html = []
        if self.best_scan is not None:
            html.append(self._scans[self.best_scan].html_table())
        else:
            html.append('<div class="noscans">No scans '\
                'from fractions with highest protein concentration.</div>')
        for sc in sorted(self._scans.values(), key = lambda sc: abs(sc.deltart)):
            if sc.in_primary and sc.scan_id != self.best_scan:
                html.append(sc.html_table())
        for sc in sorted(self._scans.values(), key = lambda sc: abs(sc.deltart)):
            if not sc.in_primary and sc.scan_id != self.best_scan:
                html.append(sc.html_table())
        for sc in sorted(self._scans.values(), key = lambda sc: abs(sc.deltart)):
            if not sc.in_primary and not sc.in_secondary:
                html.append(sc.html_table())
        html = '\n'.join(html)
        return container % ('ms2c_%u_%u' % \
            (int(self.tbl['aaa'][self.i]), self.oi), header, html)
    
    def html_table_b64(self):
        return base64.encodestring(self.html_table()).replace('\n', '')
    
    def msg(self, text):
        if self.log:
            with open(self.main.ms2log, 'a') as f:
                f.write(text)
    
    def _any_scan(self, method, **kwargs):
        for i, sc in iteritems(self._scans):
            self.msg('\t\t:: Calling method %s() on scan #%u\n' % (method, i[0]))
            if getattr(sc, method)(**kwargs):
                return True
        return False
    
    def identify(self):
        for hg in self.classes:
            self.msg('\t>>> Attempting to identify %s in all scans\n' % (hg))
            if self._any_scan('is_%s' % hg.lower()):
                self.identities.add(hg)
                self.msg('\t<<< Result: identified as %s\n' % hg)
            else:
                self.msg('\t<<< Result: not %s\n' % hg)
    
    def identify2(self, num = 1):
        
        for scanid, scan in iteritems(self._scans):
            
            for hg in self.classes2:
                
                self.msg('\t>>> Attempting to identify %s in scan %u\n' %
                         (hg, scanid[0]))
                
                identified = False
                
                if hg not in self.identities2:
                    self.identities2[hg] = []
                
                method = '%s_%s_%u' % (hg.lower(), self.mode, num)
                
                if hasattr(scan, method):
                    
                    self.identities2[hg].append(getattr(scan, method)())
                    
                    identified = any(
                        map(
                            lambda i: i['score'] >= 5,
                            self.identities2[hg]
                        )
                    )
                
                if identified:
                    self.msg('\t<<< Result: identified as %s\n' % hg)
                else:
                    self.msg('\t<<< Result: not %s\n' % hg)
            
            if hasattr(scan, 'fa_co_2'):
                del scan.fa_co_2
            
            if hasattr(scan, 'fa_list'):
                scan.fa_list = None


class MS2Scan(object):
    """
    This class represents one MS2 scan and provides methods for its analysis.
    """
    
    def __init__(self, scan, scan_id, feature):
        
        self.scan = scan
        self.scan_id = scan_id
        self.feature = feature
        self.deltart = self.feature.deltart[self.scan_id]
        self.frac_id = self.scan_id[1]
        self.frac_name = self.feature.fracsi[self.frac_id]
        
        self.ms2_file = self.feature.main.ms2files\
            [self.feature.protein][self.feature.mode][self.frac_name]
        self.in_primary = self.frac_name in \
            self.feature.main.fracs_order[self.feature.protein]['prim']
        self.in_secondary = self.frac_name in \
            self.feature.main.fracs_order[self.feature.protein]['sec']
        self.i = self.feature.i
        self.tbl = self.feature.tbl
        self.insmax = self.scan[0,2]
        self.recc = re.compile(r'.*?([0-9]{1,2}):([0-9]).*')
        self.fa = {}
        self.fa1 = {}
        self._order = None
        self.sort_by_i()
        self.fa_list = None
        self.build_fa_list()
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def print_identities(self, fname = None):
        """
        Prints identities to standard output or file.
        """
        
        if fname is None:
            sys.stdout.write(self.identities_str())
        else:
            with open(fname, 'w') as fp:
                fp.write(self.identities_str())
    
    def identities_str(self, num = 1):
        """
        Returns table of all identification attemts as string.
        """
        
        result = ['=== Scan #%u (fraction %s) ===' % (
            self.scan_id[0], self.frac_name)]
        
        for hg in self.feature.classes2:
            
            method = '%s_%s_%u' % (
                hg.lower(), self.feature.mode, num
            )
            
            if not hasattr(self, method):
                continue
            
            idd = getattr(self, method)()
            
            result.append('%s\t%u\t%s' % (
                hg,
                idd['score'],
                ', '.join(idd['fattya'])
            ))
        
        return '%s\n' % '\n'.join(result)
    
    def print_scan(self):
        """
        Prints the list of fragments as an annotated table.
        """
        
        ms1mz = self.tbl['mz'][self.i]
        header = '\tFrag. m/z\tIntensity\tIdentity%sNL mass\n'\
            '\t%s\n' % (' ' * 26, '=' * 73)
        table = '\n\t'.join(
            map(
                lambda sc:
                    '%9.4f\t%10.2f\t%s%s%9.4f' % \
                        tuple(list(sc[[1, 2, 7]]) + \
                            [' ' * (32 - len(sc[7])), ms1mz - sc[1]]),
                self.scan
            )
        )
        
        self.feature.msg('\tScan %u (fraction %s(#%u); %s %s; '\
            'intensity = %.01f (%.02f%%)):\n\n%s\t%s\n\n' % \
            (self.scan_id[0],
             self.frac_name,
             self.frac_id,
             'contains' \
                if self.feature.ifracs[self.frac_name][1] \
                else 'does not contain',
             self.feature.protein,
             self.tbl['fe'][self.i, self.frac_id] \
                 if self.frac_id < self.tbl['fe'].shape[1] else np.nan,
             (self.tbl['fe'][self.i, self.frac_id] \
                 if self.frac_id < self.tbl['fe'].shape[1] else np.nan) / \
                 np.nanmax(self.tbl['fe'][self.i, :]) * 100.0,
             header,
             table)
        )
    
    def html_table(self):
        table = '\t\t<table id="%s" class="scantbl %s">\n%s\n\t\t</table>\n'
        th = '\t\t\t\t<th>\n\t\t\t\t\t%s\n\t\t\t\t</th>\n'
        ttl = '\t\t\t<tr class="%s">\n\t\t\t\t<th colspan="4">\n\t\t\t\t\t%s'\
            '\n\t\t\t\t</th>\n\t\t\t</tr>\n'
        tr = '\t\t\t<tr class="%s">\n%s\n\t\t\t</tr>\n'
        td = '\t\t\t\t<td>\n\t\t\t\t\t%s\n\t\t\t\t</td>\n'
        ms1mz = self.tbl['mz'][self.i]
        rows = ttl % (
            'scantitle',
            'Scan %u (%s, %s; '\
            'intensity = %.01f (%.02f%%); dRT = %.03f min)' % (
                self.scan_id[0],
                self.frac_name,
                'the highest fraction' if self.in_primary \
                    else 'not the highest, but contains %s' % \
                        self.feature.protein if self.in_secondary \
                    else 'does not contain %s' % \
                        self.feature.protein,
                self.tbl['fe'][self.i, self.frac_id] \
                    if fri < self.tbl['fe'].shape[1] else np.nan,
                (self.tbl['fe'][self.i, self.frac_id] \
                    if fri < self.tbl['fe'].shape[1] else np.nan) / \
                    np.nanmax(self.tbl['fe'][self.i, :]) * 100.0,
                self.deltart
            )
        )
        rows += tr % (
            'scanhdr',
            ''.join(
                map(
                    lambda cname:
                        th % cname,
                    ['Frag m/z', 'Intensity', 'Identity', 'NL mass']
                )
            )
        )
        for rn, row in enumerate(self.scan):
            rows += tr % (
                'fragrow %s' % ('first5' if rn < 5 else 'after5'),
                ''.join([
                    td % ('%.04f' % row[1]),
                    td % ('%.02f' % row[2]),
                    td % row[7],
                    td % ('%.04f' % (ms1mz - row[1]))
                ])
        )
        return table % ('%u_%u_%u' % (
            self.tbl['i'][self.i], self.scan_id[0], self.scan_id[1]),
            'best' if self.scan_id == self.feature.best_scan \
                else 'primary' if self.in_primary \
                else 'secondary' if self.in_secondary \
                else 'noprotein',
            rows
        )
    
    def get_by_rank(self, rank = 1, min_mz = 0.0):
        this_rank = 0
        return_next = False
        prev_mz = 0.0
        intensity = ''
        ids = []
        for r in self.scan:
            if r[1] < min_mz:
                continue
            if abs(r[1] - prev_mz) > 0.0001:
                prev_mz = r[1]
                this_rank += 1
            if this_rank == rank:
                return_next = True
                intensity = '%.04f(%u)' % (r[1], r[2])
                ids.append('%s (%.03f)' % (r[7], r[1]))
            elif this_rank != rank and return_next:
                return intensity, '; '.join(ids)
        return '', ''
    
    def full_list_str(self):
        result = []
        prev_mz = self.scan[0,1]
        intensity = self.scan[0,2]
        names = set([])
        for i, r in enumerate(self.scan):
            if abs(r[1] - prev_mz) > 0.0001:
                if len(names) == 1 and  list(names)[0] == 'unknown':
                    result.append('%s (%.03f) (%u)' % ('/'.join(sorted(list(names))), r[1], intensity))
                else:
                    result.append('%s (%u)' % ('/'.join(sorted(list(names))), intensity))
                names = set([])
                intensity = r[2]
                prev_mz = r[1]
            names.add(r[7])
        result.append('%s (%u)' % ('/'.join(sorted(list(names))), intensity))
        return '; '.join(result)
    
    def most_abundant_mz(self):
        result = self.scan[0,1]
        self.feature.msg('\t\t  -- Most abundant m/z is %.03f\n' % result)
        return result
    
    def mz_match(self, mz_detected, mz):
        return abs(mz_detected - mz) <= self.feature.main.ms2_tlr
    
    def sort_by_mz(self):
        """
        Sorts the scan array by m/z increasing.
        """
        self._order = self._order[self.scan[:,1].argsort()]
        self.scan = self.scan[self.scan[:,1].argsort(),:]
    
    def sort_by_i(self, return_order = False):
        """
        Sorts the scan array by intensity decreasing.
        """
        if self._order is None:
            order = self.scan[:,2].argsort()[::-1]
            self.scan = self.scan[order,:]
            self._order = np.array(xrange(self.scan.shape[0]), dtype = np.int)
        else:
            order = self._order.argsort()
            self.scan = self.scan[order,:]
            self._order = self._order[order]
        if return_order:
            return order
    
    def mz_lookup(self, mz):
        """
        Returns the index of the closest m/z value
        detected in the scan if it is within the
        range of tolerance, otherwise None.
        """
        du = 999.0
        dl = 999.0
        self.sort_by_mz()
        ui = self.scan[:,1].searchsorted(mz)
        if ui < self.scan.shape[0]:
            du = self.scan[ui,1] - mz
        if ui > 0:
            dl = mz - self.scan[ui - 1,1]
        i = ui if du < dl else ui - 1
        i = i if self.mz_match(self.scan[i,1], mz) else None
        sort = self.sort_by_i(return_order = True)
        if i is not None:
            i = np.where(sort == i)[0][0]
        return i
    
    def has_mz(self, mz):
        result = self.mz_lookup(mz) is not None
        self.feature.msg('\t\t  -- m/z %.03f occures in this scan? -- %s\n' % \
            (mz, str(result)))
        return result
    
    def has_nl(self, nl):
        result = self.has_mz(self.feature.tbl['mz'][self.feature.i] - nl)
        self.feature.msg('\t\t  -- neutral loss of %.03f occures in '\
            'this scan? Looked up m/z %.03f - %.03f = %.03f -- %s\n' % \
            (nl, self.feature.tbl['mz'][self.feature.i], nl,
             self.feature.tbl['mz'][self.feature.i] - nl, str(result)))
        return result
    
    def most_abundant_mz_is(self, mz):
        result = self.mz_match(self.most_abundant_mz(), mz)
        self.feature.msg('\t\t  -- m/z %.03f is the most abundant? -- %s\n' % \
            (mz, str(result)))
        return result
    
    def mz_among_most_abundant(self, mz, n = 2):
        """
        Tells if an m/z is among the most aboundant `n` fragments
        in a spectrum.
        
        :param float mz: The m/z value.
        :param int n: The number of most abundant fragments considered.
        
        """
        
        result = False
        
        for i in xrange(min(n, self.scan.shape[0])):
            
            if self.mz_match(self.scan[i,1], mz):
                
                result = True
                break
        
        self.feature.msg('\t\t  -- m/z %.03f is among the %u most abundant? -- '\
            '%s\n' % (mz, n, str(result)))
        
        return result
    
    def nl_among_most_abundant(self, nl, n = 2):
        """
        Tells if a neutral loss corresponds to one of the
        most aboundant `n` fragments in a spectrum.
        
        :param float nl: The mass of the neutral loss.
        :param int n: The number of most abundant fragments considered.
        
        """
        
        result = False
        mz = self.feature.tbl['mz'][self.feature.i] - nl
        
        for i in xrange(min(n, self.scan.shape[0])):
            
            if self.mz_match(self.scan[i,1], mz):
                
                result = True
                break
        
        self.feature.msg('\t\t  -- neutral loss %.03f is among '\
            'the %u most abundant? -- '\
            '%s\n' % (nl, n, str(result)))
        
        return result
    
    def mz_percent_of_most_abundant(self, mz, percent = 80.0):
        """
        Tells if an m/z has at least certain percent of intensity
        compared to the most intensive fragment.
        
        :param float mz: The m/z value.
        :param float percent: The threshold in percent
                              of the highest intensity.
        
        """
        
        insmax = self.scan[0,2]
        result = False
        
        for frag in self.scan:
            
            if self.mz_match(frag[1], mz):
                
                result = True
                break
            
            if frag[2] < insmax * 100.0 / percent:
                result = False
                break
        
        self.feature.msg('\t\t  -- m/z %.03f has abundance at least %.01f %% of'\
            ' the highest abundance? -- %s\n' % \
            (mz, percent, str(result)))
        
        return result
    
    def fa_type_is(self, i, fa_type, sphingo = False, uns = None,
                   scan_index = True):
        """
        Tells if a fatty acid fragment is a specified type. The type
        should be a part of the string representation of the fragment,
        e.g. `-O]` for fragments with one oxygen loss.
        """
        
        ifa = None
        
        if not scan_index:
            ifa = i
            i   = self.fa_list[ifa][5]
        
        result = (
            (fa_type in self.scan[i,8] or fa_type in self.scan[i,7]) and
            (not sphingo or 'Sphingosine' in self.scan[i,7]) and
            (uns is None or ifa is None or self.fa_list[ifa][0][1] <= uns)
        )
        
        self.feature.msg('\t\t  -- Fragment #%u (%s, %s): fatty acid type '\
            'is %s?  -- %s\n' % \
                (i, self.scan[i,7], self.scan[i,8], fa_type, str(result)))
        
        return result
    
    def is_fa(self, i, sphingo = False):
        """
        Examines whether a fragment is fatty acid-like or not.
        In the labels of fatty acid fragments we always 
        """
        
        result = 'FA' in self.scan[i,7] or 'Lyso' in self.scan[i,7] or \
            (sphingo and 'Sphi' in self.scan[i,7])
        
        self.feature.msg('\t\t  -- Fragment #%u (%s): is fatty acid? '\
            '-- %s\n' % (i, self.scan[i,7], str(result)))
        
        return result
    
    def most_abundant_fa(self, fa_type, head = 1, sphingo = False):
        """
        Returns `True` if there is a fatty acid among the most abundant
        fragments and it is of the defined type; `False` if there is no
        fatty acid, or it is different type.
        
        :param str fa_type: The type of the fatty acid fragment ion.
        :param int head: The number of most abundant fragments considered.
        :param bool sphingo: Look for a sphingolipid backbone.
        """
        
        result = False
        
        for i in xrange(self.scan.shape[0]):
            
            if i == head:
                break
            
            if self.is_fa(i, sphingo = sphingo):
                result = self.fa_type_is(i, fa_type, sphingo = sphingo)
        self.feature.msg('\t\t  -- Having fatty acid %s among %u most abundant '\
            'features? -- %s\n' % (fa_type, head, str(result)))
        
        return result
    
    def fa_among_most_abundant(self, fa_type, n = 2,
                               min_mass = None, sphingo = False,
                               uns = None):
        """
        Returns `True` if there is one of the defined type of fatty acid
        fragments among the given number of most abundant fragments, and
        it has a mass greater than the given threhold.
        """
        
        self.build_fa_list()
        result = False
        
        for i, fa in enumerate(self.fa_list):
            
            if not sphingo or fa[3] and (
                    min_mass is None or
                    self.scan[fa[5],1] >= min_mass
                ):
                
                if min_mass is not None:
                    
                    self.feature.msg('\t\t\t-- Fragment #%u having mass larger '\
                        'than %.01f\n' % (i, min_mass))
                
                if self.fa_type_is(i, fa_type, sphingo = sphingo,
                                   uns = uns, scan_index = False):
                    result = True
            
            if i == n:
                break
            
            elif min_mass is not None:
                self.feature.msg('\t\t\t-- Fragment #%u having mass lower '\
                        'than %.01f\n' % (i, min_mass))
        
        self.feature.msg('\t\t  -- Having fatty acid fragment %s among %u most '\
            'abundant -- %s\n' % (fa_type, n, str(result)))
        
        return result
    
    def fa_percent_of_most_abundant(self, fa_type, percent = 80.0, sphingo = False):
        for i in xrange(self.scan.shape[0]):
            if self.is_fa(i, sphingo = sphingo):
                if self.fa_type_is(i, fa_type, sphingo = sphingo):
                    return True
            if self.scan[i,2] < self.insmax * 100.0 / percent:
                return False
        return False
    
    def mz_most_abundant_fold(self, mz, fold):
        """
        Tells if an m/z is the most abundant fragment
        and it has at least a certain
        fold higher intensity than any other fragment.
        
        :param float mz: The m/z value.
        :param float fold: The m/z must be this times higher than any other.
        """
        
        result = False
        if self.most_abundant_mz_is(mz):
            result = self.scan.shape[0] == 1 or \
                self.scan[1,2] * fold <= self.scan[0,2]
        self.feature.msg('\t\t  -- m/z %.03f is at least %u times higher than '\
            'any other? -- %s\n' % (mz, fold, str(result)))
        return result
    
    def sum_cc_is(self, cc1, cc2, cc):
        """
        Returns `True` if the sum of the 2 carbon counts and
        unsaturations is equal with the third one.
        
        :param tuple cc1: Carbon count and unsaturation 1.
        :param tuple cc2: Carbon count and unsaturation 2.
        :param str cc: Expected total carbon count and unsaturation.
        """
        
        return self.cc2str(self.sum_cc([cc1, cc2])) == cc
    
    def cer_fa_test(self, frag1, frag2):
        return \
            self.fa_type_is(frag1[5], 'CerFA(') and \
            self.fa_type_is(frag2[5], 'CerSphi-N(') and \
            frag1[4] > frag2[4] * 2
    
    def fa_combinations3(self, hg, head = None):
        """
        Finds all combinations of 3 fatty acids which match the
        total carbon count and unsaturation resulted by database
        lookups of the MS1 precursor mass.
        This can be used for example at TAG.
        
        :param str hg: The short name of the headgroup, e.g. `TAG`.
        :param int head: If `None` or `numpy.inf` all fragment ions
                         will be considered, otherwise only the first
                         most aboundant until the number `head`.
        """
        
        result = set([])
        
        if hg in self.feature.ms1fa and len(self.feature.ms1fa[hg]):
            ccs = list(self.feature.ms1fa[hg])
        else:
            return result
        
        head = np.inf if head is None else head
        
        for cc in ccs:
            
            try:
                icc = self.cc2int(cc)
            except AttributeError:
                continue
            
            for frag0 in self.fa_list:
                
                if frag0[5] >= head:
                    break
                
                cc0 = frag0[0]
                
                cc12e = '%u:%u' % tuple(map(lambda x: x[0] - x[1],
                                            zip(*[icc, cc0])))
                
                cc12s = self.fa_combinations(cc12e, head = head, by_cc = True)
                
                for cc12 in cc12s:
                    
                    cc012 = self.ccs2int(cc12) + [cc0]
                    
                    if self.sum_cc(cc012) == icc:
                        
                        result.add(self.ccs2str(cc012))
        
        return result
    
    def fa_combinations_old(self, hg, sphingo = False,
                            head = None, by_cc = False):
        """
        Finds all combinations of 2 fatty acids which match the
        total carbon count and unsaturation resulted by database
        lookups of the MS1 precursor mass.
        Alternatively a carbon count and unsaturation can be provided
        if `by_cc` is set to `True`.
        
        :param str hg: Short name of the headgroup, e.g. `PC`; or cc:unsat e.g.
                       `32:1` if `by_cc` is `True`.
        :param bool sphingo: Assume sphingolipid.
        :param int head: If `None` the total fragment list used, if a number,
                         only the most intensive fragments accordingly.
        :param bool by_cc: Use the MS1 database identification to find out
                           the possible carbon counts and unsaturations for
                           the given headgroup, or a cc:uns provided and
                           search combinations accordingly.
        
        """
        
        result = set([])
        if hg in self.feature.ms1fa and len(self.feature.ms1fa[hg]):
            ccs = list(self.feature.ms1fa[hg])
        elif by_cc:
            ccs = [hg]
        else:
            return result
        
        head = np.inf if head is None else head
        
        self.build_fa_list()
        
        for cc in ccs:
            
            for frag1 in self.fa_list:
                
                for frag2 in self.fa_list:
                    
                    result.update(
                        self.get_fa_combinations(frag1, frag2, hg,
                                                cc, sphingo, head)
                    )
        
        return result
    
    def fa_combinations_preprocess(self, regenerate = False):
        """
        Generates a lookup table for all possible combinations of two
        fatty acids.
        """
        
        if not hasattr(self, 'fa_co_2') or regenerate:
            
            self.fa_co_2 = {}
            l = self.fa_list
            
            for i, j in itertools.combinations_with_replacement(
                xrange(len(self.fa_list)), 2):
                
                key = self.sum_cc([(l[i][0][0], l[i][0][1]),
                                   (l[j][0][0], l[j][0][1])])
                
                if key not in self.fa_co_2:
                    self.fa_co_2[key] = set([])
                
                self.fa_co_2[key].add((i, j))
    
    def fa_combinations(self, hg, sphingo = False,
                               head = None, by_cc = False):
        """
        Finds all combinations of 2 fatty acids which match the
        total carbon count and unsaturation resulted by database
        lookups of the MS1 precursor mass.
        Alternatively a carbon count and unsaturation can be provided
        if `by_cc` is set to `True`.
        
        This method does the same as `fa_combinations` but works with
        a preprocessed lookup table.
        
        :param str hg: Short name of the headgroup, e.g. `PC`; or cc:unsat e.g.
                       `32:1` if `by_cc` is `True`.
        :param bool sphingo: Assume sphingolipid.
        :param int head: If `None` the total fragment list used, if a number,
                         only the most intensive fragments accordingly.
        :param bool by_cc: Use the MS1 database identification to find out
                           the possible carbon counts and unsaturations for
                           the given headgroup, or a cc:uns provided and
                           search combinations accordingly.
        
        """
        
        result = set([])
        if hg in self.feature.ms1fa and len(self.feature.ms1fa[hg]):
            ccs = list(self.feature.ms1fa[hg])
        elif by_cc:
            ccs = [hg]
        else:
            return result
        
        head = np.inf if head is None else head
        
        self.build_fa_list()
        self.fa_combinations_preprocess()
        
        for cc in ccs:
            
            icc = self.cc2int(cc)
            
            if icc in self.fa_co_2:
                
                for i, j in self.fa_co_2[icc]:
                    
                    frag1 = self.fa_list[i]
                    frag2 = self.fa_list[j]
                    
                    result.update(
                        self.get_fa_combinations(frag1, frag2, hg,
                                                 cc, sphingo, head)
                    )
        
        return result
    
    def get_fa_combinations(self, frag1, frag2, hg, cc, sphingo, head):
        """
        Processes two fatty acid fragments to decide
        if their combination is valid.
        """
        
        result = set([])
        
        if frag1[5] >= head or frag2[5] >= head:
            return result
        
        if hg == 'Cer' and not self.cer_fa_test(frag1, frag2):
            # where not the 'CerFA' is the most intensive
            # those are clearly false
            return result
        
        if frag1[0][0] is not None and frag2[0][0] is not None and \
            (frag1[1] is None or hg in frag1[1]) and \
            (frag2[1] is None or hg in frag2[1]) and \
            (not sphingo or frag1[3] or frag2[3]):
            if self.sum_cc_is(frag1[0], frag2[0], cc):
                ether_1 = 'O-' if frag1[2] else ''
                ether_2 = 'O-' if frag2[2] else ''
                fa_1 = '%s%u:%u' % (ether_1, frag1[0][0], frag1[0][1])
                fa_2 = '%s%u:%u' % (ether_2, frag2[0][0], frag2[0][1])
                if frag1[3]:
                    fa_1 = 'd%s' % fa_1
                elif frag2[3]:
                    sph = 'd%s' % fa_2
                    fa_2 = fa_1
                    fa_1 = sph
                if not frag1[3] and not frag2[3]:
                    fa = tuple(sorted([fa_1, fa_2]))
                else:
                    fa = (fa_1, fa_2)
                result.add('%s/%s' % fa)
        
        return result
    
    def matching_fa_frags_of_type(self, hg, typ, sphingo = False,
        return_details = False):
        """
        Returns carbon counts of those fragments which are of the given type
        and have complement fatty acid fragment of any type.
        
        Details is a dict with carbon counts as keys
        and fragment names as values.
        """
        result = set([])
        details = {}
        if hg in self.feature.ms1fa and len(self.feature.ms1fa[hg]):
            for cc in self.feature.ms1fa[hg]:
                self.build_fa_list()
                for frag1 in self.fa_list:
                    for frag2 in self.fa_list:
                        if frag1[0][0] is not None and \
                            frag2[0][0] is not None and \
                            (frag1[1] is None or hg in frag1[1]) and \
                            (frag2[1] is None or hg in frag2[1]) and \
                            (not sphingo or frag1[3]):
                            if self.fa_type_is(frag1[5], typ) and \
                                self.sum_cc_is(frag1[0], frag2[0], cc):
                                result.add(frag1[0])
                                if return_details:
                                    if frag1[0] not in details:
                                        details[frag1[0]] = set([])
                                    details[frag1[0]].add(self.scan[frag2[5],7])
        if return_details:
            return (result, details)
        else:
            return result
    
    def cer_missing_fa(self, cer_hg):
        """
        Infers the fatty acid carbon count and unsaturation
        by subtracting the sphingoid backbone from the total.
        This works with Cer, CerP and HexCer.
        """
        
        result = set([])
        
        cer_ccs = set([])
        
        for frag in self.scan[:5]:
            
            if 'phingo' in frag[7]:
                
                cer_ccs.add(self.get_cc(frag[7]))
        
        if cer_hg in self.feature.ms1fa:
            
            for cc in self.feature.ms1fa[cer_hg]:
                
                for cer_cc in cer_ccs:
                    
                    cc = self.get_cc(cc)
                    carb = cc[0] - cer_cc[0]
                    unsat = cc[1] - cer_cc[1]
                    
                    result.add('d%u:%u/%u:%u' % (
                        cer_cc[0], cer_cc[1], carb, unsat))
        
        return result
    
    def cer_matching_fa(self, cer_fa):
        score = 0
        if 'Cer' in self.feature.ms1fa:
            cer_cc = self.get_cc(cer_fa)
            for cc in self.feature.ms1fa['Cer']:
                cc = self.get_cc(cc)
                carb = cc[0] - cer_cc[0]
                unsat = cc[1] - cer_cc[1] + 2
                if self.frag_name_present(
                    '[FA-alkyl(C%u:%u)-H]-' % (carb, unsat)):
                    score += 1
                carb = cc[0] - cer_cc[0] - 2
                unsat = cc[1] - cer_cc[1] + 1
                if self.frag_name_present(
                    '[FA-alkyl(C%u:%u)-H]-' % (carb, unsat)):
                    score += 1
        return score
    
    def build_fa_list(self, rebuild = False):
        """
        Returns list with elements:
            carbon count, headgroups (set or None),
            esther (False) or ether (True),
            sphingosine (True) or fatty acid (False),
            fragment intensity and row index
        """
        if self.fa_list is None or rebuild:
            self.fa_list = []
            for i, frag in enumerate(self.scan):
                if frag[7] != 'unknown' and self.is_fa(i, sphingo = True):
                    cc = self.get_cc(frag[7])
                    hgs = self.get_hg(frag[7])
                    is_ether = 'alk' in frag[7]
                    is_sphingo = 'Sphi' in frag[7]
                    self.fa_list.append([cc, hgs, is_ether, is_sphingo, frag[2], i])
    
    def get_hg(self, frag_name):
        hgfrags = self.feature.main.nHgfrags \
            if self.feature.mode == 'neg' \
            else self.feature.main.pHgfrags
        return hgfrags[frag_name] \
            if frag_name in hgfrags and \
                len(hgfrags[frag_name]) \
            else None
    
    def get_cc(self, fa):
        """
        Extracts carbon count from any string, for example fatty acid names.
        Recognizes the pattern [number:number].
        E.g. from `Phosphatidylcholine (36:1)` returns the tuple `(36, 1)`.
        To convert pure cc:uns strings, use the method `cc2int` instead,
        as that one is faster.
        
        :param str fa: Any string containing carbon count and unsaturation.
        
        """
        
        m = self.recc.match(fa)
        
        if m is not None:
            return tuple(map(int, m.groups()))
        
        return (None, None)
    
    def most_abundant_fa_cc(self, fa_type = None, head = 2):
        fa_cc = []
        for i, frag in enumerate(self.scan):
            if i == head:
                break
            if self.is_fa(i) and (
                    fa_type is None or
                    self.fa_type_is(i, fa_type)
                ):
                
                cc = self.get_cc(frag[7])
                if cc[0] is not None:
                    fa_cc.append((cc, frag[2]))
        
        return fa_cc
    
    def cc2str(self, cc):
        """
        Converts carbon count and unsaturation from tuple of integers
        to string. E.g. `(18, 1)` results `18:1`.
        
        :param tuple cc: Tuple of 2 integers representing carbon count
                         and unsaturation.
        """
        
        return '%u:%u' % cc
    
    def ccs2str(self, ccs):
        """
        Converts multiple carbon counts and unsaturations from tuples
        of integers format to string. E.g. `[(18, 1), (18, 0)]` results
        `18:1/18:0`.
        
        :param list ccs: List of tuples of integers.
        """
        
        return '/'.join(map(self.cc2str, sorted(ccs)))
    
    def cc2int(self, cc):
        """
        Converts carbon count and unsaturation from string format to
        tuple of integers. E.g. `18:1` results `(18, 1)`.
        
        :param str cc: String representing carbon count and unsaturation
                       separated by colon.
        """
        
        return tuple(map(int, cc.split(':')))
    
    def ccs2int(self, ccs):
        """
        Converts a string of multiple carbon counts and unsaturations
        to a list of tuples of integers.
        
        :param str ccs: Multiple carbon counts and unsaturations in string
                        representation, e.g. `18:1/16:0`.
        
        """
        
        return list(map(self.cc2int, ccs.split('/')))
    
    def sum_cc_str(self, ccs):
        """
        Returns the sum of multiple carbon counts and unsaturations.
        Accepts string format and results string format.
        
        :param str ccs: Multiple carbon counts and unsaturations in string
                        representation, e.g. `18:1/16:0`.
        """
        
        return self.sum_cc(self.ccs2int(ccs))
    
    def sum_cc(self, ccs):
        """
        Adds numeric carbon counts and unsaturations.
        Accepts a list of tuples of integers, returns
        a tuple of integers.
        
        :param list ccs: A list with tuples of integers,
                         e.g. `[(14, 1), (16, 0)]`.
        """
        
        return (
            tuple(
                reduce(
                    lambda cu1, cu2:
                        # here `cu`: carbon count and unsaturation
                        (cu1[0] + cu2[0], cu1[1] + cu2[1]),
                    ccs
                )
            )
        )
    
    def sum_cc2(self, ccs):
        """
        Returns the total carbon count and unsaturation in tuple
        format from a list of tuples where the first element of
        the tuple is another tuple with the cc and uns, and the
        second is the intensity, which is discarded here.
        
        :param list ccs: List of the format described above. E.g.
                         `[((16, 0), 1000.0), ((18, 1), 722)]`, this
                         results `(34, 1)`.
        """
        
        return self.sum_cc(map(lambda cci: cci[0], ccs))
    
    def sum_cc2str(self, ccs):
        """
        Returns the total carbon count and unsaturation in string
        format from a list of tuples where the first element of
        the tuple is another tuple with the cc and uns, and the
        second is the intensity, which is discarded here.
        
        :param list ccs: List of the format described above. E.g.
                         `[((16, 0), 1000.0), ((18, 1), 722)]`, this
                         results `34:1`.
        """
        
        return self.cc2str(self.sum_cc2(ccs))
    
    def add_fa1(self, fa, hg):
        if hg not in self.fa1:
            self.fa1[hg] = set([])
            self.fa1[hg].add(
                tuple(
                    map(
                        lambda fai:
                            fai[0],
                        fa
                    )
                )
            )
            fastr = ', '.join(
                    map(
                        lambda fai:
                            self.cc2str(fai[0]),
                        fa
                    )
                )
            self.feature.msg('\t\t  -- Adding fatty acids %s at headgroup '\
                '%s\n' % (fastr, hg))
    
    def fa_ccs_agree_ms1(self, hg, fa_type = None, head = 2):
        fa_cc = self.most_abundant_fa_cc(fa_type = fa_type, head = head)
        if len(fa_cc) > 0:
            cc = self.sum_cc2str([fa_cc[0]] * 2)
            agr = self.fa_cc_agrees_ms1(cc, hg)
            if agr:
                self.add_fa1(fa_cc[:1], hg)
            if len(fa_cc) > 1:
                cc = self.sum_cc2str(fa_cc[:2])
                agr = self.fa_cc_agrees_ms1(cc, hg)
                if agr:
                    self.add_fa1(fa_cc[:2], hg)
        return hg in self.fa
    
    def fa_cc_agrees_ms1(self, cc, hg):
        result = False
        if hg in self.feature.ms1fa and cc in self.feature.ms1fa[hg]:
            if hg not in self.feature.fa:
                self.feature.fa[hg] = set([])
            if hg not in self.fa:
                self.fa[hg] = set([])
            self.feature.fa[hg].add(cc)
            self.fa[hg].add(cc)
            result = True
        self.feature.msg('\t\t  -- Carbon count from MS2: %s; from databases '\
            'lookup: %s -- Any of these matches: %s\n' % \
                (
                    cc,
                    str(self.feature.ms1fa[hg]) \
                        if hg in self.feature.ms1fa else '-',
                    str(result))
                )
        return result
    
    def frag_name_present(self, name):
        return name in self.scan[:,7]
    
    #### New methods
    
    def cer1p_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Ceramide-1-phosphate.
        """
        
        score = 0
        fattya = set([])
        if self.most_abundant_mz_is(78.95905658):
            score += 5
            if self.has_mz(96.96962158):
                score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def hexcer_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Hexosyl-Ceramide.
        """
        
        score = 0
        fattya = set([])
        
        if all(map(lambda mz: self.mz_among_most_abundant(mz, n = 10),
                   # these are 3 fragments found at GLTP
                   [71.0115000, 89.0220000, 101.021900])):
            
            score += 5
        
        return {'score': score, 'fattya': fattya}
    
    def hexcer_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a Hexosyl-Ceramide.
        """
        
        score = 0
        fattya = set([])
        
        hexfrags = sum(map(lambda nl: self.nl_among_most_abundant(nl, n = 15),
                           [198.073955, 180.06339, 162.052825]))
        
        if hexfrags:
            score += hexfrags + 4
        
        if score:
            
            fattya.update(self.cer_missing_fa('HexCer'))
        
        return {'score': score, 'fattya': fattya}
    
    def cer1p_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a Ceramide-1-phosphate.
        Specimen: GLTPD1 + 728.59
        """
        
        score = 0
        fattya = set([])
        
        if self.fa_among_most_abundant('-H2O-H2O+]+', n = 3, sphingo = True):
            score += 4
            
            if any(map(self.has_mz,
                       # these present at Cer too
                       # a specific difference needed!
                       [82.0651257, 107.072951, 135.104251, 149.119901])):
                score += 1
            
            fattya.update(self.cer_missing_fa('Cer1P'))
        
        return {'score': score, 'fattya': fattya}
    
    def dag_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a DAG.
        Specimen: SEC14L2 + 584.52; Enric: BNIP2 + 770.67
        """
        
        score = 0
        fattya = set([])
        
        if(self.fa_combinations('DAG', head = 10)):
            score += 4
            
            if(self.fa_combinations('DAG', head = 6)):
                
                score += 2
            
            fattya.update(self.fa_combinations('DAG'))
        
        return {'score': score, 'fattya': fattya}
    
    def dag_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a DAG.
        This method is not ready, does nothing.
        """
        score = 0
        fattya = set([])
        
        return {'score': score, 'fattya': fattya}
    
    def tag_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a TAG.
        This method is not ready, does nothing.
        """
        score = 0
        fattya = set([])
        
        return {'score': score, 'fattya': fattya}
    
    def tag_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a TAG.
        Specimen: STARD11 + 818.7187
        """
        score = 0
        fattya = set([])
        
        fattya.update(self.fa_combinations3('TAG'))
        
        if fattya:
            score += 5
        
        return {'score': score, 'fattya': fattya}
    
    def pi_pos_1(self):
        """
        Examines if a negative MS2 spectrum is Phosphatidylinositol.
        Specimen: SEC14L2 + 906.60 and 882.6
        """
        
        score = 0
        fattya = set([])
        
        fattya.update(self.fa_combinations('PI'))
        if fattya:
            score += 1
            if self.has_nl(259.021894):
                score += 4
            if self.has_nl(277.056272):
                score += 4
        
        return {'score': score, 'fattya': fattya}
    
    def ps_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a Phosphatidylserine.
        Specimen: BPI + 790.56
        """
        
        score = 0
        fattya = set([])
        
        if self.nl_among_most_abundant(185.008927, 1):
            score += 5
            
            fattya.update(self.fa_combinations('PS'))
        
        return {'score': score, 'fattya': fattya}
    
    def bmp_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum
        is a Bismonoacylglycerophosphate.
        Specimen: BPIFB2 + 792.57
        """
        
        score = 0
        fattya = set([])
        
        if self.fa_among_most_abundant('+G(', 3):
            fattya.update(self.fa_combinations('BMP'))
            if fattya:
                score += 4
        
        return {'score': score, 'fattya': fattya}
    
    def pg_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum
        is a Phosphatidylglycerol.
        """
        
        score = 0
        fattya = set([])
        
        if self.nl_among_most_abundant(189.0402, 1):
            
            score += 5
            
            fattya.update(self.fa_combinations('PG'))
            if fattya:
                score += 4
        
        return {'score': score, 'fattya': fattya}
    
    def va_pos_1(self):
        """
        Examines if a positive MS2 spectrum is vitamin A (retinol).
        Specimen: RBP1 + 269.2245 and RBP4 + 269.2245
        """
        
        score = 0
        fattya = set([])
        
        if self.mz_among_most_abundant(269.224, 3):
            score += 5
            score += sum(map(self.has_mz, [213.165, 145.1027, 157.1028]))
        
        return {'score': score, 'fattya': fattya}
    
    def bmp_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum
        is a Bismonoacylglycerophosphate.
        The result will be similar to `pg_neg1`, as in negative
        mode we do not know a way to distinguish these species.
        It only depends on the headgroup fragment 152.996, as
        this might have higher intensity at BMP.
        Specimen: GM2A - 799.54 and BPIFB2 - 773.5258 (latter might be BMP)
        """
        
        result = self.pg_neg_1()
        
        if self.mz_among_most_abundant(152.9958366, 5):
            result['score'] += 3
        else:
            result['score'] -= 3
        
        return result
    
    #### End: new methods
    
    def pe_neg_1(self):
        score = 0
        fattya = set([])
        if self.is_fa(0) and self.fa_type_is(0, '-H]-') and self.has_mz(140.0118206):
            score += 5
            fattya = self.fa_combinations('PE')
            if self.has_mz(196.0380330):
                score +=1
            fa_h_ccs = self.matching_fa_frags_of_type('PE', '-H]-')
            for fa_h_cc in fa_h_ccs:
                for fa_other in [
                    '[Lyso-PE(C%u:%u)-]-',
                    '[Lyso-PE-alkyl(C%u:%u)-H2O]-',
                    '[Lyso-PE-alkyl(C%u:%u)-]-',
                    '[FA(C%u:%u)-H-CO2]-']:
                    if self.frag_name_present(fa_other % fa_h_cc):
                        score += 1
        return {'score': score, 'fattya': fattya}
    
    def pc_neg_1(self):
        score = 0
        fattya = set([])
        if self.is_fa(0) and self.fa_type_is(0, '-H]-') and self.has_mz(168.0431206):
            score += 5
            fattya = self.fa_combinations('PC')
            fa_h_ccs = self.matching_fa_frags_of_type('PC', '-H]-')
            for fa_h_cc in fa_h_ccs:
                if self.frag_name_present('[Lyso-PC(c%u:%u)-]-' % fa_h_cc):
                    score += 1
        return {'score': score, 'fattya': fattya}
    
    def pi_neg_1(self):
        """
        Examines if a negative MS2 spectrum is Phosphatidylinositol.
        Specimen: GM2A - 835.52
        """
        
        score = 0
        fattya = set([])
        
        if self.has_mz(241.0118779) and self.has_mz(152.9958366) and \
            self.has_mz(78.95905658):
            
            score += 5
            fattya = self.fa_combinations('PI')
            for hgfrag_mz in [96.96962158, 259.0224425, 297.0380926]:
                if self.has_mz(hgfrag_mz):
                    score += 1
            fa_h_ccs = self.matching_fa_frags_of_type('PI', '-H]-')
            for fa_h_cc in fa_h_ccs:
                for fa_other in [
                    '[Lyso-PI(C%u:%u)-]-',
                    '[Lyso-PI(C%u:%u)-H2O]-]']:
                    if self.frag_name_present(fa_other % fa_h_cc):
                        score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def ps_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Phosphatidylserine.
        Specimen: ORP9 - 788.54
        """
        
        score = 0
        fattya = set([])
        
        if self.is_fa(0) and self.fa_type_is(0, '-H]-') and \
            self.has_mz(152.9958366):
            
            score += 5
            fattya = self.fa_combinations('PS')
            
            if self.has_mz(87.03202840):
                score += 1
            
            fa_h_ccs = self.matching_fa_frags_of_type('PS', '-H]-')
            
            for fa_h_cc in fa_h_ccs:
                
                for fa_other in [
                    '[Lyso-PS(C%u:%u)-]-',
                    '[Lyso-PA(C%u:%u)-]-']:
                    
                    if self.frag_name_present(fa_other % fa_h_cc):
                        score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def pg_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is Phosphatidylglycerol.
        The result will be similar to `pg_neg1`, as in negative
        mode we do not know a way to distinguish these species.
        It only depends on the headgroup fragment 152.996, as
        this might have higher intensity at BMP.
        Specimen: GM2A - 799.54 and BPIFB2 - 773.5258 (latter might be BMP)
        """
        
        score = 0
        fattya = set([])
        
        if self.is_fa(0) and self.fa_type_is(0, '-H]-') and \
            self.has_mz(152.9958366):
            
            score += 5
            
            if self.mz_among_most_abundant(152.9958366, 5):
                score -= 3
            
            fattya = self.fa_combinations('PG')
            
            if self.has_mz(171.0064016):
                score += 1
            
            fa_h_ccs = self.matching_fa_frags_of_type('PG', '-H]-')
            
            for fa_h_cc in fa_h_ccs:
                
                for fa_other in [
                    'Lyso-PG(C%u:%u)-]-',
                    'Lyso-PG(C%u:%u)-H2O]-']:
                    
                    if self.frag_name_present(fa_other % fa_h_cc):
                        score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def sm_neg_1(self):
        score = 0
        fattya = set([])
        if self.mz_among_most_abundant(168.0431206) and self.has_nl(60.02113):
            score += 5
        return {'score': score, 'fattya': fattya}
    
    def cer_neg_1(self):
        score = 0
        fattya = set([])
        if self.fa_among_most_abundant('CerFA', n = 2):
            score += 5
            fattya = self.fa_combinations('Cer', sphingo = True)
            fa_h_ccs = self.matching_fa_frags_of_type('Cer', 'CerFA(')
            for fa_h_cc in fa_h_ccs:
                for fa_other in [
                    '[CerFA-N(C%u:%u)-]-',
                    '[CerFA-C2N(C%u:%u)-]-']:
                    if self.frag_name_present(fa_other % fa_h_cc):
                        score += 1
        return {'score': score, 'fattya': fattya}
    
    def cerp_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a Ceramide-1-phosphate.
        """
        
        score = 0
        fattya = set([])
        if self.most_abundant_mz_is(78.95905658):
            score += 5
            if self.has_mz(96.96962158):
                score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def pc_pos_1(self):
        score = 0
        fattya = set([])
        
        if self.most_abundant_mz_is(184.073323) and self.has_mz(86.096425):
            score += 5
            fattya = self.fa_combinations('PC')
            if self.has_mz(104.106990):
                score += 1
            if self.has_mz(124.999822):
                score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def sm_pos_1(self):
        score = 0
        fattya = set([])
        
        if all(
            map(
                lambda mz:
                    self.has_mz(mz),
                [60.080776, 86.096425, 104.106990,
                    124.999822, 184.073323]
            )
        ):
            score += 5
            if self.has_mz(58.0651):
                score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def fa_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a fatty acid.
        Here we only check if the most abundant fragment is the
        fatty acid itself.
        """
        score = 0
        fattya = set([])
        
        self.build_fa_list()
        
        if self.is_fa(0):
            
            if 'FA' in self.feature.ms1fa:
                
                for cc in self.feature.ms1fa['FA']:
                    
                    if self.cc2int(cc)  == self.fa_list[0][0]:
                        
                        score += 5
                        fattya.add(cc)
        
        return {'score': score, 'fattya': fattya}
    
    def fa_neg_1(self):
        """
        Examines if a negative mode MS2 spectrum is a fatty acid.
        Here we only check if the most abundant fragment is the
        fatty acid itself.
        """
        
        # these are the same
        return self.fa_pos_1()
    
    def cerp_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a Ceramide-1-phosphate.
        This method is not ready, does nothing.
        """
        
        score = 0
        fattya = set([])
        
        if self.fa_among_most_abundant('-H2O-H2O+]+', n = 3, sphingo = True):
            score += 1
        
        return {'score': score, 'fattya': fattya}
    
    def pe_pos_1(self):
        score = 0
        fattya = set([])
        if self.nl_among_most_abundant(141.019097, 1):
            score += 5
            fattya = self.fa_combinations('PE')
        return {'score': score, 'fattya': fattya}
    
    def cer_pos_1(self):
        """
        Examines if a positive mode MS2 spectrum is a Ceramide.
        """
        
        score = 0
        fattya = set([])
        
        if 'Cer' not in self.feature.ms1fa:
            ms1uns = None
        else:
            # larger unsaturation than the whole molecule
            # does not make sense
            ms1uns = max(map(lambda _cc: self.cc2int(_cc)[1],
                             self.feature.ms1fa['Cer']))
        
        if self.fa_among_most_abundant('-H2O-H2O+]+', n = 10,
                                       sphingo = True, uns = ms1uns):
            score += 5
            fattya = self.fa_combinations('Cer', sphingo = True)
            sph_ccs, fa_frags = self.matching_fa_frags_of_type('Cer',
                '-H2O-H2O+]+', sphingo = True, return_details = True)
            for cc, fa_frag_names in iteritems(fa_frags):
                for fa_frag_name in fa_frag_names:
                    if '+H]+' in fa_frag_name:
                        score += 1
                    if '-O]+' in fa_frag_name:
                        score += 1
                    if 'NL' in fa_frag_name:
                        score += 1
            for sph_cc in sph_ccs:
                for fa_other in [
                    '[Sphingosine(C%u:%u)-C-H2O-H2O+]+',
                    '[Sphingosine(C%u:%u)-H2O+]+']:
                    if self.frag_name_present(fa_other % sph_cc):
                        score += 1
            if not len(
                list(
                    filter(
                        lambda mz:
                            self.has_mz(mz),
                        [58.065126, 104.106990, 124.999822, 184.073323]
                    )
                )
            ):
                score += 1
            score += len(
                list(
                    filter(
                        lambda mz:
                            self.has_mz(mz),
                        [60.0443902, 70.0651257, 82.0651257, 96.0807757,
                        107.072951, 121.088601, 135.104251, 149.119901]
                    )
                )
            )
        return {'score': score, 'fattya': fattya}
    
    def vd_pos_1(self):
        score = 0
        fattya = set([])
        return {'score': score, 'fattya': fattya}
    
    def is_pe(self):
        if self.feature.mode == 'pos':
            return self.pa_pe_ps_pg_pos('PE')
        else:
            return self.pe_pc_pg_neg('PE')
    
    def is_pc(self):
        if self.feature.mode == 'pos':
            return self.pc_pos('PC')
        else:
            return self.pe_pc_pg_neg('PC')
    
    def is_pa(self):
        if self.feature.mode == 'pos':
            return self.pa_pe_ps_pg_pos('PA')
        else:
            return self.pa_ps_neg('PA')
    
    def is_ps(self):
        if self.feature.mode == 'pos':
            return self.pa_pe_ps_pg_pos('PS')
        else:
            return self.pa_ps_neg('PS')
    
    def is_pg(self):
        if self.feature.mode == 'pos':
            return self.pa_pe_ps_pg_pos('PG')
        else:
            return self.pe_pc_pg_neg('PG')
    
    def pa_pe_ps_pg_pos(self, hg):
        return self.mz_among_most_abundant(141.0191) \
            and self.fa_among_most_abundant('-O]+', min_mass = 140.0) \
            and self.fa_ccs_agree_ms1(hg, '-O]+')
    
    def pa_ps_neg(self, hg):
        return self.has_mz(152.9958366) and self.has_mz(78.95905658) \
            and self.most_abundant_fa('-H]-') \
            and self.fa_ccs_agree_ms1(hg, '-H]-')
    
    def pe_pc_pg_neg(self, hg):
        return self.most_abundant_fa('-H]-') \
            and self.fa_ccs_agree_ms1(hg, '-H]-')
    
    def pc_pos(self, hg):
        return self.mz_most_abundant_fold(184.0733, 3) \
            and self.fa_ccs_agree_ms1(hg, head = 4)
