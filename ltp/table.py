#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `ltp` python module
#
#  Copyright (c) 2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#


import time
import imp
import openpyxl
import re

from lipyd import moldb
from lipyd import ms2
from lipyd import mgf
from lipyd import settings

resheet = re.compile(r'([A-z0-9]+)_(positive|negative)_?(best)?')
remgf   = re.compile(r'([A-Z0-9]+)_(pos|neg)_([A-Z])([0-9]{1,2})\.mgf')

class TableBase(object):
    """
    Opens an xlsx file for editing and saves a copy
    with preserving formatting.
    """
    
    def __init__(self, infile, outfile = None):
        
        self.infile  = infile
        self.outfile = outfile or '%s__%s__.xlsx' % (
            infile,
            time.strftime('%Y.%m.%d_%H.%M'),
        )
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def main(self):
        
        self.open()
        self.process_sheets()
        self.write()
    
    def open(self):
        
        self.xls_original = openpyxl.load_workbook(self.infile)
    
    def process_sheets(self):
        
        for sheet in self.itersheets():
            
            self.process_sheet()
    
    def process_sheet(self):
        """
        To be overridden by child classes.
        """
        
        pass
    
    def itersheets(self):
        
        for sheet in self.xls_original.sheetnames:
            
            self.set_sheet(name)
            
            yield self.sheet
    
    def set_sheet(self, name):
        
        self.sheet = self.xls_original[name]
    
    def insert_col(self, idx, values, title = None):
        
        idx = self.numof_cols() + 1 if idx == 'END' else idx
        
        self.sheet.insert_cols(idx)
        
        ii = 1
        
        if title is not None:
            
            ii = 2
            self.sheet[1][idx].value = title
        
        for i, val in enumerate(values):
            
            self.sheet[i + ii][idx].value = val
    
    def write(self):
        
        self.xls_original.save(self.outfile)
    
    def numof_cols(self):
        
        try:
            
            return len(next(self.sheet.rows))
        
        except StopIteration:
            
            return 0
    
    @staticmethod
    def col_letter(letter):
        
        return openpyxl.utils.column_index_from_string(letter)


class LtpTable(TableBase):
    
    def __init__(self, infile, outfile = None, mgfdir = None):
        
        TableBase.__init__(self, infile, outfile)
        
        self.prot_frac = prot_frac
        self.mgfdir = mgfdir
    
    def process_sheet(self):
        
        annot = self.resheet.search(self.sheet.title)
        
        if not annot:
            
            return
        
        self.ionmode = annot[0]
        self.protein = annot[1]
        self.best    = annot[2] is not None
        
        self.get_mgf()
        
        for idx, title, value in self.new_columns():
            
            self.insert_col(idx, values, title)
    
    def new_columns(self):
        """
        To be overridden by child classes.
        """
        
        return
        
        yield None
    
    def get_mgf(self):
        
        self.mgf_files = {}
        
        if self.mgfdir and os.path.exists(self.mgfdir):
            
            for mgfname in os.listdir(self.mgfdir):
                
                mgfdata = remgf.search(mgfname)
                
                if not mgfdata:
                    
                    continue
                
                protein, ionmode, fracrow, fraccol = mgfdata.groups()
                frac = (fracrow, int(fraccol))
                
                if (
                    ionmode == self.ionmode and
                    protein == self.protein and
                    (not self.prot_frac or frac in self.prot_frac)
                ):
                    
                    self.mgf_files[frac] = os.path.join(self.mgfdir, mgfname))
        
        if not self.mgf_files:
            
            sys.stdout.write(
                '!!! No mgf files found for `%s` %s\n' % (
                    self.protein, self.ionmode
                )
            )
    
    def open_mgf(self):
        
        self.mgf = {}
        
        for frac, mgfpath in iteritems(self.mgf_files):
            
            self.mgf[frac] = mgf.MgfReader(
                mgfpath, label = '%s%u' % frac, charge = None
            )


class MS2Identifier(LtpTable):
    
    def __init__(
            self,
            infile,
            outdir = None,
            ms1_tolerance = None,
            prot_frac = None,
            mgfdir = None,
        ):
        
        self.ms1_tolerance = ms1_tolerance
        
        infile_path = os.path.split(infile)
        basedir = infile_path[:-2] if len(infile_path) > 2 else ['.']
        
        if not outdir:
            
            outdir = '%s__%s' % (
                os.path.split(infile)[-2],
                time.strftime('%Y.%m.%d_%H.%M'),
            )
            outdir = os.path.join(basedir + [outdir])
        
        if not os.path.exists(outdir):
            
            os.mkdir(outdir)
        
        outfile = os.path.join([outdir, infile_path[-1]])
        
        LtpTable.__init__(
            self,
            infile = infile,
            outfile = outfile,
            prot_frac = prot_frac,
            mgfdir = mgfdir,
        )
    
    def new_columns(self):
        
        databases = ['all', 'swisslipids', 'lipidmaps', 'lipyd']
        
        titles = [
            'ms2_new_top',
            'ms2_new_all',
            'databases_all',
            'database_swisslipids',
            'database_lipidmaps',
            'database_lipyd',
        ]
        
        idx = {
            'ms2_new_top': self.col_letter('U'),
            'ms2_new_all': 'END',
            'database_all': 'END',
            'database_swisslipids': 'END',
            'database_lipidmaps': 'END',
            'database_lipyd': 'END',
        }
        
        data = dict((title, []) for title in titles)
        adducts = list(settings.get('ad2ex')[1][self.ionmode].keys())
        
        for add in adducts:
            
            for db in databases:
                
                title = '%s_database_%s' % (adduct, db)
                titles.append(title)
                idx[title] = 'END'
        
        if self.numof_cols() == 0:
            
            return
        
        for i, row in enumerate(self.sheet.rows):
            
            if i == 0:
                
                # header row
                continue
            
            mz = row[self.col_letter('P') - 1]
            ms1_records = moldb.adduct_lookup(
                mz, self.ionmode, tolerance = self.ms1_tolerance
            )
            rt = tuple(
                float(i.strip())
                for i in row[self.col_letter('G') - 1].split('-')
            )
        
        for title in titles:
            
            yield idx[title], title, data[title]
