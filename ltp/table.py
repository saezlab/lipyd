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

resheet = re.compile(r'([A-z0-9]+)_(positive|negative)_?(best)?')

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
        
        self.sheet.insert_cols(idx)
        
        ii = 1
        
        if title is not None:
            
            ii = 2
            self.sheet[1][idx].value = title
        
        for i, val in enumerate(values):
            
            self.sheet[i + ii][idx].value = val
    
    def write(self):
        
        self.xls_original.save(self.outfile)


class LtpTable(TableBase):
    
    def __init__(self, infile, outfile = None):
        
        TableBase.__init__(infile, outfile)
    
    def process_sheets(self):
