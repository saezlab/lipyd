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

import time
import os


class Logger(object):
    
    strftime = time.strftime
    
    def __init__(self, fname, verbosity = 0, logdir = None):
        
        self.logdir = self.get_logdir(logdir)
        self.fname  = os.path.join(self.logdir, fname)
        self.verbosity = verbosity
        self.open_logfile()
        self.msg('Welcome!')
        self.msg('Logger started, logging into `%s`.' % self.fname)
    
    def msg(self, msg = '', level = 0):
        
        if level <= self.verbosity:
            
            self.fp.write('[%s] %s\n' % (self.timestamp(), msg))
    
    def timestamp(self):
        
        return self.strftime('%Y-%m-%d %H:%M:%S')
    
    def __del__(self):
        
        self.msg('Logger shut down, logfile `%s` closed.' % self.fname)
        self.msg('Bye.')
        self.close_logfile()
    
    def get_logdir(self, dirname = None):
        
        dirname = dirname or 'lipyd_log'
        
        if not os.path.isdir(dirname):
            
            os.mkdir(dirname)
        
        return dirname
    
    def open_logfile(self):
        
        self.close_logfile()
        self.fp = open(self.fname, 'w')
    
    def close_logfile(self):
        
        if hasattr(self, 'fp') and not self.fp.closed:
            
            self.fp.close()
