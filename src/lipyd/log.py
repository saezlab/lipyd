#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2019 - EMBL
#
#  File author(s):
#  Dénes Türei (turei.denes@gmail.com)
#  Igor Bulanov
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://denes.omnipathdb.org/
#

import time
import os


def new_logger(name = 'lipyd', logdir = 'lipyd_log', verbosity = 0):
    """Returns a new logger with default settings (can be customized).

    Parameters
    ----------
    name :
         (Default value = 'lipyd')
    logdir :
         (Default value = 'lipyd_log')
    verbosity :
         (Default value = 0)

    Returns
    -------

    """
    
    return Logger(
        fname = '%s__%s.log' % (
            name,
            Logger.timestamp().replace(' ', '_').replace(':', '.'),
        ),
        verbosity = 0,
        logdir = logdir,
    )


class Logger(object):
    """ """
    
    strftime = time.strftime
    
    def __init__(self, fname, verbosity = 0, logdir = None):
        
        self.logdir = self.get_logdir(logdir)
        self.fname  = os.path.join(self.logdir, fname)
        self.verbosity = verbosity
        self.open_logfile()
        self.msg('Welcome!')
        self.msg('Logger started, logging into `%s`.' % self.fname)
    
    def msg(self, msg = '', level = 0):
        """

        Parameters
        ----------
        msg :
             (Default value = '')
        level :
             (Default value = 0)

        Returns
        -------

        """
        
        if level <= self.verbosity:
            
            self.fp.write('[%s] %s\n' % (self.timestamp(), msg))
    
    @classmethod
    def timestamp(cls):
        """ """
        
        return cls.strftime('%Y-%m-%d %H:%M:%S')
    
    def __del__(self):
        
        self.msg('Logger shut down, logfile `%s` closed.' % self.fname)
        self.msg('Bye.')
        self.close_logfile()
    
    def get_logdir(self, dirname = None):
        """

        Parameters
        ----------
        dirname :
             (Default value = None)

        Returns
        -------

        """
        
        dirname = dirname or 'lipyd_log'
        
        if not os.path.isdir(dirname):
            
            os.mkdir(dirname)
        
        return dirname
    
    def open_logfile(self):
        """ """
        
        self.close_logfile()
        self.fp = open(self.fname, 'w')
    
    def close_logfile(self):
        """ """
        
        if hasattr(self, 'fp') and not self.fp.closed:
            
            self.fp.close()
    
    def flush(self):
        """ """
        
        if hasattr(self, 'fp') and not self.fp.closed:
            
            self.fp.flush()
