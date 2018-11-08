#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `lipyd` python module
#
#  Copyright (c) 2015-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import random
import os
import sys

import lipyd.log as log


class Session(object):
    
    def __init__(self, label = None, log_verbosity = 0):
        
        self.label = label or self.gen_session_id()
        self.log_verbosity = log_verbosity
        self.start_logger()
        self.log.msg('Session `%s` started.' % self.label)
    
    @staticmethod
    def gen_session_id(length = 5):
        
        abc = '0123456789abcdefghijklmnopqrstuvwxyz'
        return ''.join(random.choice(abc) for i in range(length))
    
    def start_logger(self):
        
        self.logfile = 'lipyd-%s.log' % self.label
        self.log = log.Logger(self.logfile, verbosity = self.log_verbosity)
    
    def __del__(self):
        
        self.log.msg('Session `%s` finished.' % self.label)


def get_session():
    
    mod = sys.modules[__name__]
    
    if not hasattr(mod, 'session'):
        
        new_session()
    
    return sys.modules[__name__].session


def get_log():
    
    return get_session().log


def new_session(label = None, log_verbosity = 0):
    
    mod = sys.modules[__name__]
    
    if hasattr(mod, 'session'):
        
        delattr(mod, 'session')
    
    setattr(mod, 'session', Session(label, log_verbosity))
