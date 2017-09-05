#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `emese` python module
#
#  Copyright (c) 2014-2017 - EMBL
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from __future__ import print_function

import sys
import imp
import time

import tqdm

__all__ = ['Progress', 'OldProgress']

class Progress(object):
    
    """
    Before I had my custom progressbar here.
    Now it is a wrapper around the great progressbar `tqdm`.
    Old implementation moved to `OldProgress` class.
    """
    
    def __init__(self, total = None, name = "Progress",
             interval = None, percent = True, status = 'initializing',
             done = 0, init = True, unit = 'it'):
        
        self.name = name
        self.interval = max(int(total / 100), 1) if interval is None else interval
        self.total = total
        self.done = done
        self.status = status
        self.unit = unit
        self.start_time = time.time()
        self.min_update_interval = 0.1
        self.last_printed_value = 0
        
        if init:
            self.init_tqdm()
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def init_tqdm(self):
        self.tqdm = tqdm.tqdm(total = self.total,
                              desc = '%s: %s' % (self.name, self.status),
                              unit_scale = True,
                              unit = self.unit)
        self.last_updated = time.time()
    
    def step(self, step = 1, msg = None, status = 'busy', force = False):
        """
        Updates the progressbar by the desired number of steps.
        
        :param int step: Number of steps or items.
        """
        self.done += step
        
        if force or (self.done % self.interval < 1.0 and \
            time.time() - self.last_updated > self.min_update_interval):
            
            self.set_status(status)
            
            this_update = max(0, self.done - self.last_printed_value)
            
            if this_update == 0:
                self.tqdm.refresh()
                self.tqdm.fp.flush()
            else:
                self.tqdm.update(int(this_update))
                
            self.last_printed_value = self.done
            self.last_updated = time.time()
    
    def terminate(self, status = 'finished'):
        """
        Terminates the progressbar and destroys the tqdm object.
        """
        self.step(self.total - self.done, force = True, status = status)
        self.tqdm.close()
    
    def set_total(self, total):
        """
        Changes the total value of the progress bar.
        """
        self.total = total
        self.tqdm.total = total
        self.step(0)
    
    def set_done(self, done):
        """
        Sets the position of the progress bar.
        """
        self.done = done
        self.tqdm.n = self.done
        self.tqdm.last_print_n = self.done
        self.step(0)
    
    def set_status(self, status):
        """
        Changes the prefix of the progressbar.
        """
        if status != self.status:
            self.status = status
            self.tqdm.set_description(self.get_desc())
            self.tqdm.refresh()
    
    def get_desc(self):
        """
        Returns a formatted string of the description, consisted of
        the name and the status. The name supposed something constant
        within the life of the progressbar, while the status is there
        to give information about the current stage of the task.
        """
        return '%s%s%s%s' % (' ' * 8,
                             self.name,
                             ' -- ' if len(self.name) else '',
                             self.status)


class OldProgress(object):
    
    def __init__(self, total = None, name = "Progress",
        interval = 3000, percent = True, status = 'initializing'):
        self.status = status
        self.name = name
        self.interval = interval
        self.set_total(total)
        self.done = 0
        self.percent = percent
        self.in_ipython()
        self.last_updated = time.time()
        self.start_time = time.time()
        self.min_update_interval = 0.0 \
            if not self.ipython or self.ipython == 'terminal' \
            else 1.0
        sys.stdout.write("\r"+" "*150)
        if self.percent:
            sys.stdout.write("\r\t:: %s: %s 0.00%%" % \
                (self.name, '%s,'%self.status))
        else:
            self.total = int(self.total)
            sys.stdout.write("\r\t:: %s: %s 0/%u" % \
                (self.name, '%s,'%self.status, self.total))
        sys.stdout.flush()
    
    def step(self, step = 1, msg = None, status = 'working on it', force = False):
        self.status = status
        self.done += step
        if force or (self.done % self.interval < 1.0 and \
            time.time() - self.last_updated > self.min_update_interval):
            sys.stdout.write("\r"+" "*150)
            if self.percent:
                sys.stdout.write(
                    "\r\t:: %s: %s %.2f%% %s" % (
                        self.name,
                        '' if self.status is None else'%s,'%self.status,
                        float(self.done)/float(self.total)*100.0,
                        '' if msg is None else '[%s]'%msg))
            else:
                sys.stdout.write(
                    "\r\t:: %s: %s %u/%u %s" % (
                        self.name, '%s,'%self.status, self.done, self.total,
                        '' if msg is None else '[%s]'%msg))
            sys.stdout.flush()
            self.last_updated = time.time()
    
    def set_total(self, total):
        self.total = total or 9999999999.0
    
    def set_done(self, done):
        self.done = done
    
    def set_status(self, status):
        #if self.ipython == 'notebook':
        #    status = '%s; progress indicator disabled in IPython notebook.' % status
        self.step(step = 0, status = status, force = True)
    
    def terminate(self, status = 'finished'):
        self.finished_time = time.time()
        self.seconds_elapsed =  self.finished_time - self.start_time
        m, s = divmod(self.seconds_elapsed, 60)
        h, m = divmod(m, 60)
        self.time_elapsed = "%d:%02d:%02d" % (h, m, s)
        sys.stdout.write("\r"+" "*150)
        self.status = status
        if self.percent:
            sys.stdout.write("\r\t:: %s: %s 100.0%% [%s elapsed]" % \
                (self.name, '%s,'%self.status, self.time_elapsed))
        else:
            sys.stdout.write("\r\t:: %s: %s %u/%u [%s elapsed]" % (\
                self.name,
                '%s,'%self.status,
                self.total,
                self.total,
                self.time_elapsed)
            )
        sys.stdout.write("\n")
        sys.stdout.flush()
        self.last_updated = time.time()
    
    def in_ipython(self):
        self.ipython = False
        if 'ipykernel' in sys.modules:
            self.ipython = 'notebook'
        elif 'Ipython' in sys.modules:
            self.ipython = 'terminal'
