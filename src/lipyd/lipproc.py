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

import collections


Chain = namedtuple('Chain', ['c', 'u', 't', 'p'])
# by default type id `FA` the prefix is empty tuple
Chain.__new__.__defaults__ = ('FA', (,))

LipidLabel = namedtuple('LipidLabel', ['cls', 'db_id', 'db', 'names'])
# names are empty tuple by default
LipidLabel.__new__.defaults__ = ((,),)
