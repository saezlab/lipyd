#!/usr/bin/env python

#
# This file is part of the `lipyd` Python module
#
# Copyright 2024
# Heidelberg University Hospital
#
# File author(s): Saez Lab (omnipathdb@gmail.com)
#
# Distributed under the GPLv3 license
# See the file `LICENSE` or read a copy at
# https://www.gnu.org/licenses/gpl-3.0.txt
#

"""
Python module for lipidomics LC MS/MS data processing.
"""

__all__ = [
    '__version__',
    '__author__',
]

from ._metadata import __author__, __version__

import pprint
import lipyd.session as session
import lipyd.pprint_namedtuple as pprint_namedtuple

_session = session.get_session()
_session.log.msg('This is %s module version %s' % (__name__, __version__))

pprint.PrettyPrinter = pprint_namedtuple.PrettyPrinter
