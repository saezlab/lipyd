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


import lipyd._version as _version
__version__ = _version.__version__

import pprint
import lipyd.session as session
import lipyd.pprint_namedtuple as pprint_namedtuple

_session = session.get_session()
_session.log.msg('This is %s module version %s' % (__name__, __version__))

pprint.PrettyPrinter = pprint_namedtuple.PrettyPrinter
