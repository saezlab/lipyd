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


import lipyd._version as _version
__version__ = _version.__version__

import lipyd.session as session

_session = session.get_session()
_session.log.msg('This is %s module version %s' % (__name__, __version__))

from lipyd.main import *
