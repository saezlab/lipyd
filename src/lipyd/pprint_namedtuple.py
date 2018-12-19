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

# from https://stackoverflow.com/a/43823671/854988

import sys
from io import StringIO
import pprint


class PrettyPrinter(pprint.PrettyPrinter):
    """ """
    
    def format_namedtuple(self, object, stream, indent, allowance, context, level):
        """

        Parameters
        ----------
        object :
            
        stream :
            
        indent :
            
        allowance :
            
        context :
            
        level :
            

        Returns
        -------

        """
        # Code almost equal to _format_dict, see pprint code
        write = stream.write
        write(object.__class__.__name__ + '(')
        object_dict = object._asdict()
        length = len(object_dict)
        multiline = False
        if length:
            # We first try to print inline, and if it is too large then we print it on multiple lines
            inline_stream = StringIO()
            self.format_namedtuple_items(object_dict.items(), inline_stream, indent, allowance + 1, context, level, inline=True)
            max_width = self._width - indent - allowance
            if len(inline_stream.getvalue()) > max_width:
                self.format_namedtuple_items(object_dict.items(), stream, indent, allowance + 1, context, level, inline=False)
                multiline = True
            else:
                stream.write(inline_stream.getvalue())
        
        if multiline:
            write('\n' + ' ' * (level * 2 - 2) + ')')
        else:
            write(')')

    def format_namedtuple_items(self, items, stream, indent, allowance, context, level, inline=False):
        """

        Parameters
        ----------
        items :
            
        stream :
            
        indent :
            
        allowance :
            
        context :
            
        level :
            
        inline :
             (Default value = False)

        Returns
        -------

        """
        # Code almost equal to _format_dict_items, see pprint code
        #print(level, file = sys.stdout)
        indent += self._indent_per_level
        write = stream.write
        last_index = len(items) - 1
        if inline:
            delimnl = ', '
        else:
            delimnl = ',\n' + '  ' * level
            write('\n' + '  ' * level)
        for i, (key, ent) in enumerate(items):
            last = i == last_index
            write(key + ' = ')
            self._format(ent, stream, indent,
                         allowance if last else 1,
                         context, level + 1)
            if not last:
                write(delimnl)
    
    def _pprint_tuple(self, object, stream, indent, allowance, context, level):
        """

        Parameters
        ----------
        object :
            
        stream :
            
        indent :
            
        allowance :
            
        context :
            
        level :
            

        Returns
        -------

        """
        
        stream.write('(\n' + ' ' * (level * 2))
        endchar = ',\n)' if len(object) == 1 else '\n' + ' ' * (level  * 2) + ')'
        self._format_items(
            object, stream, level,
            allowance + (len(endchar) if '\n' not in endchar else 0),
            context, level)
        stream.write(endchar)
    
    _dispatch = {}
    _dispatch[tuple.__repr__] = _pprint_tuple
    
    def _format(self, object, stream, indent, allowance, context, level):
        """

        Parameters
        ----------
        object :
            
        stream :
            
        indent :
            
        allowance :
            
        context :
            
        level :
            

        Returns
        -------

        """
        # We dynamically add the types of our namedtuple and namedtuple like 
        # classes to the _dispatch object of pprint that maps classes to
        # formatting methods
        # We use a simple criteria (_asdict method) that allows us to use the
        # same formatting on other classes but a more precise one is possible
        
        if hasattr(object, '_asdict') and type(object).__repr__ not in self._dispatch:
            self._dispatch[type(object).__repr__] = PrettyPrinter.format_namedtuple
        super()._format(object, stream, indent, allowance, context, level)
