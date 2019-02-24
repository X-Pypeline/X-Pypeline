# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017-)
#
# This file is part of the XPypeline python package.
#
# hveto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hveto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hveto.  If not, see <http://www.gnu.org/licenses/>.

# ---- Import standard modules to the python path.

from scipy import sparse

class csc_sparse_map(sparse.csc_matrix):
    _metadata_slots = ('energy', 'tindex', 'findex',
                       'yindex', 'xindex', 'labels', 'name')
    def __init__(self, matrix, **kwargs):
        self.yindex = kwargs.pop('yindex', None)
        self.xindex = kwargs.pop('xindex', None)
        self.tindex = kwargs.pop('tindex', None)
        self.findex = kwargs.pop('findex', None)
        self.energy = kwargs.pop('energy', None)
        self.labels = kwargs.pop('labels', None)
        self.labels = kwargs.pop('name', None)
        super(csc_sparse_map, self).__init__(matrix, **kwargs)

    def _repr_helper(self, print_):
        if print_ is repr:
            opstr = '='
        else:
            opstr = ': '

        # get prefix and suffix
        prefix = '{}('.format(type(self).__name__)
        suffix = ')'
        if print_ is repr:
            prefix = '<{}'.format(prefix)
            suffix += '>'

        indent = ' ' * len(prefix)

        # format value
        nnz = self.getnnz()
        arrstr = ("<%dx%d sparse matrix of type '%s'\n"
                 "\twith %d stored elements in %s format>" % \
                 (self.shape + (self.dtype.type, nnz, "Compressed Sparse Column")))

        # format unit
        metadata = [('unit', 'dimensionless')]

        # format other metadata
        try:
            attrs = self._print_slots
        except AttributeError:
            attrs = self._metadata_slots
        for key in attrs:
            try:
                val = getattr(self, key)
            except (AttributeError, KeyError):
                val = None
            thisindent = indent + ' ' * (len(key) + len(opstr))
            metadata.append((
                key.lstrip('_'),
                print_(val).replace('\n', '\n{}'.format(thisindent)),
            ))
        metadata = (',\n{}'.format(indent)).join(
            '{0}{1}{2}'.format(key, opstr, value) for key, value in metadata)

        return "{0}{1}\n{2}{3}{4}".format(
            prefix, arrstr, indent, metadata, suffix)

    def __repr__(self):
        """Return a representation of this object

        This just represents each of the metadata objects appropriately
        after the core data array
        """
        return self._repr_helper(repr)

    def __str__(self):
        """Return a printable string format representation of this object

        This just prints each of the metadata objects appropriately
        after the core data array
        """
        return self._repr_helper(str)
