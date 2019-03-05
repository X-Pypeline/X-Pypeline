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
import numpy

class csc_sparse_map(sparse.csc_matrix):
    _metadata_slots = ('energy', 'tindex', 'findex',
                       'yindex', 'xindex', 'dx', 'dy', 'x0', 'y0',
                       'pixel_labels', 'name', 'map_type',
                       'phi', 'theta',)
    def __init__(self, matrix, **kwargs):
        self.yindex = kwargs.pop('yindex', None)
        self.xindex = kwargs.pop('xindex', None)
        self.dy = kwargs.pop('dy', None)
        self.dx = kwargs.pop('dx', None)
        self.x0 = kwargs.pop('x0', None)
        self.y0 = kwargs.pop('y0', None)
        self.tindex = kwargs.pop('tindex', None)
        self.findex = kwargs.pop('findex', None)
        self.energy = kwargs.pop('energy', None)
        self.pixel_labels = kwargs.pop('pixel_labels', None)
        self.name = kwargs.pop('name', None)
        self.phi = kwargs.pop('phi', None)
        self.theta = kwargs.pop('theta', None)
        self.map_type = kwargs.pop('map_type', 'excess_energy')
        super(csc_sparse_map, self).__init__(matrix, **kwargs)

    def __metadata_finalize__(self, obj, force=False):
        # apply metadata from obj to self if creating a new object
        for attr in self._metadata_slots:
            _attr = '%s' % attr  # use private attribute (not property)
            # if attribute is unset, default it to None, then update
            # from obj if desired
            try:
                getattr(self, _attr)
            except AttributeError:
                update = True
            else:
                update = force

            if update:
                try:
                    val = getattr(obj, _attr)
                except AttributeError:
                    continue
                else:
                    setattr(self, _attr, val)

    def __add__(self, other):
        obj = super(csc_sparse_map, self).__add__(other) 
        obj.__metadata_finalize__(self, force=True)
        if (obj.energy is not None) and (self.energy is not None):
            obj.energy = numpy.add(self.energy, other.energy)
        else:
            obj.energy = None
        return obj

    def __abs__(self):
        obj = super(csc_sparse_map, self).__abs__()
        obj.__metadata_finalize__(self, force=True)
        if obj.energy is not None:
            obj.energy = numpy.abs(obj.energy)
        return obj

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
