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
import numpy
import pandas

from xpipeline.cluster import superclustertriggers

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['XCluster']

_default_columns = ['min_time_of_cluster',
                    'weighted_center_time', 'max_time_of_cluster',
                    'min_frequency_of_cluster',
                    'weighted_center_frequency',
                    'max_frequency_of_cluster',
                    'number_of_pixels', 'energy_of_cluster']

_default_dtype = numpy.float64

class XCluster(pandas.DataFrame):
    """Initiate an `XCluster` class

    This class is designed to take the output
    from `clusterproperties.clusterproperities_wrapper`
    and provide wrapper methods
    for calculating a variety of statistics
    for the cluster of pixels

    Returns:
        `XCluster`
    """
    @classmethod
    def nearest_neighbor(cls, cluster_array, labelled_map, **kwargs):
        columns = kwargs.pop('columns', _default_columns)
        dtype = kwargs.pop('dtype', _default_dtype)
        tab = cls(cluster_array, columns=columns, dtype=dtype)
        tab['dt'] = tab.max_time_of_cluster - tab.min_time_of_cluster
        tab['df'] = tab.max_frequency_of_cluster - tab.min_frequency_of_cluster
        return tab

    def supercluster(self, statistic_column='energy_of_cluster'):
        bounding_box_columns = ['min_time_of_cluster', 'min_frequency_of_cluster',
                                 'dt', 'df']
        bounding_box_columns.append(statistic_column)
        bound_box = self[bounding_box_columns].values
        mask = superclustertriggers.fastsupercluster_wrapper(len(bound_box), bound_box).astype(bool)
        return self.loc[mask]
