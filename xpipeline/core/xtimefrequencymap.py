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
from collections import OrderedDict
from functools import reduce

from gwpy.spectrogram import Spectrogram
from gwpy.timeseries import TimeSeries

from xpipeline.cluster import nearestneighbor
from xpipeline.cluster import clusterproperties
from ..cluster.cluster import XCluster
from .xfrequencyseries import XFrequencySeriesDict
from .sparse import csc_sparse_map
from filelock import Timeout, FileLock

import operator
import numpy
import h5py
import tables

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['csc_XSparseTimeFrequencyMap',
           'XSparseTimeFrequencyMapDict'
           'XTimeFrequencyMapDict',
           'residual_time_shift', 'XTimeFrequencyMap']

class XTimeFrequencyMapDict(OrderedDict):
    def abs(self):
        """Take the absolute value of all maps in dict

           Returns:
               `XTimeFrequencyMapDict`:
                   power_map of all Fourier Grams in Dict
        """
        return XTimeFrequencyMapDict({k: v.abs() for k,v in self.items()})

    def to_coherent(self):
        """Sum all maps in the dict

           Returns:
               `XTimeFrequencyMap`:
                   A coherent TF-Map
        """
        return reduce(operator.add, self.values())

    def to_dominant_polarization_frame(self, projected_asds):
        """Project Tfmap to an antenna frame give a dict of asds

           Parameters:
               projected_asds : `dict`
                   key-wise dict of antenna pattern name
                   and values `XFrequencySeriesDict`

           Returns:
               `OrderedDict`:
                   A key-wise dict of antenna pattern name
                   and value `XTimeFrequencyMapDict`
        """
        projected_time_frequency_maps = OrderedDict()
        for pattern, asds in projected_asds.items():
            projected_time_frequency_maps[pattern] = XTimeFrequencyMapDict()
            for det, asd in asds.items():
                mask = numpy.in1d(asd.xindex, self[det].yindex)
                projected_time_frequency_maps[pattern][det] = self[det] * asd[mask]
        return projected_time_frequency_maps

    def blackout_pixels(self, blackpixel_percentile, **kwargs):
        """Set pixels below certain energy level to zero

        Parameters:

            blackpixel_percentile : `dict`, `int`
                either a `dict` of (channel, `int`) pairs for key-wise
                significant pixel calc, or a single int to use as the
                threshold of all items.

        Returns:
            `dict`: key-wise pair of channel : freq,time indices
        """
        if not isinstance(blackpixel_percentile, dict):
            blackpixel_percentile = dict((c, blackpixel_percentile)
                                         for c in self)

        return XSparseTimeFrequencyMapDict({k: self[k].blackout_pixels(v, **kwargs)
                                     for k,v in blackpixel_percentile.items()})

    def to_sparse(self, tindex, findex, **kwargs):
        """

        Parameters:

            tindex : `dict`,
                a `dict` of (channel, array) pairs for key-wise
                of the time indices to extract for the
                creation of a dict of sparse tfmaps from
                a dict of full tfmaps

            findex : `dict`,
                a `dict` of (channel, array) pairs for key-wise
                of the frequency indices to extract for the
                creation of a dict of sparse tfmaps from
                a dict of full tfmaps

        Returns:
            `XSparseTimeFrequencyMapDict`:
                Sparse matrices using supplied t-f indices
                and extracting the data at those points
                from the full tfmap.
        """
        if not isinstance(tindex, dict):
            raise ValueError("Must be a dict")

        if not isinstance(findex, dict):
            raise ValueError("Must be a dict")

        return XSparseTimeFrequencyMapDict(
            {k : v.to_sparse(tindex=tindex[k], findex=findex[k], **kwargs)
            for k, v in self.items()})

    def plot(self, label='key', **kwargs):
        """Plot the data for this `XTimeFrequencyMapDict`.

        Parameters
        ----------
        label : `str`, optional

            labelling system to use, or fixed label for all elements
            Special values include

            - ``'key'``: use the key of the `XTimeFrequencyMapDict`,
            - ``'name'``: use the :attr:`~XTimeSeries.name` of each element

            If anything else, that fixed label will be used for all lines.

        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        from matplotlib import pyplot
        from gwpy.plot import Plot

        figargs = dict()
        for key in ['figsize', 'dpi']:
            if key in kwargs:
                figargs[key] = kwargs.pop(key)

        nrows = len(self.keys())
        plot_, axes = pyplot.subplots(nrows=nrows, ncols=1, sharex=True,
                                      subplot_kw={'xscale': 'auto-gps'},
                                      FigureClass=Plot, **figargs)
        iaxis = 0
        for lab, imap in self.items():
            if label.lower() == 'name':
                lab = tfmap.name
            elif label.lower() != 'key':
                lab = label
            axes[iaxis].imshow(imap)
            axes[iaxis].set_title(str(lab).replace('_','-'))
            if iaxis != nrows-1:
                axes[iaxis].set_xlabel(' ')
            else:
                break
            iaxis += 1
        return plot_

    def gaussianity(self):
        """Calculate the gaussianity of all maps in this dict

           (.abs().percentile(99)/ .abs().median(0))**2
           Returns:
               `XFrequencySeriesDict`:
                   The gaussianity measure of each frequency bin
        """
        # make sure we are dealing with energy map not fftgram
        return XFrequencySeriesDict({k: v.gaussianity()
                                     for k,v in self.items()}
                                    )


class XTimeFrequencyMap(Spectrogram):
    def blackout_pixels(self, blackpixel_percentile, **kwargs):
        """Set pixels below certain energy level to zero

           Parameters:

           blackpixel_percentile : `int`
               the x-percentile loudest tile all pixels
               below which will be set to energy of 0
        """
        energy_threshold  = numpy.percentile(self, blackpixel_percentile,
                                             interpolation='midpoint')
        tf_indices = numpy.nonzero(self.value >= energy_threshold)
        time = tf_indices[0]
        freq = tf_indices[1]
        data = self.value[time, freq]
        return csc_XSparseTimeFrequencyMap((data, (time, freq)),
                                           shape=self.shape,
                                           tindex=time,
                                           findex=freq,
                                           energy=data,
                                           xindex=self.xindex,
                                           yindex=self.yindex,
                                           dx=self.xindex[1] - self.xindex[0],
                                           x0=self.xindex[0],
                                           y0=self.yindex[0],
                                           dy=self.yindex[1] - self.yindex[0],
                                           name=self.name, **kwargs)

    def gaussianity(self):
        """Calculate the gaussianity of this map

           (.abs().percentile(99)/ .abs().median(0))**2
           Returns:
               `gwpy.frequency.FrequencySeries`:
                   The gaussianity measure of each frequency bin
        """
        self_ = self.abs()
        return self_.percentile(99) / self_.median(0)

    def phaseshift(self, delta):
        """Phase shift this `spectrogram` by ``delta``

        This modifies the spectrogram in-place.

        Parameters
        ----------
        delta : `float`, `~astropy.units.Quantity`, `str`
            The amount by which to shift (in seconds if `float`), give
            a negative value to shift backwards in time
        """
        frequency_shift = residual_time_shift(delta,
                                              self.frequencies.to_value())
        return self * frequency_shift

    def to_dominant_polarization_frame(self, dpf_asd):
        return self * dpf_asd

    def to_sparse(self, tindex, findex, **kwargs):
        """

        Parameters:

            tindex : `array`,
                an array of the time indices to extract for the
                creation of a dict of sparse tfmaps from
                a dict of full tfmaps

            findex : `dict`,
                an array of the frequency indices to extract for the
                creation of a dict of sparse tfmaps from
                a dict of full tfmaps

        Returns:
            `XSparseTimeFrequencyMapDict`:
                Sparse matrices using supplied t-f indices
                and extracting the data at those points
                from the full tfmap.
        """
        shape = self.shape
        data = self.value[tindex, findex]
        return csc_XSparseTimeFrequencyMap((data, (tindex, findex)),
                                           shape=shape, yindex=self.yindex,
                                           xindex=self.xindex, tindex=tindex,
                                           findex=findex, energy=data,
                                           dx=self.xindex[1] - self.xindex[0],
                                           dy=self.yindex[1] - self.yindex[0],
                                           x0=self.xindex[0],
                                           y0=self.yindex[0],
                                           name=self.name, **kwargs)

class XSparseTimeFrequencyMapDict(OrderedDict):
    def to_coherent(self):
        """Sum all maps in the dict

           Returns:
               `XTimeFrequencyMap`:
                   A coherent TF-Map
        """
        return reduce(operator.add, self.values())

    def abs(self):
        """Take the absolute value of all maps in dict

           Returns:
               `XTimeFrequencyMapDict`:
                   power_map of all Fourier Grams in Dict
        """
        return XSparseTimeFrequencyMapDict({k: csc_XSparseTimeFrequencyMap(numpy.abs(v), yindex=v.yindex,
                                           xindex=v.xindex, tindex=v.tindex,
                                           findex=v.findex, energy=numpy.abs(v.energy),
                                           dx=v.xindex[1] - v.xindex[0],
                                           dy=v.yindex[1] - v.yindex[0],
                                           x0=v.xindex[0],
                                           y0=v.yindex[0],
                                           name=v.name,
                                           phi=v.phi, theta=v.theta,) for k,v in self.items()})

    def to_xtimefrequencymapdict(self):
        """Convert dict fo sparse matrix to `XTimeFrequencyMapDict`
        """
        maps = {k : XTimeFrequencyMap(v.toarray(), xindex=v.xindex,
                                      yindex=v.yindex)
                for k, v in self.items()}
        return XTimeFrequencyMapDict(maps)

    def label(self, connectivity=8):
        """Convert dict fo sparse matrix to `XTimeFrequencyMapDict`
        """
        tfmap = list(self.values())[0]
        pixels = numpy.vstack([tfmap.tindex, tfmap.findex])
        coord_dim_array = (tfmap.xindex.size, tfmap.yindex.size)
        npixels = pixels.shape[1]
        labelled_map = nearestneighbor.fastlabel_wrapper(pixels + 1, coord_dim_array,
                                                         connectivity, npixels).astype(int)

        for k, v in self.items():
            v.pixel_labels = labelled_map
        return


    def cluster(self, method='nearest_neighbors', **kwargs):
        """Convert dict fo sparse matrix to `XTimeFrequencyMapDict`
        """
        if method=='nearest_neighbors':
            connectivity = kwargs.pop('connectivity', 8)
            total_energy = 0
            for k, v in self.items():
                total_energy += v.energy
            pixels = numpy.vstack([v.tindex, v.findex])
            coord_dim_array = (v.xindex.size, v.yindex.size)

            npixels = pixels.shape[1]

            labelled_map = nearestneighbor.fastlabel_wrapper(pixels + 1, coord_dim_array, connectivity, npixels).astype(int)

            dim_array = numpy.array([total_energy.shape[0], 1, 2.0])

            cluster_array = clusterproperties.clusterproperities_wrapper(dim_array, labelled_map, total_energy, pixels[0,:] + 1, pixels[1,:] + 1, True).T

            cluster_array[:, 0:3] = cluster_array[:, 0:3] * (v.xindex[2] - v.xindex[1])  + v.xindex[0]

            cluster_array[:, 3:6] = cluster_array[:, 3:6] * (v.yindex[2] - v.yindex[1])  + v.yindex[0]

            return XCluster.nearest_neighbor(cluster_array, labelled_map)
        else:
            raise ValueError('Clustering method undefined')

    def to_cnn(self, dim=2):
        """Sum all maps in the dict

           Returns:
               `XTimeFrequencyMap`:
                   A coherent TF-Map
        """
        return

    def plot(self, label='key', **kwargs):
        """Plot the data for this `XTimeFrequencyMapDict`.

        Parameters
        ----------
        label : `str`, optional

            labelling system to use, or fixed label for all elements
            Special values include

            - ``'key'``: use the key of the `XTimeFrequencyMapDict`,
            - ``'name'``: use the :attr:`~XTimeSeries.name` of each element

            If anything else, that fixed label will be used for all lines.

        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        return self.to_xtimefrequencymapdict().plot(label='key', **kwargs)

    def write(self, *args, **kwargs):
        """Plot the data for this `XTimeFrequencyMapDict`.

        Parameters
        ----------
        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        for k, v in self.items():
            v.write(*args, **kwargs)
        return

class csc_XSparseTimeFrequencyMap(csc_sparse_map):
    def write(self, filename, path, table_description=None,**kwargs):
        """Plot the data for this `XTimeFrequencyMapDict`.

        Parameters
        ----------
        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        if table_description is None:
            table_description = {
            'dx' : tables.Float64Col(),
            'dy' : tables.Float64Col(),
            'x0' : tables.Float64Col(),
            'y0' : tables.Float64Col(),
            'shape' : tables.Int64Col(2),
            'phi' : tables.Float64Col(),
            'theta' : tables.Float64Col(),
            'map_type' : tables.StringCol(20),
            'ifo' : tables.StringCol(2),
            'x' : tables.Float64Col(shape=self.tindex.size),
            'y' : tables.Float64Col(shape=self.findex.size),
            'energy' : tables.ComplexCol(itemsize=16, shape=self.energy.size),
            }

        lock = FileLock(filename + '.lock',)
        with lock:
            h5file = tables.open_file('{0}'.format(filename), mode="a", title="Sparse Time Frequency Maps")
            # check if table already exists in this path
            try:
                table = list(h5file.walk_nodes(path, "Table"))[0]
            except:
                # If these groups do not exist then create a new table inside of path
                table = h5file.create_table(path, 'event', table_description, "maps", createparents=True)

            event = table.row
            x = numpy.zeros(event['x'].size)
            x[:self.tindex.size] = self.tindex

            y = numpy.zeros(event['y'].size)
            y[:self.findex.size] = self.findex

            energy = numpy.zeros(event['energy'].size, dtype='complex')
            energy[:self.energy.size] = self.energy

            event['dx'] = self.dx.value
            event['dy'] = self.dy.value
            event['x0'] = self.x0.value
            event['y0'] = self.y0.value
            event['phi'] = self.phi
            event['theta'] = self.theta
            event['shape'] = self.shape
            event['map_type'] = self.map_type
            event['ifo'] = self.name
            event['x'] = x
            event['y'] = y
            event['energy'] = energy
            event.append()
            table.flush()
            h5file.close()
        return

    @classmethod
    def read(cls, row):
        """This can table a row of a PyTable and return a sparse tf map.

        Parameters
        ----------
        row (`tables.Table.row`):
            The is a row of a `tables.Table` generated by
            xpipeline-analysis

        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        dx = row['dx']
        dy = row['dy']
        x0 = row['x0']
        y0 = row['y0']
        phi = row['phi']
        theta = row['theta']
        shape = row['shape']
        map_type = row['map_type']
        ifo = row['ifo']
        x = row['x']
        y = row['y']
        energy = row['energy']
        tf_map = cls((energy, (x, y)), shape=shape, dx=dx, dy=dy,
                    x0=x0, y0=y0, map_type=map_type, name=ifo, phi=phi, theta=theta)
        tf_map.tindex = x[:tf_map.size-1].astype(int)
        tf_map.findex = y[:tf_map.size-1].astype(int)
        tf_map.energy = energy[:tf_map.size-1]
        tf_map.xindex = numpy.arange(x0, x0 + dx*shape[0], dx)
        tf_map.yindex = numpy.arange(y0, y0 + dy*shape[1], dy)
        return tf_map

    def label(self, connectivity=8):
        """Convert dict fo sparse matrix to `XTimeFrequencyMapDict`
        """
        pixels = numpy.vstack([self.tindex, self.findex])
        coord_dim_array = (self.xindex.size, self.yindex.size)
        npixels = pixels.shape[1]
        labelled_map = nearestneighbor.fastlabel_wrapper(pixels + 1, coord_dim_array,
                                                         connectivity, npixels).astype(int)

        self.pixel_labels = labelled_map
        return

    def plot(self, **kwargs):
        """Plot the data for this `XTimeFrequencyMapDict`.

        Parameters
        ----------
        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        return XTimeFrequencyMap(self.toarray(), xindex=self.xindex,
                                 yindex=self.yindex).plot()

def residual_time_shift(seconds, frequencies):
    # define sqrt of -1
    sqrt_of_neg_1 = numpy.sqrt(numpy.array([-1],dtype=complex))
    residual_time_shift_phases = numpy.exp(sqrt_of_neg_1 * 2 * numpy.pi *
                                           frequencies * seconds)

    return residual_time_shift_phases
