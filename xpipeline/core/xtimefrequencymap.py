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
from gwpy.signal.fft.ui import seconds_to_samples

from xpipeline.cluster import nearestneighbor
from xpipeline.cluster import clusterproperties
from ..cluster.cluster import XCluster
from .xfrequencyseries import XFrequencySeriesDict
from .sparse import csc_sparse_map

import operator
import numpy
import h5py

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

    def blackout_pixels(self, blackpixel_percentile):
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

        return XSparseTimeFrequencyMapDict({k: self[k].blackout_pixels(v)
                                     for k,v in blackpixel_percentile.items()})

    def circular_time_slide(self, npixels_to_shift):
        """slide all maps in dict by specified number of pixels

        Parameters:

            npixels_to_shift : `dict`, `int`
                either a `dict` of (channel, `int`) pairs for key-wise
                internal time slide calc. Since it does not make sense
                to slide all the tfmaps the same number of pixels
                you must specify a dict with all detectors and pixels to slide.

        Returns:
            `dict`: key-wise pair of channel
        """
        if not isinstance(npixels_to_shift, dict):
            raise ValueError("Must be a dict")

        return XTimeFrequencyMapDict({key :
                self[key].circular_time_slide(npixels_to_shift=npix_to_shift)
                for key, npix_to_shift in npixels_to_shift.items()})

    def to_sparse(self, tindex, findex):
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
            {k : v.to_sparse(tindex=tindex[k], findex=findex[k])
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
    def blackout_pixels(self, blackpixel_percentile):
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
                                           yindex=self.yindex)

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
        sqrt_of_neg_1 = numpy.sqrt(numpy.array([-1], dtype=complex))
        frequency_shift = residual_time_shift(delta,
                                              self.frequencies.to_value())
        return self * frequency_shift

    def circular_time_slide(self, npixels_to_shift):
        """Slide the TF pixels of this map

        This should move the appropriate number of time bins
        such that slide represents a slide in seconds.

        Parameters
        ----------
        seconds : `int`,
            How many seconds we are sliding the map

        sample_frequency : `float`
            what is the sample frequency of the data

        offset : `float`
            what offset was used to make this spectrogram

        Returns:
            `XTimeFrequencyMap`:
                A time frequency map slide by the appropriate number
                of seconds
        """
        if not npixels_to_shift:
            return self

        idx = numpy.searchsorted(self.xindex.value,
                                 numpy.roll(self.xindex.value,
                                            npixels_to_shift
                                ))

        return self[idx]

    def to_dominant_polarization_frame(self, dpf_asd):
        return self * dpf_asd

    def to_sparse(self, tindex, findex):
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
                                           findex=findex, energy=data)

class XSparseTimeFrequencyMapDict(OrderedDict):
    def to_coherent(self):
        """Sum all maps in the dict

           Returns:
               `XTimeFrequencyMap`:
                   A coherent TF-Map
        """
        return reduce(operator.add, self.values())

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
            v.labels = labelled_map
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

            cluster_array = clusterproperties.clusterproperities_wrapper(dim_array, labelled_map, total_energy, pixels[0,:] + 1, pixels[1,:] + 1).T

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

class csc_XSparseTimeFrequencyMap(csc_sparse_map):
    def write(self, filename, path):
        """Plot the data for this `XTimeFrequencyMapDict`.

        Parameters
        ----------
        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        f = h5py.File(filename,'w')
        g = f.create_group(path)
        g.create_dataset('value', data=self.energy)
        g.create_dataset('tindex', data=self.tindex)
        g.create_dataset('findex', data=self.findex)
        g.attrs['shape'] = self.shape
        import pdb
        pdb.set_trace()

    @classmethod
    def read(cls, filename, path):
        """Plot the data for this `XTimeFrequencyMapDict`.

        Parameters
        ----------
        **kwargs
            all other keyword arguments are passed to the plotter as
            appropriate
        """
        f = h5py.File(filename,'r')
        g = f[path]
        import pdb
        pdb.set_trace()
        return cls((g['value'], (g2['tindex'], g2['tindex'])),
                    shape=g.attrs['shape'], tindex=tindex, findex=findex,
                    energy=g['value'])

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
