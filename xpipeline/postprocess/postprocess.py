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
from xpipeline.core import XSparseTimeFrequencyMapDict, csc_XSparseTimeFrequencyMap
from ..cluster.cluster import XCluster

_default_columns = ['min_time_of_cluster',
                    'weighted_center_time', 'max_time_of_cluster',
                    'min_frequency_of_cluster',
                    'weighted_center_frequency',
                    'max_frequency_of_cluster',
                    'number_of_pixels',]

def extract_clusters_from_dict(maps, statistic_column='standard_energy', connectivity_list=[8,24,48,80],):
    all_clusters = XCluster()
    for fft_length in maps.keys():
        for skypostion, imap in maps[fft_length].items():
            # Obtain all sparse time frequency maps for this fftlength and sky position
            all_energies = []
            all_columns = _default_columns.copy()
            all_columns.append('standard_energy')
            projected_asd_magnitude_squared =[]
            bayesian_statistics = []
            for k,v in imap.items():
                all_energies.append(v.to_coherent().power2(2).energy)
                all_columns.append('coherent_' + k)
                all_energies.append(v.power2().to_coherent().energy)
                all_columns.append('incoherent_' + k)
                try:
                    # if possible will will calculate all the bayesian statistics for each pixel
                    if k == 'f_plus':
                        projected_asd_magnitude_squared.append(v.projected_asd_magnitude_squared.value)
                        bayesian_statistics.append('loghbayesian')
                    elif k == 'f_cross':
                        projected_asd_magnitude_squared.append(v.projected_asd_magnitude_squared.value)
                    elif k == 'f_right':
                        projected_asd_magnitude_squared.append(v.projected_asd_magnitude_squared.value)
                        bayesian_statistics.append('loghbayesiancirc')
                    elif k =='f_left':
                        projected_asd_magnitude_squared.append(v.projected_asd_magnitude_squared.value)
                    else:
                        pass
                except:
                    pass
            all_columns.extend(bayesian_statistics)
            if projected_asd_magnitude_squared:
                projected_asd_magnitude_squared =  numpy.asarray(projected_asd_magnitude_squared).flatten(order='F')

            # Just assign the energy attribute of the last sparse maps to be
            # all the coherent and incoherent energies and get all cluster properties
            # with one call to the cluster method. (this method calls fastsparseclusterprop.cpp,
            # which does the heavy lifting)
            tmp_sparse_map = list(v.values())[0]
            all_energies = numpy.asarray(all_energies)
            all_energies = numpy.vstack((all_energies[0] + all_energies[2], all_energies))
            tmp_sparse_map.energy = all_energies
            for connectivity in connectivity_list:
                clusters = tmp_sparse_map.cluster(columns=all_columns, connectivity=connectivity, projected_asd_magnitude_squared=projected_asd_magnitude_squared)
                clusters['connectivity'] = connectivity

                # append the cluster to other clusters from same sky locations
                all_clusters = all_clusters.append(clusters)

    all_clusters = all_clusters.groupby(['connectivity']).apply(lambda x, statistic_column: XCluster(x).supercluster(statistic_column=statistic_column), statistic_column=statistic_column).reset_index(drop=True)

    return all_clusters

def extract_clusters_from_table(table, event_type, **kwargs):
    all_clusters = pandas.DataFrame()
    for fft_length in set(table.cols.dx):
        for phi, theta in set(zip(table.cols.phi, table.cols.theta)):
            # Obtain all sparse time frequency maps for this fftlength and sky position
            sparse_maps = [csc_XSparseTimeFrequencyMap.read(row)
                            for row in table.where("""(dx == {0}) & (phi == {1}) & (theta == {2})""".format(fft_length, phi, theta,))]

            # Reformat the above list of projected sparse time frequency maps into
            # a nested dictionary of key : {key1 :value}} where
            # key=projection (i.e. 'f_plus') key1 is detectors
            projected_sparse_maps = XSparseTimeFrequencyMapDict({imap.map_type : XSparseTimeFrequencyMapDict() for imap in sparse_maps})
            for imap in sparse_maps:
                projected_sparse_maps[imap.map_type][imap.name] = imap

            all_energies = []
            all_columns = _default_columns.copy()
            for k,v in projected_sparse_maps.items():
                all_energies.append(v.to_coherent().power2(2).energy)
                all_columns.append('coherent_' + k.decode("utf-8"))
                all_energies.append(v.power2().to_coherent().energy)
                all_columns.append('incoherent_' + k.decode("utf-8"))

            # Just assign the energy attribute of the last sparse maps to be
            # all the coherent and incoherent energies and get all clsuter properities
            tmp_sparse_map = list(v.values())[0]
            tmp_sparse_map.energy = numpy.asarray(all_energies)
            clusters = tmp_sparse_map.cluster(columns=all_columns, **kwargs)

        # append the cluster to other clusters from same sky locations
        all_clusters = all_clusters.append(clusters)

    # super cluster over sky locations and ffitlengths for this event
    all_clusters = all_clusters.supercluster(statistic_column='coherent_f_plus')

    # extract event info
    trigger_info = table._v_pathname.split('/')
    if event_type in ['background', 'onsource']:
        event_info = list(filter(lambda x: 'event' in x, trigger_info))[0]
        internal_time_slide = int(list(filter(lambda x: 'internal_slide' in x, trigger_info))[0].split('_')[-1])
        all_clusters['event'] = event_info
        all_clusters['internal_time_slide'] = internal_time_slide
    else:
        event_info = list(filter(lambda x: 'event' in x, trigger_info))[0]
        waveform_info = list(filter(lambda x: 'waveform' in x, trigger_info))[0]
        injection_scale = float(list(filter(lambda x: 'injection_scale' in x, trigger_info))[0].split('_')[-1].replace('d','.'))
        injection_number = int(list(filter(lambda x: 'injection_number' in x, trigger_info))[0].split('_')[-1])
        all_clusters['event'] = event_info
        all_clusters['waveform'] = waveform_info
        all_clusters['injection_scale'] = injection_scale
        all_clusters['injection_number'] = injection_number

    return all_clusters

def extract_cnn_maps_from_table(table,):
    for fft_length in set(table.cols.dx):
        if fft_length == 0.25:
            continue
        for phi, theta in set(zip(table.cols.phi, table.cols.theta)):
            # Obtain all sparse time frequency maps for this fftlength and sky position
            sparse_maps = [csc_XSparseTimeFrequencyMap.read(row)
                            for row in table.where("""(dx == {0}) & (phi == {1}) & (theta == {2})""".format(fft_length, phi, theta,))]

            # Reformat the above list of projected sparse time frequency maps into
            # a nested dictionary of key : {key1 :value}} where
            # key=projection (i.e. 'f_plus') key1 is detectors
            projected_sparse_maps = XSparseTimeFrequencyMapDict({imap.map_type : XSparseTimeFrequencyMapDict() for imap in sparse_maps})
            for imap in sparse_maps:
                projected_sparse_maps[imap.map_type][imap.name] = imap

            # Determine the size of the input to the CNN
            input_size = list(sparse_maps[0].shape)
            input_size.append( len(list(projected_sparse_maps.keys()))*2)

            # initialize input to all zeros
            data = numpy.zeros(input_size)

            i = 0
            for k,v in projected_sparse_maps.items():
                data[:,:,i] = v.to_coherent().power2(2).todense()
                data[:,:,i+1] = v.power2().to_coherent().todense()
                i +=2
    return data

def extract_energy_maps_from_table(table, **kwargs):

    injection_time = kwargs.pop('injection_time', None)
    window = kwargs.pop('window', None)
    for fft_length in set(table.cols.dx):
        if fft_length == 0.25:
            continue
        for phi, theta in set(zip(table.cols.phi, table.cols.theta)):
            # Obtain all sparse time frequency maps for this fftlength and sky position
            sparse_maps = [csc_XSparseTimeFrequencyMap.read(row)
                            for row in table.where("""(dx == {0}) & (phi == {1}) & (theta == {2})""".format(fft_length, phi, theta,))]

            # assume you will use all pixels as samples unless this is in an injection in whcih casse
            # we only want to provide weight to pixels within some window of the injection
            mask = numpy.ones(sparse_maps[0].tindex.size)
            if injection_time is not None:
                injection_time_indices = numpy.searchsorted(sparse_maps[0].xindex, [injection_time-window, injection_time+window])
                mask[numpy.searchsorted(injection_time_indices,sparse_maps[0].tindex) != 1] = 0
            # Reformat the above list of projected sparse time frequency maps into
            # a nested dictionary of key : {key1 :value}} where
            # key=projection (i.e. 'f_plus') key1 is detectors
            projected_sparse_maps = XSparseTimeFrequencyMapDict({imap.map_type : XSparseTimeFrequencyMapDict() for imap in sparse_maps})
            for imap in sparse_maps:
                projected_sparse_maps[imap.map_type][imap.name] = imap

            # number of rows
            num_rows = len(list(projected_sparse_maps.keys()))*2 + 1
            # Determine the size of the input to the CNN
            input_size = [num_rows]
            input_size.extend(list(sparse_maps[0].energy.shape))

            # initialize input to all zeros
            data = numpy.zeros(input_size)
            data[num_rows-1, :] = mask

            i = 0
            for k,v in projected_sparse_maps.items():
                data[i, :] = v.to_coherent().power2(2).energy
                data[i+1, :] = v.power2().to_coherent().energy
                i +=2
    return data

def undo_slides(x,block_time=256):
    # (blocktime + (cluster_min_time_hanford + external slide livingston) -
    # (event_time + extenral slide livingsto) - half blocktime - internal slide) + 
    # (cluster_min_time_hanford + external slide livingston)
    # min time of cluster hanford
    
    h_min_time = x['min_time_of_cluster']
    external_slide = int(x['event'].split('_')[-1])
    internal_slide = x['internal_time_slide']
    event_time = int(x['event'].split('_')[1])
    return event_time+numpy.mod(h_min_time-event_time-180,block_time) + external_slide
