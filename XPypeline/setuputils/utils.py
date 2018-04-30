# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
#
# This file is part of XPypeline.
#
# GWpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWpy.  If not, see <http://www.gnu.org/licenses/>.

"""This module contains utility functions for setUpJobs 
"""
from gwpy.segments import Segment, DataQualityDict, SegmentList, DataQualityFlag
import operator

def query_veto_definer_file(ifo, start, end, cp):
    # Obtain segments that are analysis ready
    analysis_ready = DataQualityFlag.query('{0}:DMT-ANALYSIS_READY:1'.format(ifo), start, end)
    # Define the start and stop time of analysis ready
    analysis_ready_start = analysis_ready.active.extent()[0]
    analysis_ready_end = analysis_ready.active.extent()[1]
    # Query for vetos found during this time.
    vdf = DataQualityDict.from_veto_definer_file(cp.get('segfind', 'veto-file'), analysis_ready_start, analysis_ready_end)
    # Populate only for analysis ready segments.
    vdf.populate(segments=analysis_ready.active)

    return vdf


def filter_for_cat_type(vdf, ifo, cat):
    return reduce(operator.or_, [f.active for f in vdf.values() if f.ifo == ifo and f.category in cat], SegmentList())
   

def validate_segments(ifos, start, end, cp, trigger_time=None):
    """determine analysis ready segments during requested analysis time
    Parameters
    ----------
    ifos : `str`
        list of ifos used in X-Pipeline analysis

    start : `float`, :class:`~gwpy.time.LIGOTimeGPS`

    end : `float`, `~gwpy.time.LIGOTimeGPS`

    cp : `object` ConfigParser object

    trigger_time

    Returns
    -------
    DataQualityDict : ~gwpy.segments.DataQualityDict`
    """
    analysis_seg_files = []

    # If simulating noise skip all segment, veto 
    # and network validation checks. The on source and off source is
    # simulated and therefore has no DQ issues.
    if cp.has_option('parameters','makeSimulatedNoise'):
        for ifo in ifos:
            print 'Making simulated noise, creating temp segment file'
            f = open('segments_{0}.txt'.format(ifo),'w')
            f.write('0 {0} {1} {2}'.format(start, end ,end - start))
            f.close()
            analysis_seg_files.append('segments_{0}.txt'.format(ifo))
    else:
        for ifo in ifos:
            if cp.has_option(ifo,'segment-list'):
                if not os.path.isfile(cp.get(ifo,'segment-list')):
                    raise ValueError('Please uncomment the '
                                     'the segment file in ini file '
                                     'as it does not exist. '
                                     'If you want to use a '
                                     'supplied file please provide one.')
                else:
                    analysis_seg_files.append(cp.get(ifo,'segment-list'))
            else:
                # Query for veto definer file
                vdf = query_veto_definer_file(ifo, start, end, cp) 

                # Filter for cat1 vetos
                segs = filter_for_cat_type(vdf, ifo, [1])

                # ---- Write out cat1 veto to text file.
                filename_cat1 = "input/" + ifo +  "-veto-cat1.txt"
                segs.write(filename_cat1)

                # Compute analysis ready segments in order to
                # subtract out cat1 vetos
                analysis_ready = DataQualityFlag.query('{0}:DMT-ANALYSIS_READY:1'.format(ifo), start, end)

                # Subtract cat 1 veto from analysis_ready
                analysis_ready_minus_cat1 = analysis_ready.active - segs
                # Save new segment list to file
                filename_analysis_ready_minus_cat1 = "input/" + ifo + "_science_cat1.txt"
                analysis_ready_minus_cat1.write(filename_analysis_ready_minus_cat1)
                analysis_seg_files.append(filename_analysis_ready_minus_cat1)

    return analysis_seg_files

def validate_network():

    return


def validate_vetos(ifos, start, end, cp):
    """determine vetos during requested analysis time
    Parameters
    ----------
    ifos : `str`
	list of ifos used in X-Pipeline analysis

    start : `float`, :class:`~gwpy.time.LIGOTimeGPS`

    end : `float`, `~gwpy.time.LIGOTimeGPS`

    cp : `object` ConfigParser object

    Returns
    -------
    DataQualityDict : ~gwpy.segments.DataQualityDict`
    """
    veto_seg_files = []

    # If simulating noise skip all segment, veto 
    # and network validation checks. The on source and off source is
    # simulated and therefore has no DQ issues.
    if cp.has_option('parameters','makeSimulatedNoise'):
        for ifo in ifos:
            veto_seg_files.append("None")
    else:
        for ifo in ifos:
            if cp.has_option(ifo,'veto-list'):
                if not os.path.isfile(cp.get(ifo,'veto-list')):
                    raise ValueError('Please uncomment the '
                                     'the veto file in ini file '
                                     'as it does not exist. '
                                     'If you want to use a '
                                     'supplied file please provide one.')
                else:
                    veto_seg_files.append(cp.get(ifo,'veto-list'))
            else:
                # Query for veto definer file
                vdf = query_veto_definer_file(ifo, start, end, cp)

                # Filter for cat24 vetos
                cat = [2, 4]
                for iCat in cat:
                    segs = filter_for_cat_type(vdf, ifo, [iCat])
                    filename = "input/" + ifo +  "-veto-cat{0}.txt".format(iCat)
                    segs.write(filename)

                segs = filter_for_cat_type(vdf, ifo, cat)
                filename = "input/" + ifo + "_cat24veto.txt"
                segs.write(filename)
                veto_seg_files.append(filename)
                
    return veto_seg_files


def make_chunks(segment, length, overlap=0):
    """ 
    Divides the science segment into chunks of length seconds overlapped by 
    overlap seconds. 
    when generating chunks 
    """
    chunks = []
    for iseg in segment:
        segtemp = SegmentList([iseg])
        time_left = abs(segtemp)
        start = segtemp.extent()[0]
        increment = length - overlap
        while time_left >= length:
            end = start + length
            chunks.append(Segment(start, end))
            start += increment
            time_left -= increment

    return SegmentList(chunks)
