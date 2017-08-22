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
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWpy.  If not, see <http://www.gnu.org/licenses/>.

"""This module provides the basis for selecting a valid network of detectors to use for an analysis based on some requiremnets for the quality of the on-source data 
"""
import numpy as np

def main(detectorList,
   centerTimesFile,
   timeOffsetsFile,
   cat1_segment_file,
   cat24_segment_file,
   output_file,
   verbosity, windowBeginOffset, windowEndOffset, transientTime, windowLength):
    """
      detectors_considered        cell array of ifos to be considered, e.g.,
                                  ['H1' 'H2' 'L1'], size (1 * nIfo).
      centerTimes                 vector of segment center times, size (nJobs * 1).
      timeOffsets                 array of segment time offsets, 
                                  size(nJobs * nIfo).
      cat1_seglist                structure containing cat1 DQ segments (i.e., 
                                  livetimes). 
                                  cat1_seglist(nIfo).gpsStart()
                                  cat1_seglist(nIfo).duration()
      cat2_seglist                Similar to cat1_seglist but containing cat2 
                                  veto segments (deadtime).
      verbosity                   string used to control verbosity of output,
                                  either 'verbose' or 'quiet'. 
      output_file                  String. Name of file to which we will write
                                   the detector networks. 
      windowBeginOffset            String. Offset [sec] between start of on-source
                                   window and trigger time, e.g., -600. 
      windowEndOffset              String. Offset [sec] between end of on-source 
                                   window and trigger time, e.g., 60.
      transientTime                String. Length of time [sec] discarded at 
                                   beginning and end of each analysis block, e.g., 4.
      windowLength                 String. 2*transientTime+jobsPerWindow*(blockTime-2*transientTime)

     XNETWORKSELECTION is a wrapper function for calling XAPPLYNETWORKTESTS.
     It is intended to be compiled, allowing it to be called from the command
     line or from python scripts such as grb.py.

     XNETWORKSELECTION takes in lists of event times and on-source windows,
     plus lists of science mode livetime and veto deadtime segment files. It  
     uses XAPPLYNETWORKTESTS to determine which detectors satisfy the data
     quality criteria (defined in XAPPLYNETWORKTESTS) to be used for the
     analysis of each event. The list of acceptable detectors for each event
     is written to the file specified by the input argument output_file, with
     one line per event, of the form 'H1L1V1', or 'H1H2', etc. 


      Our criteria for including an ifo in our network is as follows:w
      i) No deadtime incurred by category 1 flags in an interval [-128,+128]s 
         about the GRB trigger time. 
         X-pipeline requires 256s in order to estimate the PSD, category 1 
         flags raised within these 256s will lead to inaccurate/incorrect PSD 
         estimation.
     ii) Less than 9s (5#) deadtime incurred by category 2 flags in the onsource, 
         i.e., an interval [-120,+60]s about the GRB trigger time. 
         We tune our glitch rejection cuts using calculation of the 90# UL on 
         hrss amplitude. If we lose much more than 5# deadtime in our (dummy) 
         onsource we may not be able to accurately estimate the hrss amplitude 
         leading to 90# detection efficiency.
    iii) No deadtime incurred by category 2 flags in an interval [-5,+1]s about 
         the GRB trigger time. Many GRB models predict similar arrival times of 
         EM and GW emission from a GRB. We therefore require that inclusion of 
         an ifo does not kill times close to the GRB trigger time.

    """

    # ---- Number of jobs equals number of centerTimes pass_inited in.
    nJobs = len(centerTimes)

    # ---- Initially assume all ifos FAIL our criteria for all of
    #      the jobs.
    # ---- Note that unlike other pass_init flags this has one element
    #      per job (rather than one per trigger etc.).
    pass_init = np.zeros(nJobs,len(detectorList))

    ###############################################################################
    #                             Apply network test.
    ###############################################################################

    # ---- Loop over jobs, each job represents a choice
    #      of segment centerTimes and timeOffsets.
    for iJob in range(0,nJobs):

       # ---- Loop over ifos we are considering.
       for iIfo in range(0,len(detectorList)):

          #########################################################################
          #                     Unslide centerTimes.
          #########################################################################

          # ---- Before applying vetoSeg cuts we must unslide the triggers.
          unslidCenterTimes[iJob,iIfo] = centerTimes[iJob] + timeOffsets[iJob,iIfo]

          #########################################################################
          #  Check that we have  cat1 flags raised in [-128,+128]s interval about 
          #                         unslid centerTimes.
          #########################################################################

          # ---- Define our [-128,+128]s interval about the (unslid) centerTimes.
          goodBefore = -windowBeginOffset+transientTime
          duration   = windowLength

          # ---- Find intersections between cat1 segments and the [-128,+128]s 
          #      interval about unslid centerTimes.
          # ---- Passing the shorter list to Coincidence2 first
          #      speeds it up. 

          coincOut=Coincidence2(unslidCenterTimes(iJob,iIfo)-goodBefore,duration,
             cat1_seglist(iIfo).gpsStart,
             cat1_seglist(iIfo).duration)

          # ---- We require an info to have cat1 flags (i.e., be in science mode)
          #      for the full duration of the [-128,+128]s.
          if ~isempty(coincOut):
             # ---- If our interval is completely contained within a cat1 segment
             #      we should only have one intersection and it should have the 
             #      duration of our interval.
             if (size(coincOut,1) == 1) and (coincOut(1,2) == duration):
                pass_init[iJob,iIfo] = 1

                if (verboseFlag == 1):
                   print >> sys.stdout, '#s has good cat 1 DQ for job #d \n'.format(detectors_considered[iIfo], iJob )

          #########################################################################
          #    Check for cat2 flags raised in [-5,+1]s interval about unslid 
          #                        centerTimes.
          #########################################################################

          # ---- Define our [-5,+1]s interval about the (unslid) centerTimes.
          goodBefore = 5
          duration   = 6

          # ---- Find intersections between cat2 segments and the [-5,+1]s 
          #      interval about unslid centerTimes.
          # ---- Passing the shorter list to Coincidence2 first
          #      speeds it up. 
          coincOut=Coincidence2(unslidCenterTimes(iJob,iIfo)-goodBefore,duration,
             cat2_seglist(iIfo).gpsStart,
             cat2_seglist(iIfo).duration)

          # ---- If there were any coincidences between our interval and the
          #      vetoSegs we must discard this job. We don't count
          #      coincidences with zero duration (edge overlap)
          if not isempty(coincOut) and any(coincOut[:,2]>0):
             # ---- Setting pass_init for this job to zero.
             pass_init[iJob,iIfo] = 0

             if (verboseFlag == 1):
                print >> sys.stdout, 'Time killed in [-#f,+#f]s interval about center time \n'.format(goodBefore,duration-goodBefore)
                print >> sys.stdout,'startTime     stopTime\n'
                print >> sys.stdout,'--------------------------\n'
                for iCoinc in range(0,len(coincOut[:,1])):
                   print >> sys.stdout,'#9.2f  #9.2f \n'.format(
                      coincOut(iCoinc,1),
                      coincOut(iCoinc,1) + coincOut(iCoinc,2))
                print >> sys.stdout, 'Setting pass_init to zero \n'

          #########################################################################
          #     Find deadtime in [-120,+60]s interval about unslid centerTime. 
          #########################################################################

          # ---- Define our [-120,+60]s interval about the (unslid) centerTimes.
          goodBefore = -windowBeginOffset
          duration   = windowEndOffset-windowBeginOffset

          # ---- Allowed Threshold on deadtime in on source window is 5#
          killedThresh = 0.05*duration

          # ---- Find intersections between cat2 veto segs and our [-120,+60]s
          #      on-source region.
          # ---- Passing the shorter list to Coincidence2 first
          #      speeds it up 
          coincOut=Coincidence2(unslidCenterTimes(iJob,iIfo)-goodBefore,duration,
             cat2_seglist(iIfo).gpsStart,
             cat2_seglist(iIfo).duration)

          # ---- Record times killed by segments. 
          if isempty(coincOut):
             #killedTimes[iJob,iIfo] = [0,0]
             # ---- No times killed by cat2 veto segments.
             liveTimes[iJob,iIfo] = [unslidCenterTimes(iJob,iIfo)-goodBefore,duration]
          else:
             #killedTimes[iJob,iIfo] = [coincOut(:,1),coincOut(:,2)]
             # ---- Some times killed by cat2 veto segments.
             #      If any one ifo has deadtime >= killedThresh
             #      we will discard it straight away.
             # ---- Sum up duration of intersections to find deadtime.
             totDeadTime_cat2[iJob,iIfo] = sum(coincOut[:,2])

             if totDeadTime_cat2[iJob,iIfo] >= killedThresh:
                pass_init[iJob,iIfo] = 0

             liveTimes[iJob,iIfo] = ComplementSegmentList(coincOut[:,1],coincOut[:,2],
                unslidCenterTimes(iJob,iIfo)-goodBefore,
                unslidCenterTimes(iJob,iIfo)-goodBefore+duration)

          # ---- Resliding liveTimes in order to measure how much of
          #      deadtime in our interval.
          #      To do this we must SUBTRACT timeOffsets.
          reslidLiveTimes[iJob,iIfo] = [liveTimes[iJob,iIfo][:,1] -
             timeOffsets(iJob,iIfo),
             liveTimes[iJob,iIfo][:,2]]


       # ---- For current job, find intersection of reslid liveTimes from all ifos
       #      that are still in our network.

       intersectReslidLiveTimes[iJob] = [0, Inf]
       for iIfo in range(0, len(detectors_considered)):
          # ---- We only care about the ifo's livetime if it is still in
          #      our network. 
          if pass_init(iJob,iIfo):
             coincOut = Coincidence2(intersectReslidLiveTimes[iJob][:,1],
                                     intersectReslidLiveTimes[iJob][:,2],
                                     reslidLiveTimes[iJob,iIfo][:,1],
                                     reslidLiveTimes[iJob,iIfo][:,2])

             if isempty(coincOut):
                intersectReslidLiveTimes[iJob] = [0,0]
             else:
                intersectReslidLiveTimes[iJob] = [coincOut[:,1], coincOut[:,2]]


       # ---- For current job, find time killed by cat2 veto flags.
       finalKilledTimes[iJob] = ComplementSegmentList(
          intersectReslidLiveTimes[iJob][:,1],
          intersectReslidLiveTimes[iJob][:,2],
          centerTimes[iJob]-goodBefore,centerTimes[iJob]-goodBefore+duration)
       deadtime[iJob] = sum(finalKilledTimes[iJob][:,2])

       ############################################################################
       #         Discard jobs with >=9s deadtime in [-120,+60]s interval 
       #                     about centerTime. 
       ############################################################################

       if (verboseFlag == 1):
          print >> sys.stdout, 'Total time killed in [' +  num2str(windowBeginOffset) + ',' + num2str(windowEndOffset) + ']s interval about centre time: #9.2f \n'.format(deadtime(iJob))
          print >> sys.stdout,'startTime     stopTime\n'
          print >> sys.stdout,'--------------------------\n'

          for iCoinc in range(0, len(finalKilledTimes[iJob][:,1])):
             print >> sys.stdout,'#9.2f  #9.2f \n'.format(
                finalKilledTimes[iJob][iCoinc,1],
                finalKilledTimes[iJob][iCoinc,1] +
                finalKilledTimes[iJob][iCoinc,2])

       # ---- Check total duration of killed times
       if deadtime(iJob) >= killedThresh:

          warning(['Killing network due to cat2 flags in on-source. '
             'We should investigate this network further'])

          networkCell[iJob] = ['X']

          for iIfo in range(0, len(detectors_considered)):
            if pass_init(iJob,iIfo):
              networkCell[iJob] = [networkCell[iJob], detectors_considered[iIfo]]
            pass_init[iJob,iIfo] = 0


          if (verboseFlag == 1):
             print >> sys.stdout,'More than #ds killed, setting pass_init to zero \n'.fromat(killedThresh)
       else:
          networkCell[iJob] = []

       #############################################################################
       #                     Construct network cell array.
       #############################################################################

       # ---- Loop over ifos we are considering.
       for iIfo in range(0,len(detectors_considered)):
          if pass_init[iJob,iIfo]:
             networkCell[iJob] = [networkCell[iJob], detectors_considered[iIfo]]
