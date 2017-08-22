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

"""This module provides function to help write param files needed for xdetection 
"""


def write_onsource(filename, outputType, frameCacheAll, skyPositionList, likelihoodTypeStr, parameters, onSourceEndOffset, cp):

    # ---- Parameters file for on-source analysis.
    fParam=open(filename, 'w')
    # ---- First write framecache file, channel file, event file, and sky position.
    fParam.write('channelFileName:input/channels.txt' + '\n')
    fParam.write('frameCacheFile:' + frameCacheAll + '\n')
    fParam.write('eventFileName:input/event_on_source.txt' + '\n')
    fParam.write('skyPositionList:' + skyPositionList + '\n')
    fParam.write('skyCoordinateSystem:earthfixed' + '\n')
    fParam.write('likelihoodtype:' + likelihoodTypeStr + '\n')
    # ---- Now write all of the other parameters from the parameters section.
    #      We ignore the likelihoodType_* lines since this is handled above.
    for i in range(0,len(parameters)) :
        if not(parameters[i].startswith("likelihoodtype")):
            value = cp.get('parameters',parameters[i])
            if parameters[i] == "onsourceendoffset":
                fParam.write('onsourceendoffset:' + str(onSourceEndOffset) + '\n')
            elif parameters[i] == "circtimeslidestep":
                continue
            else:
                fParam.write(parameters[i] + ':' + value + '\n')
    if outputType == 'seedless' :
        fParam.write('seedlessparams:input/seedless_onsource.txt\n')
    fParam.close()


def write_ul_source(filename, outputType, frameCacheAll, skyPositionList, likelihoodTypeStr, parameters, onSourceEndOffset, cp):

    # ---- Parameters file for ul-source analysis.
    fParam=open(filename, 'w')
    # ---- First write framecache file, channel file, event file, and sky position.
    fParam.write('channelFileName:input/channels.txt' + '\n')
    fParam.write('frameCacheFile:' + frameCacheAll + '\n')
    fParam.write('eventFileName:input/event_on_source.txt' + '\n')
    fParam.write('skyPositionList:' + skyPositionList + '\n')
    fParam.write('skyCoordinateSystem:earthfixed' + '\n')
    fParam.write('likelihoodtype:' + likelihoodTypeStr + '\n')
    # ---- Now write all of the other parameters from the parameters section.
    #      We ignore the likelihoodType_* lines since this is handled above.
    for i in range(0,len(parameters)) :
        if not(parameters[i].startswith("likelihoodtype") or parameters[i].endswith("ecimateRate")):
            value = cp.get('parameters',parameters[i])
            if parameters[i] == "onsourceendoffset":
                fParam.write('onsourceendoffset:' + str(onSourceEndOffset) + '\n')
            elif parameters[i] == "circtimeslidestep":
                continue
            else:
                fParam.write(parameters[i] + ':' + value + '\n')
    fParam.write('decimateRate:100\n')
    fParam.write('superDecimateRate:100\n')
    if outputType == 'seedless' :
        fParam.write('seedlessparams:input/seedless_ul.txt\n')
    fParam.close()


def write_offsource(filename, outputType, frameCacheAll, skyPositionList, likelihoodTypeStr, parameters, onSourceEndOffset, cp):
    # ---- Parameters file for off-source analysis.
    fParam=open(filename, 'w')
    # ---- First write framecache file, channel file, event file, and sky position.
    fParam.write('channelFileName:input/channels.txt' + '\n')
    fParam.write('frameCacheFile:' + frameCacheAll + '\n')
    fParam.write('eventFileName:input/event_off_source.txt' + '\n')
    fParam.write('skyPositionList:' + skyPositionList + '\n')
    fParam.write('skyCoordinateSystem:earthfixed' + '\n')
    fParam.write('likelihoodtype:' + likelihoodTypeStr + '\n')
    # ---- Now write all of the other parameters from the parameters section.
    #      We ignore the likelihoodType_* lines since this is handled above.
    for i in range(0,len(parameters)) :
        if not(parameters[i].startswith("likelihoodtype")):
            value = cp.get('parameters',parameters[i])
            if parameters[i] == "onsourceendoffset":
                fParam.write('onsourceendoffset:' + str(onSourceEndOffset) + '\n')
            else:
                fParam.write(parameters[i] + ':' + value + '\n')
    if outputType == 'seedless' :
        fParam.write('seedlessparams:input/seedless_offsource.txt\n')
    fParam.close()


def write_waveforms(filename, outputType, frameCacheAll, skyPositionList, likelihoodTypeStr, parameters, onSourceEndOffset, cp, injectionScalesList, disableFastInjections, waveform):
              fParam=open(filename, 'w')
              # ---- First write framecache file, channel file, event file, and sky position.
              fParam.write('channelFileName:input/channels.txt' + '\n')
              fParam.write('frameCacheFile:' + frameCacheAll + '\n')
              fParam.write('eventFileName:input/event_inj_source.txt' + '\n')
              fParam.write('catalogdirectory:input/' + '\n')
              fParam.write('skyPositionList:' + skyPositionList + '\n')
              fParam.write('skyCoordinateSystem:earthfixed' + '\n')
              fParam.write('likelihoodtype:' + likelihoodTypeStr + '\n')
              # ---- Now write all of the other parameters from the parameters section.
              #      We ignore the likelihoodType_* lines since this is handled above.
              for i in range(0,len(parameters)) :
                  if not(parameters[i].startswith("likelihoodtype")) :
                      value = cp.get('parameters',parameters[i])
                      if parameters[i] == "outputtype"  and value == "clusters" and not(disableFastInjections):
                          fParam.write('outputtype:injectionclusters\n')
                      elif parameters[i] == "onsourceendoffset":
                          fParam.write('onsourceendoffset:' + str(onSourceEndOffset) + '\n')
                      elif parameters[i] == "circtimeslidestep":
                          continue
                      else:
                          fParam.write(parameters[i] + ':' + value + '\n')

              # ---- Write simulations info.  The specified injection file
              #      will be written later.
              fParam.write('injectionFileName:input/injection_' + waveform + '.txt' + '\n')
              fParam.write('injectionScale:' + str(injectionScalesList[:]) + '\n')
              # fParam.write('catalogdirectory:' + XPIPELINE_ROOT  + '/waveforms\n')
              if outputType == 'seedless' :
                  fParam.write('seedlessparams:input/seedless_simulation.txt\n')
              fParam.close()

def write_mdc(filename, outputType, frameCacheAll, skyPositionList, likelihoodTypeStr, parameters, onSourceEndOffset, cp, injectionScalesList, disableFastInjections, waveform):
    # ---- Read [mdc] sets and injection scales.  If mdc_sets is empty or specifies
    #      unknown MDC sets then the script will have already exited when trying to
    #      write the mdcchannel file above.
    mdc_setsList = cp.get('mdc','mdc_sets')
    mdc_sets = mdc_setsList.split(',')
    # ---- Write one parameters file for each (mdc set, injection scale) pair.
    for iMDC in mdc_sets :
        if cp.has_section(iMDC) & cp.has_option(iMDC,'injectionScales') :
            # ---- This check lets you specify different injection scales for
            #      each MDC set.
            injectionScalesList = cp.get(iMDC,'injectionScales')
            injectionScales = injectionScalesList.split(',')
        else:
            # ---- Otherwise, use the injection scales specified in
            #      the [injection] section.
            injectionScalesList = cp.get('injection','injectionScales')
            injectionScales = injectionScalesList.split(',')

        # ---- Write a separate parameter file for each injection scale.
        scale_counter = 0
        for injectionScale in injectionScales :
            f=open(filename + iMDC + "_" + str(scale_counter) + ".txt", 'w')
            # ---- First write framecache file, channel file, event file, and sky position.
            f.write('channelFileName:input/channels.txt' + '\n')
            f.write('frameCacheFile:' + frameCacheAll + '\n')
            f.write('eventFileName:input/event_'+ iMDC + '.txt' + '\n')
            f.write('skyPositionList:' + skyPositionList + '\n')
            f.write('skyCoordinateSystem:earthfixed' + '\n')
            f.write('likelihoodtype:' + likelihoodTypeStr + '\n')
            # ---- Now write all of the other parameters from the parameters section.
            #      We ignore the likelihoodType_* lines since this is handled above.
            for i in range(0,len(parameters)) :
                if not(parameters[i].startswith("likelihoodtype")):
                    value = cp.get('parameters',parameters[i])
                    if parameters[i] == "outputtype"  and value == "clusters" and not(disableFastInjections):
                        f.write('outputtype:injectionclusters\n')
                    elif parameters[i] == "onsourceendoffset":
                        f.write('onsourceendoffset:' + str(onSourceEndOffset) + '\n')
                    elif parameters[i] == "circtimeslidestep":
                        continue
                    else:
                        f.write(parameters[i] + ':' + value + '\n')
            # ---- Write mdc info.
            f.write('mdcChannelFileName:input/channels_' + iMDC + '.txt' + '\n')
            f.write('injectionFileName:input/injection_' + iMDC + '.txt' + '\n')
            f.write('injectionScale:' + str(injectionScale) + '\n')
            if outputType == 'seedless' :
                f.write('seedlessparams:input/seedless_simulation.txt\n')
            f.close()
            scale_counter = scale_counter + 1

    # ---- Status message.
    print >> sys.stdout, "Writing Matlab-formatted MDC channel files "\
        "and framecache files... "

    # ---- Get list of MDC sets to process.
    mdc_setsList = cp.get('mdc','mdc_sets')
    mdc_sets = mdc_setsList.split(',')
    # print >> sys.stdout, "mdc_sets:", mdc_sets

    # ---- We will create a list of all the log files for each mdc set
    mdc_log_files = []

    # ---- Make MDC channels file for each set.
    for setIdx in range(len(mdc_sets)) :
        set = mdc_sets[setIdx]
        # print >> sys.stdout, "set:", set
        if cp.has_section(set) :
            # ---- Read channel parameters for this mdc set.
            mdcChannelListLine = cp.get(set,'channelList')
            mdcChannelList = mdcChannelListLine.split(',')
            mdcFrameTypeListLine = cp.get(set,'frameTypeList')
            mdcFrameTypeList = mdcFrameTypeListLine.split(',')
            numberOfChannels = cp.get(set,'numberOfChannels')
            # ---- Keep only info for detectors requested for the analysis.
            mdcChannel = []
            mdcFrameType = []
            for jj in range(0,len(keepIndex)):
                mdcChannel.append(mdcChannelList[keepIndex[jj]])
                mdcFrameType.append(mdcFrameTypeList[keepIndex[jj]])

            # ---- For each detector, write the
            #      corresponding channel name and frame type to a file.
            f=open('input/channels_' + set + '.txt', 'w')
            for i in range(0,len(detector)) :
                # ---- if we are doing a MOCK analysis we must not add on
                #      grb_name to mdcFrameType.
                if not(mockAnalysis):
                    # ---- mdcFrameTypes have format: H1_SGC554Q8d9_ON_GRB060211B
                    # ---- from params file mdcFrameTypes of format: H1_SGC554Q8d9_ON_
                    #      we append last underscore if missing and append grb_name.
                    if not(mdcFrameType[i].endswith('_')):
                        mdcFrameType[i] = mdcFrameType[i] + '_'
                    mdcFrameType[i] = mdcFrameType[i] + grb_name

                f.write(detector[i] + ':' + mdcChannel[i] + ' ' + mdcFrameType[i] + '\n')

            f.close()
            # ---- end loop over ifo

            # -------------------------------------------------------------------------
            #                         Find MDC frames.
            # -------------------------------------------------------------------------

            # ---- Check to see if mdc frame cache has been specified in config file.
            try:
                mdcFrameCache = cp.get(set,'frameCacheFile')
            except:
                print >> sys.stdout, "Warning: No frameCacheFile file specified " \
                    "in [" + set + "] section of configuration file."
                print >> sys.stdout, "         A frameCache for the ifo data file " \
                    "will be generated automatically."
                mdcFrameCache = None

            if mdcFrameCache:
                # ---- Check that the frame cache file specified actually exists.
                if not os.path.isfile(mdcFrameCache):
                    print >> sys.stderr,"Error: non existant framecache: ", mdcFrameCache
                    sys.exit(1)

                # ---- If the specified mdc frame cache exists then concat it
                #      other frame caches.
                command = 'cat ' + mdcFrameCache + '  >> ' + frameCacheAll
                os.system(command)

            # ---- If the mdcFrameCache was not specified we generate it here.
            #      Note that we look for frames using the virtual detector name
            #      to ensure the MDC injection data is consistent with the 
            #      detector being simulated.
            else:
                for i in range(0,len(detector)) :
                    # ---- generate frame cache file for MDCs
                    if not mdcFrameCache:
                        # ---- Status message.
                        print >> sys.stdout, "Writing MDC framecache file for " \
                            "mdcset: " + set + ", ifo: " + detector[i] + " ..."
                        # ---- Clear away any pre-existing mdcframecache files.
                        os.system('rm -f mdcframecache_temp.txt')

                        # ---- Construct dataFind command.
                        dataFindCommand = ' '.join([datafind_exec,
                        "--server", datafind_server,
                        "--observatory",detector[i][0],
                        "--type",mdcFrameType[i],
                        "--gps-start-time", str(start_time),
                        "--gps-end-time",str(end_time),
                        "--url-type file",
                        "--lal-cache",
                        " > mdclalcache.txt"])
                        # ---- Issue dataFind command.
                        print "calling dataFind:", dataFindCommand
                        os.system(dataFindCommand)
                        print "... finished call to dataFind."

                        # ---- Convert lalframecache file to readframedata format.
                        print "calling convertlalcache:"
                        os.system('convertlalcache.pl mdclalcache.txt mdcframecache_temp.txt')
                        os.system('cat mdcframecache_temp.txt >> ' + frameCacheAll)
                        print "... finished call to convertlalcache."
                        # ---- Clean up.
                        os.system('rm -f mdcframecache_temp.txt mdclalcache.txt')

                        # ---- Status message.
                        print >> sys.stdout, "... finished writing MDC framecache file."
                        print >> sys.stdout

            #-------------------------------------------------------------------------
            #   Create a list of the mdc log files and create an mdc segment list.
            #-------------------------------------------------------------------------

            if not(mdc_path):
                print >> sys.stderr, "Error: mdc_path must be specified if we are using mdcs."
                print >> sys.stderr, "Use --mdc-path to specify it."
                sys.exit(1)

            # ---- check dir names end in '/'
            if not(mdc_path.endswith('/')):
                mdc_path = mdc_path + '/'
            mdc_grbdir = grb_name + '/'

            # ---- mdcFrameType[i] will have name in format H1_SGC100Q8d9_ON_GRB051105
            #      logs are in dir structure:
            #      /data/node5/eharstad/ExtTrigMDC/GRB051105/SGC1000Q8d9_ON/logs/
            #      ExtTrigMDC-SGC1000Q8d9_ON_GRB051105-815207086-256-Log.txt

            # ---- Remove <ifo>_ and _GRB* part from mdcFrameType.
            if mockAnalysis:
                strList0    = mdcFrameType[i].rsplit('-')
                mdcType     = strList0[1]
            else:
                strList0    = mdcFrameType[i].rsplit('_')
                mdcType     = strList0[1] + '_' + strList0[2]

            mdc_log_glob = mdc_path + grb_name + '/'  + mdcType +  '/logs/ExtTrigMDC-' + mdcType + '*.txt'
            print >> sys.stdout, ''
            print >> sys.stdout, 'Getting mdc log files from ' + mdc_log_glob
            mdc_log_files.append(glob.glob(mdc_log_glob))

            # ---- Report error if no log files found.
            if not(len(mdc_log_files[setIdx])):
                print >> sys.stderr, 'Error: No ExtTrigMDC*.txt log files in ' + mdc_log_glob
                sys.exit(1)
            # ---- If these are on source mdcs check we only have one log file.
            if mdcType.count('ON'):
                if len(mdc_log_files[setIdx]) > 1:
                    print >> sys.stderr, 'Error: We only expect one ON source ExtTrigMDC*.txt log file in ' + mdc_log_glob
                    sys.exit(1)

            # ---- Extract start time of each mdc block from its log file name
            #      use this to create a segment list which we later use to
            #      figure out which mdc blocks have good data quality flags.
            mdc_segs = []
            for fileIdx in range(0,len(mdc_log_files[setIdx])):
                # ---- We need to extract start times for each log file.
                # ---- We are in a loop over mdc sets.
                strList = mdc_log_files[setIdx][fileIdx].rsplit('-')

                # ---- Choose element of strList counting back from end of strList in case
                #      dir names contain a hyphen which would throw off our counting.
                mdc_start_time = int(strList[len(strList)-3])
                mdc_block_time = int(strList[len(strList)-2])
                mdc_segs.append([mdc_start_time,mdc_block_time])

            # ---- Sort our segments on mdc_start_time.
            mdc_segs_sorted=sorted(mdc_segs, key=operator.itemgetter(0))

            # ---- Write mdc segtment list.
            fmdcseg=open('input/segment_' + set + '.txt', 'w')
            for segIdx in range(0,len(mdc_segs_sorted)):
                mdc_start_time = mdc_segs_sorted[segIdx][0]
                mdc_block_time = mdc_segs_sorted[segIdx][1]
                mdc_end_time   = mdc_start_time + mdc_block_time
                time_range_string = str(segIdx) + ' ' + \
                    str(mdc_start_time) + ' ' + str(mdc_end_time)  + \
                    ' ' + str(mdc_block_time) + '\n'
                fmdcseg.write(time_range_string)
            fmdcseg.close()

        else:
            print >> sys.stdout, "Error: MDC set ", set, \
                " is not defined in the parameters file.  Exiting."
            print >> sys.stdout
            sys.exit(1)

    # ---- Status message.
    print >> sys.stdout, "... finished writing MDC channel and frame files.   "
    print >> sys.stdout

    # ---- Status message.
    print >> sys.stdout, "Writing MDC event files ...          "

    # ---- Get list of MDC sets to process.
    mdc_setsList = cp.get('mdc','mdc_sets')
    mdc_sets = mdc_setsList.split(',')

    for set in mdc_sets:
        # ---- Read mdc segments into a "ScienceData" object.
        mdc_segment = pipeline.ScienceData()
        mdc_segment.read('input/segment_' + set + '.txt', blockTime )

        # ---- Now get the intersection of all of the detector segment lists with this
        #      on-source list. Note that the segment lists mst be sorted by start time
        #      for intersection to work properly
        mdc_coincidence_segment = copy.deepcopy(mdc_segment)
        for det in full_segment_list:
            mdc_coincidence_segment.intersection(det)

        # ---- Do some checking... check we have some MDCs in coinc...

        # ---- At this point, the ScienceData object mdc_segment contains the
        #      mdc segments.  Write this to the mdc event file.
        f=open('input/event_' + set + '.txt', 'w')
        # ---- We also rewrite our segment files so they only contain the segs
        #      we are going to analyse
        fmdcseg=open('input/segment_' + set + '.txt','w')
        # ---- loop over coinc segs
        for i in range(mdc_coincidence_segment.__len__()):
            duration = mdc_coincidence_segment.__getitem__(i).end() - \
                mdc_coincidence_segment.__getitem__(i).start()
            if (duration == blockTime):
                time_range_string = str(mdc_coincidence_segment.__getitem__(i).start() \
                    + blockTime / 2) + '\n'
                f.write(time_range_string)
                time_range_string = str(i) + ' ' + \
                    str(mdc_coincidence_segment.__getitem__(i).start()) \
                    + ' ' + str(mdc_coincidence_segment.__getitem__(i).end())  \
                    + ' ' + str(mdc_coincidence_segment.__getitem__(i).end() \
                    -mdc_coincidence_segment.__getitem__(i).start()) + '\n'
                fmdcseg.write(time_range_string)
            elif (duration > blockTime):
                # ---- this should not be possible
                print >> sys.stderr,"Error: something has gone wrong in creating MDC segment list "
                sys.exit(1)

        f.close()
        fmdcseg.close()

    # ---- Status message.
    print >> sys.stdout, "... finished writing MDC event files."
    print >> sys.stdout

