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
