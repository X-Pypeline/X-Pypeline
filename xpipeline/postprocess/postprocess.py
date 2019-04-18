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

def xapplyratiocuts(triggers, ePlusIndex,eCrossIndex,eNullIndex,
                    iPlusIndex,iCrossIndex,iNullIndex,
                    vetoPlusRange,vetoCrossRange,vetoNullRange,
                    FAR_Tuning,typeOfCutPlus,typeOfCutCross,
                    detection_statistic_column_name):
    
    if float(FAR_Tuning).is_integer():
        FAR_Tuning = (FAR_Tuning/100)

    # Calculate all Ratio Values
    plusRatioEoverI = numpy.log(triggers[ePlusIndex] / triggers[iPlusIndex])
    crossRatioEoverI = numpy.log(triggers[eCrossIndex] / triggers[iCrossIndex])
    plusRatioIoverE = numpy.log(triggers[iPlusIndex] / triggers[ePlusIndex])
    crossRatioIoverE = numpy.log(triggers[iCrossIndex] / triggers[eCrossIndex])

    if iNullIndex == 0:
        triggers['nullenergy'] = numpy.ones(plusRatioEoverI.size) 
        nullRatioIoverE = triggers['nullenergy']  
    else:
        nullRatioIoverE = numpy.log(triggers[iNullIndex] / triggers[eNullIndex])

    # Depending on wether you are requesting a two sided or one sided cut.
    # construct specific ratio and cut arrays.
    if typeOfCutPlus == 'twosided':

        ratioCutsPlus = numpy.log(vetoPlusRange)
        ratioCutsCross = numpy.log(vetoCrossRange)
        ratioCutsNull = numpy.log(vetoNullRange)
        ratioArrayPlus = [plusRatioEoverI,plusRatioIoverE]
        ratioArrayCross = [crossRatioEoverI,crossRatioIoverE]
        ratioArrayNull = nullRatioIoverE
        sidedCut = 2

    else:

        ratioCutsPlus = numpy.log(abs(vetoPlusRange))
        ratioCutsCross = numpy.log(abs(vetoCrossRange))
        ratioCutsNull = numpy.log(abs(vetoNullRange))
        if typeOfCutPlus == 'EoverI':
            ratioArrayPlus = plusRatioEoverI
        else:
            ratioArrayPlus = plusRatioIoverE

        if typeOfCutCross == 'EoverI':
            ratioArrayCross = crossRatioEoverI
        else:
            ratioArrayCross = crossRatioIoverE

        ratioArrayNull = nullRatioIoverE
        sidedCut = 1

    # Determine how many indepdent trials there are
    background_trials = triggers.groupby(['event','internal_time_slide'])
    number_of_trials = len(background_trials)
    indexFARBackground = int(numpy.ceil(number_of_trials* FAR_Tuning))
    loudestBackgroundRatioJob = []
    trial_name = []
    

    for key, item in background_trials:
        triggers_from_this_trial = background_trials.get_group(key).index
        
        detection_statistic = triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy()
        
        ratioArrayPlusJobNumber = ratioArrayPlus.loc[triggers_from_this_trial].to_numpy()
        ratioArrayCrossJobNumber = ratioArrayCross.loc[triggers_from_this_trial].to_numpy()
        ratioArrayNullJobNumber = ratioArrayNull.loc[triggers_from_this_trial].to_numpy()

        # Reshape Calculated Ratios and Cut Values such that they are the same size
        # Specifically we will make it so that the shape is
        # (# of cuts to apply by number of triggers), this means that to find
        # that means that ratioArray we need to repeat the triggers by # of cuts
        # and for veto* we need to repeat the threshold to be applied by how many triggers said threshold needs
        # to be applied to.

        ratioArrayTempPlus = numpy.repeat(numpy.atleast_2d(ratioArrayPlusJobNumber),
                                          len(ratioCutsPlus),
                                          axis=0)
        ratioArrayTempCross = numpy.repeat(numpy.atleast_2d(ratioArrayCrossJobNumber),
                                           len(ratioCutsCross),
                                           axis=0)
        ratioArrayTempNull = numpy.repeat(numpy.atleast_2d(ratioArrayNullJobNumber),
                                          len(ratioCutsNull),
                                          axis=0)

        vetoPlusRep = numpy.kron(numpy.atleast_2d(ratioCutsPlus).T,
                                      numpy.ones((sidedCut, len(ratioArrayPlusJobNumber))))
        vetoCrossRep = numpy.kron(numpy.atleast_2d(ratioCutsCross).T,
                                      numpy.ones((sidedCut, len(ratioArrayCrossJobNumber),)))
        vetoNullRep = numpy.kron(numpy.atleast_2d(ratioCutsNull).T,
                                     numpy.ones((1, len(ratioArrayNullJobNumber),)))

        # Determine what clusters passed all the Ratio cuts
        ratioPassCut = ((ratioArrayTempPlus  > vetoPlusRep) &
                        (ratioArrayTempCross  > vetoCrossRep) & 
                        (ratioArrayTempNull > vetoNullRep))

        # Find Surviving Offsource
        backGroundArray = ratioPassCut*numpy.repeat(numpy.atleast_2d(detection_statistic),
                                                   len(ratioCutsPlus),
                                                   axis=0)
        
        # Find loudest surviving trigger for this trial for the grid of applied cuts
        # to do this we find the max over the columns
        #(remember each row represented an applied a cut)
        loudestBackgroundRatioJob.append(numpy.max(backGroundArray,1))
        trial_name.append(key)

    loudestBackgroundRatioJob = numpy.vstack(loudestBackgroundRatioJob)
    loudestBackgroundRatioJob.sort(axis=0)

    return loudestBackgroundRatioJob[indexFARBackground,:]

def xapplyratiocutsinjections(triggers, ePlusIndex,eCrossIndex,eNullIndex,
                              iPlusIndex,iCrossIndex,iNullIndex,
                              vetoPlusRange,vetoCrossRange,vetoNullRange,
                              loudestBackgroundRatioJob,typeOfCutPlus,typeOfCutCross,
                              detection_statistic_column_name,):

    # Calculate all Ratio Values
    plusRatioEoverI = numpy.log(triggers[ePlusIndex] / triggers[iPlusIndex])
    crossRatioEoverI = numpy.log(triggers[eCrossIndex] / triggers[iCrossIndex])
    plusRatioIoverE = numpy.log(triggers[iPlusIndex] / triggers[ePlusIndex])
    crossRatioIoverE = numpy.log(triggers[iCrossIndex] / triggers[eCrossIndex])

    if iNullIndex == 0:
        triggers['nullenergy'] = numpy.ones(plusRatioEoverI.size) 
        nullRatioIoverE = triggers['nullenergy']  
    else:
        nullRatioIoverE = numpy.log(triggers[iNullIndex] / triggers[eNullIndex])

    # Depending on wether you are requesting a two sided or one sided cut.
    # construct specific ratio and cut arrays.
    if typeOfCutPlus == 'twosided':

        ratioCutsPlus = numpy.log(vetoPlusRange)
        ratioCutsCross = numpy.log(vetoCrossRange)
        ratioCutsNull = numpy.log(vetoNullRange)
        ratioArrayPlus = [plusRatioEoverI,plusRatioIoverE]
        ratioArrayCross = [crossRatioEoverI,crossRatioIoverE]
        ratioArrayNull = nullRatioIoverE
        sidedCut = 2

    else:

        ratioCutsPlus = numpy.log(abs(vetoPlusRange))
        ratioCutsCross = numpy.log(abs(vetoCrossRange))
        ratioCutsNull = numpy.log(abs(vetoNullRange))
        if typeOfCutPlus == 'EoverI':
            ratioArrayPlus = plusRatioEoverI
        else:
            ratioArrayPlus = plusRatioIoverE

        if typeOfCutCross == 'EoverI':
            ratioArrayCross = crossRatioEoverI
        else:
            ratioArrayCross = crossRatioIoverE

        ratioArrayNull = nullRatioIoverE
        sidedCut = 1
    
    injection_trials = triggers.groupby(['injection_scale','injection_number'])
    all_ratio_boolean = []
    all_indices = []

    for key, item in injection_trials:
        triggers_from_this_trial = injection_trials.get_group(key).index
        
        detection_statistic = triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy()
        ratioArrayPlusJobNumber = ratioArrayPlus.loc[triggers_from_this_trial].to_numpy()
        ratioArrayCrossJobNumber = ratioArrayCross.loc[triggers_from_this_trial].to_numpy()
        ratioArrayNullJobNumber = ratioArrayNull.loc[triggers_from_this_trial].to_numpy()

        # Reshape Calculated Ratios and Cut Values such that they are the same size
        # Specifically we will make it so that the shape is
        # (# of cuts to apply by number of triggers), this means that to find
        # that means that ratioArray we need to repeat the triggers by # of cuts
        # and for veto* we need to repeat the threshold to be applied by how many triggers said threshold needs
        # to be applied to.

        ratioArrayTempPlus = numpy.repeat(numpy.atleast_2d(ratioArrayPlusJobNumber),
                                          len(ratioCutsPlus),
                                          axis=0)
        ratioArrayTempCross = numpy.repeat(numpy.atleast_2d(ratioArrayCrossJobNumber),
                                           len(ratioCutsCross),
                                           axis=0)
        ratioArrayTempNull = numpy.repeat(numpy.atleast_2d(ratioArrayNullJobNumber),
                                          len(ratioCutsNull),
                                          axis=0)
        detection_statistic_tmp = numpy.repeat(numpy.atleast_2d(detection_statistic),
                                              len(ratioCutsNull),
                                              axis=0)

        vetoPlusRep = numpy.kron(numpy.atleast_2d(ratioCutsPlus).T,
                                      numpy.ones((sidedCut, len(ratioArrayPlusJobNumber))))
        
        vetoCrossRep = numpy.kron(numpy.atleast_2d(ratioCutsCross).T,
                                      numpy.ones((sidedCut, len(ratioArrayCrossJobNumber),)))
        
        vetoNullRep = numpy.kron(numpy.atleast_2d(ratioCutsNull).T,
                                     numpy.ones((1, len(ratioArrayNullJobNumber),)))
        
        loudest_background_job = numpy.kron(numpy.atleast_2d(loudestBackgroundRatioJob).T,
                                            numpy.ones((1, len(detection_statistic),)))

        # Determine what clusters passed all the Ratio cuts and are louder that the loudest
        # surviving background trigger fwhen using that same ratio cut
        ratioPassCut = ((ratioArrayTempPlus  > vetoPlusRep) &
                        (ratioArrayTempCross  > vetoCrossRep) & 
                        (ratioArrayTempNull > vetoNullRep) &
                        (detection_statistic_tmp > loudest_background_job)
                       )
        
        all_ratio_boolean.append(ratioPassCut)
        all_indices.extend(triggers_from_this_trial)

    # temporarly create a pandas dataframe indicating whether a given triggers passed
    # one of the cut combinations and the loudest background trigger from that cut combination
    pass_ratio_cuts = pandas.DataFrame(numpy.hstack(all_ratio_boolean).T, index=all_indices)

    return triggers.join(pass_ratio_cuts)

def xapplyalphacuts(triggers,
                    ePlusIndex,eCrossIndex,eNullIndex,
                    iPlusIndex,iCrossIndex,iNullIndex,
                    vetoPlusRange,vetoCrossRange,vetoNullRange,
                    vetoPlusRange2,vetoCrossRange2,vetoNullRange2,
                    FAR,typeOfCutPlus,typeOfCutCross,detection_statistic_column_name):

    if float(FAR).is_integer():
        FAR = (FAR/100)

    # Calculate all ratio values
    denominator_plus = (triggers[ePlusIndex] + triggers[iPlusIndex])**0.8
    denominator_cross = (triggers[eCrossIndex] + triggers[iCrossIndex])**0.8
    
    # Calculate the so called alpha values
    plusAlphaEoverI  = 2*(triggers[ePlusIndex] - triggers[iPlusIndex])/denominator_plus
    crossAlphaEoverI = 2*(triggers[eCrossIndex] - triggers[iCrossIndex])/denominator_cross
    plusAlphaIoverE  = 2*(triggers[iPlusIndex] - triggers[ePlusIndex])/denominator_plus
    crossAlphaIoverE = 2*(triggers[iCrossIndex] - triggers[eCrossIndex])/denominator_cross
    
    # Calculate all Ratio value
    plusRatioEoverI = numpy.log(triggers[ePlusIndex] / triggers[iPlusIndex])
    crossRatioEoverI = numpy.log(triggers[eCrossIndex] / triggers[iCrossIndex])
    plusRatioIoverE = numpy.log(triggers[iPlusIndex] / triggers[ePlusIndex])
    crossRatioIoverE = numpy.log(triggers[iCrossIndex] / triggers[eCrossIndex])

    if iNullIndex == 0:
        triggers['nullenergy'] = numpy.ones(plusRatioEoverI.size)
        nullAlphaIoverE = triggers['nullenergy']
        nullRatioIoverE = triggers['nullenergy']
    else:
        nullAlphaIoverE = 2*(triggers[iNullIndex]  - triggers[eNullIndex])/(triggers[iNullIndex] + triggers[eNullIndex])**0.8
        nullRatioIoverE = numpy.log(triggers[iNullIndex] / triggers[eNullIndex])


    # Depending on wether you are requesting a two sided or one sided cut.
    # construct specific ratio and cut arrays.
    if typeOfCutPlus == 'twosided':

        alphaCutsPlus = vetoPlusRange2
        alphaCutsCross = vetoCrossRange2
        alphaCutsNull = abs(vetoNullRange2)
        alphaArrayPlus = abs(plusAlphaEoverI)+1
        alphaArrayCross = abs(crossAlphaEoverI)+1
        alphaArrayNull = nullAlphaIoverE+1
        sidedCut = 2
        
        ratioCutsPlus = numpy.log(vetoPlusRange)
        ratioCutsCross = numpy.log(vetoCrossRange)
        ratioCutsNull = numpy.log(vetoNullRange)
        ratioArrayPlus = [plusRatioEoverI,plusRatioIoverE]
        ratioArrayCross = [crossRatioEoverI,crossRatioIoverE]
        ratioArrayNull = nullRatioIoverE
        sidedCut = 2

    else:

        alphaCutsPlus = abs(vetoPlusRange2)
        alphaCutsCross = abs(vetoCrossRange2)
        alphaCutsNull = abs(vetoNullRange2)

        if typeOfCutPlus == 'EoverI':
            alphaArrayPlus = plusAlphaEoverI + 1
        else:
            alphaArrayPlus = plusAlphaIoverE + 1

        if typeOfCutCross == 'EoverI':
            alphaArrayCross = crossAlphaEoverI + 1
        else:
            alphaArrayCross = crossAlphaIoverE + 1
   
        alphaArrayNull = nullAlphaIoverE + 1
        sidedCut = 1
        
        ratioCutsPlus = numpy.log(abs(vetoPlusRange))
        ratioCutsCross = numpy.log(abs(vetoCrossRange))
        ratioCutsNull = numpy.log(abs(vetoNullRange))
        if typeOfCutPlus == 'EoverI':
            ratioArrayPlus = plusRatioEoverI
        else:
            ratioArrayPlus = plusRatioIoverE

        if typeOfCutCross == 'EoverI':
            ratioArrayCross = crossRatioEoverI
        else:
            ratioArrayCross = crossRatioIoverE

        ratioArrayNull = nullRatioIoverE
        sidedCut = 1


    # Determine how many indepdent trials there are
    background_trials = triggers.groupby(['event','internal_time_slide'])
    number_of_trials = len(background_trials)
    indexFARBackground = int(numpy.ceil(number_of_trials* FAR))
    loduest_background_that_survives_cut = []
    trial_name = []
    

    for key, item in background_trials:
        triggers_from_this_trial = background_trials.get_group(key).index
        
        # Extract detction statisc for triggers associated with this event
        detection_statistic = triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy()
        
        # Get all the ratio values for triggers associated with this event
        ratioArrayPlusJobNumber = ratioArrayPlus.loc[triggers_from_this_trial].to_numpy()
        ratioArrayCrossJobNumber = ratioArrayCross.loc[triggers_from_this_trial].to_numpy()
        ratioArrayNullJobNumber = ratioArrayNull.loc[triggers_from_this_trial].to_numpy()

        # Reshape Calculated Ratios and Cut Values such that they are the same size
        # Specifically we will make it so that the shape is
        # (# of cuts to apply by number of triggers), this means that to find
        # that means that ratioArray we need to repeat the triggers by # of cuts
        # and for veto* we need to repeat the threshold to be applied by how many triggers said threshold needs
        # to be applied to.

        ratioArrayTempPlus = numpy.repeat(numpy.atleast_2d(ratioArrayPlusJobNumber),
                                          len(ratioCutsPlus),
                                          axis=0)
        ratioArrayTempCross = numpy.repeat(numpy.atleast_2d(ratioArrayCrossJobNumber),
                                           len(ratioCutsCross),
                                           axis=0)
        ratioArrayTempNull = numpy.repeat(numpy.atleast_2d(ratioArrayNullJobNumber),
                                          len(ratioCutsNull),
                                          axis=0)

        vetoPlusRep = numpy.kron(numpy.atleast_2d(ratioCutsPlus).T,
                                      numpy.ones((len(alphaCutsPlus), len(ratioArrayPlusJobNumber))))
        vetoCrossRep = numpy.kron(numpy.atleast_2d(ratioCutsCross).T,
                                       numpy.ones((len(alphaCutsCross), len(ratioArrayCrossJobNumber),)))
        vetoNullRep = numpy.kron(numpy.atleast_2d(ratioCutsNull).T,
                                      numpy.ones((len(alphaCutsNull), len(ratioArrayNullJobNumber),)))

        
        # Get all the relevant alpha ratios for triggers associated with this event
        alphaArrayPlusJobNumber = alphaArrayPlus.loc[triggers_from_this_trial].to_numpy()
        alphaArrayCrossJobNumber = alphaArrayCross.loc[triggers_from_this_trial].to_numpy()
        alphaArrayNullJobNumber = alphaArrayNull.loc[triggers_from_this_trial].to_numpy()

        # Reshape Calculated Ratios and Cut Values such that they are the same size
        # Specifically we will make it so that the shape is
        # (# of cuts to apply by number of triggers), this means that to find
        # that means that ratioArray we need to repeat the triggers by # of cuts
        # and for veto* we need to repeat the threshold to be applied by how many triggers said threshold needs
        # to be applied to.

        alphaArrayTempPlus = numpy.repeat(numpy.atleast_2d(alphaArrayPlusJobNumber),
                                          len(alphaCutsPlus),
                                          axis=0)
        alphaArrayTempCross = numpy.repeat(numpy.atleast_2d(alphaArrayCrossJobNumber),
                                           len(alphaCutsCross),
                                           axis=0)
        alphaArrayTempNull = numpy.repeat(numpy.atleast_2d(alphaArrayNullJobNumber),
                                          len(alphaCutsNull),
                                          axis=0)

        vetoPlusRange2Rep = numpy.kron(numpy.atleast_2d(alphaCutsPlus).T,
                                      numpy.ones((1, len(alphaArrayPlusJobNumber))))
        vetoCrossRange2Rep = numpy.kron(numpy.atleast_2d(alphaCutsCross).T,
                                      numpy.ones((1, len(alphaArrayCrossJobNumber),)))
        vetoNullRange2Rep = numpy.kron(numpy.atleast_2d(alphaCutsNull).T,
                                     numpy.ones((1, len(alphaArrayNullJobNumber),)))

        
        # Take absolute value of vetoRange2Rep.
        vetoPlusRange2Rep = abs(vetoPlusRange2Rep)
        vetoCrossRange2Rep = abs(vetoCrossRange2Rep)
        vetoNullRange2Rep = abs(vetoNullRange2Rep)
        
        vetoPlusRange2Rep[vetoPlusRange2Rep==0] = -numpy.inf
        vetoCrossRange2Rep[vetoCrossRange2Rep==0] = -numpy.inf
        vetoNullRange2Rep[vetoNullRange2Rep==0] = -numpy.inf

        # Determine what clusters passed all the Ratio cuts
        coherent_cuts = (
                        (alphaArrayTempPlus  >= vetoPlusRange2Rep) &
                        (alphaArrayTempCross >= vetoCrossRange2Rep) &
                        (alphaArrayTempNull  >= vetoNullRange2Rep) &
                        (ratioArrayTempPlus  > vetoPlusRep) &
                        (ratioArrayTempCross  > vetoCrossRep) & 
                        (ratioArrayTempNull > vetoNullRep)
                       )

        # Find loudest for each column
        # Find Surviving Offsource
        backGroundArray = coherent_cuts*numpy.repeat(numpy.atleast_2d(detection_statistic),
                                                          len(alphaCutsPlus),
                                                          axis=0)
        # Find loudest surviving trigger for this trial for the grid of applied cuts
        # to do this we find the max over the columns
        #(remember each row represented an applied a cut)
        loduest_background_that_survives_cut.append(numpy.max(backGroundArray,1))
        trial_name.append(key)

    loduest_background_that_survives_cut = numpy.vstack(loduest_background_that_survives_cut)
    loduest_background_that_survives_cut.sort(axis=0)

    return loduest_background_that_survives_cut[indexFARBackground,:]

def xapplyalphacutsinjection(triggers,
                    ePlusIndex,eCrossIndex,eNullIndex,
                    iPlusIndex,iCrossIndex,iNullIndex,
                    vetoPlusRange,vetoCrossRange,vetoNullRange,
                    vetoPlusRange2,vetoCrossRange2,vetoNullRange2,
                    loudestBackgroundAlpha,typeOfCutPlus,typeOfCutCross,detection_statistic_column_name):

    # Calculate all ratio values
    denominator_plus = (triggers[ePlusIndex] + triggers[iPlusIndex])**0.8
    denominator_cross = (triggers[eCrossIndex] + triggers[iCrossIndex])**0.8
    
    # Calculate the so called alpha values
    plusAlphaEoverI  = 2*(triggers[ePlusIndex] - triggers[iPlusIndex])/denominator_plus
    crossAlphaEoverI = 2*(triggers[eCrossIndex] - triggers[iCrossIndex])/denominator_cross
    plusAlphaIoverE  = 2*(triggers[iPlusIndex] - triggers[ePlusIndex])/denominator_plus
    crossAlphaIoverE = 2*(triggers[iCrossIndex] - triggers[eCrossIndex])/denominator_cross
    
    # Calculate all Ratio value
    plusRatioEoverI = numpy.log(triggers[ePlusIndex] / triggers[iPlusIndex])
    crossRatioEoverI = numpy.log(triggers[eCrossIndex] / triggers[iCrossIndex])
    plusRatioIoverE = numpy.log(triggers[iPlusIndex] / triggers[ePlusIndex])
    crossRatioIoverE = numpy.log(triggers[iCrossIndex] / triggers[eCrossIndex])

    if iNullIndex == 0:
        triggers['nullenergy'] = numpy.ones(plusRatioEoverI.size)
        nullAlphaIoverE = triggers['nullenergy']
        nullRatioIoverE = triggers['nullenergy']
    else:
        nullAlphaIoverE = 2*(triggers[iNullIndex]  - triggers[eNullIndex])/(triggers[iNullIndex] + triggers[eNullIndex])**0.8
        nullRatioIoverE = numpy.log(triggers[iNullIndex] / triggers[eNullIndex])


    # Depending on wether you are requesting a two sided or one sided cut.
    # construct specific ratio and cut arrays.
    if typeOfCutPlus == 'twosided':

        alphaCutsPlus = vetoPlusRange2
        alphaCutsCross = vetoCrossRange2
        alphaCutsNull = abs(vetoNullRange2)
        alphaArrayPlus = abs(plusAlphaEoverI)+1
        alphaArrayCross = abs(crossAlphaEoverI)+1
        alphaArrayNull = nullAlphaIoverE+1
        sidedCut = 2
        
        ratioCutsPlus = numpy.log(vetoPlusRange)
        ratioCutsCross = numpy.log(vetoCrossRange)
        ratioCutsNull = numpy.log(vetoNullRange)
        ratioArrayPlus = [plusRatioEoverI,plusRatioIoverE]
        ratioArrayCross = [crossRatioEoverI,crossRatioIoverE]
        ratioArrayNull = nullRatioIoverE
        sidedCut = 2

    else:

        alphaCutsPlus = abs(vetoPlusRange2)
        alphaCutsCross = abs(vetoCrossRange2)
        alphaCutsNull = abs(vetoNullRange2)

        if typeOfCutPlus == 'EoverI':
            alphaArrayPlus = plusAlphaEoverI + 1
        else:
            alphaArrayPlus = plusAlphaIoverE + 1

        if typeOfCutCross == 'EoverI':
            alphaArrayCross = crossAlphaEoverI + 1
        else:
            alphaArrayCross = crossAlphaIoverE + 1
   
        alphaArrayNull = nullAlphaIoverE + 1
        sidedCut = 1
        
        ratioCutsPlus = numpy.log(abs(vetoPlusRange))
        ratioCutsCross = numpy.log(abs(vetoCrossRange))
        ratioCutsNull = numpy.log(abs(vetoNullRange))
        if typeOfCutPlus == 'EoverI':
            ratioArrayPlus = plusRatioEoverI
        else:
            ratioArrayPlus = plusRatioIoverE

        if typeOfCutCross == 'EoverI':
            ratioArrayCross = crossRatioEoverI
        else:
            ratioArrayCross = crossRatioIoverE

        ratioArrayNull = nullRatioIoverE
        sidedCut = 1


    # Determine how many indepdent trials there are
    injection_trials = triggers.groupby(['injection_scale','injection_number'])
    all_alpha_boolean = []
    all_indices = []

    for key, item in injection_trials:
        triggers_from_this_trial = injection_trials.get_group(key).index
        
        # Extract detction statisc for triggers associated with this event
        detection_statistic = triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy()
        
        # Get all the ratio values for triggers associated with this event
        ratioArrayPlusJobNumber = ratioArrayPlus.loc[triggers_from_this_trial].to_numpy()
        ratioArrayCrossJobNumber = ratioArrayCross.loc[triggers_from_this_trial].to_numpy()
        ratioArrayNullJobNumber = ratioArrayNull.loc[triggers_from_this_trial].to_numpy()

        # Reshape Calculated Ratios and Cut Values such that they are the same size
        # Specifically we will make it so that the shape is
        # (# of cuts to apply by number of triggers), this means that to find
        # that means that ratioArray we need to repeat the triggers by # of cuts
        # and for veto* we need to repeat the threshold to be applied by how many triggers said threshold needs
        # to be applied to.

        ratioArrayTempPlus = numpy.repeat(numpy.atleast_2d(ratioArrayPlusJobNumber),
                                          len(ratioCutsPlus),
                                          axis=0)
        ratioArrayTempCross = numpy.repeat(numpy.atleast_2d(ratioArrayCrossJobNumber),
                                           len(ratioCutsCross),
                                           axis=0)
        ratioArrayTempNull = numpy.repeat(numpy.atleast_2d(ratioArrayNullJobNumber),
                                          len(ratioCutsNull),
                                          axis=0)

        vetoPlusRep = numpy.kron(numpy.atleast_2d(ratioCutsPlus).T,
                                      numpy.ones((len(alphaCutsPlus), len(ratioArrayPlusJobNumber))))
        vetoCrossRep = numpy.kron(numpy.atleast_2d(ratioCutsCross).T,
                                       numpy.ones((len(alphaCutsCross), len(ratioArrayCrossJobNumber),)))
        vetoNullRep = numpy.kron(numpy.atleast_2d(ratioCutsNull).T,
                                      numpy.ones((len(alphaCutsNull), len(ratioArrayNullJobNumber),)))

        
        # Get all the relevant alpha ratios for triggers associated with this event
        alphaArrayPlusJobNumber = alphaArrayPlus.loc[triggers_from_this_trial].to_numpy()
        alphaArrayCrossJobNumber = alphaArrayCross.loc[triggers_from_this_trial].to_numpy()
        alphaArrayNullJobNumber = alphaArrayNull.loc[triggers_from_this_trial].to_numpy()

        # Reshape Calculated Ratios and Cut Values such that they are the same size
        # Specifically we will make it so that the shape is
        # (# of cuts to apply by number of triggers), this means that to find
        # that means that ratioArray we need to repeat the triggers by # of cuts
        # and for veto* we need to repeat the threshold to be applied by how many triggers said threshold needs
        # to be applied to.

        alphaArrayTempPlus = numpy.repeat(numpy.atleast_2d(alphaArrayPlusJobNumber),
                                          len(alphaCutsPlus),
                                          axis=0)
        alphaArrayTempCross = numpy.repeat(numpy.atleast_2d(alphaArrayCrossJobNumber),
                                           len(alphaCutsCross),
                                           axis=0)
        alphaArrayTempNull = numpy.repeat(numpy.atleast_2d(alphaArrayNullJobNumber),
                                          len(alphaCutsNull),
                                          axis=0)

        vetoPlusRange2Rep = numpy.kron(numpy.atleast_2d(alphaCutsPlus).T,
                                      numpy.ones((1, len(alphaArrayPlusJobNumber))))
        vetoCrossRange2Rep = numpy.kron(numpy.atleast_2d(alphaCutsCross).T,
                                      numpy.ones((1, len(alphaArrayCrossJobNumber),)))
        vetoNullRange2Rep = numpy.kron(numpy.atleast_2d(alphaCutsNull).T,
                                     numpy.ones((1, len(alphaArrayNullJobNumber),)))

        
        # Take absolute value of vetoRange2Rep.
        vetoPlusRange2Rep = abs(vetoPlusRange2Rep)
        vetoCrossRange2Rep = abs(vetoCrossRange2Rep)
        vetoNullRange2Rep = abs(vetoNullRange2Rep)
        
        vetoPlusRange2Rep[vetoPlusRange2Rep==0] = -numpy.inf
        vetoCrossRange2Rep[vetoCrossRange2Rep==0] = -numpy.inf
        vetoNullRange2Rep[vetoNullRange2Rep==0] = -numpy.inf
        
        # Now repeat the background
        detection_statistic_tmp = numpy.repeat(numpy.atleast_2d(detection_statistic),
                                              len(alphaCutsNull),
                                              axis=0)
        
        loudest_background_job = numpy.kron(numpy.atleast_2d(loudestBackgroundAlpha).T,
                                            numpy.ones((1, len(detection_statistic),)))

        # Determine what clusters passed all the Ratio cuts
        coherent_cuts = (
                        (alphaArrayTempPlus >= vetoPlusRange2Rep) &
                        (alphaArrayTempCross >= vetoCrossRange2Rep) &
                        (alphaArrayTempNull >= vetoNullRange2Rep) &
                        (ratioArrayTempPlus > vetoPlusRep) &
                        (ratioArrayTempCross > vetoCrossRep) & 
                        (ratioArrayTempNull > vetoNullRep) &
                        (detection_statistic_tmp > loudest_background_job)
                       )

        all_alpha_boolean.append(coherent_cuts)
        all_indices.extend(triggers_from_this_trial)

    # temporarly create a pandas dataframe indicating whether a given triggers passed
    # one of the cut combinations and the loudest background trigger from that cut combination
    pass_alpha_cuts = pandas.DataFrame(numpy.hstack(all_alpha_boolean).T, index=all_indices)

    return triggers.join(pass_alpha_cuts)
