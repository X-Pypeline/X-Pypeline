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

def prep_triggers_for_ratio(triggers, ePlusIndex,eCrossIndex,eNullIndex,
                            iPlusIndex,iCrossIndex,iNullIndex,
                            vetoPlusRange,vetoCrossRange,vetoNullRange,
                            typeOfCutPlus,typeOfCutCross,):
    # Depending on wether you are requesting a two sided or one sided cut.
    # construct specific ratio and cut arrays.

    # Calculate all Ratio Values
    plusRatioEoverI = numpy.log10(triggers[ePlusIndex] / triggers[iPlusIndex])
    crossRatioEoverI = numpy.log10(triggers[eCrossIndex] / triggers[iCrossIndex])
    plusRatioIoverE = numpy.log10(triggers[iPlusIndex] / triggers[ePlusIndex])
    crossRatioIoverE = numpy.log10(triggers[iCrossIndex] / triggers[eCrossIndex])

    if iNullIndex == 0:
        triggers['nullenergy'] = numpy.ones(plusRatioEoverI.size)
        nullRatioIoverE = triggers['nullenergy']
    else:
        nullRatioIoverE = numpy.log10(triggers[iNullIndex] / triggers[eNullIndex])

    ratioArrayNull = nullRatioIoverE

    if typeOfCutPlus == 'twosided':

        ratioCutsPlus = numpy.log10(vetoPlusRange)
        ratioCutsCross = numpy.log10(vetoCrossRange)
        ratioCutsNull = numpy.log10(vetoNullRange)

        ratioCutsPlus_repeated = numpy.atleast_2d(numpy.repeat(ratioCutsPlus,2)).T
        ratioCutsCross_repeated = numpy.atleast_2d(numpy.repeat(ratioCutsCross,2)).T
        ratioCutsNull_repeated = numpy.atleast_2d(ratioCutsNull).T

        ratioArrayPlus = pandas.concat((plusRatioEoverI, plusRatioIoverE), join='outer', axis=1)
        ratioArrayCross = pandas.concat((crossRatioEoverI, crossRatioIoverE), join='outer', axis=1)

    else:

        ratioCutsPlus = numpy.log10(abs(vetoPlusRange))
        ratioCutsCross = numpy.log10(abs(vetoCrossRange))
        ratioCutsNull = numpy.log10(abs(vetoNullRange))

        ratioCutsPlus_repeated = numpy.atleast_2d(ratioCutsPlus).T
        ratioCutsCross_repeated = numpy.atleast_2d(ratioCutsCross).T
        ratioCutsNull_repeated = numpy.atleast_2d(ratioCutsNull).T

        if typeOfCutPlus == 'EoverI':
            ratioArrayPlus = plusRatioEoverI
        else:
            ratioArrayPlus = plusRatioIoverE

        if typeOfCutCross == 'EoverI':
            ratioArrayCross = crossRatioEoverI
        else:
            ratioArrayCross = crossRatioIoverE

    return (ratioArrayPlus, ratioArrayCross, ratioArrayNull, ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated,
           ratioCutsPlus, ratioCutsCross, ratioCutsNull)

def apply_ratio_cuts_this_trial(triggers,triggers_from_this_trial,typeOfCutPlus,
                                ratioArrayPlus,ratioArrayCross,ratioArrayNull,
                                ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated,
                                ratioCutsPlus, ratioCutsCross, ratioCutsNull,
                                detection_statistic_column_name,
                                applying_to_injections,loudestBackgroundRatioJob=None):

    detection_statistic = numpy.atleast_2d(triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy())

    ratioArrayPlusJobNumber = numpy.atleast_2d(ratioArrayPlus.loc[triggers_from_this_trial].to_numpy().T)
    ratioArrayCrossJobNumber = numpy.atleast_2d(ratioArrayCross.loc[triggers_from_this_trial].to_numpy().T)
    ratioArrayNullJobNumber = numpy.atleast_2d(ratioArrayNull.loc[triggers_from_this_trial].to_numpy())

    # Reshape Calculated Ratios and Cut Values such that they are the same size
    # Specifically we will make it so that the shape is
    # (# of cuts to apply by number of triggers), this means that to find
    # that means that ratioArray we need to repeat the triggers by # of cuts
    # and for veto* we need to repeat the threshold to be applied by how many triggers said threshold needs
    # to be applied to.

    ratioArrayTempPlus = numpy.repeat(ratioArrayPlusJobNumber,
                                      len(ratioCutsPlus),
                                      axis=0)
    ratioArrayTempCross = numpy.repeat(ratioArrayCrossJobNumber,
                                       len(ratioCutsCross),
                                       axis=0)
    ratioArrayTempNull = numpy.repeat(ratioArrayNullJobNumber,
                                      len(ratioCutsNull),
                                      axis=0)

    vetoPlusRep = numpy.kron(ratioCutsPlus_repeated, numpy.ones((1, ratioArrayPlusJobNumber.shape[1])))
    vetoCrossRep = numpy.kron(ratioCutsCross_repeated, numpy.ones((1, ratioArrayCrossJobNumber.shape[1],)))
    vetoNullRep = numpy.kron(ratioCutsNull_repeated, numpy.ones((1, ratioArrayNullJobNumber.shape[1],)))

    if applying_to_injections:
        detection_statistic_tmp = numpy.repeat(detection_statistic,
                                              len(ratioCutsNull),
                                              axis=0)
        loudest_background_job = numpy.kron(numpy.atleast_2d(loudestBackgroundRatioJob).T,
                                            numpy.ones((1, detection_statistic.shape[1],)))

    # Determine what clusters passed all the Ratio cuts
    if applying_to_injections:
        if typeOfCutPlus == 'twosided':
            ratioPassCut = (((ratioArrayTempPlus  >= vetoPlusRep)[0:len(ratioCutsPlus)] | (ratioArrayTempPlus  >= vetoPlusRep)[len(ratioCutsPlus)::]) &
            ((ratioArrayTempCross  >= vetoCrossRep)[0:len(ratioCutsCross)] | (ratioArrayTempCross  >= vetoCrossRep)[len(ratioCutsCross)::]) &
            (ratioArrayTempNull >= vetoNullRep) &
             (detection_statistic_tmp > loudest_background_job))
        else:
            ratioPassCut = ((ratioArrayTempPlus  >= vetoPlusRep) &
                            (ratioArrayTempCross  >= vetoCrossRep) &
                            (ratioArrayTempNull >= vetoNullRep) &
                            (detection_statistic_tmp > loudest_background_job))
    else:
        if typeOfCutPlus == 'twosided':
            ratioPassCut = (((ratioArrayTempPlus  >= vetoPlusRep)[0:len(ratioCutsPlus)] | (ratioArrayTempPlus  >= vetoPlusRep)[len(ratioCutsPlus)::]) &
            ((ratioArrayTempCross  >= vetoCrossRep)[0:len(ratioCutsCross)] | (ratioArrayTempCross  >= vetoCrossRep)[len(ratioCutsCross)::]) &
            (ratioArrayTempNull >= vetoNullRep))
        else:
            ratioPassCut = ((ratioArrayTempPlus  >= vetoPlusRep) &
                            (ratioArrayTempCross  >= vetoCrossRep) &
                            (ratioArrayTempNull >= vetoNullRep))
    return ratioPassCut

def xapplyratiocuts(triggers, ePlusIndex,eCrossIndex,eNullIndex,
                    iPlusIndex,iCrossIndex,iNullIndex,
                    vetoPlusRange,vetoCrossRange,vetoNullRange,
                    FAR_Tuning,typeOfCutPlus,typeOfCutCross,
                    detection_statistic_column_name):

    if float(FAR_Tuning).is_integer():
        FAR_Tuning = (FAR_Tuning/100)

    # Depending on what was asked prep triggers to apply ratio array
    ratioArrayPlus, ratioArrayCross, ratioArrayNull, ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated, ratioCutsPlus, ratioCutsCross, ratioCutsNull = \
        prep_triggers_for_ratio(triggers, ePlusIndex,eCrossIndex,eNullIndex,
                                iPlusIndex,iCrossIndex,iNullIndex,
                                vetoPlusRange,vetoCrossRange,vetoNullRange,
                                typeOfCutPlus,typeOfCutCross,)

    # Determine how many indepdent trials there are
    background_trials = triggers.groupby(['event','internal_time_slide'])
    number_of_trials = len(background_trials)
    indexFARBackground = int(numpy.ceil(number_of_trials* FAR_Tuning))
    loudestBackgroundRatioJob = []
    trial_name = []

    for key, item in background_trials:
        triggers_from_this_trial = background_trials.get_group(key).index

        detection_statistic = numpy.atleast_2d(triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy())

        ratioPassCut = apply_ratio_cuts_this_trial(triggers,triggers_from_this_trial,typeOfCutPlus,
                                ratioArrayPlus,ratioArrayCross,ratioArrayNull,
                                ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated,
                                ratioCutsPlus, ratioCutsCross, ratioCutsNull,
                                detection_statistic_column_name,
                                False,loudestBackgroundRatioJob=None)

        # Find Surviving Offsource
        backGroundArray = ratioPassCut*numpy.repeat(detection_statistic,
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

    # Depending on what was asked prep triggers to apply ratio array
    ratioArrayPlus, ratioArrayCross, ratioArrayNull, ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated, ratioCutsPlus, ratioCutsCross, ratioCutsNull = \
        prep_triggers_for_ratio(triggers, ePlusIndex,eCrossIndex,eNullIndex,
                                iPlusIndex,iCrossIndex,iNullIndex,
                                vetoPlusRange,vetoCrossRange,vetoNullRange,
                                typeOfCutPlus,typeOfCutCross,)

    injection_trials = triggers.groupby(['injection_scale','injection_number'])
    all_ratio_boolean = []
    all_indices = []

    for key, item in injection_trials:
        triggers_from_this_trial = injection_trials.get_group(key).index
        ratioPassCut = apply_ratio_cuts_this_trial(triggers,triggers_from_this_trial,typeOfCutPlus,
                                ratioArrayPlus,ratioArrayCross,ratioArrayNull,
                                ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated,
                                ratioCutsPlus, ratioCutsCross, ratioCutsNull,
                                detection_statistic_column_name,
                                True,loudestBackgroundRatioJob=loudestBackgroundRatioJob)

        all_ratio_boolean.append(ratioPassCut)
        all_indices.extend(triggers_from_this_trial)

    # temporarly create a pandas dataframe indicating whether a given triggers passed
    # one of the cut combinations and the loudest background trigger from that cut combination
    pass_ratio_cuts = pandas.DataFrame(numpy.hstack(all_ratio_boolean).T, index=all_indices)

    return triggers.join(pass_ratio_cuts)

def prep_triggers_for_alpha_cuts(triggers,
                                ePlusIndex,eCrossIndex,eNullIndex,
                                iPlusIndex,iCrossIndex,iNullIndex,
                                vetoPlusRange2,vetoCrossRange2,vetoNullRange2,
                                typeOfCutPlus,typeOfCutCross,):
    # Calculate all ratio values
    denominator_plus = (triggers[ePlusIndex] + triggers[iPlusIndex])**0.8
    denominator_cross = (triggers[eCrossIndex] + triggers[iCrossIndex])**0.8

    # Calculate the so called alpha values
    plusAlphaEoverI  = 2*(triggers[ePlusIndex] - triggers[iPlusIndex])/denominator_plus
    crossAlphaEoverI = 2*(triggers[eCrossIndex] - triggers[iCrossIndex])/denominator_cross
    plusAlphaIoverE  = 2*(triggers[iPlusIndex] - triggers[ePlusIndex])/denominator_plus
    crossAlphaIoverE = 2*(triggers[iCrossIndex] - triggers[eCrossIndex])/denominator_cross

    if iNullIndex == 0:
        triggers['nullenergy'] = numpy.ones(plusAlphaEoverI.size)
        nullAlphaIoverE = triggers['nullenergy']
    else:
        nullAlphaIoverE = 2*(triggers[iNullIndex]  - triggers[eNullIndex])/(triggers[iNullIndex] + triggers[eNullIndex])**0.8

    if typeOfCutPlus == 'twosided':
        alphaCutsPlus = vetoPlusRange2
        alphaCutsCross = vetoCrossRange2
        alphaCutsNull = abs(vetoNullRange2)
        alphaArrayPlus = abs(plusAlphaEoverI)+1
        alphaArrayCross = abs(crossAlphaEoverI)+1
        alphaArrayNull = nullAlphaIoverE+1

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
    return alphaArrayPlus, alphaArrayCross, alphaArrayNull, alphaCutsPlus, alphaCutsCross, alphaCutsNull

def apply_alpha_cuts_to_trial(triggers,triggers_from_this_trial,
                            alphaArrayPlus,alphaArrayCross,alphaArrayNull,
                            alphaCutsPlus, alphaCutsCross, alphaCutsNull,
                            detection_statistic_column_name,
                            applying_to_injections,loudestBackgroundAlpha=None):

    detection_statistic = numpy.atleast_2d(triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy())
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

    if applying_to_injections:
        # Now repeat the background
        detection_statistic_tmp = numpy.repeat(numpy.atleast_2d(detection_statistic),
                                              len(alphaCutsNull),
                                              axis=0)

        loudest_background_job = numpy.kron(numpy.atleast_2d(loudestBackgroundAlpha).T,
                                            numpy.ones((1,  detection_statistic.shape[1],)))


    # Take absolute value of vetoRange2Rep.
    vetoPlusRange2Rep = abs(vetoPlusRange2Rep)
    vetoCrossRange2Rep = abs(vetoCrossRange2Rep)
    vetoNullRange2Rep = abs(vetoNullRange2Rep)

    vetoPlusRange2Rep[vetoPlusRange2Rep==0] = -numpy.inf
    vetoCrossRange2Rep[vetoCrossRange2Rep==0] = -numpy.inf
    vetoNullRange2Rep[vetoNullRange2Rep==0] = -numpy.inf

    if applying_to_injections:
        coherent_cuts = (
                        (alphaArrayTempPlus >= vetoPlusRange2Rep) &
                        (alphaArrayTempCross >= vetoCrossRange2Rep) &
                        (alphaArrayTempNull >= vetoNullRange2Rep) &
                        (detection_statistic_tmp > loudest_background_job)
                        )
    else:
        coherent_cuts = (
                        (alphaArrayTempPlus >= vetoPlusRange2Rep) &
                        (alphaArrayTempCross >= vetoCrossRange2Rep) &
                        (alphaArrayTempNull >= vetoNullRange2Rep)
                        )
    return coherent_cuts

def xapplyalphacuts(triggers,
                    ePlusIndex,eCrossIndex,eNullIndex,
                    iPlusIndex,iCrossIndex,iNullIndex,
                    vetoPlusRange,vetoCrossRange,vetoNullRange,
                    vetoPlusRange2,vetoCrossRange2,vetoNullRange2,
                    FAR,typeOfCutPlus,typeOfCutCross,detection_statistic_column_name):

    if float(FAR).is_integer():
        FAR = (FAR/100)

    alphaArrayPlus, alphaArrayCross, alphaArrayNull, alphaCutsPlus, alphaCutsCross, alphaCutsNull = \
        prep_triggers_for_alpha_cuts(triggers,
                                ePlusIndex,eCrossIndex,eNullIndex,
                                iPlusIndex,iCrossIndex,iNullIndex,
                                vetoPlusRange2,vetoCrossRange2,vetoNullRange2,
                                typeOfCutPlus,typeOfCutCross,)

    # Depending on what was asked prep triggers to apply ratio array
    ratioArrayPlus, ratioArrayCross, ratioArrayNull, ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated, ratioCutsPlus, ratioCutsCross, ratioCutsNull = \
        prep_triggers_for_ratio(triggers, ePlusIndex,eCrossIndex,eNullIndex,
                                iPlusIndex,iCrossIndex,iNullIndex,
                                vetoPlusRange,vetoCrossRange,vetoNullRange,
                                typeOfCutPlus,typeOfCutCross,)

    # Depending on wether you are requesting a two sided or one sided cut.
    # construct specific ratio and cut arrays.

    # Determine how many indepdent trials there are
    background_trials = triggers.groupby(['event','internal_time_slide'])
    number_of_trials = len(background_trials)
    indexFARBackground = int(numpy.ceil(number_of_trials* FAR))
    loudest_background_that_survives_cut = []
    trial_name = []


    for key, item in background_trials:
        triggers_from_this_trial = background_trials.get_group(key).index

        # Extract detction statisc for triggers associated with this event
        detection_statistic = numpy.atleast_2d(triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy())

        ratioPassCut = apply_ratio_cuts_this_trial(triggers,triggers_from_this_trial,typeOfCutPlus,
                                                ratioArrayPlus,ratioArrayCross,ratioArrayNull,
                                                ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated,
                                                ratioCutsPlus, ratioCutsCross, ratioCutsNull,
                                                detection_statistic_column_name,
                                                False,loudestBackgroundRatioJob=None)

        alpha_cuts = apply_alpha_cuts_to_trial(triggers,triggers_from_this_trial,
                            alphaArrayPlus,alphaArrayCross,alphaArrayNull,
                            alphaCutsPlus, alphaCutsCross, alphaCutsNull,
                            detection_statistic_column_name,
                            False,loudestBackgroundAlpha=None)

        coherent_cuts = numpy.repeat(ratioPassCut,len(alphaCutsNull),0) & alpha_cuts


        # Find loudest for each column
        # Find Surviving Offsource
        backGroundArray = coherent_cuts*numpy.repeat(numpy.atleast_2d(detection_statistic),
                                                          len(alphaCutsPlus),
                                                          axis=0)
        # Find loudest surviving trigger for this trial for the grid of applied cuts
        # to do this we find the max over the columns
        #(remember each row represented an applied a cut)
        loudest_background_that_survives_cut.append(numpy.max(backGroundArray,1))
        trial_name.append(key)

    loudest_background_that_survives_cut = numpy.vstack(loudest_background_that_survives_cut)
    sorted_indices = loudest_background_that_survives_cut.argsort(axis=0)[:,0]
    loudest_background_that_survives_cut.sort(axis=0)

    return loudest_background_that_survives_cut[indexFARBackground,:], numpy.array(trial_name)[sorted_indices], loudest_background_that_survives_cut

def xapplyalphacutsinjection(triggers,
                    ePlusIndex,eCrossIndex,eNullIndex,
                    iPlusIndex,iCrossIndex,iNullIndex,
                    vetoPlusRange,vetoCrossRange,vetoNullRange,
                    vetoPlusRange2,vetoCrossRange2,vetoNullRange2,
                    loudestBackgroundAlpha,typeOfCutPlus,typeOfCutCross,detection_statistic_column_name):

    alphaArrayPlus, alphaArrayCross, alphaArrayNull, alphaCutsPlus, alphaCutsCross, alphaCutsNull = \
        prep_triggers_for_alpha_cuts(triggers,
                                ePlusIndex,eCrossIndex,eNullIndex,
                                iPlusIndex,iCrossIndex,iNullIndex,
                                vetoPlusRange2,vetoCrossRange2,vetoNullRange2,
                                typeOfCutPlus,typeOfCutCross,)

    # Depending on what was asked prep triggers to apply ratio array
    ratioArrayPlus, ratioArrayCross, ratioArrayNull, ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated, ratioCutsPlus, ratioCutsCross, ratioCutsNull = \
        prep_triggers_for_ratio(triggers, ePlusIndex,eCrossIndex,eNullIndex,
                                iPlusIndex,iCrossIndex,iNullIndex,
                                vetoPlusRange,vetoCrossRange,vetoNullRange,
                                typeOfCutPlus,typeOfCutCross,)

    # Determine how many indepdent trials there are
    injection_trials = triggers.groupby(['injection_scale','injection_number'])
    all_alpha_boolean = []
    all_indices = []

    for key, item in injection_trials:
        triggers_from_this_trial = injection_trials.get_group(key).index

        # Extract detction statisc for triggers associated with this event
        detection_statistic = triggers[detection_statistic_column_name].loc[triggers_from_this_trial].to_numpy()

        ratioPassCut = apply_ratio_cuts_this_trial(triggers,triggers_from_this_trial,typeOfCutPlus,
                                                ratioArrayPlus,ratioArrayCross,ratioArrayNull,
                                                ratioCutsPlus_repeated, ratioCutsCross_repeated, ratioCutsNull_repeated,
                                                ratioCutsPlus, ratioCutsCross, ratioCutsNull,
                                                detection_statistic_column_name,
                                                True,loudestBackgroundRatioJob=[0])

        alpha_cuts = apply_alpha_cuts_to_trial(triggers,triggers_from_this_trial,
                            alphaArrayPlus,alphaArrayCross,alphaArrayNull,
                            alphaCutsPlus, alphaCutsCross, alphaCutsNull,
                            detection_statistic_column_name,
                            True,loudestBackgroundAlpha=loudestBackgroundAlpha)

        coherent_cuts = numpy.repeat(ratioPassCut,len(alphaCutsNull),0) & alpha_cuts

        all_alpha_boolean.append(coherent_cuts)
        all_indices.extend(triggers_from_this_trial)

    # temporarly create a pandas dataframe indicating whether a given triggers passed
    # one of the cut combinations and the loudest background trigger from that cut combination
    pass_alpha_cuts = pandas.DataFrame(numpy.hstack(all_alpha_boolean).T, index=all_indices)

    return triggers.join(pass_alpha_cuts)
