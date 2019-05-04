import pandas
import os
import random
import numpy

likelihood_columns = ['standard_energy','coherent_f_plus','incoherent_f_plus', 'coherent_f_cross', 'incoherent_f_cross','circenergy', 'circnullenergy', 'circinc', 'circnullinc',]

def window_cut(injection_triggers_df, injection_info_dir='input', duration=0.1):

    # Perform Window Cut
    groups = injection_triggers_df.groupby(['waveform', 'injection_number']).groups
    survive_window_cut_indices = []
    for key, indices in groups.items():
        with open(os.path.join(injection_info_dir, 'injection_' + key[0].split('_')[-1] + '.txt'), 'r') as f:
            tmp = f.readlines()
        injection_info = tmp[key[1]].split(' ')
        injection_time = float(injection_info[0]) + 10e-10*float(injection_info[1])
        window = [injection_time - duration, injection_time + duration]
        triggers_this_trial = injection_triggers_df.loc[indices]
        triggers_this_trial = triggers_this_trial.loc[((triggers_this_trial.min_time_of_cluster > window[0]) &
                                                       (triggers_this_trial.min_time_of_cluster < window[1])) |
                                                     ((triggers_this_trial.max_time_of_cluster > window[0]) &
                                                       (triggers_this_trial.max_time_of_cluster < window[1])) |
                                                     ((triggers_this_trial.min_time_of_cluster < window[0]) &
                                                       (triggers_this_trial.max_time_of_cluster > window[1])) |
                                                      ((triggers_this_trial.min_time_of_cluster > window[0]) &
                                                       (triggers_this_trial.max_time_of_cluster < window[1]))
                                                     ]
        survive_window_cut_indices.extend(triggers_this_trial.index.tolist())

    injection_triggers_df = injection_triggers_df.loc[survive_window_cut_indices]

    return injection_triggers_df

def prep_data(injection_triggers_df, background_triggers_df, injection_info_dir='input', duration=0.1):

    # Perform the window cut
    injection_triggers_df = window_cut(injection_triggers_df, injection_info_dir=injection_info_dir, duration=duration)

    # calculate some likelihoods
    injection_triggers_df['standard_energy'] = injection_triggers_df['coherent_f_plus'] + injection_triggers_df['coherent_f_cross']
    injection_triggers_df['circenergy'] = injection_triggers_df[['coherent_f_right', 'coherent_f_left']].max(1)
    injection_triggers_df['circnullenergy'] = injection_triggers_df[['coherent_f_right_null', 'coherent_f_left_null']].min(1)
    injection_triggers_df['circinc'] = injection_triggers_df[['incoherent_f_left','incoherent_f_right']].max(1)
    injection_triggers_df['circnullinc'] = injection_triggers_df[['incoherent_f_left_null','incoherent_f_right_null']].max(1)

    background_triggers_df['standard_energy'] = background_triggers_df['coherent_f_plus'] + background_triggers_df['coherent_f_cross']
    background_triggers_df['circenergy'] = background_triggers_df[['coherent_f_right', 'coherent_f_left']].max(1)
    background_triggers_df['circnullenergy'] = background_triggers_df[['coherent_f_right_null', 'coherent_f_left_null']].min(1)
    background_triggers_df['circinc'] = background_triggers_df[['incoherent_f_left','incoherent_f_right']].max(1)
    background_triggers_df['circnullinc'] = background_triggers_df[['incoherent_f_left_null','incoherent_f_right_null']].max(1)

    # Log all likelihoods
    injection_triggers_df[likelihood_columns] = numpy.log(injection_triggers_df[likelihood_columns])
    background_triggers_df[likelihood_columns] = numpy.log(background_triggers_df[likelihood_columns])


    # Calculate all ratio values
    injection_triggers_df['circ_ratio_e_over_i'] = injection_triggers_df['circenergy'] - injection_triggers_df['circinc']
    injection_triggers_df['circ_ratio_i_over_e'] = injection_triggers_df['circinc'] - injection_triggers_df['circenergy']
    injection_triggers_df['circnull_ratio_e_over_i'] = injection_triggers_df['circnullenergy'] - injection_triggers_df['circnullinc']
    injection_triggers_df['circnull_ratio_i_over_e'] = injection_triggers_df['circnullinc'] - injection_triggers_df['circnullenergy']

    injection_triggers_df['plus_ratio_e_over_i'] = injection_triggers_df['coherent_f_plus'] - injection_triggers_df['incoherent_f_plus']
    injection_triggers_df['plus_ratio_i_over_e'] = injection_triggers_df['incoherent_f_plus'] - injection_triggers_df['coherent_f_plus']
    injection_triggers_df['cross_ratio_e_over_i'] = injection_triggers_df['coherent_f_cross'] - injection_triggers_df['incoherent_f_cross']
    injection_triggers_df['cross_ratio_i_over_e'] = injection_triggers_df['incoherent_f_cross'] - injection_triggers_df['coherent_f_cross']

    # Calculate all alpha values
    denominator_plus = (injection_triggers_df['coherent_f_plus'] + injection_triggers_df['incoherent_f_plus'])**0.8
    denominator_cross = (injection_triggers_df['coherent_f_cross'] + injection_triggers_df['incoherent_f_cross'])**0.8
    denominator_circenergy = (injection_triggers_df['circenergy'] + injection_triggers_df['circinc'])**0.8
    denominator_circnull = (injection_triggers_df['circnullenergy'] + injection_triggers_df['circnullinc'])**0.8

    injection_triggers_df['plus_alpha_e_over_i'] = numpy.abs(2*(injection_triggers_df['coherent_f_plus'] - injection_triggers_df['incoherent_f_plus'])/denominator_plus) + 1
    injection_triggers_df['plus_alpha_i_over_e'] = numpy.abs(2*(injection_triggers_df['incoherent_f_plus'] - injection_triggers_df['coherent_f_plus'])/denominator_plus) + 1
    injection_triggers_df['cross_alpha_e_over_i'] = numpy.abs(2*(injection_triggers_df['coherent_f_cross'] - injection_triggers_df['incoherent_f_cross'])/denominator_cross) + 1
    injection_triggers_df['cross_alpha_i_over_e'] = numpy.abs(2*(injection_triggers_df['incoherent_f_cross'] - injection_triggers_df['coherent_f_cross'])/denominator_cross) + 1

    injection_triggers_df['circ_alpha_e_over_i'] = (2*(injection_triggers_df['circenergy'] - injection_triggers_df['circinc'])/denominator_circenergy) + 1
    injection_triggers_df['circ_alpha_i_over_e'] = (2*(injection_triggers_df['circinc'] - injection_triggers_df['circenergy'])/denominator_circenergy) + 1
    injection_triggers_df['circnull_alpha_e_over_i'] = (2*(injection_triggers_df['circnullenergy'] - injection_triggers_df['circnullinc'])/denominator_circnull) + 1
    injection_triggers_df['circnull_alpha_i_over_e'] = (2*(injection_triggers_df['circnullinc'] - injection_triggers_df['circnullenergy'])/denominator_circnull) + 1

    # Calculate all ratio values
    background_triggers_df['circ_ratio_e_over_i'] = background_triggers_df['circenergy'] - background_triggers_df['circinc']
    background_triggers_df['circ_ratio_i_over_e'] = background_triggers_df['circinc'] - background_triggers_df['circenergy']
    background_triggers_df['circnull_ratio_e_over_i'] = background_triggers_df['circnullenergy'] - background_triggers_df['circnullinc']
    background_triggers_df['circnull_ratio_i_over_e'] = background_triggers_df['circnullinc'] - background_triggers_df['circnullenergy']

    background_triggers_df['plus_ratio_e_over_i'] = background_triggers_df['coherent_f_plus'] - background_triggers_df['incoherent_f_plus']
    background_triggers_df['plus_ratio_i_over_e'] = background_triggers_df['incoherent_f_plus'] - background_triggers_df['coherent_f_plus']
    background_triggers_df['cross_ratio_e_over_i'] = background_triggers_df['coherent_f_cross'] - background_triggers_df['incoherent_f_cross']
    background_triggers_df['cross_ratio_i_over_e'] = background_triggers_df['incoherent_f_cross'] - background_triggers_df['coherent_f_cross']

    # Calculate all alpha values
    denominator_plus = (background_triggers_df['coherent_f_plus'] + background_triggers_df['incoherent_f_plus'])**0.8
    denominator_cross = (background_triggers_df['coherent_f_cross'] + background_triggers_df['incoherent_f_cross'])**0.8
    denominator_circenergy = (background_triggers_df['circenergy'] + background_triggers_df['circinc'])**0.8
    denominator_circnull = (background_triggers_df['circnullenergy'] + background_triggers_df['circnullinc'])**0.8

    background_triggers_df['plus_alpha_e_over_i'] = numpy.abs(2*(background_triggers_df['coherent_f_plus'] - background_triggers_df['incoherent_f_plus'])/denominator_plus) + 1
    background_triggers_df['plus_alpha_i_over_e'] = numpy.abs(2*(background_triggers_df['incoherent_f_plus'] - background_triggers_df['coherent_f_plus'])/denominator_plus) + 1
    background_triggers_df['cross_alpha_e_over_i'] = numpy.abs(2*(background_triggers_df['coherent_f_cross'] - background_triggers_df['incoherent_f_cross'])/denominator_cross) + 1
    background_triggers_df['cross_alpha_i_over_e'] = numpy.abs(2*(background_triggers_df['incoherent_f_cross'] - background_triggers_df['coherent_f_cross'])/denominator_cross) + 1

    background_triggers_df['circ_alpha_e_over_i'] = (2*(background_triggers_df['circenergy'] - background_triggers_df['circinc'])/denominator_circenergy) + 1
    background_triggers_df['circ_alpha_i_over_e'] = (2*(background_triggers_df['circinc'] - background_triggers_df['circenergy'])/denominator_circenergy) + 1
    background_triggers_df['circnull_alpha_e_over_i'] = (2*(background_triggers_df['circnullenergy'] - background_triggers_df['circnullinc'])/denominator_circnull) + 1
    background_triggers_df['circnull_alpha_i_over_e'] = (2*(background_triggers_df['circnullinc'] - background_triggers_df['circnullenergy'])/denominator_circnull) + 1


    return injection_triggers_df, background_triggers_df
