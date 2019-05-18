import numpy
import itertools
import random

def choose_background_injection_training(f, injection_type='onsource_injection', randomseed=1986,):
    """We need to select half of background and half of injection events for training
    """
    # set set for reproducibility
    numpy.random.seed(randomseed)
    random.seed(randomseed)

    try:
        for idx, node in enumerate(f.list_nodes('/')):
            if '/' + injection_type == node._v_pathname:
                break

        # First we figure out how many waveforms we injected
        injection_events = [i for i in f.list_nodes('/')[idx]._v_children.keys()]

        waveforms = [waveform for waveform in f.get_node('/{0}/{1}'.format(injection_type, injection_events[0]))._v_children.keys()]

        inj_scales = [inj_scale for inj_scale in f.get_node('/{0}/{1}/{2}'.format(injection_type, injection_events[0], waveforms[0]))._v_children.keys()]

        all_training_injection_events = []
        all_validation_injection_events = []
        for injection_event in injection_events:
            injections = numpy.asarray([injection for injection in f.get_node('/{0}/{1}/{2}/{3}'.format(injection_type, injection_event, waveforms[0], inj_scales[0]))._v_children.keys()])

            number_of_injections = injections.size
            training_injection_events_idx = numpy.array(random.sample(range(number_of_injections), max(int(0.5*number_of_injections),1)))
            validation_injection_events_idx = numpy.setxor1d(numpy.indices(numpy.arange(len(injections)).shape), training_injection_events_idx)
            
            list_of_training_injection_paths = [['/' + injection_type], [injection_event], waveforms,
                                                inj_scales, injections[training_injection_events_idx].tolist()]

            training_injection_events = list(map(lambda x: '/'.join(x), itertools.product(*list_of_training_injection_paths)))
            all_training_injection_events.extend(training_injection_events)

            list_of_validation_injection_paths = [['/' + injection_type], [injection_event], waveforms, inj_scales, injections[validation_injection_events_idx].tolist()]

            validation_injection_events = list(map(lambda x: '/'.join(x), itertools.product(*list_of_validation_injection_paths)))
            all_validation_injection_events.extend(validation_injection_events)
    except:
        print('This file does not have any injections will return no training or testing injection events')
        all_training_injection_events = []
        all_validation_injection_events = []

    try:
        # We need to select have of background and half of injection events for training
        background_events = numpy.asarray([group for group in f.walk_groups('/background') if 'internal_slide' in group._v_name])
        number_of_background_events = background_events.size
        # randomly select half of these events for training
        training_events = numpy.array(random.sample(range(number_of_background_events), int(0.5*number_of_background_events)))
        # Use the other indices for testing
        validation_events = numpy.setxor1d(numpy.indices(numpy.arange(number_of_background_events).shape), training_events)
        # get training event table names
        training_background_events = background_events[training_events]
        # get testing
        validation_background_events = background_events[validation_events]
    except:
        print('This file does not have any background maps, will return no training or testing background events')
        training_background_events = []
        validation_background_events = []

    return training_background_events, validation_background_events, all_training_injection_events, all_validation_injection_events
