import numpy
import itertools

def choose_background_injection_training(f, randomseed=1986):
    """We need to select half of background and half of injection events for training
    """
    # set set for reproducibility
    numpy.random.seed(randomseed)

    # First we figure out how many waveforms we injected
    injection_events = [i for i in f.root.injection._v_children.keys()]

    waveforms = [waveform for waveform in f.get_node('/injection/{0}'.format(injection_events[0]))._v_children.keys()]

    inj_scales = [inj_scale for inj_scale in f.get_node('/injection/{0}/{1}'.format(injection_events[0], waveforms[0]))._v_children.keys()]

    injections = numpy.asarray([injection for injection in f.get_node('/injection/{0}/{1}/{2}'.format(injection_events[0], waveforms[0], inj_scales[0]))._v_children.keys()])

    training_injection_events_idx = numpy.random.randint(0, injections.size, size=int(0.5*injections.size))

    list_of_training_injection_paths = [['/injection'], injection_events, waveforms,
                                        inj_scales, injections[training_injection_events_idx].tolist()]

    training_injection_events = list(map(lambda x: '/'.join(x), itertools.product(*list_of_training_injection_paths)))

    list_of_validation_injection_paths = [['/injection'], injection_events, waveforms, inj_scales, injections[~training_injection_events_idx].tolist()]

    validation_injection_events = list(map(lambda x: '/'.join(x), itertools.product(*list_of_validation_injection_paths)))

    # We need to select have of background and half of injection events for training
    background_events = numpy.asarray([group for group in f.walk_groups('/background') if 'internal_slide' in group._v_name])
    training_events = numpy.random.randint(0, background_events.size, size=int(0.5*background_events.size))
    training_background_events = background_events[training_events]
    validation_background_events = background_events[~training_events]

    return training_background_events, validation_background_events, training_injection_events, validation_injection_events
