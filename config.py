import numpy as np


def IQ_imbalance(g, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    N = 1 / ((1 - g ** 2) * (2 * c ** 2 - 1))
    return [float(N * x) for x in [(1 - g) * c, (1 + g) * s,
                                   (1 - g) * s, (1 + g) * c]]


def catch_sech_shape(amp, total_time, kappa, lambd):
    if total_time < 16:
        pad_time = 16 - total_time
    else:
        pad_time = (4 - total_time % 4) % 4
    t = np.arange(total_time) - total_time // 2
    tau = kappa / 2. * t
    pulse = np.sqrt(lambd / 2. / (np.exp(lambd * tau) + 1. - lambd / 2.))
    pulse = pulse * (1 + lambd / 2. * np.tanh(lambd * tau / 2.))
    return list(amp * pulse) + [0.0] * pad_time


def sech_shape(amp, total_time, n_sigma, padding=True):
    pad_time = 0
    if padding:
        if total_time < 16:
            pad_time = 16 - total_time
        else:
            pad_time = (4 - total_time % 4) % 4
    t = np.array(range(-total_time // 2, total_time // 2, 1))
    sigma = total_time / n_sigma
    return list(amp / np.cosh(np.sqrt(np.pi / 2) * t / sigma)) + \
        [0.0] * pad_time

def ramsey_shape(amp1, amp2, time1, time2, time_wait, n_sigma):
    padding = (4 - (time1 + time2 + time_wait) % 4) % 4
    return sech_shape(amp1, time1, n_sigma, padding=False) + [0.0] * time_wait\
        + sech_shape(amp2, time2, n_sigma, padding=False) + [0.0] * padding


def ramsey_drag_shape(drag_amp1, drag_amp2, time1, time2, time_wait, n_sigma):

    padding = (4 - (time1 + time2 + time_wait) % 4) % 4
    return drag_sech(drag_amp1, time1, n_sigma=n_sigma) + [0.0] * time_wait\
        + drag_sech(drag_amp2, time2, n_sigma=n_sigma) + [0.0] * padding

def drag_sech(amp, total_time, n_sigma=4):
    t = np.array(range(-total_time // 2, total_time // 2, 1))
    sigma = total_time / n_sigma
    return list(-amp * np.tanh(np.sqrt(np.pi / 2) * t / sigma) /
                np.cosh(np.sqrt(np.pi / 2) * t / sigma))


def padded_time(t):
    return max(t + (4 - t % 4) % 4, 16)


parameters = {
    'sech_swap_amp': 0.1,
    'disp_amp': 0.0558,
    'phase_mod2_onoff': 1,
    'phase_mod4_onoff': 1,
    'correction_angle': -1.44 + 0.25,
    'buffer_sech_amp': 0.1,

    'delay_catch': 15,
}

optimized_parameters = {

    'offset_memory': 0.017,
    'offset_buffer_I': -0.0037,
    'offset_buffer_Q': -0.0016,
    'offset_pump_I': -0.0047,
    'offset_pump_Q': -0.0084,
    'offset_ampli': 0.,
    'offset_qubit_I': -0.0013,
    'offset_qubit_Q': -0.0038,
    'drag_amp': 0.05,

    'qubit_if': 200e6,

    'wait_mod2': 140,
    'wait_mod4': 64,
    'buffer_I': 1,
    'unconditional_pi_2_length': 16,
    'unconditional_pi_2_amp': 0.2267,
    'unconditional_pi_length': 32,
    'unconditional_pi_amp': 0.2,

    'buffer_sech_length': 300,
    'theta': 0,

    'y_scaling': 1,

    'disp_length': 52,
    'sech_catch_amp': 0.165,
    'sech_catch_length': 300,

    'lambd': 0.532,
    'kappa': 0.02 * 2 * np.pi,
    'ki': 0.001 * 2 * np.pi,
    'meas_length': 252,

    'sech_swap_length': 240,

    'memory_if': 120e6,
    'pump_if': 170e6,
    'buffer_if': 50e6,

    'short_ro_length': 252,
    'short_ro_amp': 0.16,

    'clear_ro_amp_up': 0.16,
    'clear_ro_amp_steady': 0,
    'clear_ro_amp_down': -0.115,
    'clear_ro_amp_down_steady': 0.055,
    'clear_ro_length_up': 252,
    'clear_ro_length_steady': 0,
    'clear_ro_length_down': 100,
    'clear_ro_length_down_steady': 52,

    'ro_if': 51e6,

    'weights_path': 'weights.npy',
}


def get_parameters():
    return parameters


def get_config(params):
    params.update(optimized_parameters)

    offset_qubit_I = float(params['offset_qubit_I'])
    offset_qubit_Q = float(params['offset_qubit_Q'])
    offset_memory = float(params['offset_memory'])
    offset_buffer_I = float(params['offset_buffer_I'])
    offset_buffer_Q = float(params['offset_buffer_Q'])
    offset_pump_I = float(params['offset_pump_I'])
    offset_pump_Q = float(params['offset_pump_Q'])
    offset_ampli = float(params['offset_ampli'])

    drag_amp = float(params['drag_amp'])
    y_scaling = float(params['y_scaling'])
    theta = float(params['theta'])

    qubit_if = float(params['qubit_if'])
    buffer_if = float(params['buffer_if'])
    pump_if = float(params['pump_if'])
    memory_if = float(params['memory_if'])

    ro_if = float(params['ro_if'])

    clear_ro_length_up = int(params['clear_ro_length_up'])
    clear_ro_amp_up = float(params['clear_ro_amp_up'])
    clear_ro_length_steady = int(params['clear_ro_length_steady'])
    clear_ro_amp_steady = float(params['clear_ro_amp_steady'])

    clear_ro_length_down = int(params['clear_ro_length_down'])
    clear_ro_amp_down = float(params['clear_ro_amp_down'])
    clear_ro_length_down_steady = int(params['clear_ro_length_down_steady'])
    clear_ro_amp_down_steady = float(params['clear_ro_amp_down_steady'])

    short_ro_length = int(params['short_ro_length'])
    short_ro_amp = float(params['short_ro_amp'])

    unconditional_pi_length = int(params['unconditional_pi_length'])
    unconditional_pi_amp = float(params['unconditional_pi_amp'])
    unconditional_pi_2_length = int(params['unconditional_pi_2_length'])
    unconditional_pi_2_amp = float(params['unconditional_pi_2_amp'])

    meas_length = int(params['meas_length'])

    buffer_I = float(params['buffer_I'])

    wait_mod2 = int(params['wait_mod2'])
    wait_mod4 = int(params['wait_mod4'])
    phase_mod2_onoff = int(params['phase_mod2_onoff'])
    phase_mod4_onoff = int(params['phase_mod4_onoff'])

    correction_angle = float(params['correction_angle'])
    if params['weights_path'] != '':
        weights = np.load(params['weights_path'])
    else:
        weights = np.array([1.0] * (padded_time(meas_length) // 4))
    weights1 = (np.cos(correction_angle) * weights).tolist()
    weights2 = (np.sin(correction_angle) * weights).tolist()

    sech_catch_amp = float(params['sech_catch_amp'])
    sech_catch_length = int(params['sech_catch_length'])
    sech_swap_amp = float(params['sech_swap_amp'])
    sech_swap_length = int(params['sech_swap_length'])

    buffer_sech_amp = float(params['buffer_sech_amp'])
    buffer_sech_length = int(params['buffer_sech_length'])

    delay_catch = int(params['delay_catch'])
    padding_catch = (4 - delay_catch % 4) % 4

    lambd = float(params['lambd'])
    kappa = float(params['kappa'])
    ki = float(params['ki'])

    config1 = {
        'version': 1,

        'controllers': {
            'con1': {
                'type': 'opx1',
                'analog_outputs': {
                    1: {'offset': offset_qubit_I},
                    2: {'offset': offset_qubit_Q},
                    4: {'offset': offset_memory},
                    5: {'offset': offset_buffer_I},
                    6: {'offset': offset_buffer_Q},
                    7: {'offset': offset_pump_I},
                    8: {'offset': offset_pump_Q},
                },
                'digital_outputs': {
                    2: {}
                },
                'analog_inputs': {
                    1: {'offset': 0.1945},
                    2: {'offset': 0.196},
                }
            }
        },

        'elements': {
            'qubit_short': {
                'mixInputs': {
                    'I': ('con1', 1),
                    'Q': ('con1', 2),
                    'mixer': 'mixer_qubit_short',
                    'lo_frequency': 4.4792e9
                },
                'digitalInputs': {
                    'port': {
                        'port': ('con1', 2),
                        'buffer': 0,
                        'delay': 130
                    }
                },
                'intermediate_frequency': qubit_if,
                'operations': {
                    'unconditional_pi': 'unconditional_pi_pulse',
                    'mod_2': 'mod_2_pulse',
                    'mod_4': 'mod_4_pulse',
                    'mod_4_y': 'mod_4_y_pulse',
                }
            },
            'readout': {
                'singleInput': {
                    'port': ('con1', 3)
                },
                'intermediate_frequency': ro_if,
                'operations': {
                    'meas_clear': 'meas_clear_pulse',
                    'clear_ro': 'clear_ro_pulse',
                },
                'outputs': {
                    'out1': ('con1', 2)
                },
                'digitalInputs': {
                    'port': {
                        'port': ('con1', 2),
                        'buffer': 0,
                        'delay': 130
                    }
                },
                'time_of_flight': 56,
                'smearing': 0,
            },
            'pump': {
                'mixInputs': {
                    'I': ('con1', 7),
                    'Q': ('con1', 8),
                    'mixer': 'mixer_pump',
                    'lo_frequency': 10.222e9 - 3.627e9 + 170e6
                },
                'digitalInputs': {
                    'port': {
                        'port': ('con1', 2),
                        'buffer': 0,
                        'delay': 130
                    }
                },
                'intermediate_frequency': pump_if,
                'operations': {
                    'sech_catch': 'sech_catch_pulse',
                    'sech_swap': 'sech_swap_pulse',
                }
            },
            'buffer': {
                'mixInputs': {
                    'I': ('con1', 5),
                    'Q': ('con1', 6),
                    'mixer': 'mixer_buffer',
                    'lo_frequency': 10.222e9 + 50e6
                },
                'digitalInputs': {
                    'port': {
                        'port': ('con1', 2),
                        'buffer': 0,
                        'delay': 130
                    }
                },
                'intermediate_frequency': buffer_if,
                'operations': {
                    'sech_kick': 'buffer_sech_pulse',
                },
                'outputs': {
                    'out1': ('con1', 1)
                },
                'time_of_flight': 48,
                'smearing': 0
            },
        },
        'pulses': {
            'mod_2_pulse': {
                'operation': 'control',
                'length': padded_time(wait_mod2 + 2 * unconditional_pi_2_length),
                'waveforms': {
                    'I': 'mod2_I_wf',
                    'Q': 'mod2_Q_wf'
                },
                'digital_marker': 'trigger'
            },
            'mod_4_pulse': {
                'operation': 'control',
                'length': padded_time(wait_mod4 + 2 * unconditional_pi_2_length),
                'waveforms': {
                    'I': 'mod4_I_wf',
                    'Q': 'mod4_Q_wf'
                },
                'digital_marker': 'trigger'
            },
            'mod_4_y_pulse': {
                'operation': 'control',
                'length': padded_time(wait_mod4 + 2 * unconditional_pi_2_length),
                'waveforms': {
                    'I': 'mod4_y_I_wf',
                    'Q': 'mod4_y_Q_wf'
                },
                'digital_marker': 'trigger'
            },
            'unconditional_pi_pulse': {
                'operation': 'control',
                'length': padded_time(unconditional_pi_length),
                'waveforms': {
                    'I': 'unconditional_pi_wf',
                    'Q': 'zero_wf'
                },
                'digital_marker': 'trigger'
            },
            'clear_ro_pulse': {
                'operation': 'measurement',
                'length': padded_time(clear_ro_length_up + clear_ro_length_steady),
                'waveforms': {
                    'single': 'clear_ro_wf'
                },
                'integration_weights': {
                },
                'digital_marker': 'ON'
            },
            'meas_clear_pulse': {
                'operation': 'measurement',
                'length': padded_time(meas_length),
                'waveforms': {
                    'single': 'meas_clear_wf'
                },
                'integration_weights': {
                    'optiW1': 'opti_clear_W1',
                    'optiW2': 'opti_clear_W2',
                },
                'digital_marker': 'ON'
            },
            'sech_catch_pulse': {
                'operation': 'control',
                'length': padded_time(sech_catch_length),
                'waveforms': {
                    'I': 'sech_catch_wf',
                    'Q': 'zero_wf'
                },
                'digital_marker': 'trigger'
            },
            'sech_swap_pulse': {
                'operation': 'control',
                'length': padded_time(sech_swap_length),
                'waveforms': {
                    'I': 'sech_swap_wf',
                    'Q': 'zero_wf'
                },
                'digital_marker': 'trigger'
            },
            'buffer_sech_pulse': {
                'operation': 'control',
                'length': padded_time(buffer_sech_length) + delay_catch + padding_catch,
                'waveforms': {
                    'I': 'buffer_sech_I_wf',
                    'Q': 'buffer_sech_Q_wf'
                },
                'digital_marker': 'trigger'
            },
        },
        'waveforms': {
            'short_ro_wf': {
                'type': 'constant',
                'sample': short_ro_amp
            },
            'clear_ro_wf': {
                'type': 'arbitrary',
                'samples': [clear_ro_amp_up] * clear_ro_length_up\
                + [clear_ro_amp_steady] * clear_ro_length_steady\
            },
            'meas_clear_wf': {
                'type': 'arbitrary',
                'samples': [clear_ro_amp_down] * clear_ro_length_down\
                + [clear_ro_amp_down_steady] * clear_ro_length_down_steady\
                + [0] * (padded_time(meas_length) - clear_ro_length_down - clear_ro_length_down_steady)
            },
            'zero_wf': {
                'type': 'constant',
                'sample': 0.0
            },
            'mod2_I_wf': {
                'type': 'arbitrary',
                'samples': ramsey_shape(unconditional_pi_2_amp,
                                        (-1.)**int(phase_mod2_onoff) * unconditional_pi_2_amp,
                                        unconditional_pi_2_length, unconditional_pi_2_length, wait_mod2, 4)
            },
            'mod2_Q_wf': {
                'type': 'arbitrary',
                'samples': ramsey_drag_shape(drag_amp,
                                             (-1.)**int(phase_mod2_onoff) * drag_amp,
                                             unconditional_pi_2_length, unconditional_pi_2_length, wait_mod2, 4)
                #
            },
            'mod4_I_wf': {
                'type': 'arbitrary',
                'samples': list(np.array(ramsey_shape(unconditional_pi_2_amp,
                                                      (-1)**phase_mod4_onoff * unconditional_pi_2_amp * (1 - phase_mod2_onoff),
                                                      unconditional_pi_2_length, unconditional_pi_2_length, wait_mod4, 4))\
                                + np.array(ramsey_drag_shape(0,
                                                             -(-1)**phase_mod4_onoff * drag_amp * (phase_mod2_onoff),
                                                             unconditional_pi_2_length, unconditional_pi_2_length, wait_mod4, 4)))
            },
            'mod4_Q_wf': {
                'type': 'arbitrary',
                'samples': list(np.array(ramsey_shape(0,
                                                      (-1)**phase_mod4_onoff * unconditional_pi_2_amp * phase_mod2_onoff,
                                                      unconditional_pi_2_length, unconditional_pi_2_length, wait_mod4, 4))\
                                + np.array(ramsey_drag_shape(drag_amp,
                                                             (-1)**phase_mod4_onoff * drag_amp * (1 - phase_mod2_onoff),
                                                             unconditional_pi_2_length, unconditional_pi_2_length, wait_mod4, 4)))
            },

            # mod4_y can be on x depending on phase_mod2_onoff...
            'mod4_y_I_wf': {
                'type': 'arbitrary',
                'samples': list(np.array(ramsey_shape(unconditional_pi_2_amp,
                                                      (-1)**phase_mod4_onoff * unconditional_pi_2_amp * phase_mod2_onoff,
                                                      unconditional_pi_2_length, unconditional_pi_2_length, wait_mod4, 4))\
                                + np.array(ramsey_drag_shape(-drag_amp,
                                                             -(-1)**phase_mod4_onoff * drag_amp * (1 - phase_mod2_onoff),
                                                             unconditional_pi_2_length, unconditional_pi_2_length, wait_mod4, 4)))
            },
            'mod4_y_Q_wf': {
                'type': 'arbitrary',
                'samples': list(np.array(ramsey_shape(0, (-1)**phase_mod4_onoff * unconditional_pi_2_amp * (1 - phase_mod2_onoff),
                                                      unconditional_pi_2_length, unconditional_pi_2_length, wait_mod4, 4))\
                                + np.array(ramsey_drag_shape(drag_amp,
                                                             (-1)**phase_mod4_onoff * drag_amp * (phase_mod2_onoff),
                                                             unconditional_pi_2_length, unconditional_pi_2_length, wait_mod4, 4)))
            },
            'unconditional_pi_wf': {
                'type': 'arbitrary',
                'samples': sech_shape(unconditional_pi_amp, unconditional_pi_length, 4)
            },
            'sech_catch_wf': {
                'type': 'arbitrary',
                'samples': catch_sech_shape(sech_catch_amp, sech_catch_length, kappa, lambd)
            },
            'sech_swap_wf': {
                'type': 'arbitrary',
                'samples': sech_shape(sech_swap_amp, sech_swap_length, 4)
            },
            'buffer_sech_wf': {
                'type': 'arbitrary',
                'samples': [0.0] * delay_catch + sech_shape(buffer_sech_amp, buffer_sech_length, 4) \
                + [0.0] * padding_catch
            },
            'buffer_sech_I_wf': {
                'type': 'arbitrary',
                'samples': list((buffer_I) * np.array([0.0] * delay_catch + sech_shape(buffer_sech_amp, buffer_sech_length, 4) \
                                                      + [0.0] * padding_catch))
            },
            'buffer_sech_Q_wf': {
                'type': 'arbitrary',
                'samples': list((1 - buffer_I) * np.array([0.0] * delay_catch + sech_shape(buffer_sech_amp, buffer_sech_length, 4) \
                                                          + [0.0] * padding_catch))
            },

        },
        'digital_waveforms': {
            'trigger': {
                'samples': [(1, 20), (0, 0)]
            },
            'ON': {
                'samples': [(1, 0)]
            }
        },
        'integration_weights': {
            'optiW1': {
                'cosine': weights1,
                'sine': [-i for i in weights2]
            },
            'optiW2': {
                'cosine': weights2,
                'sine': weights1
            },
        },

        'mixers': {
            'mixer_qubit_short': [{
                'intermediate_frequency': qubit_if,
                'correction': IQ_imbalance(0.255, np.pi * 0.011),
                'lo_frequency': 4.4792e9,
            }],
            'mixer_pump': [{
                'intermediate_frequency': pump_if,
                'correction': IQ_imbalance(0.005, np.pi * 0.014),
                'lo_frequency': 10.222e9 - 3.627e9 + 170e6,
            }],
            'mixer_buffer': [{
                'intermediate_frequency': buffer_if,
                'correction': IQ_imbalance(-0.001, -np.pi * 0.054),
                'lo_frequency': 10.222e9 + 50e6,
            }],
        }
    }
    return config1
