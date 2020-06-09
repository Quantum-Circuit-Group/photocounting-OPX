from qm.qua import *
import numpy as np


averaging = 10000
wait_min = 0
wait_max = 30000
n_points = 51
wait_list = np.linspace(
    wait_min // 4,
    wait_max // 4,
    n_points,
    dtype=int).tolist()
delay = 20 // 4

with program() as envelop:
    N = declare(int)
    w = declare(int)
    I1 = declare(fixed)
    Q1 = declare(fixed)
    I2 = declare(fixed)
    Q2 = declare(fixed)

    with for_each_(w, wait_list):
        save(w, 'wait')

        with for_(N, 0, N < averaging, N + 1):
            align('buffer', 'pump', 'qubit_long')
            wait(20, 'buffer')
            play('arbitrary', 'buffer')
            wait(delay, 'pump')
            play('square_swap', 'pump')

            align('pump', 'qubit_long')
            wait(w, 'qubit_long')
            play('selective_pi', 'qubit_long')

            align('qubit_long', 'readout')
            play('clear_ro', 'readout')
            measure('meas_clear', 'readout', None,
                    ('optiW1', I1),
                    ('optiW2', Q1))

            save(I1, 'I1')
            save(Q1, 'Q1')
            wait(12500, 'qubit_long')

            align('buffer', 'pump', 'qubit_long')
            wait(20, 'buffer')
            play('arbitrary', 'buffer')
            wait(delay, 'pump')
            play('square_swap', 'pump')

            align('pump', 'qubit_long')
            wait(w, 'qubit_long')
            wait(200, 'qubit_long')

            align('qubit_long', 'readout')
            play('clear_ro', 'readout')
            measure('meas_clear', 'readout', None,
                    ('optiW1', I2),
                    ('optiW2', Q2))

            save(I2, 'I2')
            save(Q2, 'Q2')
            wait(12500, 'qubit_long')
