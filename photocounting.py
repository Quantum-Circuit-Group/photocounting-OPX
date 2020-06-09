from qm.qua import *

averaging = 20000
threshold = 1e-5

with program() as catch_photocount_and_release:
    N = declare(int)
    I0 = declare(fixed)
    I1 = declare(fixed)
    Q1 = declare(fixed)
    I2 = declare(fixed)
    Q2 = declare(fixed)

    with for_(N, 0, N < averaging, N + 1):
        play('clear_ro', 'readout')
        measure('meas_clear', 'readout', None,
                ('optiW1', I0))

        # As long as we're not certain to be in |g>
        with for_(cond=I0 < threshold):
            play('unconditional_pi', 'qubit_short')
            align('readout', 'qubit_short')
            play('clear_ro', 'readout')
            measure('meas_clear', 'readout', None,
                    ('optiW1', I0))

        align('buffer', 'pump', 'readout')
        play('sech_kick', 'buffer')
        play('sech_catch', 'pump')

        align('pump', 'qubit_short')
        play('mod_2', 'qubit_short')
        align('qubit_short', 'readout')

        play('clear_ro', 'readout')
        measure('meas_clear', 'readout', None,
                ('optiW1', I1),
                ('optiW2', Q1))

        # If we are in |e> (0 mod 2)
        with if_(I1 < 0):
            play('mod_4', 'qubit_short')
        with else_():
            play('mod_4_y', 'qubit_short')
        align('qubit_short', 'readout')

        play('clear_ro', 'readout')
        measure('meas_clear', 'readout', None,
                ('optiW1', I2),
                ('optiW2', Q2))

        align('pump', 'readout')
        play('sech_swap', 'pump')

        save(I1, 'I1')
        save(Q1, 'Q1')
        save(I2, 'I2')
        save(Q2, 'Q2')
