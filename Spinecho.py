import numpy as np

import matplotlib.pyplot as plt

import pdb



from math import pi

import pypulseq as pp

from pypulseq import Sequence

from pypulseq.calc_duration import calc_duration

from pypulseq.make_adc import make_adc

from pypulseq.make_delay import make_delay

from pypulseq.make_sinc_pulse import make_sinc_pulse

from pypulseq.make_block_pulse import make_block_pulse

from pypulseq.make_trapezoid import make_trapezoid

from pypulseq.opts import Opts



# User-defined parameters

TR = 1  # Repetition time (s)

TE = 10e-3  # Echo time (s)

p90=0.3e-3

p180=0.1e-3

flip90 = 90 * pi / 180  # Flip angle for 90-degree pulse

flip180 = 180 * pi / 180  # Flip angle for 180-degree pulse

slice_thickness = 3e-3  # Slice thickness (m)

FOV = 20e-3  # Field of view (m)

Nx = 128  # Readout points

Ny = 128  # Phase encoding points

TBW = 7  # Time-bandwidth product

BW = 50e3  # Bandwidth (Hz)

pre_dur = 0.3e-3  # Prephasing duration (s)

phase_dur = 0.3e-3  # Phase encoding duration (s)

rewinder_length = 0.3e-3  # Rewinder duration (s)

n_slices = 3  # Number of slices

slice_gap = 5e-3  # Gap between slices (m)



# Delays and amplitudes from gr-mri

p90_delay=0.1e-3

gr_ramp=p90_delay

p180_delay = 5.2e-3

readout_delay = 8.968e-3

readout_duration = 2.56e-3

gx_dur=[0.5e-3,2.76e-3]

gy_step = 0.005927755480086448

gy_dur=0.34e-3

phase_grad_delay = 0.5e-3

gy_start_amp = -0.37937635072553266*1e3

gy_amp=gy_start_amp+gy_step

slice_grad_dead_time = 4.9e-3

slice_grad_amps = [0.7466666666666667*1e3, 0.5656565656565659*1e3]

gz_dur = 0.5e-3



# System limits

system = Opts(

    max_grad=320,

    grad_unit='mT/m',

    max_slew=2000,

    slew_unit='T/m/s',

)

seq = Sequence(system)



# Define RF pulses

rf90,gz90,__ = make_sinc_pulse(

    flip_angle=flip90,

    system=system,

    duration=300e-6,

    slice_thickness=slice_thickness,

    apodization=0.5,

    time_bw_product=TBW,

    return_gz=True,

)

gz90.amplitude = slice_grad_amps[0]

rf180 = make_block_pulse(

    flip_angle=flip180,

    system=system,

    duration=100e-6

)



# Define gradients

delta_k = 1 / FOV

k_width = Nx * delta_k

gx = make_trapezoid(channel='x', system=system,amplitude=240,duration=2.76e-3,max_slew=240/gr_ramp,flat_time=readout_duration)

adc = make_adc(num_samples=Nx, duration=gx.flat_time, delay=gx.rise_time)

phase_areas = (np.arange(Ny) - (Ny / 2)) * delta_k



# Prephasing gradients

gx_pre = make_trapezoid(channel='x', system=system,flat_time=pre_dur,amplitude= -804.84,max_slew=804.84/gr_ramp,duration=gy_dur)

gy_pre = make_trapezoid(channel='y', system=system,flat_time=phase_dur,amplitude= gy_amp,max_slew=-gy_amp*5/gr_ramp,duration=phase_grad_delay)

# Spoiler gradient

gz_spoil = make_trapezoid(channel='z', system=system, area=gz90.area * 4, duration=gz_dur,amplitude=slice_grad_amps[1],max_slew=slice_grad_amps[1]/(gr_ramp),flat_time=0.3e-3)

# Delays

delay1 = TE / 2 + calc_duration(rf90) / 2 - calc_duration(rf180)/2 + gr_ramp

t0=delay1-calc_duration(gy_pre)-calc_duration(gz90)

dead=4.9e-3

t1=dead-t0-calc_duration(rf180)-calc_duration(gy_pre)

delay3 = 8.368e-3 -calc_duration(gz90)-t0-t1-calc_duration(gz_spoil)-calc_duration(rf180)-calc_duration(gy_pre)

delay_TR = TR - TE - calc_duration(rf90) / 2 - calc_duration(gx) / 2 - calc_duration(gz_spoil)

delay_TR = make_delay(delay_TR)

# Slice selection

delta_z = n_slices * slice_gap

z_positions = np.linspace(-delta_z / 2, delta_z / 2, n_slices)



for slice_idx, z in enumerate(z_positions):

    rf90.freq_offset = gz90.amplitude * z

    rf180.freq_offset = gz90.amplitude * z



    for pe_idx, phase_area in enumerate(phase_areas):

        # Excitation

        # Define sequence blocks



        gy_pre.amplitude=gy_pre.amplitude+gy_step

        seq.add_block(rf90, gz90)

        seq.add_block(gy_pre)

        seq.add_block(make_delay(t0))

        #delay1=delay1-pp.calc_duration(gy_pre)

        seq.add_block(rf180)

        seq.add_block(make_delay(t1))

        seq.add_block(gz_spoil)

        seq.add_block(make_delay(delay3))

         # Readout

        seq.add_block(gx_pre)

        seq.add_block(gx,adc)

        #seq.add_block(make_delay(delay2))





        # Repetition delay

        seq.add_block(delay_TR)



# Plot and output sequence

seq.plot(time_range=[0.0e-3,12e-3], time_disp='ms')



print(gz90.rise_time)
