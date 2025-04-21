!pip install PyPulseq
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

TR = 0.75  # Repetition time (s)

TE = 1.94e-3  # Echo time (s)

p90=0.3e-3

flip90 = 90 * pi / 180  # Flip angle for 90-degree pulse

flip180 = 180 * pi / 180  # Flip angle for 180-degree pulse

slice_thickness = 3e-3  # Slice thickness (m)

FOV = 20e-3  # Field of view (m)

Nx = 128  # Readout points

Ny = 128  # Phase encoding points

TBW = 8  # Time-bandwidth product

BW = 100000  # Bandwidth (Hz)

pre_dur = 0.5e-3  # Prephasing duration (s)

phase_dur = 0.5e-3  # Phase encoding duration (s)

rewinder_length = 0.5e-3  # Rewinder duration (s)

n_slices = 3  # Number of slices

slice_gap = 5e-3  # Gap between slices (m)



# Delays and amplitudes from gr-mri

p90_delay=0.1e-3
gr_ramp=p90_delay
readout_delay = 1.21e-3
readout_duration = 1.28e-3
gx_dur=[0.7e-3,1.48e-3]
gy_step = 0.0036301758366420883*1000
gy_dur=0.54e-3
phase_grad_delay = 0.5e-3
gy_start_amp = -232.331
gy_amp=gy_start_amp+gy_step
slice_grad_amps = [640,-322.2148]
gz_dur = [0.5e-3,0.7e-3]



# System limits

system = Opts(

    max_grad=900,
    grad_unit='mT/m',
    max_slew=2100,
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




# Define gradients

delta_k = 1 / FOV
k_width = Nx * delta_k
gx = make_trapezoid(channel='x', system=system,amplitude=480,duration=1.48e-3,max_slew=480/gr_ramp,flat_time=readout_duration)
adc = make_adc(num_samples=Nx, duration=gx.flat_time, delay=gx.rise_time)
phase_areas = (np.arange(Ny) - (Ny / 2)) * delta_k



# Prephasing gradients

gx_pre = make_trapezoid(channel='x', system=system,flat_time=pre_dur,amplitude= -554.1,max_slew=554.1/gr_ramp,duration=gy_dur)
gy_pre = make_trapezoid(channel='y', system=system,flat_time=phase_dur,amplitude= -232.331,max_slew=232.331*25/gr_ramp,duration=phase_grad_delay)

# Spoiler gradient
gz_spoil = make_trapezoid(channel='z', system=system, area=gz90.area * 4, duration=0.7e-3,amplitude=slice_grad_amps[1],max_slew=-slice_grad_amps[1]/(gr_ramp),flat_time=0.5e-3)

# Delays


# Slice selection
delta_z = n_slices * slice_gap
z_positions = np.linspace(-delta_z / 2, delta_z / 2, n_slices)



for slice_idx, z in enumerate(z_positions):

    rf90.freq_offset = gz90.amplitude * z

    for pe_idx, phase_area in enumerate(phase_areas):

        # Excitation

        # Define sequence blocks



        gy_pre.amplitude=gy_pre.amplitude+gy_step
        seq.add_block(rf90, gz90)
        seq.add_block(gy_pre,gz_spoil,pp.make_delay(0.3e-3),gx_pre)
        #seq.add_block(gz_spoil)


         # Readout

        #seq.add_block(make_delay(delay))

        seq.add_block(gx,adc)





        # Repetition des
        seq.add_block(TR)



# Plot and output sequence

seq.plot(time_range=[0,4e-3], time_disp='ms')



#print(gz90.rise_time)
print(gx.flat_time*1000)
