import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pynlo_util.pynlo_util as pyutil
import numpy as np
import matplotlib.pyplot as plt
import pynlo
from scipy.signal import hilbert
from matplotlib.colors import LogNorm
import json

c = 299792458.0 # m/s

# Input pulse parameters
fwhm_in = 0.120 # ps
center_wl = 3500 # nm
energy_in = 120e-06 # J
gdd = -160 * 1e-06 # ps^2
tod = -2900 * 1e-09 # ps^3

# Grid parameters
time_span = 10.0 # ps
prop_steps = 80
time_steps = 2**13

# Si parameters
n2_Si = 10e-19 # m^2 / W
beta2_Si = 441.71 # ps^2 / km
beta3_Si = 0.85 # ps^3 / km
beta4_Si = -0.01 # ps^4 / km

# YAG parameters
n2_YAG = 0.5e-19 # m^2 / W
beta2_YAG = -645.59         # (ps^2/km)
beta3_YAG = 4.55            # (ps^3/km)
beta4_YAG = -0.04

# CaF2 paramaters
n2_CAF2 = 0.2e-19
beta2_CAF2 = -179.86
beta3_CAF2 = 1.11
beta4_CAF2 = -0.01
radius_final = 0.0035

# Simulation parameters
lengths = [0.002, 0.002, 0.001, 0.004]
peak_intensities = [7e13, 1000e13, 50e13]
raman = True
steep = True
alpha = 0.001

# start simulation
input_pulse = pynlo.light.DerivedPulses.GaussianPulse(1, fwhm_in, center_wl, time_window_ps = time_span,
                                                GDD = gdd, TOD = tod, NPTS = time_steps, frep_MHz = 100, power_is_avg = False)
input_pulse.set_epp(energy_in)

radius_1 = pyutil.calc_radius(energy_in, fwhm_in * 1e-12, peak_intensities[0])
print('radius 1 (m): ' +  str(radius_1))
gamma_1 = pyutil.eff_nonlin(n2_Si, center_wl * 1e-09, radius_1)
print('gamma 1 (W/m): ' + str(gamma_1))
fiber_1 = pynlo.media.fibers.fiber.FiberInstance()
fiber_1.generate_fiber(lengths[0], center_wl_nm=center_wl, betas=(beta2_Si, beta3_Si, beta4_Si),
                       gamma_W_m=gamma_1, gvd_units="ps^n/km", gain=-alpha)

evol = pynlo.interactions.FourWaveMixing.SSFM.SSFM(local_error=0.001, USE_SIMPLE_RAMAN=True,
                                                    disable_Raman=np.logical_not(raman),
                                                    disable_self_steepening=np.logical_not(steep))
y1, AW1, AT1, output_1 = evol.propagate(pulse_in=input_pulse, fiber=fiber_1, n_steps=prop_steps)

output_power_1 = np.abs(output_1.AT)**2
output_duration_1 = pyutil.fwhm(input_pulse.T_ps, output_power_1) * 1e-12 # (s)

radius_2 = pyutil.calc_radius(output_1.calc_epp(), output_duration_1, peak_intensities[1])
print('radius 2 (m): ' + str(radius_2))
gamma_2 = pyutil.eff_nonlin(n2_YAG, center_wl * 1e-09, radius_2)
print('gamma 2 (W/m): ' + str(gamma_2))
fiber_2 = pynlo.media.fibers.fiber.FiberInstance()
fiber_2.generate_fiber(lengths[1], center_wl_nm=center_wl, betas=(beta2_YAG, beta3_YAG, beta4_YAG),
                       gamma_W_m=gamma_2, gvd_units="ps^n/km", gain=-alpha)

y2, AW2, AT2, output_2 = evol.propagate(pulse_in=output_1, fiber=fiber_2, n_steps=prop_steps)

output_power_2 = np.abs(output_2.AT)**2
output_duration_2 = pyutil.fwhm(input_pulse.T_ps, output_power_2) * 1e-12 # (s)

radius_3 = pyutil.calc_radius(output_2.calc_epp(), output_duration_2, peak_intensities[2])
print('radius 3 (m): ' + str(radius_3))
gamma_3 = pyutil.eff_nonlin(n2_Si, center_wl * 1e-09, radius_3)
print('gamma 3 (W/m): ' + str(gamma_3))
fiber_3 = pynlo.media.fibers.fiber.FiberInstance()
fiber_3.generate_fiber(lengths[2], center_wl_nm=center_wl, betas=(beta2_Si, beta3_Si, beta4_Si),
                       gamma_W_m=gamma_3, gvd_units="ps^n/km", gain=-alpha)

y3, AW3, AT3, output_3 = evol.propagate(pulse_in=output_2, fiber=fiber_3, n_steps=prop_steps)

gamma_4 = pyutil.eff_nonlin(n2_CAF2, center_wl * 1e-09, radius_final)
print('radius 4 (m): ' + str(radius_final))
print('gamma 4 (W/m): ' + str(gamma_4))
fiber_4 = pynlo.media.fibers.fiber.FiberInstance()
fiber_4.generate_fiber(lengths[3], center_wl_nm=center_wl, betas=(beta2_CAF2, beta3_CAF2, beta4_CAF2),
                       gamma_W_m=gamma_4, gvd_units="ps^n/km", gain=-alpha)

y4, AW4, AT4, output_4 = evol.propagate(pulse_in=output_3, fiber=fiber_4, n_steps=prop_steps)

# Calculate final pulse duration
output_power_4 = np.abs(output_4.AT)**2
output_duration_4 = pyutil.fwhm(input_pulse.T_ps, output_power_4) # (ps)
print("Final Pulse Durations (ps): " + str(output_duration_4))

# Create Figure
F = input_pulse.W_mks / (2 * np.pi) * 1e-12
zW0 = np.abs(np.transpose(AW1)[:, (F > 0)])
zW = np.abs(np.transpose(AW4)[:, (F > 0)])
plt.figure(figsize=(10,5))
plt.plot( (c / (np.array(F[F > 0])[::-1] * 1e12)) * 1e6,
         np.array(zW[-1])[::-1]/np.max(zW[-1]),
         color='r', label="Output Spectrum")
plt.plot( (c / (np.array(F[F > 0])[::-1] * 1e12)) * 1e6,
         np.array(zW0[0])[::-1]/np.max(zW0[0]),
         color='b', label="Input Spectrum")
plt.yscale('log')
plt.xlabel("Wavelength (um)")
plt.ylabel("Normalized Intensity (a.u.)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()