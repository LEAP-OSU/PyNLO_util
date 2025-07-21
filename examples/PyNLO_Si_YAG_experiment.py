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

# Input Pulse Parameters (before first Si plate)
################################################
FWHM_DURATION = 0.120   # (ps)
CENTER_WL = 3500        # (nm)
PULSE_ENERGY = 120e-06    # (J)
GDD = 0.0              # (ps^2)
TOD = 0.0         # (ps^3)

# Grid Parameters
#################
TIME_SPAN = 10.0       # (ps)
PROP_STEPS = 50     # This is how many steps the Z distance is divided into
TIME_STEPS = 2**13  # This is how many temporal points there are

# First Fiber Parameters (Si, 2mm)
########################
FIBER_LENGTH_1 = 0.002    # (m)
KERR_COEFF_Si = 10e-19     # (m^2/W)
ALPHA = 0.0             # attenuation const (dB/cm)
BETA2_Si = 441.71         # (ps^2/km)
BETA3_Si = 0.85            # (ps^3/km)
BETA4_Si = -0.01           # (ps^4/km)
radius = pyutil.calc_radius(PULSE_ENERGY, FWHM_DURATION * 1e-12, 7e13)
print(radius)
GAMMA_1 = pyutil.eff_nonlin(KERR_COEFF_Si, CENTER_WL * 1e-9, radius)   # (1/Wm)
print("Si (1) eff nonlinearity: " + str(GAMMA_1) + " (1/Wm)")
RAMAN = True
STEEP = True

# Second Fiber Parameters (YAG, 2mm)
####################################
FIBER_LENGTH_2 = 0.002    # (m)
KERR_COEFF_YAG = 0.5e-19     # (m^2/W)
BETA2_YAG = -645.59         # (ps^2/km)
BETA3_YAG = 4.55            # (ps^3/km)
BETA4_YAG = -0.04
PEAK_INTENSITY_2 = 1000 * 1e13  #(W/m^2)
# Gamma calculated after first propagation

# Third Fiber Parameters (Si, 1mm)
##################################
FIBER_LENGTH_3 = 0.001      # (m)
PEAK_INTENSITY_3 = 50 * 1e13    #(W/m^2)
# Gamma calculated after second propagation

# Fourth FIber Parameters (CaF2)
################################
FIBER_LENGTH_4 = 0.004
KERR_COEFF_CAF2 = 0.2e-19
PEAK_INTENSITY_4 = PEAK_INTENSITY_3 / 100
BETA2_CAF2 = -179.86
BETA3_CAF2 = 1.11
BETA4_CAF2 = -0.01



# Create Pulse
##############
input_pulse = pynlo.light.DerivedPulses.GaussianPulse(1, FWHM_DURATION, CENTER_WL, time_window_ps = TIME_SPAN,
                                                GDD = GDD, TOD = TOD, NPTS = TIME_STEPS, frep_MHz = 100, power_is_avg = False)
input_pulse.set_epp(PULSE_ENERGY)

# Create first fiber
####################
fiber_1 = pynlo.media.fibers.fiber.FiberInstance()
fiber_1.generate_fiber(FIBER_LENGTH_1, center_wl_nm=CENTER_WL, betas=(BETA2_Si, BETA3_Si, BETA4_Si),
                        gamma_W_m=GAMMA_1, gvd_units="ps^n/km", gain=-ALPHA)

# Initiate propagation
######################
evol = pynlo.interactions.FourWaveMixing.SSFM.SSFM(local_error=0.001, USE_SIMPLE_RAMAN=True,
                                                    disable_Raman=np.logical_not(RAMAN),
                                                    disable_self_steepening=np.logical_not(STEEP))
y1, AW1, AT1, output_1 = evol.propagate(pulse_in=input_pulse, fiber=fiber_1, n_steps=PROP_STEPS)

# Time axis
time_ps = input_pulse.T_ps 

# Calculate new pulse duration
output_power_1 = np.abs(output_1.AT)**2
output_duration_1 = pyutil.fwhm(time_ps, output_power_1) * 1e-12 # (s)

# Calculate new radius and effective nonlinearity using peak intensity from literature
output_energy_1 = output_1.calc_epp() # (J)
output_radius_1 = pyutil.calc_radius(output_energy_1, output_duration_1, PEAK_INTENSITY_2)
print("Beam radius after first Si plate: " + str(output_radius_1) + " (m)")
gamma_2 = pyutil.eff_nonlin(KERR_COEFF_YAG, CENTER_WL * 1e-9, output_radius_1)
print("YAG eff nonlinearity: " + str(gamma_2) + " (1/Wm)")

# Create Second Fiber
fiber_2 = pynlo.media.fibers.fiber.FiberInstance()
fiber_2.generate_fiber(FIBER_LENGTH_2, center_wl_nm=CENTER_WL, betas=(BETA2_YAG, BETA3_YAG, BETA4_YAG),
                        gamma_W_m=gamma_2, gvd_units="ps^n/km", gain=-ALPHA)

# Initiate propagation
y2, AW2, AT2, output_2 = evol.propagate(pulse_in=output_1, fiber=fiber_2, n_steps=PROP_STEPS)

# Calculate new pulse duration
output_power_2 = np.abs(output_2.AT)**2
output_duration_2 = pyutil.fwhm(time_ps, output_power_2) * 1e-12 # (s)

# Calculate new radius and effective nonlinearity using peak intensity from literature
output_energy_2 = output_2.calc_epp() # (J)
output_radius_2 = pyutil.calc_radius(output_energy_2, output_duration_2, PEAK_INTENSITY_3)
print("Beam radius after YAG plate: " + str(output_radius_2) + " (m)")
gamma_3 = pyutil.eff_nonlin(KERR_COEFF_Si, CENTER_WL * 1e-9, output_radius_2)
print("Si (2) eff nonlinearity: " + str(gamma_3) + " (1/Wm)")

# Create Third Fiber
fiber_3 = pynlo.media.fibers.fiber.FiberInstance()
fiber_3.generate_fiber(FIBER_LENGTH_3, center_wl_nm=CENTER_WL, betas=(BETA2_Si, BETA3_Si, BETA4_Si),
                        gamma_W_m=gamma_3, gvd_units="ps^n/km", gain=-ALPHA)

# Initiate propagation
y3, AW3, AT3, output_3 = evol.propagate(pulse_in=output_2, fiber=fiber_3, n_steps=PROP_STEPS)

# Calculate new pulse duration
output_power_3 = np.abs(output_3.AT)**2
output_duration_3 = pyutil.fwhm(time_ps, output_power_3) * 1e-12 # (s)

# Calculate new radius and effective nonlinearity using peak intensity from literature
output_energy_3 = output_3.calc_epp() # (J)
output_radius_3 = pyutil.calc_radius(output_energy_3, output_duration_3, PEAK_INTENSITY_4)
print("Beam radius after Si (2) plate: " + str(output_radius_3) + " (m)")
gamma_4 = pyutil.eff_nonlin(KERR_COEFF_CAF2, CENTER_WL * 1e-9, output_radius_3)
print("CaF2 effective nonlinearity: " + str(gamma_4) + " (1/Wm)")

# Create fourth fiber
fiber_4 = pynlo.media.fibers.fiber.FiberInstance()
fiber_4.generate_fiber(FIBER_LENGTH_4, center_wl_nm=CENTER_WL, betas=(BETA2_CAF2, BETA3_CAF2, BETA4_CAF2),
                        gamma_W_m=gamma_4, gvd_units="ps^n/km", gain=-ALPHA)

# Initiate propagation
y4, AW4, AT4, output_4 = evol.propagate(pulse_in=output_3, fiber=fiber_4, n_steps=PROP_STEPS)

# Calculate final pulse duration
output_power_4 = np.abs(output_4.AT)**2
output_duration_4 = pyutil.fwhm(time_ps, output_power_4) # (ps)
print("Final Pulse Durations: " + str(output_duration_4))

# Create Figure
F = input_pulse.W_mks / (2 * np.pi) * 1e-12
zW0 = np.abs(np.transpose(AW1)[:, (F > 0)])
zW = np.abs(np.transpose(AW4)[:, (F > 0)]) # absolute intensity of E, freq domain
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

