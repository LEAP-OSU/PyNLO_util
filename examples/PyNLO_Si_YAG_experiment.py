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
BETA2_Si = -441.71         # (ps^2/km)
BETA3_Si = 0.85            # (ps^3/km)
BETA4_Si = -0.01           # (ps^4/km)
GAMMA_1 = pyutil.eff_nonlin(KERR_COEFF_Si, CENTER_WL * 1e-9, 2923.85 * 1e-6)   # (1/Wm)
print("Si (1) eff nonlinearity: " + str(GAMMA_1) + " (1/Wm)")
RAMAN = True
STEEP = True

# Second Fiber Parameters (YAG, 2mm)
####################################
FIBER_LENGTH_2 = 0.002    # (m)
KERR_COEFF_YAG = 0.5e-19     # (m^2/W)
BETA2_YAG = -408.15         # (ps^2/km)
BETA3_YAG = 2.54            # (ps^3/km)
BETA4_YAG = -0.02
PEAK_INTENSITY_2 = 1000 * 1e13  #(W/m^2)
# Gamma calculated after first propagation

# Third Fiber Parameters (Si, 1mm)
##################################
FIBER_LENGTH_3 = 0.001      # (m)
PEAK_INTENSITY_3 = 50 * 1e13    #(W/m^2)
# Gamma calculated after second propagation


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
y, AW1, AT, pulse_2 = evol.propagate(pulse_in=input_pulse, fiber=fiber_1, n_steps=PROP_STEPS)

# Calculate YAG gamma value
###########################
time_ps = input_pulse.T_ps
power_2 = np.abs(pulse_2.AT)**2
pulse_duration_2 = pyutil.fwhm(time_ps, power_2) * 1e-12 # (s)
pulse_energy_2 = pulse_2.calc_epp() # (J)
radius_2 = pyutil.calc_radius(pulse_energy_2, pulse_duration_2, PEAK_INTENSITY_2)
print("Beam radius after first Si plate: " + str(radius_2) + " (m)")
gamma_2 = pyutil.eff_nonlin(KERR_COEFF_YAG, CENTER_WL * 1e-9, radius_2)
print("YAG eff nonlinearity: " + str(gamma_2) + " (1/Wm)")

# Create Second Fiber
#####################
fiber_2 = pynlo.media.fibers.fiber.FiberInstance()
fiber_2.generate_fiber(FIBER_LENGTH_2, center_wl_nm=CENTER_WL, betas=(BETA2_YAG, BETA3_YAG, BETA4_YAG),
                        gamma_W_m=gamma_2, gvd_units="ps^n/km", gain=-ALPHA)

# Initiate propagation
######################
y, AW, AT, pulse_3 = evol.propagate(pulse_in=pulse_2, fiber=fiber_2, n_steps=PROP_STEPS)

# Calculate Second Si gamma value
#################################
power_3 = np.abs(pulse_3.AT)**2
pulse_duration_3 = pyutil.fwhm(time_ps, power_3) * 1e-12 # (s)
pulse_energy_3 = pulse_3.calc_epp() # (J)
radius_3 = pyutil.calc_radius(pulse_energy_3, pulse_duration_3, PEAK_INTENSITY_3)
print("Beam radius after YAG plate: " + str(radius_3) + " (m)")
gamma_3 = pyutil.eff_nonlin(KERR_COEFF_Si, CENTER_WL * 1e-9, radius_3)
print("Si (2) eff nonlinearity: " + str(gamma_3) + " (1/Wm)")

# Create Third Fiber
####################
fiber_3 = pynlo.media.fibers.fiber.FiberInstance()
fiber_3.generate_fiber(FIBER_LENGTH_3, center_wl_nm=CENTER_WL, betas=(BETA2_Si, BETA3_Si, BETA4_Si),
                        gamma_W_m=gamma_3, gvd_units="ps^n/km", gain=-ALPHA)

# Initiate propagation
######################
y, AW, AT, pulse_4 = evol.propagate(pulse_in=pulse_3, fiber=fiber_3, n_steps=PROP_STEPS)

# Calculate final pulse duration 
#############################################
power_4 = np.abs(pulse_4.AT)**2
pulse_duration_4 = pyutil.fwhm(time_ps, power_4) # (ps)
print("final pulse duration: " + str(pulse_duration_4) + " (ps)")



# Create Figure
F = input_pulse.W_mks / (2 * np.pi) * 1e-12
zW0 = np.abs(np.transpose(AW1)[:, (F > 0)])
zW = np.abs(np.transpose(AW)[:, (F > 0)]) # absolute intensity of E, freq domain
zT = np.abs(np.transpose(AT)) # absolute intensity of E, time domain

plt.figure(figsize=(10,5))
plt.plot(F[F > 0], zW[-1]/np.max(zW[-1]), color='r', label="Output Spectrum")
plt.plot(F[F > 0], zW0[0]/np.max(zW0[0]), color='b', label="Input Spectrum")
plt.xlabel("Freq (THz)")
plt.ylabel("Intensity (a.u.)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

