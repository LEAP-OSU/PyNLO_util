import json
import matplotlib.pyplot as plt
import numpy as np
c = 299792458


# Compare Spectrums
with open("./data/lwe_spectrums.json", "r") as json_file:
    lwe_spectrums = json.load(json_file)
with open("./data/pynlo_spectrums_3uj_YAG.json", "r") as json_file:
    pynlo_spectrums_3uj = json.load(json_file)
# with open("./data/pynlo_spectrums_149uj_YAG.json", "r") as json_file:
#     pynlo_spectrums_149uj = json.load(json_file)
    

plt.figure(figsize=(10,5))
plt.plot( (c / (np.array(pynlo_spectrums_3uj["freq"])[::-1] * 1e12)) * 1e6,
         np.array(pynlo_spectrums_3uj["seed"])[::-1]/np.max(pynlo_spectrums_3uj["seed"]),
         color='black', label="PyNLO seed", linewidth=2)
plt.plot( (c / (np.array(lwe_spectrums["lwe_freq"])[::-1])) * 1e6,
         np.array(lwe_spectrums["seed"])[::-1]/np.max(lwe_spectrums["seed"]),
         color="orange", label="LWE seed", linewidth=2)
plt.plot( (c / (np.array(pynlo_spectrums_3uj["freq"])[::-1] * 1e12)) * 1e6,
         np.array(pynlo_spectrums_3uj["0 mm"])[::-1]/np.max(pynlo_spectrums_3uj["0 mm"]),
         color='blue', label="PyNLO, at focus", linewidth=2)
plt.plot( (c / (np.array(lwe_spectrums["lwe_freq"])[::-1])) * 1e6,
         np.array(lwe_spectrums["NL = 2e-20"])[::-1]/np.max(lwe_spectrums["NL = 2e-20"]),
         color="red", label="LWE, at focus", linewidth=2)
# plt.plot(np.array(lwe_spectrums["lwe_freq"]) * 1e-12, np.array(lwe_spectrums["NL = 4e-20"])/np.max(lwe_spectrums["NL = 4e-20"]), color="blue", label="LWE, NL = 4e-20, 0mm")
# plt.plot(np.array(lwe_spectrums["lwe_freq"]) * 1e-12, np.array(lwe_spectrums["NL = 6e-20"])/np.max(lwe_spectrums["NL = 6e-20"]), color="orange", label="LWE, NL = 6e-20, 0mm")
plt.xlabel("Wavelength (um)", fontsize=16)
plt.ylabel("Intensity (arb. units)", fontsize=16)
plt.title("Pulse Spectrum", fontsize=16)
plt.legend(fontsize=16)
plt.grid(True, alpha=0.3)
plt.show()


# Compare Pulse Durations
with open("./data/lwe_pulse_durations.json", "r") as json_file:
    lwe_pulse_durations = json.load(json_file)
with open("./data/pynlo_pulse_durations_3uj_YAG.json", "r") as json_file:
    pynlo_pulse_durations = json.load(json_file)

x_axis_lwe = [-6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 
              0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
x_axis_hemmer = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
x_axis_pynlo = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

plt.figure(figsize=(10,5))
plt.plot(x_axis_lwe, lwe_pulse_durations["NL = 2e-20"], color="red", label="LWE", linewidth=2)
plt.plot(x_axis_hemmer, lwe_pulse_durations["hemmer"], color="blue", label="Hemmer et al.", linewidth=2)
plt.plot(x_axis_pynlo, np.array(pynlo_pulse_durations["pulse durations"]) * 1e3, color="magenta", label="PyNLO", linewidth=2)
plt.xlabel("YAG Plate Position (mm)", fontsize=16)
plt.ylabel("Pulse Duration (fs)", fontsize=16)
plt.title("Pulse Duration vs. YAG Position", fontsize=16)
plt.legend(fontsize=16)
plt.grid(True, alpha=0.3)
plt.show()