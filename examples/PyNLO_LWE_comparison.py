import json
import matplotlib.pyplot as plt
import numpy as np
c = 299792458


# Compare Spectrums
with open("./data/lwe_spectrums.json", "r") as json_file:
    lwe_spectrums = json.load(json_file)
with open("./data/pynlo_spectrums_3uj_YAG.json", "r") as json_file:
    pynlo_spectrums_3uj = json.load(json_file)
with open("./data/pynlo_spectrums_149uj_YAG.json", "r") as json_file:
    pynlo_spectrums_149uj = json.load(json_file)
    

plt.figure(figsize=(10,5))
plt.plot( (c / (np.array(pynlo_spectrums_3uj["freq"])[::-1] * 1e12)) * 1e6,
         np.array(pynlo_spectrums_3uj["seed"])[::-1]/np.max(pynlo_spectrums_3uj["seed"]),
         color='black', label="PyNLO seed 3 uJ")
plt.plot( (c / (np.array(lwe_spectrums["lwe_freq"])[::-1])) * 1e6,
         np.array(lwe_spectrums["seed"])[::-1]/np.max(lwe_spectrums["seed"]),
         color="grey", label="LWE seed")
plt.plot( (c / (np.array(pynlo_spectrums_3uj["freq"])[::-1] * 1e12)) * 1e6,
         np.array(pynlo_spectrums_3uj["0 mm"])[::-1]/np.max(pynlo_spectrums_3uj["0 mm"]),
         color='magenta', label="PyNLO, 0 mm, 3 uJ", linestyle="--")
plt.plot( (c / (np.array(lwe_spectrums["lwe_freq"])[::-1])) * 1e6,
         np.array(lwe_spectrums["NL = 2e-20"])[::-1]/np.max(lwe_spectrums["NL = 2e-20"]),
         color="red", label="LWE, NL = 2e-20, 0mm", linestyle="--")
# plt.plot(np.array(lwe_spectrums["lwe_freq"]) * 1e-12, np.array(lwe_spectrums["NL = 4e-20"])/np.max(lwe_spectrums["NL = 4e-20"]), color="blue", label="LWE, NL = 4e-20, 0mm")
# plt.plot(np.array(lwe_spectrums["lwe_freq"]) * 1e-12, np.array(lwe_spectrums["NL = 6e-20"])/np.max(lwe_spectrums["NL = 6e-20"]), color="orange", label="LWE, NL = 6e-20, 0mm")
plt.xlabel("Wavelength (um)")
plt.ylabel("Intensity (arb. units)")
plt.title("Pulse Spectrum")
plt.legend()
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
x_axis_pynlo = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

plt.figure(figsize=(10,5))
plt.plot(x_axis_lwe, lwe_pulse_durations["NL = 2e-20"], color="red", label="LWE, NL = 2e-20")
plt.plot(x_axis_hemmer, lwe_pulse_durations["hemmer"], color="blue", label="Hemmer et al.")
plt.plot(x_axis_pynlo, np.array(pynlo_pulse_durations["pulse durations"]) * 1e3, color="magenta", label="PyNLO")
plt.xlabel("YAG Plate Position (mm)")
plt.ylabel("Pulse Duration (fs)")
plt.title("Pulse Duration vs. YAG Position")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()