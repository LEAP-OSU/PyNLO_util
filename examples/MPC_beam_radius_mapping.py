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

waist = 0.000062
rayleigh = 0.0005
max_z = 0.1
precision = 0.001
z_axis, radius = pyutil.map_beam_radius(waist, rayleigh, max_z, precision, verbose=True)