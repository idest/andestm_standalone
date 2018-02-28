import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb

# --- Colormap from a list ---

def get_diff_cmap(bins):
    colors = [(1, 0, 0), (1, 1, 1), (0, 0, 1)]  # Red -> White -> Blue
    n_bins = bins  # Discretizes the interpolation into bins
    cmap_name = 'diff'
    diff_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
    return diff_cmap

# --- Custom colormap

custom_cmap_name = 'diff_custom'

cdict = {'red':   ((0.0, 1.00, 1.00),
                   (0.1, 1.00, 1.00),
                   (0.2, 1.00, 1.00),
                   (0.3, 1.00, 1.00),
                   (0.4, 1.00, 1.00),
                   (0.6, 1.00, 0.75),
                   (0.7, 0.75, 0.50),
                   (0.8, 0.50, 0.25),
                   (0.9, 0.25, 0.00),
                   (1.0, 0.00, 0.00)),

         'green': ((0.0, 0.00, 0.00),
                   (0.1, 0.00, 0.25),
                   (0.2, 0.25, 0.50),
                   (0.3, 0.50, 0.75),
                   (0.4, 0.75, 1.00),
                   (0.6, 1.00, 0.75),
                   (0.7, 0.75, 0.50),
                   (0.8, 0.50, 0.25),
                   (0.9, 0.25, 0.00),
                   (1.0, 0.00, 0.00)),

         'blue':  ((0.0, 0.00, 0.00),
                   (0.1, 0.00, 0.25),
                   (0.2, 0.25, 0.50),
                   (0.3, 0.50, 0.75),
                   (0.4, 0.75, 1.00),
                   (0.6, 1.00, 1.00),
                   (0.7, 1.00, 1.00),
                   (0.8, 1.00, 1.00),
                   (0.9, 1.00, 1.00),
                   (1.0, 1.00, 1.00)),

        'alpha':  ((0.0, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (1.0, 1.0, 1.0))
         }

diff_cmap_custom = LinearSegmentedColormap(custom_cmap_name, cdict)
