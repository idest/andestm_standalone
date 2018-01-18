import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --- Colormaps from a list ---

colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # R -> G -> B
n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
cmap_name = 'diff'


cdict = {'red':   ((0.0, 0.0, 0.0),
                  #(0.25, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                  #(0.75, 1.0, 1.0),
                   (1.0, 0.0, 0.0)),

         'green': ((0.0, 0.0, 0.0),
                  #(0.25, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                  #(0.75, 0.0, 0.0),
                   (1.0, 1.0, 1.0)),

         'blue':  ((0.0, 1.0, 1.0),
                  #(0.25, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                  #(0.75, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

        'alpha':  ((0.0, 1.0, 1.0),
                  #(0.25, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                  #(0.75, 1.0, 1.0),
                   (1.0, 1.0, 1.0))
         }

diff_cmap = LinearSegmentedColormap(cmap_name, cdict)
