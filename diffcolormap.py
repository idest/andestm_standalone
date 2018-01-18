import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --- Colormaps from a list ---

colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # R -> G -> B
n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
cmap_name = 'diff'


cdict = {'red':   ((0.0, 0.00, 0.00),
                   (0.2, 0.00, 0.50),
                   (0.4, 0.50, 1.00),
                   (0.6, 1.00, 1.00),
                   (0.8, 1.00, 1.00),
                   (1.0, 1.00, 0.00)),

         'green': ((0.0, 0.0, 0.0),
                   (0.2, 0.0, 0.5),
                   (0.4, 0.5, 1.0),
                   (0.6, 1.0, 0.5),
                   (0.8, 0.5, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.00, 1.00),
                   (0.2, 1.00, 1.00),
                   (0.4, 1.00, 1.00),
                   (0.6, 1.00, 0.50),
                   (0.8, 0.50, 0.00),
                   (1.0, 0.00, 0.00)),

        'alpha':  ((0.0, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (1.0, 1.0, 1.0))
         }

diff_cmap = LinearSegmentedColormap(cmap_name, cdict)
