import sys
sys.path.insert(0, 'src/')
import numpy as np
import pandas as pd
from setup import data_setup, input_setup, exec_setup
from compute import compute
from scipy.interpolate import RegularGridInterpolator
from tm_detachment import mdetachment, maptmdet
from structural_det import sdetachment, mapsmdet
from det_stats import map_diff, plot_ccr, map_int, plot_hist
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
from scipy.stats import pearsonr


# NaN helper function
def nan_helper(z):
    return np.isnan(z), lambda w: w.nonzero()[0]


# TM MODEL
# Load data and TM model
smax, smin = 100, 50 # MPa
smoother = -20 # Smooth nan values
gm_data, areas, trench_age, rhe_data, coast = data_setup()
t_input, m_input = input_setup()


def termomecanico(t_input, m_input):
    model = compute(gm_data, areas, trench_age, rhe_data, coast, t_input, m_input)
    cs = model.mm.get_coordinate_system()
    gm = model.mm.get_geometric_model()
    grid = cs.get_3D_grid()[2]
    crust = grid.crop(top_z=gm.get_topo(), bottom_z=gm.get_moho())
    yset, ysec = model.mm.get_yse()
    ysed = model.mm.get_ductile_yield_strength()
    # ysed = np.ma.masked_where(np.isnan(crust), ysed)
    return cs, ysec, gm, ysed


# Running TM model and obtaining detachment
cs, ysec, gm, ysed = termomecanico(t_input, m_input)
lon = cs.get_x_axis()
lat = cs.get_y_axis()
det, detdat = mdetachment(cs, ysed, gm, smax, smin)
export = detdat.to_csv(r'/home/julvelillo/Documentos/detdat.csv')


# Mapping detachment
mmap = maptmdet(lon, lat, det)
# mmap.figure.savefig('det_results/tm_map{}_{}.eps'.format(smax, smin),
#                      format='eps', dpi=500)
mmap.savefig('det_results/tm_map{}_{}'.format(smax, smin))
plt.close()


# SM MODEL
# Runing SM model and obtaining detachment
slon, slat, depth = sdetachment()


# Mapping detachment
smap = mapsmdet(slon, slat, depth)
# smap.figure.savefig('det_results/sm_map.eps', format='eps', dpi=500)
smap.savefig('det_results/sm_map')
plt.close()


# STATISTICS
# Interpolating tm model to sm model points
def interpolate_det(z, x, y, i, j):
    # z = np.ma.masked_invalid(z)
    # z[np.isnan(z)] = smoother
    # df = pd.DataFrame(z)
    # export = df.to_csv(r'det_results/z.csv')
    # interpolate = RectBivariateSpline(x, y[::-1], z[:, ::-1])
    # idet = interpolate.ev(i, j)
    # nans, k = nan_helper(z)
    interpolate = RegularGridInterpolator((x, y[::-1]), z[:, ::-1], method='linear')
    points = i, j
    idet = interpolate(points)
    return idet


idet = interpolate_det(det, lon, lat, slon, slat)


# Mapping Interpolation
intd = map_int(idet, slon, slat)
# intd.figure.savefig('int_map{}_{}.eps'.format(smax, smin), format='eps', dpi=500)
intd.savefig('det_results/int_map{}_{}'.format(smax, smin))
plt.close()

# Correlation Coefficient
diff = idet - depth
ccr, p = pearsonr(idet, depth)


# Plotting Correlation Coefficient
plotccr = plot_ccr(idet, depth, ccr)
# plt.savefig('det_results/ccr{}_{}.eps'.format(smax, smin), format='eps', dpi=500)
plotccr.savefig('det_results/ccr{}_{}'.format(smax, smin))
plt.close()


# Mapping Difference
mapdiff = map_diff(diff, slon, slat)
# mapdiff.figure.savefig('det_results/diff_map{}_{}.eps'.format(smax, smin),
#                         format='eps', dpi=500)
mapdiff.savefig('det_results/diff_map{}_{}'.format(smax, smin))
plt.close()


#Plotting hist and PDF
hist = plot_hist(diff)
# hist.savefig('det_results/hist{}_{}.eps'.format(smax, smin), format=eps, dpi=500)
hist.savefig('det_results/hist{}_{}'.format(smax, smin))
plt.close()
