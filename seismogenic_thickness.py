import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from src.compute import Data, CoordinateSystem, SpatialArray2D
from termomecanico import termomecanico
from src.setup import data_setup, input_setup
from src.plot import heatmap_map


data = Data(*data_setup(), *input_setup()) 
cs = CoordinateSystem(data.get_cs_data(), 0.2, 1)
x_axis = cs.get_x_axis()
y_axis = cs.get_y_axis()

# Earthquakes CSN
eqs_csn = pd.read_excel('data/earthquakes/CSN_2000_2018_C_SB30.xlsx')
eqs_csn_ISA = eqs_csn[(eqs_csn['ISA'] == True)].copy()
eqs_csn_ISA_OSB = eqs_csn[(eqs_csn['ISA'] == True) & (eqs_csn['OSB'] == True)].copy()
# Earthquakes USGS
#eqs_usgs = pd.read_excel('data/earthquakes/USGS_1900_2018_C_SB30.xlsx')
#eqs_usgs_ISA = eqs_usgs[(eqs_usgs['ISA'] == True)].copy()
#eqs_usgs_ISA_OSB = eqs_usgs[(eqs_usgs['ISA'] == True) & (eqs_usgs['OSB'] == True)].copy()
## Earthquakes CSN & USGS
#eqs = pd.concat([eqs_usgs, eqs_csn], ignore_index=True)
#eqs_ISA = eqs[(eqs['ISA'] == True)].copy()
#eqs_ISA_OSB = eqs[(eqs['ISA'] == True) & (eqs['OSB'] == True)].copy()


#eqs_usgs = pd.read_excel('data/earthquakes/USGS_1900_2018_C_SB30.xlsx')
lat_bins = np.append(y_axis[::-1]-0.1, y_axis[0]+0.1)
lat_labels = y_axis[::-1]
lon_bins = np.append(x_axis-0.1, x_axis[-1]+0.1)
lon_labels = x_axis

def get_seismogenic_thickness(eqs):
    lat_bin = pd.cut(eqs['latitude'].values, lat_bins, labels=lat_labels)
    lon_bin = pd.cut(eqs['longitude'].values, lon_bins, labels=lon_labels)
    eqs['latitude_bin'] = lat_bin
    eqs['longitude_bin'] = lon_bin

    seismogenic_ths = np.zeros((len(x_axis), len(y_axis)))*np.nan
    for lon_idx,lon in enumerate(x_axis):
        for lat_idx,lat in enumerate(y_axis):
            current_bin_eqs = eqs[
                    (eqs['longitude_bin'] == lon) &
                    (eqs['latitude_bin'] == lat)]
            seismogenic_ths[lon_idx, lat_idx] = current_bin_eqs['depth'].max()
    return SpatialArray2D(seismogenic_ths, cs)

seismogenic_thickness_csn = get_seismogenic_thickness(eqs_csn)
seismogenic_thickness_csn_ISA = get_seismogenic_thickness(eqs_csn_ISA) 
seismogenic_thickness_csn_ISA_OSB = get_seismogenic_thickness(eqs_csn_ISA_OSB)
#seismogenic_thickness_usgs = get_seismogenic_thickness(eqs_usgs)
#seismogenic_thickness_usgs_ISA = get_seismogenic_thickness(eqs_usgs_ISA) 
#seismogenic_thickness_usgs_ISA_OSB = get_seismogenic_thickness(eqs_usgs_ISA_OSB)
#seismogenic_thickness = get_seismogenic_thickness(eqs)
#seismogenic_thickness_ISA = get_seismogenic_thickness(eqs_ISA) 
#seismogenic_thickness_ISA_OSB = get_seismogenic_thickness(eqs_ISA_OSB)
#fig = plt.figure(figsize=(12,18))
#gs = gridspec.GridSpec(3,3)
#ax1 = fig.add_subplot(gs[0,0])
#heatmap_map(
#    seismogenic_thickness_csn, ax=ax1,
#    title='Eqs. CSN ({:d})'.format(len(eqs_csn.index)),
#    colormap='viridis', cbar_limits=[0,400])
#ax2 = fig.add_subplot(gs[0,1])
#heatmap_map(
#    seismogenic_thickness_csn_ISA, ax=ax2,
#    title='Eqs. CSN - ISA ({:d})'.format(len(eqs_csn_ISA.index)),
#    colormap='viridis', cbar_limits=[0,400])
#ax3 = fig.add_subplot(gs[0,2])
#heatmap_map(
#    seismogenic_thickness_csn_ISA_OSB, ax=ax3,
#    title='Eqs. CSN - ISA - OSB (30 km.) ({:d})'.format(len(eqs_csn_ISA_OSB.index)),
#    colormap='viridis', cbar_limits=[0,400])
#ax4 = fig.add_subplot(gs[1,0])
#heatmap_map(
#    seismogenic_thickness_usgs, ax=ax4,
#    title='Eqs. USGS ({:d})'.format(len(eqs_usgs.index)),
#    colormap='viridis', cbar_limits=[0,400])
#ax5 = fig.add_subplot(gs[1,1])
#heatmap_map(
#    seismogenic_thickness_usgs_ISA, ax=ax5,
#    title='Eqs. USGS - ISA ({:d})'.format(len(eqs_usgs_ISA.index)),
#    colormap='viridis', cbar_limits=[0,400])
#ax6 = fig.add_subplot(gs[1,2])
#heatmap_map(
#    seismogenic_thickness_usgs_ISA_OSB, ax=ax6,
#    title='Eqs. USGS - ISA - OSB (30 km.) ({:d})'.format(len(eqs_usgs_ISA_OSB.index)),
#    colormap='viridis', cbar_limits=[0,400])
#ax7 = fig.add_subplot(gs[2,0])
#heatmap_map(
#    seismogenic_thickness, ax=ax7,
#    title='Eqs. USGS & CSN ({:d})'.format(len(eqs.index)),
#    colormap='viridis', cbar_limits=[0,400])
#ax8 = fig.add_subplot(gs[2,1])
#heatmap_map(
#    seismogenic_thickness_ISA, ax=ax8,
#    title='Eqs. USGS & CSN - ISA ({:d})'.format(len(eqs_ISA.index)),
#    colormap='viridis', cbar_limits=[0,400])
#ax9 = fig.add_subplot(gs[2,2])
#heatmap_map(
#    seismogenic_thickness_ISA_OSB, ax=ax9,
#    title='Eqs. USGS & CSN - ISA - OSB (30 km.) ({:d})'.format(len(eqs_ISA_OSB.index)),
#    colormap='viridis', cbar_limits=[0,400])
model, _, _ = termomecanico(*input_setup())
s1 = seismogenic_thickness_csn[~np.isnan(seismogenic_thickness_csn)]
e1 = model.mm.get_eet()[~np.isnan(seismogenic_thickness_csn)]
s2 = seismogenic_thickness_csn_ISA[~np.isnan(seismogenic_thickness_csn_ISA)]
e2 = model.mm.get_eet()[~np.isnan(seismogenic_thickness_csn_ISA)]
s3 = seismogenic_thickness_csn_ISA_OSB[~np.isnan(seismogenic_thickness_csn_ISA_OSB)]
e3 = model.mm.get_eet()[~np.isnan(seismogenic_thickness_csn_ISA_OSB)]
fig = plt.figure(figsize=(12,4))
gs = gridspec.GridSpec(1,3)
ax1 = fig.add_subplot(gs[0,0])
ax1.scatter(s1,e1)
ax2 = fig.add_subplot(gs[0,1])
ax2.scatter(s2,e2)
ax3 = fig.add_subplot(gs[0,2])
ax3.scatter(s3,e3)
plt.tight_layout()
plt.show()


#print('d')
