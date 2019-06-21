import sys
sys.path.insert(0, 'src/')
import numpy as np
import numpy.ma as ma
import pandas as pd
from setup import data_setup, input_setup, exec_setup
from compute import compute
from src.stats import evaluate_model
from src.datos_q import shf_data
from src.colormaps import (jet_white_r, get_diff_cmap,
                           get_elevation_diff_cmap, jet_white_r_2)
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


#tm model data
gm_data, areas, trench_age, rhe_data, coast = data_setup()
t_input, m_input = input_setup()


#running tm model
def termomecanico(t_input, m_input):
    model = compute(gm_data, areas, trench_age, rhe_data, coast, t_input, m_input)
    cs = model.mm.get_coordinate_system()
    gm = model.mm.get_geometric_model()
    return cs, gm, model


#surface heat flow, tm output and estimators
cs, gm, model = termomecanico(t_input, m_input)
shf = model.tm.get_surface_heat_flow(format='positive milliwatts')
lon = cs.get_x_axis()
lat = cs.get_y_axis()
estimators, df = evaluate_model(shf, shf_data, return_dataframe=True)
diff = df['diffs']


#multimap
def multimap(lon, lat, shf, diff):
        #laoding data and variables
        xloc = np.arange(-80, -60+5.0, 5.0)
        yloc = np.arange(-45, -10+5.0, 5.0)
        mlon, mlat = np.meshgrid(lon, lat)
        diff_step = 5
        divisions = np.arange(-80, 80+diff_step, diff_step)
        ticks = np.arange(-80, 80+diff_step, 10)
        bins = len(divisions) - 1
        diff_cmap = get_diff_cmap(bins)
        norm = MidPointNorm(midpoint=0, vmin=-80, vmax=80)
        #plotting figure
        fig = plt.figure()
        heat_flow = shf
        hf_masked = np.ma.masked_invalid(heat_flow)
        ax = plt.axes(projection=ccrs.PlateCarree())
        coastline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
        border = cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m')
        plt.title('Thermal model and residual values')
        ax.background_patch.set_fill(False)
        ax.add_feature(border, facecolor='None', edgecolor='gray', alpha=0.7)
        ax.add_feature(coastline, facecolor='None', edgecolor='k')
        cmap = plt.cm.afmhot
        hfcf = ax.contourf(mlon, mlat, heat_flow.T, 24, transform=ccrs.PlateCarree(),
                           cmap=cmap, vmin=0, vmax=120)
        #fix the white borders at each contourf
        for c in hfcf.collections:
            c.set_edgecolor("face")
        cbar = plt.colorbar(hfcf)
        cbar.set_label('Heat Flow [W/m²]', rotation=90)
        gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.3)
        gl.xlabels_top=False
        gl.ylabels_right=False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlocator = mticker.FixedLocator(xloc)
        gl.ylocator = mticker.FixedLocator(yloc)
        ax.set_extent([-80, -60, -45, -10])
        return fig


hfmap = multimap(lon, lat, shf, diff)
hfmap.savefig('map_final/multi_map.pdf', transparent='True')
plt.close()
