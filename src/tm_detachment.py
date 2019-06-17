import numpy as np
import numpy.ma as ma
import pandas as pd
from setup import data_setup, input_setup, exec_setup
from compute import compute
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# NaN helper function
def nan_helper(z):
    return np.isnan(z), lambda w: w.nonzero()[0]


# Obtaining detachment
def mdetachment(cs, ysed, gm, smax, smin):
    depth = cs.get_3D_grid()[2]
    depthc = depth.crop(top_z=gm.get_topo(), bottom_z=gm.get_moho())
    det = np.empty(np.shape(depthc))
    index = np.where(np.logical_and(ysed<=smax, ysed>=smin))
    det[index] = depthc[index]
    det = np.nanmax(det, axis=2)
    # nans, k = nan_helper(det) #Selecting nan values
    # i, j = np.indices(det.shape) #Indices of det array
    # det[nans] = griddata((i[~nans], j[~nans]), det[~nans], (i[nans], j[nans]), method='cubic') #Interpolate nan values
    detdat = pd.DataFrame(det)
    return det, detdat


# Mapping detachment
def maptmdet(lon, lat, det):
    # mdet = ma.masked_invalid(det)
    xloc = np.arange(-76.0, -65+3.0, 3.0)
    yloc = np.arange(-45, -32+2.0, 2.0)
    mlon, mlat = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    coastline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
    border = cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m')
    faults = cfeature.ShapelyFeature(shpreader.Reader('data/shp/FALLAS/FALLAS.shp').geometries(),
                                     ccrs.PlateCarree())
    lofz = cfeature.ShapelyFeature(shpreader.Reader('data/shp/FALLAS/lofz.shp').geometries(),
                                   ccrs.PlateCarree())
    plt.title('TM Detachment')
    ax.add_feature(border, facecolor='None', edgecolor='gray', alpha=0.7)
    ax.add_feature(coastline, facecolor='None', edgecolor='k')
    ax.add_feature(faults, facecolor='None', edgecolor='k', linewidth=0.5)
    ax.add_feature(lofz, facecolor='None', edgecolor='b', linewidth=0.5)
    # det = ax.contourf(mlon, mlat, det.T, transform=ccrs.PlateCarree(), cmap='rainbow', vmin=-35, vmax=5)
    cmap = plt.cm.rainbow
    cmap.set_bad('w',1.0)
    det = ax.imshow(det.T, origin='upper', transform=ccrs.PlateCarree(),
                    extent=(np.amin(mlon), np.amax(lon), np.amin(lat), np.amax(lat)),
                    cmap='rainbow', interpolation='bicubic', vmin = -35.0, vmax = -5.0)
    # det = ax.scatter(mlon, mlat, c=det.T, transform=ccrs.PlateCarree(),
    #                  s=20, marker='o', cmap='rainbow')
    cbar = plt.colorbar(det)
    cbar.set_label('Depth [Km]', rotation=90)
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.8)
    gl.xlabels_top=False
    gl.ylabels_right=False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(xloc)
    gl.ylocator = mticker.FixedLocator(yloc)
    ax.set_extent([-75, -68, -45, -32])
    return fig
