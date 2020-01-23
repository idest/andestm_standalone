import numpy as np
import numpy.ma as ma
import pandas as pd
import geopandas as gpd
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
    icd = gm.get_icd()
    depthc = depth.crop(top_z=gm.get_topo(), bottom_z=gm.get_moho())
    det = np.empty(np.shape(depthc))
    index = np.where(np.logical_and(ysed<=smax, ysed>=smin))
    det[index] = depthc[index]
    det = np.nanmax(det, axis=2)
    # nans, k = nan_helper(det) #Selecting nan values
    # i, j = np.indices(det.shape) #Indices of det array
    # det[nans] = griddata((i[~nans], j[~nans]), det[~nans], (i[nans], j[nans]), method='cubic') #Interpolate nan values
    # det = np.ma.masked_where(det >= -1, det)
    detdat = pd.DataFrame(det)
    return det, detdat


# Mapping detachment
def maptmdet(lon, lat, det, csn_e, usgs_e, icd, topo, labslab, moho, map):
    mdet = ma.masked_invalid(det)
    xloc = np.arange(-76.0, -65+3.0, 3.0)
    yloc = np.arange(-44, -33+2.0, 2.0)
    usgs_e = usgs_e[usgs_e['depth'] >= -15.0]
    csn_e = csn_e[csn_e['depth'] >= -15.0]
    uex = pd.DataFrame(usgs_e.loc[:, 'longitude']).to_numpy()
    uey = pd.DataFrame(usgs_e.loc[:, 'latitude']).to_numpy()
    cex = pd.DataFrame(csn_e.loc[:, 'longitude']).to_numpy()
    cey = pd.DataFrame(csn_e.loc[:, 'latitude']).to_numpy()
    mlon, mlat = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    coastline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
    border = cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m')
    rfaults = cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/reverse/reverse.shp').geometries(),
                                      ccrs.PlateCarree())
    sfaults = cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/sinestral/sinestral.shp').geometries(),
                                      ccrs.PlateCarree())
    nfaults = cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/normal/normal.shp').geometries(),
                                      ccrs.PlateCarree())
    dfaults =  cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/dextral/dextral.shp').geometries(),
                                      ccrs.PlateCarree())
    ufaults = cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/undefined/undefined.shp').geometries(),
                                      ccrs.PlateCarree())
    # tz = cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/trust_zone/trust_zone.shp').geometries(),
    #                                   ccrs.PlateCarree())
    plt.title('TM Rigid-Ductile Transition')
    ax.add_feature(border, facecolor='None', edgecolor='gray', alpha=0.7)
    ax.add_feature(coastline, facecolor='None', edgecolor='gray')
    if map == 'faults':
        reverse = ax.add_feature(rfaults, facecolor='None', edgecolor='k', linewidth=1.0)
        normal = ax.add_feature(nfaults, facecolor='None', edgecolor='m', linewidth=1.0)
        dextral = ax.add_feature(dfaults, facecolor='None', edgecolor='r', linewidth=1.0)
        sinestral = ax.add_feature(sfaults, facecolor='None', edgecolor='g', linewidth=1.0)
        undefined = ax.add_feature(ufaults, facecolor='None', linestyle='--', edgecolor='k', linewidth=1.0)
    elif map == 'earthquakes':
        ax.scatter(uex, uey, transform=ccrs.PlateCarree(), s=0.5, color='k', alpha=0.5)
    elif map == 'moho':
        cs = ax.contour(mlon, mlat, moho.T, levels=15, colors='black', transform=ccrs.PlateCarree(),
                        linewidths=1.0, linestyles='solid')
        ax.clabel(cs, cs.levels, inline_spacing=1, fontsize=8, manual=True, fmt='%1.1f',
                  use_clabeltext=True)
    elif map == 'icd':
        cs = ax.contour(mlon, mlat, icd.T, levels=8, colors='black', transform=ccrs.PlateCarree(),
                        linewidths=1.0, linestyles='solid')
        ax.clabel(cs, cs.levels, inline_spacing=1, fontsize=8, manual=True, fmt='%1.1f',
                  use_clabeltext=True)
    elif map == 'labslab':
        cs = ax.contour(mlon, mlat, labslab.T, levels=15, colors='black', transform=ccrs.PlateCarree(),
                        linewidths=1.0, linestyles='solid')
        ax.clabel(cs, cs.levels, inline_spacing=1, fontsize=8, manual=True, fmt='%1.1f',
                  use_clabeltext=True)
    cmap = plt.cm.RdYlBu_r
    cmap.set_bad('w',1.0)
    # det = ax.imshow(det.T, origin='upper', transform=ccrs.PlateCarree(),
    #                 extent=(np.amin(mlon), np.amax(lon), np.amin(lat), np.amax(lat)),
    #                 cmap=cmap, interpolation='bicubic', vmin = -35.0, vmax = 0.)
    det = ax.contourf(mlon, mlat, det.T, 56, transform=ccrs.PlateCarree(),
                      cmap=cmap, vmin=-35., vmax=0.)
    #fix fow white borders in contourf
    for c in det.collections:
        c.set_edgecolor("face")
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(det)
    m.set_clim(-35., 0.)
    cbar = plt.colorbar(m, ticks=np.arange(-35., 1., 5.),
                        boundaries=np.linspace(-35., 0, 28.))
    cbar.set_ticklabels(np.arange(-35., 1., 5.))
    # cbar = plt.colorbar(det)
    cbar.set_label('Depth [Km]', rotation=90)
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.8)
    gl.xlabels_top=False
    gl.ylabels_right=False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(xloc)
    gl.ylocator = mticker.FixedLocator(yloc)
    ax.set_extent([-75, -68, -44, -33])
    return fig
