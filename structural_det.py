import pandas as pd
from pandas import DataFrame
from pyproj import Proj
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def sdetachment():
    # Loading data
    data = pd.read_csv('data/structural_detachment.csv', header=None)
    depth = pd.DataFrame(data.loc[:,2]/1000).to_numpy() #To km depth
    x = pd.DataFrame(data.loc[:,0]).to_numpy()
    y = pd.DataFrame(data.loc[:,1]).to_numpy()
    # depth = depth.reshape(depth.shape[0], y.shape)


    # Transforming x,y to lon,lat
    myProj = Proj(init='epsg:24878', proj='utm', zone=18, south=True)
    lon, lat = myProj(x, y, inverse=True)
    return lon, lat, depth


# Creating map
def mapsmdet(lon, lat, slon, slat, depth, icd, topo, labslab, moho, map):
    # lon, lat, depth = np.meshgrid(lon, lat, depth, sparse = True)
    xloc = np.arange(-76.0, -65+3.0, 3.0)
    yloc = np.arange(-44, -33+2.0, 2.0)
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
    tz = cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/trust_zone/trust_zone.shp').geometries(),
                                      ccrs.PlateCarree())
    plt.title('SM Detachment')
    ax.add_feature(border, facecolor='None', edgecolor='gray', alpha=0.7)
    ax.add_feature(coastline, facecolor='None', edgecolor='k')
    if map == 'faults':
        reverse = ax.add_feature(rfaults, facecolor='None', edgecolor='k', linewidth=1.0)
        normal = ax.add_feature(nfaults, facecolor='None', edgecolor='m', linewidth=1.0)
        dextral = ax.add_feature(dfaults, facecolor='None', edgecolor='r', linewidth=1.0)
        sinestral = ax.add_feature(sfaults, facecolor='None', edgecolor='g', linewidth=1.0)
        undefined = ax.add_feature(ufaults, facecolor='None', linestyle='--', edgecolor='k', linewidth=1.0)
        ax.add_feature(tz, linestyle='--', facecolor='None', edgecolor='k', linewidth=1.0, alpha=0.5)
    elif map == 'moho':
        cs = ax.contour(mlon, mlat, moho.T, levels=15, colors='black', transform=ccrs.PlateCarree(),
                        linewidths=1.0, linestyles='solid')
        ax.clabel(cs, cs.levels, inline_spacing=1, fontsize=8, manual=True, fmt='%1.1f',
                  use_clabeltext=True)
    elif map == 'icd':
        cs = ax.contour(mlon, mlat, icd.T, levels=6, colors='black', transform=ccrs.PlateCarree(),
                        linewidths=1.0, linestyles='solid')
        ax.clabel(cs, cs.levels, inline_spacing=1, fontsize=8, manual=True, fmt='%1.1f',
                  use_clabeltext=True)
    elif map == 'labslab':
        cs = ax.contour(mlon, mlat, labslab.T, levels=15, colors='black', transform=ccrs.PlateCarree(),
                        linewidths=1.0, linestyles='solid')
        ax.clabel(cs, cs.levels, inline_spacing=1, fontsize=8, manual=True, fmt='%1.1f',
                  use_clabeltext=True)
    det = ax.scatter(slon, slat, c=depth, transform=ccrs.PlateCarree(), s=0.5, marker='s', cmap='RdYlBu_r',
                     vmin=-35.0, vmax=0.)
    cbar = plt.colorbar(det)
    cbar.set_label('Depth [Km]', rotation=90)
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.8)
    ax.background_patch.set_facecolor('silver')
    gl.xlabels_top=False
    gl.ylabels_right=False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(xloc)
    gl.ylocator = mticker.FixedLocator(yloc)
    ax.set_extent([-75, -68, -44, -33])
    return fig
