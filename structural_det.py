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
def mapsmdet(lon, lat, depth):
    # lon, lat, depth = np.meshgrid(lon, lat, depth, sparse = True)
    xloc = np.arange(-76.0, -65+3.0, 3.0)
    yloc = np.arange(-45, -32+2.0, 2.0)
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    coastline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
    border = cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m')
    faults = cfeature.ShapelyFeature(shpreader.Reader('data/shp/FALLAS/FALLAS.shp').geometries(),
                                     ccrs.PlateCarree())
    lofz = cfeature.ShapelyFeature(shpreader.Reader('data/shp/FALLAS/lofz.shp').geometries(),
                                   ccrs.PlateCarree())
    plt.title('SM Detachment')
    ax.add_feature(border, facecolor='None', edgecolor='gray', alpha=0.7)
    ax.add_feature(coastline, facecolor='None', edgecolor='k')
    f = ax.add_feature(faults, facecolor='None', edgecolor='k', linewidth=0.5)
    l = ax.add_feature(lofz, facecolor='None', edgecolor='b', linewidth=0.5)
    # det = ax.contourf(lon, lat, depth, transform=ccrs.PlateCarree(), cmap='rainbow', vmin=-35, vmax=5)
    # det = ax.imshow(depth, origin='right', transform=ccrs.PlateCarree(),
    #                 extent=(np.amin(lon), np.amax(lon), np.amin(lat), np.amax(lat)),
    #                 cmap='rainbow', interpolation='bicubic')
    det = ax.scatter(lon, lat, c=depth, transform=ccrs.PlateCarree(), s=0.5, cmap='rainbow',
                     vmin=-35.0, vmax=-5.0)
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
