import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats import norm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# Mapping diff
def map_diff(diff, lon, lat):
    # mdiff = ma.masked_invalid(diff)
    xloc = np.arange(-76.0, -65+3.0, 3.0)
    yloc = np.arange(-44, -33+2.0, 2.0)
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
    plt.title('Residual')
    ax.add_feature(border, facecolor='None', edgecolor='gray', alpha=0.7)
    ax.add_feature(coastline, facecolor='None', edgecolor='k')
    reverse = ax.add_feature(rfaults, facecolor='None', edgecolor='k', linewidth=1.0)
    normal = ax.add_feature(nfaults, facecolor='None', edgecolor='m', linewidth=1.0)
    dextral = ax.add_feature(dfaults, facecolor='None', edgecolor='k', linewidth=1.0)
    sinestral = ax.add_feature(sfaults, facecolor='None', edgecolor='g', linewidth=1.0)
    undefined = ax.add_feature(ufaults, facecolor='None', linestyle='--', edgecolor='k', linewidth=1.0)
    ax.add_feature(tz, linestyle='--', facecolor='None', edgecolor='k', linewidth=1.0, alpha=0.5)
    ax.background_patch.set_facecolor('silver')
    # mapdiff = ax.contourf(lon, lat, mdiff, transform=ccrs.PlateCarree(), cmap='rainbow', vmin=-35, vmax=5)
    # mapdiff = ax.imshow(diff, origin='upper', transform=ccrs.PlateCarree(),
    #                 extent=(np.amin(lon), np.amax(lon), np.amin(lat), np.amax(lat)),
    #                 cmap='rainbow', interpolation='bicubic')
    cmap = colors.ListedColormap(['navy', 'darkblue', 'blue', 'skyblue', 'lightblue', 'white',
                                  'mistyrose', 'lightcoral', 'red', 'darkred', 'maroon'])
    norm = colors.BoundaryNorm([-15, -12.5, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 12.5, 15], cmap.N)
    mapdiff = ax.scatter(lon, lat, c=diff, transform=ccrs.PlateCarree(),
                         s=0.5, marker='o', cmap=cmap, norm=norm)
    cbar = plt.colorbar(mapdiff)
    cbar.set_label('Residual Depth [Km]', rotation=90)
    cbar.set_ticks([-15, -12.5, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 12.5, 15])
    cbar.set_ticklabels([-15, -12.5, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 12.5, 15])
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.8)
    gl.xlabels_top=False
    gl.ylabels_right=False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(xloc)
    gl.ylocator = mticker.FixedLocator(yloc)
    ax.set_extent([-75, -68, -44, -33])
    return fig


# Mapping int
def map_int(idet, lon, lat):
    xloc = np.arange(-76.0, -65+3.0, 3.0)
    yloc = np.arange(-44, -33+2.0, 2.0)
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    coastline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
    border = cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m')
    afaults = cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/active_structure/active_faults.shp').geometries(),
                                      ccrs.PlateCarree())
    ufaults = cfeature.ShapelyFeature(shpreader.Reader('/home/julvelillo/Documentos/andestm_standalone/src/data/shp/unactive_faults/unactive_faults.shp').geometries(),
                                      ccrs.PlateCarree())
    plt.title('Interpolation')
    ax.add_feature(border, facecolor='None', edgecolor='gray', alpha=0.7)
    ax.add_feature(coastline, facecolor='None', edgecolor='k')
    ax.add_feature(afaults, facecolor='None', edgecolor='g', linewidth=1.0)
    ax.add_feature(ufaults, facecolor='None', edgecolor='k', linewidth=1.0)
    ax.background_patch.set_facecolor('silver')
    int = ax.scatter(lon, lat, c=idet, transform=ccrs.PlateCarree(),
                         s=0.5, marker='s', cmap='RdYlBu_r', vmax=0., vmin=-35.)
    cbar = plt.colorbar(int)
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


# Plotting CCR
def plot_ccr(x, y, ccr):
    fig = plt.figure(figsize=(8,8))
    plt.fill_between([0,-40], [-12.5,-52.5], -40, color='maroon')
    plt.fill_between([0,-40], [-10.0,-50.0], [-12.5,-52.5], color='darkred')
    plt.fill_between([0,-40], [-7.5,-47.5], [-10.0,-50.0], color='red')
    plt.fill_between([0,-40], [-5.0,-45.0], [-7.5,-47.5], color='lightcoral')
    plt.fill_between([0,-40], [-2.5,-42.5], [-5.0,-45.0], color='mistyrose')
    plt.fill_between([0,-40], [2.5,-37.5], [5.0,-35.0], color='lightblue')
    plt.fill_between([0,-40], [5.0,-35.0], [7.5,-32.5], color='skyblue')
    plt.fill_between([0,-40], [7.5,-32.5], [10.0,-30.0], color='blue')
    plt.fill_between([0,-40], [10.0,-30.0], [12.5,-27.5], color='darkblue')
    plt.fill_between([0,-40], [12.5,-27.5], 0, color='navy')
    plt.scatter(x, y, marker='.', s=0.5, c='k', alpha=0.3)
    plt.plot([40,-40], [40,-40], color ='black', linewidth=1)
    plt.xlim(-5, -40.0)
    plt.ylim(-5, -40.0)
    plt.xlabel('TMM [km]')
    plt.ylabel('SM [km]')
    plt.text(-10.0, -35.0, 'Overestimated', rotation=45, color='white', weight='semibold')
    plt.text(-28, -20.0, 'Underestimated', rotation=45, color='white', weight='semibold')
    plt.text(-25.0, -6.0, 'Correlation = {}'.format(ccr), color='white')
    return fig


# Plotting hist
def plot_hist(x):
    mu, std = norm.fit(x)
    color = ['navy', 'darkblue', 'blue', 'skyblue', 'lightblue', 'white',
             'mistyrose', 'lightcoral', 'red', 'darkred', 'maroon']
    bins = [-15, -12.5, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 12.5, 15]
    fig = plt.figure(figsize=(8,8))
    n, bins, patches = plt.hist(x, bins, density=True, ec='k')
    for i in range(len(patches)):
        patches[i].set_facecolor(color[i])
    wmin, wmax = plt.xlim()
    w = np.linspace(wmin, wmax, 100)
    p = norm.pdf(w, mu, std)
    plt.plot(w, p, c='k', lw=2, label='PDF')
    plt.legend(loc='best')
    plt.title('Residual Histogram')
    plt.xlabel('Residual [km]')
    plt.ylabel('Density')
    plt.figtext(0.2, 0.8, 'std={:.3f}'.format(std))
    plt.figtext(0.2, 0.7, 'mu={:.3f}'.format(mu))
    plt.xticks([-15, -12.5, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 12.5, 15], rotation='vertical')
    plt.grid(alpha=0.5)
    return fig
