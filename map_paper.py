import sys
sys.path.insert(0, 'src/')
import numpy as np
import numpy.ma as ma
import math
from scipy.stats import norm
import pandas as pd
from setup import data_setup, input_setup, exec_setup
from compute import compute
from src.stats import evaluate_model
from src.datos_q import shf_data
from src.utils import MidPointNorm, round_to_1, get_magnitude
from src.colormaps import (jet_white_r, get_diff_cmap,
                           get_elevation_diff_cmap, jet_white_r_2)
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import colors
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt


# tm model data
gm_data, areas, trench_age, rhe_data, coast = data_setup()
t_input, m_input = input_setup()
exec_input, direTer, direMec = exec_setup()
names = ['lon', 'lat', 'slab_lab', 'moho', 'icd', 'topo', 'topo_r']
d_model = pd.read_csv('data/Modelo2.dat', sep="\s+", names=names)


# running tm model
def termomecanico(t_input, m_input):
    model = compute(gm_data, areas, trench_age, rhe_data, coast, t_input, m_input)
    cs = model.mm.get_coordinate_system()
    gm = model.mm.get_geometric_model()
    return cs, gm, model


# surface heat flow, tm output, estimators and map variables
cs, gm, model = termomecanico(t_input, m_input)
shf = model.tm.get_surface_heat_flow(format='positive milliwatts')
shf = pd.DataFrame(shf)
shf.to_csv('Output/{}shfmodel.csv' .format(exec_input.temcaso))
qdata = pd.read_csv('datos_Q/datos_q_def.csv', sep='\t')
q = qdata['Heat Flow [mW/m2]']
lon = cs.get_x_axis()
lat = cs.get_y_axis()
weight = [False,True]


# # heat flow map
# def qmap(data):
#     xloc = np.arange(-80, -60+5.0, 5.0)
#     yloc = np.arange(-45, -10+5.0, 5.0)
#     lon = data['Longitude']
#     lat = data['Latitude']
#     # im_ex = (-85.079166667,-58.887500000,-58.183333333,-8.437500000)
#     im_ex = (-82.8250, -52.1917, -47.7583, -8.3417)
#     # im_ex = (-80, -60, -45, -10)
#     im = plt.imread('/home/julvelillo/Documentos/general data/raster/neic etopo/basetif.png')
#     map = plt.axes(projection=ccrs.PlateCarree())
#     coastline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
#     border = cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m')
#     plt.title('Heat Flow Map')
#     map.add_feature(border, facecolor='None', edgecolor='k', linewidth=0.1, alpha=0.7)
#     map.add_feature(coastline, facecolor='None', edgecolor='k', linewidth=0.3)
#     map.imshow(im, interpolation='nearest', origin='upper', extent=im_ex, transform=ccrs.PlateCarree())
#     cmap = plt.cm.afmhot
#     # iteration for each data type
#     for i in np.unique(data['Data Type']):
#         datai = data.loc[data['Data Type'] == i]
#         marker = ['o','^','p','s']
#         label = ['ODP Borehole', 'Land Borehole', 'Geochemestry', 'Marine Geophysics']
#         q = datai['Heat Flow [mW/m2]']
#         lon, lat = datai['Longitude'], datai['Latitude']
#         for w in datai['Error [mW/m2]']:
#             if math.isnan(w) == False:
#                 scat = map.scatter(lon, lat, c=q, cmap=cmap, edgecolor='k', linewidth=0.4,
#                                     label=label[i-1], marker=marker[i-1], vmin=0, vmax=100,
#                                     transform=ccrs.PlateCarree())
#             else:
#                 scat = map.scatter(lon, lat, c=q, cmap=cmap, edgecolor='k', linewidth=0,
#                                     label=label[i-1], marker=marker[i-1], vmin=0, vmax=100,
#                                     transform=ccrs.PlateCarree())
#     cbar = plt.colorbar(scat, ticks=np.arange(0, 100, 10), boundaries=np.linspace(0, 100, 1000))
#     cbar.set_label('Heat Flow [W/m²]', rotation=90)
#     # # odp borehole
#     # data1 = data.loc[data['Data Type'] == 1]
#     # q1 = data1['Heat Flow [mW/m2]']
#     # lon1, lat1 = data1['Longitude'], data1['Latitude']
#     # scat1 = map.scatter(lon1, lat1, c=q1, cmap=cmap, edgecolor='k', linewidth=0.2,
#     #                     label='ODP Borehole', marker='o',
#     #                     transform=ccrs.PlateCarree())
#     # # land borehole
#     # data2 = data.loc[data['Data Type'] == 2]
#     # q2 = data2['Heat Flow [mW/m2]']
#     # lon2, lat2 = data2['Longitude'], data2['Latitude']
#     # scat2 = map.scatter(lon2, lat2, c=q2, cmap=cmap, edgecolor='k', linewidth=0.2,
#     #                     label='Land Borehole', marker='^',
#     #                     transform=ccrs.PlateCarree())
#     # # geochemestry
#     # data3 = data.loc[data['Data Type'] == 3]
#     # q3 = data3['Heat Flow [mW/m2]']
#     # lon3, lat3 = data3['Longitude'], data3['Latitude']
#     # scat3 = map.scatter(lon3, lat3, c=q3, cmap=cmap, edgecolor='k', linewidth=0.2,
#     #                     label='Geochemestry', marker='p',
#     #                     transform=ccrs.PlateCarree())
#     # # marine geophysics
#     # data4 = data.loc[data['Data Type'] == 4]
#     # q4 = data4['Heat Flow [mW/m2]']
#     # lon4, lat4 = data4['Longitude'], data4['Latitude']
#     # scat4 = map.scatter(lon4, lat4, c=q4, cmap=cmap, edgecolor='k', linewidth=0.2,
#     #                     label='Marine Geophysics', marker='s',
#     #                     transform=ccrs.PlateCarree())
#     gl = map.gridlines(draw_labels=True, linewidth=0.1, color='gray', alpha=0.1)
#     gl.xlabels_top=False
#     gl.ylabels_right=False
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlocator = mticker.FixedLocator(xloc)
#     gl.ylocator = mticker.FixedLocator(yloc)
#     map.set_extent([-80, -60, -45, -10], crs=ccrs.PlateCarree())
#     plt.savefig('map_final/heatflow_map.eps', format='eps', dpi=1000)
#     # plt.show()
#     plt.close()


# # multimap
# def multimap(lon, lat, shf, df_noerror, df_error, st_error, st_noerror):
#     #laoding data and variables
#     xloc = np.arange(-80, -60+5.0, 5.0)
#     yloc = np.arange(-45, -10+5.0, 5.0)
#     mlon, mlat = np.meshgrid(lon, lat)
#     diff_noerror = df_noerror['diffs']
#     diff_error = df_error['diffs']
#     sx_noerror = df_noerror['lons']
#     sy_noerror = df_noerror['lats']
#     sx_error = df_error['lons']
#     sy_error = df_error['lats']
#     dt_noerror = df_noerror['data_types']
#     dt_error = df_error['data_types']
#     sigmas_noerror = st_noerror['sigmas']
#     sigmas_error = st_error['sigmas']
#     diff_max_noerror = np.nanmax(diff_noerror)
#     diff_min_noerror = np.nanmin(diff_noerror)
#     diff_limit_noerror = np.nanmax([abs(diff_max_noerror),
#                                     abs(diff_min_noerror)])
#     diff_limit_noerror = round_to_1(diff_limit_noerror, 'ceil')
#     diff_max_error = np.nanmax(diff_error)
#     diff_min_error = np.nanmin(diff_error)
#     diff_limit_error = np.nanmax([abs(diff_max_error),
#                                     abs(diff_min_error)])
#     diff_limit_error = round_to_1(diff_limit_error, 'ceil')
#     #diff_step = 10**get_magnitude(diff_limit)
#     diff_step = 5.0
#     divisions = np.arange(-70, 70+1, diff_step)
#     # divisions = np.arange(-diff_limit, diff_limit, diff_step)
#     ticks = np.arange(-70, 70+1, 10)
#     # ticks = np.arange(-diff_limit, diff_limit+diff_step, 2*diff_step)
#     dbins = len(divisions) - 1
#     diff_cmap = get_diff_cmap(dbins)
#     ccolors = plt.get_cmap(diff_cmap)(np.arange(dbins, dtype='int'))
#     ccolors[math.floor((dbins-1)/2)] = [1., 1., 1., 1.]
#     ccolors[math.ceil((dbins-1)/2)] = [1., 1., 1., 1.]
#     diff_cmap = colors.ListedColormap(ccolors)
#     # diff_cmap = cm.RdBu
#     dnorm = MidPointNorm(midpoint=0, vmin=-70, vmax=70)
#     # dnorm = colors.BoundaryNorm([-70, -65, -60, -55, -50, -45, -40, -35, -30, -25,
#     #                             -20, -15, -10, -5, 5, 10, 15, 20, 25, 30, 35, 40,
#     #                              45, 50, 55, 60, 65, 70], diff_cmap)
#     heat_flow = shf
#     hf_masked = np.ma.masked_invalid(heat_flow)
#
#
#     #plotting figure
#     fig = plt.figure()
#     gs = gridspec.GridSpec(1,2)
#     map = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
#     coastline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m')
#     border = cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '10m')
#     plt.title('Thermal model')
#     map.background_patch.set_fill(False)
#     map.add_feature(border, facecolor='None', edgecolor='gray', linewidth=0.3, alpha=0.7)
#     map.add_feature(coastline, facecolor='None', edgecolor='k', linewidth=0.5)
#
#
#     #plotting heat flow
#     cmap = plt.cm.afmhot
#     hfcf = map.contourf(mlon, mlat, heat_flow.T, 24, transform=ccrs.PlateCarree(),
#                         cmap=cmap, vmin=0, vmax=100)
#     m = plt.cm.ScalarMappable(cmap=cm.afmhot)
#     m.set_array(heat_flow)
#     m.set_clim(0., 100.)
#     cbar = plt.colorbar(m, ticks=np.arange(0, 100, 10), boundaries=np.linspace(0, 100, 1000))
#     # cbar = plt.colorbar(hfcf)
#     cbar.set_ticklabels(np.arange(0, 100, 10))
#     cbar.set_label('Heat Flow [W/m²]', rotation=90)
#     #fix the white borders at each contourf
#     for c in hfcf.collections:
#         c.set_edgecolor("face")
#
#
#     #plotting diff scatter
#     diff_noerror = np.ma.masked_invalid(diff_noerror)
#     diff_error = np.ma.masked_invalid(diff_error)
#
#     # data type 1
#     df1 = df_noerror.loc[df_noerror['data_types'] == 1.0]
#     diff1 = df1['diffs']
#     sx1, sy1 = df1['lons'], df1['lats']
#     for de1 in df1['data_errors']:
#         if math.isnan(de1) == False:
#             dscatter1 = map.scatter(sx1, sy1, c=diff1, cmap=diff_cmap, marker='o',
#                                     norm=dnorm, label='ODP Borehole', edgecolors='k',
#                                     zorder=10, linewidth=1.0, transform=ccrs.PlateCarree())
#         else:
#             dscatter1 = map.scatter(sx1, sy1, c=diff1, cmap=diff_cmap, marker='o',
#                                     norm=dnorm, label='ODP Borehole', edgecolors='k',
#                                     zorder=10, linewidth=0, transform=ccrs.PlateCarree())
#
#     # data type 2
#     df2 = df_noerror.loc[df_noerror['data_types'] == 2.0]
#     diff2 = df2['diffs']
#     sx2, sy2 = df2['lons'], df2['lats']
#     for de2 in df2['data_errors']:
#         if math.isnan(de2) == False:
#             dscatter2 = map.scatter(sx2, sy2, c=diff2, cmap=diff_cmap, marker='^',
#                                     norm=dnorm, label='ODP Borehole', edgecolors='k',
#                                     zorder=10, linewidth=1.0, transform=ccrs.PlateCarree())
#         else:
#             dscatter2 = map.scatter(sx2, sy2, c=diff2, cmap=diff_cmap, marker='^',
#                                     norm=dnorm, label='ODP Borehole', edgecolors='k',
#                                     zorder=10, linewidth=0, transform=ccrs.PlateCarree())
#
#     # data type 3
#     df3 = df_noerror.loc[df_noerror['data_types'] == 3.0]
#     diff3 = df3['diffs']
#     sx3, sy3 = df3['lons'], df3['lats']
#     for de3 in df3['data_errors']:
#         if math.isnan(de3) == False:
#             dscatter3 = map.scatter(sx3, sy3, c=diff3, cmap=diff_cmap, marker='p',
#                                     norm=dnorm, label='ODP Borehole', edgecolors='k',
#                                     zorder=10, linewidth=1.0, transform=ccrs.PlateCarree())
#         else:
#             dscatter3 = map.scatter(sx3, sy3, c=diff3, cmap=diff_cmap, marker='p',
#                                     norm=dnorm, label='ODP Borehole', edgecolors='k',
#                                     zorder=10, linewidth=0, transform=ccrs.PlateCarree())
#
#     #data type 4
#     df4 = df_noerror.loc[df_noerror['data_types'] == 4.0]
#     diff4 = df4['diffs']
#     sx4, sy4 = df4['lons'], df4['lats']
#     for de4 in df4['data_errors']:
#         if math.isnan(de4) == False:
#             dscatter4 = map.scatter(sx4, sy4, c=diff4, cmap=diff_cmap, marker='s',
#                                     norm=dnorm, label='ODP Borehole', edgecolors='k',
#                                     zorder=10, linewidth=1.0, transform=ccrs.PlateCarree())
#         else:
#             dscatter4 = map.scatter(sx4, sy4, c=diff4, cmap=diff_cmap, marker='s',
#                                     norm=dnorm, label='ODP Borehole', edgecolors='k',
#                                     zorder=10, linewidth=0, transform=ccrs.PlateCarree())
#     gl = map.gridlines(draw_labels=True, linewidth=0.1, color='gray', alpha=0.3)
#     gl.xlabels_top=False
#     gl.ylabels_right=False
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlocator = mticker.FixedLocator(xloc)
#     gl.ylocator = mticker.FixedLocator(yloc)
#     map.set_extent([-80, -60, -45, -10])
#     map.set_aspect('auto', adjustable=None)
#
#
#     #plotting histogram
#     ax2 = fig.add_subplot(gs[0, 1])
#     N_ne, bins_ne, patches_ne = ax2.hist(diff_noerror, bins=divisions,
#                                          orientation='horizontal', ec='grey',
#                                          linewidth=0.2, alpha=0.5)
#     N_e, bins_e, patches_e = ax2.hist(diff_error, bins=divisions,
#                                          orientation='horizontal', ec='k',
#                                          linewidth=0.2, alpha=1.0)
#     ax2.set_yticks(ticks)
#     ax2.set_ylim([-70, 70])
#     ax2.yaxis.tick_right()
#     hcmap = cm.RdBu
#     bnorm_ne = Normalize(bins_ne.min(), bins_ne.max())
#     for bin_ne, patch_ne in zip(bins_ne, patches_ne):
#         color = diff_cmap(bnorm_ne(bin_ne))
#         # color = hcmap(dnorm(bin))
#         patch_ne.set_facecolor(color)
#     bnorm_e = Normalize(bins_e.min(), bins_e.max())
#     for bin_e, patch_e in zip(bins_e, patches_e):
#         color = diff_cmap(bnorm_e(bin_e))
#         patch_e.set_facecolor(color)
#     # non error stats
#     ax2.axhline(y=sigmas_noerror.n_1_sigma, color='purple', linewidth=0.5, alpha=0.2)
#     ax2.text(25,sigmas_noerror.n_1_sigma-0.6*diff_step,r'-$\sigma$',fontsize=6)
#     ax2.axhline(y=sigmas_noerror.p_1_sigma, color='purple', linewidth=0.5, alpha=0.2)
#     ax2.text(25,sigmas_noerror.p_1_sigma+0.1*diff_step,r'+$\sigma$',fontsize=6)
#     ax2.axhline(y=sigmas_noerror.n_2_sigma, color='purple', linewidth=0.5, alpha=0.2)
#     ax2.text(25,sigmas_noerror.n_2_sigma-0.65*diff_step,r'-2$\sigma$',fontsize=6)
#     ax2.axhline(y=sigmas_noerror.p_2_sigma, color='purple', linewidth=0.5, alpha=0.2)
#     ax2.text(25,sigmas_noerror.p_2_sigma+0.1*diff_step,r'+2$\sigma$',fontsize=6)
#     # error stats
#     ax2.axhline(y=sigmas_error.n_1_sigma, color='magenta', linewidth=0.5, alpha=0.2)
#     ax2.text(15,sigmas_error.n_1_sigma-0.6*diff_step,r'-$\sigma$',fontsize=6)
#     ax2.axhline(y=sigmas_error.p_1_sigma, color='magenta', linewidth=0.5, alpha=0.2)
#     ax2.text(15,sigmas_error.p_1_sigma+0.1*diff_step,r'+$\sigma$',fontsize=6)
#     ax2.axhline(y=sigmas_error.n_2_sigma, color='magenta', linewidth=0.5, alpha=0.2)
#     ax2.text(15,sigmas_error.n_2_sigma-0.65*diff_step,r'-2$\sigma$',fontsize=6)
#     ax2.axhline(y=sigmas_error.p_2_sigma, color='magenta', linewidth=0.5, alpha=0.2)
#     ax2.text(15,sigmas_error.p_2_sigma+0.1*diff_step,r'+2$\sigma$',fontsize=6)
#     ax2.set_aspect(1.2)
#     ax2.set_ylabel('Residual Values [W/m²]')
#     ax2.yaxis.set_label_position('right')
#     # plt.title('Residual Values [W/m²]')
#     gs.set_width_ratios([map, ax2])
#     fig.savefig('map_final/multi_map.pdf', transparent='True', format='pdf')
#     plt.close()


# # iterating for each weight possibility
# for w in weight:
#     if w == False:
#         st_noerror, df_noerror = evaluate_model(shf, shf_data, weigh_errors=w, return_dataframe=True)
#     elif w == True:
#         st_error, df_error = evaluate_model(shf, shf_data, weigh_errors=w, return_dataframe=True)


# print(df_noerror['diffs'].describe())
# ax = df_noerror['diffs'].hist(bins=50)
# plt.show()
# hfmap = multimap(lon, lat, shf, df_noerror, df_error, st_noerror, st_error)
# qmap = qmap(qdata)
