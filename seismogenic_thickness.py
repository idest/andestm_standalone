import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from src.compute import Data, CoordinateSystem, SpatialArray2D
from termomecanico import termomecanico
from src.setup import data_setup, input_setup, exec_setup
from src.plot import heatmap_map, diff_map
from src.colormaps import jet_white_r, get_elevation_diff_cmap
from src.utils import makedir_from_filename


data = Data(*data_setup(), *input_setup())
cs = CoordinateSystem(data.get_cs_data(), 0.2, 1)
x_axis = cs.get_x_axis()
y_axis = cs.get_y_axis()

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
#################### READ AND WRITE FILE WITH LON LAT BINS ####################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################
#-----------------------------------------------------------------------------#
########################## Define Lot Lat Bins ################################
#-----------------------------------------------------------------------------#
lat_bins = np.append(y_axis[::-1]-0.1, y_axis[0]+0.1)
lat_labels = y_axis[::-1]
lon_bins = np.append(x_axis-0.1, x_axis[-1]+0.1)
lon_labels = x_axis
def get_lon_lat_bins(eqs, lat_bins, lon_bins, lat_labels, lon_labels):
    lat_bin = pd.cut(eqs['latitude'].values, lat_bins, labels=lat_labels)
    lon_bin = pd.cut(eqs['longitude'].values, lon_bins, labels=lon_labels)
    eqs['latitude_bin'] = lat_bin
    eqs['longitude_bin'] = lon_bin
    return eqs
#-----------------------------------------------------------------------------#
########################## Earthquakes CSN ####################################
#-----------------------------------------------------------------------------#
#eqs_csn = pd.read_excel('data/earthquakes/CSN_2000_2018_C_SB30.xlsx')
#eqs_csn = get_lon_lat_bins(eqs_csn, lat_bins, lon_bins, lat_labels, lon_labels)
#writer = pd.ExcelWriter('data/earthquakes/CSN_2000_2018_C_SB30_BINS.xlsx')
#eqs_csn.to_excel(writer, 'Sheet1')
#writer.save()
#-----------------------------------------------------------------------------#
######################### Earthquakes USGS ####################################
#-----------------------------------------------------------------------------#
#eqs_usgs = pd.read_excel('data/earthquakes/USGS_1900_2018_C_SB30.xlsx')
#eqs_usgs = get_lon_lat_bins(eqs_usgs, lat_bins, lon_bins, lat_labels, lon_labels)
#writer = pd.ExcelWriter('data/earthquakes/USGS_1900_2018_C_SB30_BINS.xlsx')
#eqs_csn.to_excel(writer, 'Sheet1')
#writer.save()
###############################################################################


print('Reading and processing seismogenic thickness')
###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
######################### GET SEISMOGENIC THICKNESS ###########################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################
#-----------------------------------------------------------------------------#
######################## Define Seismogenic Thickness #########################
#-----------------------------------------------------------------------------#
def get_seismogenic_thickness(eqs, ignore_outliers=False):
    seismogenic_ths = np.zeros((len(x_axis), len(y_axis)))*np.nan
    total_per_bin = np.zeros(len(x_axis)*len(y_axis))
    i = 0
    for lon_idx,lon in enumerate(x_axis):
        for lat_idx,lat in enumerate(y_axis):
            #current_bin_eqs = eqs[
            #        (eqs['longitude_bin'] == lon) &
            #        (eqs['latitude_bin'] == lat)]
            # current_bin_earthquakes
            cbe = eqs[
                    (np.isclose(eqs['longitude_bin'], lon)) &
                    (np.isclose(eqs['latitude_bin'], lat))]
            current_total = len(cbe.index)
            total_per_bin[i] = current_total
            if current_total > 1 and ignore_outliers:
                cbe = cbe[
                    np.abs(cbe['depth'] - cbe['depth'].mean()) <= (1*cbe['depth'].std())
                ]
                current_total_without_outliers = len(cbe.index)
                total_per_bin[i] = current_total_without_outliers
                if current_total != current_total_without_outliers:
                    print('total', current_total)
                    print('total without outliers', current_total_without_outliers)
            seismogenic_ths[lon_idx, lat_idx] = cbe['depth'].max()
            i += 1
    num_eqs = sum(total_per_bin)
    return [SpatialArray2D(seismogenic_ths, cs), num_eqs]
#-----------------------------------------------------------------------------#
########################## Earthquakes CSN ####################################
#-----------------------------------------------------------------------------#
eqs_csn = pd.read_excel('data/earthquakes/CSN_2000_2018_C_SB30_BINS.xlsx',
   usecols=[7,8,9,10,22,23,24,25])
#print(eqs_csn['mag'].min())
#print(eqs_csn['mag'].max())
#eqs_csn = eqs_csn[(eqs_csn['mag'] >= 3.0)]
########## ALL
st_csn = get_seismogenic_thickness(eqs_csn)

############ ISA
#eqs_csn_ISA = eqs_csn[(eqs_csn['ISA'] == True)].copy()
#st_csn_ISA = get_seismogenic_thickness(eqs_csn_ISA)
############ ISA <50 km
#eqs_csn_ISA_LT50 = eqs_csn[
#    (eqs_csn['ISA'] == True) & (eqs_csn['depth'] <= 50.)].copy()
#st_csn_ISA_LT50 = get_seismogenic_thickness(eqs_csn_ISA_LT50)
############ ISA Without Outliers
#st_csn_ISA_WO = get_seismogenic_thickness(eqs_csn_ISA, ignore_outliers=True)
############ ISA Without Outliers < 50 km.
#st_csn_ISA_WO_LT50 = get_seismogenic_thickness(eqs_csn_ISA_LT50, ignore_outliers=True)

############ ISA & OSB
eqs_csn_ISA_OSB = eqs_csn[
    (eqs_csn['ISA'] == True) & (eqs_csn['OSB'] == True)].copy()
st_csn_ISA_OSB = get_seismogenic_thickness(eqs_csn_ISA_OSB)
############ ISA & OSB < 50 km
eqs_csn_ISA_OSB_LT50 = eqs_csn[
    (eqs_csn['ISA'] == True)& (eqs_csn['OSB'] == True)
    & (eqs_csn['depth'] <= 50.)].copy()
st_csn_ISA_OSB_LT50 = get_seismogenic_thickness(eqs_csn_ISA_OSB_LT50)
############ ISA & OSB Without Outliers
st_csn_ISA_OSB_WO = get_seismogenic_thickness(eqs_csn_ISA_OSB, ignore_outliers=True)
############ ISA & OSB Without Outliers
st_csn_ISA_OSB_WO_LT50 = get_seismogenic_thickness(eqs_csn_ISA_OSB_LT50, ignore_outliers=True)

#-----------------------------------------------------------------------------#
######################### Earthquakes USGS ####################################
#-----------------------------------------------------------------------------#
#eqs_usgs = pd.read_excel('data/earthquakes/USGS_1900_2018_C_SB30_BINS.xlsx',
#   usecols=[2,3,4,5,24,25,26,27])
############ ALL
#st_usgs = get_seismogenic_thickness(eqs_usgs)
############ ISA
#eqs_usgs_ISA = eqs_usgs[(eqs_usgs['ISA'] == True)].copy()
#st_usgs_ISA = get_seismogenic_thickness(eqs_usgs_ISA)
############ ISA & OSB
#eqs_usgs_ISA_OSB = eqs_usgs[
#    (eqs_usgs['ISA'] == True) & (eqs_usgs['OSB'] == True)].copy()
#st_usgs_ISA_OSB = get_seismogenic_thickness(eqs_usgs_ISA_OSB)
#-----------------------------------------------------------------------------#
###################### Earthquakes USGS & CSN #################################
#-----------------------------------------------------------------------------#
#eqs_csn = pd.read_excel('data/earthquakes/CSN_2000_2018_C_SB30_BINS.xlsx',
#   usecols=[7,8,9,10,22,23,24,25])
#eqs_usgs = pd.read_excel('data/earthquakes/USGS_1900_2018_C_SB30_BINS.xlsx',
#   usecols=[2,3,4,5,24,25,26,27])
#eqs = pd.concat([eqs_usgs, eqs_csn], ignore_index=True)
############ ALL
#st = get_seismogenic_thickness(eqs)
############ ISA
#eqs_ISA = eqs[(eqs['ISA'] == True)].copy()
############ ISA & OSB
#eqs_ISA_OSB = eqs[(eqs['ISA'] == True) & (eqs['OSB'] == True)].copy()
###############################################################################


###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
######################## DEFINE PLOTTING FUNCTIONS ############################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################
def plot_seismogenic_thickness(
        st_list, titles_list=None, cbar_limits=None, filename=None):
    fig = plt.figure(figsize=(len(st_list)*4,6))
    gs = gridspec.GridSpec(1,len(st_list))
    axs = [fig.add_subplot(gs[0,n]) for n in range(len(st_list))]
    for idx, st in enumerate(st_list):
        quantity_str = ' ({:.1f})'.format(st[1])
        if titles_list:
            title = titles_list[idx] + quantity_str
        else:
            title = str(idx)
        heatmap_map(
            st[0], ax=axs[idx], title=title, colormap=jet_white_r,
            cbar_limits=[0,100], draw_land=False)
    makedir_from_filename(filename)
    plt.savefig(filename + '.png', transparent=True, dpi=900)
    plt.close()

def plot_seismogenic_thickness_vs_eet_scatter(
        st_list, eet_list, titles_list=None, lims_list=None,
        figsize_tuple=None, filename=None):
    if figsize_tuple:
        fig = plt.figure(figsize=figsize_tuple)
    else:
        fig = plt.figure(figsize=(len(st_list)*4,10))
    gs = gridspec.GridSpec(1,len(st_list))
    axs_scatter = [fig.add_subplot(gs[0,n]) for n in range(len(st_list))]
    #gs = gridspec.GridSpec(5,len(st_list))
    #axs_scatter = [fig.add_subplot(gs[0:2,n]) for n in range(len(st_list))]
    #axs_diff = [fig.add_subplot(gs[2:,n]) for n in range(len(st_list))]
    for idx, (st, eet) in enumerate(zip(st_list, eet_list)):
        s = st[0][~np.isnan(st[0])]
        e = eet[~np.isnan(st[0])]
        quantity_str = ' ({:.1f})'.format(st[1])
        if titles_list:
            title = titles_list[idx] + quantity_str
        else:
            title = str(idx)
        axs_scatter[idx].scatter(s, e, edgecolors='black')
        axs_scatter[idx].set_title(title)
        if lims_list:
            axs_scatter[idx].set_xlim(0, lims_list[idx][0])
            axs_scatter[idx].set_ylim(0, lims_list[idx][1])
        #diff = eet - st[0]
        #ax_diff = axs_diff[idx]
        #heatmap_map(diff, colormap=get_elevation_diff_cmap(100), ax=ax_diff,
        #    cbar_limits=[-35,35])
        #ax_diff.set_title(title + '_diff')
    plt.tight_layout()
    makedir_from_filename(filename)
    plt.savefig(filename + '.png', transparent=True)
    plt.close()

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
######################## GET MODEL EET FOR COMPARISON #########################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

print('Running tasks')

model = termomecanico(*input_setup())
eet_m = model.mm.get_eet()
#eet = eet_m
#eet = np.loadtxt('Output/13_Termal/01_Mecanico/Exploracion/Teq_vs_Tef/lc_uc/Archivos/prom.txt')
#eet = np.loadtxt('Output/13_Termal/01_Mecanico/Exploracion/Teq_vs_Tef/lc_uc/Archivos/Quartzite__Clpx-nita_W__Ol-dunite__Antigorite/Quartzite__Clpx-nita_W__Ol-dunite__Antigorite.txt')
#eet = SpatialArray2D(eet, eet_m.cs)
#eet_T07 = np.loadtxt('data/Te_invertido/Interpolados/Te_Tassara.txt')
#eet_T07 = SpatialArray2D(eet_T07, cs).mask_irrelevant_eet()

#heatmap_map(eet, colormap=jet_white_r, cbar_limits=[0,100], title='eet')
#heatmap_map(eet2, colormap=jet_white_r, cbar_limits=[0,100], title='eet2')

exec_input, direTer, direMec = exec_setup()
save_dir =  'Output/EETs/Inputs/'
save_dir_encoded = os.fsencode(save_dir)

st_csn_ISA = st_csn_ISA_OSB
st_csn_ISA_WO_LT50 = st_csn_ISA_OSB_WO_LT50

st_csn_ISA[0] = st_csn_ISA[0].mask_irrelevant_eet()
st_csn_ISA_WO_LT50[0] = st_csn_ISA_WO_LT50[0].mask_irrelevant_eet()

for file in os.listdir(save_dir_encoded):
    filename = os.fsdecode(file)
    eet = np.loadtxt('Output/EETs/Inputs/' + filename)
    eet = SpatialArray2D(eet, eet_m.cs)

    plot_seismogenic_thickness(
        [st_csn_ISA,
         st_csn_ISA_WO_LT50],
        ['CSN_ISA', 'CSN_ISA_WO_LT50'],
        filename='Output/EETs/Outputs/st_maps/' + filename +  '_st_map' )

    plot_seismogenic_thickness_vs_eet_scatter(
        [st_csn_ISA],
        [eet],
        ['CSN_ISA'],
        lims_list=[[125,50]],
        figsize_tuple=[12,5],
        filename='Output/EETs/Outputs/scatter_ISA/' + filename + '_scatter_ISA')

    plot_seismogenic_thickness_vs_eet_scatter(
        [st_csn_ISA_WO_LT50],
        [eet],
        ['CSN_ISA_WO_LT50'],
        lims_list=[[50,50]],
        figsize_tuple=[5,5],
        filename='Output/EETs/Outputs/scatter_ISA_WO_LT50/' + filename + '_scatter_ISA_WO_LT50')

#st_csn_ISA[0] = st_csn_ISA[0].mask_irrelevant_eet()
#st_csn_ISA_LT50[0] = st_csn_ISA_LT50[0].mask_irrelevant_eet()
#st_csn_ISA_WO[0] = st_csn_ISA_WO[0].mask_irrelevant_eet()
#st_csn_ISA_WO_LT50[0] = st_csn_ISA_WO_LT50[0].mask_irrelevant_eet()
#
#plot_seismogenic_thickness(
#    [st_csn_ISA,
#     st_csn_ISA_LT50,
#     st_csn_ISA_WO,
#     st_csn_ISA_WO_LT50],
#    ['CSN_ISA', 'CSN_ISA_LT50', 'CSN_ISA_WO', 'CSN_ISA_WO_LT50'],
#    [len(eqs_csn_ISA.index),
#     len(eqs_csn_ISA_LT50.index),
#     len(eqs_csn_ISA.index),
#     len(eqs_csn_ISA_LT50.index)])
#
#plot_seismogenic_thickness_vs_eet_scatter(
#    [st_csn_ISA,
#     st_csn_ISA_LT50,
#     st_csn_ISA_WO,
#     st_csn_ISA_WO_LT50],
#    #[eet_T07, eet_T07, eet_T07, eet_07],
#    [eet, eet, eet, eet],
#    ['CSN_ISA', 'CSN_ISA_LT50', 'CSN_ISA_WO', 'CSN_ISA_WO_LT50'],
#    [len(eqs_csn_ISA.index),
#     len(eqs_csn_ISA_LT50.index),
#     len(eqs_csn_ISA.index),
#     len(eqs_csn_ISA_LT50.index)])

#st_csn_ISA_OSB[0] = st_csn_ISA_OSB[0].mask_irrelevant_eet()
#st_csn_ISA_OSB_LT50[0] = st_csn_ISA_OSB_LT50[0].mask_irrelevant_eet()
#st_csn_ISA_OSB_WO[0] = st_csn_ISA_OSB_WO[0].mask_irrelevant_eet()
#st_csn_ISA_OSB_WO_LT50[0] = st_csn_ISA_OSB_WO_LT50[0].mask_irrelevant_eet()
#
#plot_seismogenic_thickness(
#    [st_csn_ISA_OSB,
#     st_csn_ISA_OSB_LT50,
#     st_csn_ISA_OSB_WO,
#     st_csn_ISA_OSB_WO_LT50],
#    ['CSN_ISA_OSB', 'CSN_ISA_OSB_LT50', 'CSN_ISA_OSB_WO', 'CSN_ISA_OSB_WO_LT50'],
#    [len(eqs_csn_ISA_OSB.index),
#     len(eqs_csn_ISA_OSB_LT50.index),
#     len(eqs_csn_ISA_OSB.index),
#     len(eqs_csn_ISA_OSB_LT50.index)])
#
#plot_seismogenic_thickness_vs_eet_scatter(
#    [st_csn_ISA_OSB,
#     st_csn_ISA_OSB_LT50,
#     st_csn_ISA_OSB_WO,
#     st_csn_ISA_OSB_WO_LT50],
#    #[eet_T07, eet_T07, eet_T07, eet_T07],
#    [eet, eet, eet, eet],
#    ['CSN_ISA_OSB', 'CSN_ISA_OSB_LT50', 'CSN_ISA_OSB_WO', 'CSN_ISA_OSB_WO_LT50'],
#    [len(eqs_csn_ISA_OSB.index),
#     len(eqs_csn_ISA_OSB_LT50.index),
#     len(eqs_csn_ISA_OSB.index),
#     len(eqs_csn_ISA_OSB_LT50.index)])
#plt.show()
#model = termomecanico(*input_setup())
