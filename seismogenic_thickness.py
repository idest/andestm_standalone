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
def get_seismogenic_thickness(eqs):
    seismogenic_ths = np.zeros((len(x_axis), len(y_axis)))*np.nan
    for lon_idx,lon in enumerate(x_axis):
        for lat_idx,lat in enumerate(y_axis):
            #current_bin_eqs = eqs[
            #        (eqs['longitude_bin'] == lon) &
            #        (eqs['latitude_bin'] == lat)]
            current_bin_eqs = eqs[
                    (np.isclose(eqs['longitude_bin'], lon)) &
                    (np.isclose(eqs['latitude_bin'], lat))]
            seismogenic_ths[lon_idx, lat_idx] = current_bin_eqs['depth'].max()
    return SpatialArray2D(seismogenic_ths, cs)
#-----------------------------------------------------------------------------#
########################## Earthquakes CSN ####################################
#-----------------------------------------------------------------------------#
eqs_csn = pd.read_excel('data/earthquakes/CSN_2000_2018_C_SB30_BINS.xlsx',
   usecols=[7,8,9,10,22,23,24,25])
########### ALL
st_csn = get_seismogenic_thickness(eqs_csn)
############ ISA
eqs_csn_ISA = eqs_csn[(eqs_csn['ISA'] == True)].copy()
st_csn_ISA = get_seismogenic_thickness(eqs_csn_ISA) 
############ ISA & OSB
eqs_csn_ISA_OSB = eqs_csn[
    (eqs_csn['ISA'] == True) & (eqs_csn['OSB'] == True)].copy()
st_csn_ISA_OSB = get_seismogenic_thickness(eqs_csn_ISA_OSB)
############ ISA <50 km
eqs_csn_ISA_LT50 = eqs_csn[
    (eqs_csn['ISA'] == True) & (eqs_csn['depth'] <= 50.)].copy()
st_csn_ISA_LT50 = get_seismogenic_thickness(eqs_csn_ISA_LT50) 
############ ISA & OSB < 50 km
eqs_csn_ISA_OSB_LT50 = eqs_csn[
    (eqs_csn['ISA'] == True)& (eqs_csn['OSB'] == True)
    & (eqs_csn['depth'] <= 50.)].copy()
st_csn_ISA_OSB_LT50 = get_seismogenic_thickness(eqs_csn_ISA_OSB_LT50)
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
        st_list, titles_list=None, quantity_list=None, cbar_limits=None):
    fig = plt.figure(figsize=(len(st_list)*4,6))
    gs = gridspec.GridSpec(1,len(st_list))
    axs = [fig.add_subplot(gs[0,n]) for n in range(len(st_list))]
    for idx, st in enumerate(st_list):
        quantity_str = ''
        if quantity_list:
            quantity_str = ' ({:d})'.format(quantity_list[idx])
        if titles_list:
            title = titles_list[idx] + quantity_str
        else:
            title = str(idx)
        heatmap_map(
            st, ax=axs[idx], title=title, colormap='viridis',
            cbar_limits=cbar_limits)
    plt.tight_layout()
    #plt.show()
def plot_seismogenic_thickness_vs_eet_scatter(
        st_list, eet_list, titles_list=None, quantity_list=None):
    fig = plt.figure(figsize=(len(st_list)*4,4))
    gs = gridspec.GridSpec(1,len(st_list))
    axs = [fig.add_subplot(gs[0,n]) for n in range(len(st_list))]
    for idx, (st, eet) in enumerate(zip(st_list, eet_list)):
        s = st[~np.isnan(st)]
        e = eet[~np.isnan(st)]
        quantity_str = ''
        if quantity_list:
            quantity_str = ' ({:d})'.format(quantity_list[idx])
        if titles_list:
            title = titles_list[idx] + quantity_str
        else:
            title = str(idx)
        axs[idx].scatter(s,e)
        axs[idx].set_title(title)
    plt.tight_layout()
    #plt.show()

###############################################################################
#/////////////////////////////////////////////////////////////////////////////#
######################## GET MODEL EET FOR COMPARISON #########################
#/////////////////////////////////////////////////////////////////////////////#
###############################################################################

print('Running tasks')
plot_seismogenic_thickness(
    [st_csn_ISA,
     st_csn_ISA_OSB,
     st_csn_ISA_LT50,
     st_csn_ISA_OSB_LT50],
    ['CSN_ISA', 'CSN_ISA_OSB', 'CSN_ISA_LT50', 'CSN_ISA_OSB_LT50'],
    [len(eqs_csn_ISA.index),
     len(eqs_csn_ISA_OSB.index),
     len(eqs_csn_ISA_LT50.index),
     len(eqs_csn_ISA_OSB_LT50.index)])

#model, _, _ = termomecanico(*input_setup())
#eet = model.mm.get_eet()
eet_T07 = np.loadtxt('data/Te_invertido/Interpolados/Te_Tassara.txt')
eet_T07 = SpatialArray2D(eet_T07, cs).mask_irrelevant_eet()

plot_seismogenic_thickness_vs_eet_scatter(
    [st_csn_ISA,
     st_csn_ISA_OSB,
     st_csn_ISA_LT50,
     st_csn_ISA_OSB_LT50],
    [eet_T07, eet_T07, eet_T07, eet_T07],
    ['CSN_ISA', 'CSN_ISA_OSB', 'CSN_ISA_LT50', 'CSN_ISA_OSB_LT50'],
    [len(eqs_csn_ISA.index),
     len(eqs_csn_ISA_OSB.index),
     len(eqs_csn_ISA_LT50.index),
     len(eqs_csn_ISA_OSB_LT50.index)])
plt.show()
#model, _, _ = termomecanico(*input_setup())
