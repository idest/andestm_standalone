import setup
from datos_q import *
from rmse import interpolated_shf_bsi
from calc_rmse import calc_rmse
from calc_rmse import calc_rmse_error
import numpy as np
np.set_printoptions(threshold=np.nan)
from plot import map_q_surface_2
from plot import map_surface_heat_flow
from plot import plot_diffs
from utils import DotDict
import os
import sys
import resource

rmse_model,datos_rmse_error_shf = calc_rmse_error(datos_q_shf,
                                                  interpolated_shf_bsi,
                                                  datos_q_min_shf,
                                                  datos_q_max_shf, error)

map_surface_heat_flow(x_axis, y_axis, tmc, direTer, datos_q,
                      surface_heat_flow=surface_heat_flow,
                      topo=False, rmse=rmse_model,
                      datos_rmse_error=datos_rmse_error_shf)

#Mapa Diffs Max
#rmse_max = calc_rmse(datos_q_max_shf, interpolated_shf_bsi)
#map_q_surface_2(x_axis, y_axis, tmc, direTer, data_q=datos_q_max,
#                data_cmap='diff', interpolated_heat_flow=interpolated_shf_bsi,
#                topo=False, name='Max_Diff_Map',rmse=rmse_max)

#Mapa Diffs Min
#rmse_min = calc_rmse(datos_q_min_shf, interpolated_shf_bsi)
#map_q_surface_2(x_axis, y_axis, tmc, direTer, data_q=datos_q_min,
#                data_cmap='diff', interpolated_heat_flow=interpolated_shf_bsi,
#                topo=False, name='Min_Diff_Map',rmse=rmse_min)

#Mapa Diffs Prom
rmse_prom = calc_rmse(datos_q_shf, interpolated_shf_bsi)
map_q_surface_2(x_axis, y_axis, tmc, direTer, data_q=datos_q, data_cmap='diff',
                interpolated_heat_flow=interpolated_shf_bsi, topo=False,
                name='Prom_Diff_Map',rmse=rmse_prom)

#Mapa Diff RMSE error
rmse_error,data_salida = calc_rmse_error(datos_q_shf, interpolated_shf_bsi,
                                         datos_q_min_shf,
                                         datos_q_max_shf, error)

map_q_surface_2(x_axis, y_axis, tmc, direTer, data_q=datos_q, data_cmap='diff',
                interpolated_heat_flow=interpolated_shf_bsi, topo=False,
                name='Prom_Diff_RMSE_corrected', rmse=rmse_error,
                datos_rmse_error=data_salida)

#Plotter Diffs
plot_diffs(interpolated_shf_bsi,datos_q,error)


#detachment = plot.get_detachment(CS,GM,MM)

