import numpy as np
import pandas as pd

shf_data_file = np.loadtxt('datos_Q/QsObs_definitivo.txt', comments='#')

#print('hola')
#np.savetxt('shf_dt', shf_data_file)

# Limitar datos al area de estudio
shf_data_file = shf_data_file[shf_data_file[:,0] > -80.]
shf_data_file = shf_data_file[shf_data_file[:,0] < -60.]
shf_data_file = shf_data_file[shf_data_file[:,1] > -45.]
shf_data_file = shf_data_file[shf_data_file[:,1] < -10.]

# Usar solo datos con valores de surface heat flow < 130 mW/m2
shf_data_file = shf_data_file[shf_data_file[:,2] <= 130]

# Formatear datos a W/m2
#shf_data_file[:,2] = shf_data_file[:,2] * 1.e-3
#shf_data_file[:,3] = shf_data_file[:,3] * 1.e-3
#np.savetxt('datos_Q/QsObsFormatted.txt', shf_data_file)

# All data
shf_data_lons = shf_data_file[:,0]
shf_data_lats = shf_data_file[:,1]
shf_data_values = shf_data_file[:,2]
shf_data_errors = abs(shf_data_file[:,3])
shf_data_refs = shf_data_file[:,4]
shf_data_types = shf_data_file[:,5]

shf_data_coords = np.append(shf_data_lons[:, np.newaxis],
                            shf_data_lats[:, np.newaxis], axis=1)
shf_data_max = shf_data_values + shf_data_errors
shf_data_min = shf_data_values - shf_data_errors

shf_data = pd.DataFrame({
    'lons': shf_data_lons,
    'lats': shf_data_lats,
    'data_values': shf_data_values,
    'data_errors': shf_data_errors,
    'data_types': shf_data_types,
    'data_refs': shf_data_refs
})

#shf_df_ne = shf_df.drop(columns=['data_error'])
#shf_df_e = shf_df.dropna(subset=['data_error'])

## Marine Geophysics
#type_1 = shf_data_table[:,5]==1
#shf_data_x_1 = shf_data_x[np.where(type_1)]
#shf_data_y_1 = shf_data_y[np.where(type_1)]
#shf_data_1 = shf_data[np.where(type_1)]
#shf_data_types_1 = shf_data_types[np.where(type_1)]
#shf_data_error_1 = shf_data_error[np.where(type_1)]
#shf_data_max_1 = shf_data_1  + shf_data_error_1
#shf_data_min_1 = shf_data_1 - shf_data_error_1
#
## Geochemical
#type_2 = shf_data_table[:,5]==2
#shf_data_x_2 = shf_data_x[np.where(type_2)]
#shf_data_y_2 = shf_data_y[np.where(type_2)]
#shf_data_2 = shf_data[np.where(type_2)]
#shf_data_types_2 = shf_data_types[np.where(type_2)]
#shf_data_error_2 = shf_data_error[np.where(type_2)]
#shf_data_max_2 = shf_data_2  + shf_data_error_2
#shf_data_min_2 = shf_data_2 - shf_data_error_2
#
## Land Borehole
#type_3 = shf_data_table[:,5]==3
#shf_data_x_3 = shf_data_x[np.where(type_3)]
#shf_data_y_3 = shf_data_y[np.where(type_3)]
#shf_data_3 = -shf_data[np.where(type_3)]
#shf_data_types_3 = shf_data_types[np.where(type_3)]
#shf_data_error_3 = shf_data_error[np.where(type_3)]
#shf_data_max_3 = shf_data_3  + shf_data_error_3
#shf_data_min_3 = shf_data_3 - shf_data_error_3
#
## ODP Borehole
#type_4 = shf_data_table[:,5]==4
#shf_data_x_4 = shf_data_x[np.where(type_4)]
#shf_data_y_4 = shf_data_y[np.where(type_4)]
#shf_data_4 = -shf_data[np.where(type_4)]
#shf_data_types_4 = shf_data_types[np.where(type_4)]
#shf_data_error_4 = shf_data_error[np.where(type_4)]
#shf_data_max_4 = shf_data_4  + shf_data_error_4
#shf_data_min_4 = shf_data_4 - shf_data_error_4
