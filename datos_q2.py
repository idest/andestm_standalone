shf_data = np.loadtxt('datos_Q/QsObs.txt', comments='#')

# Limitar datos al area de estudio
shf_data = shf_data[shf_data[:,0] > -80.]
shf_data = shf_data[shf_data[:,0] < -60.]
shf_data = shf_data[shf_data[:,1] > -45.]
shf_data = shf_data[shf_data[:,1] < -10.]

# Usar solo datos con valores de surface heat flow < 120 mW/m2
shf_data = shf_data[shf_data[:,2] <= 120]

# Error de los datos
error = shf_data[:,-1]

# Datos + error
shf_data_max = shf_data.copy()
shf_data_max[:,2] = shf_data_max[:,2] + error
shf_data_max_formatted = -shf_data_max[:,2]*1e-3

# Datos - error
shf_data_min = shf_data.copy()
shf_data_min[:,2] = shf_data_min[:,2] - error
shf_data_min_formatted = -shf_data_min[:,2]*1e-3

# All data
shf_data_x = shf_data[:,0]
shf_data_y = shf_data[:,1]
shf_data_formatted = -shf_data[:,2]*1e-3
# Marine Geophysics
shf_data_x_1 = shf_data[:,0][np.where(shf_data[:,-2]==1)]
shf_data_y_1 = shf_data[:,1][np.where(shf_data[:,-2]==1)]
shf_data_formatted_1 = -shf_data[:,2][np.where(shf_data[:,-2]==1)]*1.e-3
# Geochemical
shf_data_x_2 = shf_data[:,0][np.where(shf_data[:,-2]==2)]
shf_data_y_2 = shf_data[:,1][np.where(shf_data[:,-2]==2)]
shf_data_formatted_2 = -shf_data[:,2][np.where(shf_data[:,-2]==2)]*1.e-3
# Land Borehole
shf_data_x_3 = shf_data[:,0][np.where(shf_data[:,-2]==3)]
shf_data_y_3 = shf_data[:,1][np.where(shf_data[:,-2]==3)]
shf_data_formatted_3 = -shf_data[:,2][np.where(shf_data[:,-2]==3)]*1.e-3
# ODP Borehole
shf_data_x_4 = shf_data[:,0][np.where(shf_data[:,-2]==4)]
shf_data_y_4 = shf_data[:,1][np.where(shf_data[:,-2]==4)]
shf_data_formatted_4 = -shf_data[:,2][np.where(shf_data[:,-2]==4)]*1.e-3
