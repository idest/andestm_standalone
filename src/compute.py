import math
import numpy as np
from dotmap import DotMap
from box import Box
import inspect
#TODO: change all references to DotMap with Box
np.seterr(divide='ignore', invalid='ignore')

class Data(object):

    @staticmethod
    def __get_max(values):
        max_v = np.nanmax(values)
        return max_v

    @staticmethod
    def __get_min(values):
        min_v = np.nanmin(values)
        return min_v

    def __init__(self, initial_data, slab_lab_areas, trench_age, rheologic_data,
                 coast, t_input, m_input):

        self.coordinate_data = self.__set_coordinate_data(initial_data)
        self.geometric_data = self.__set_geometric_data(initial_data)
        self.slab_lab_areas = np.asarray(slab_lab_areas, dtype=bool)
        self.coast = np.asarray(coast)
        self.trench_age = trench_age
        self.rheologic_data = rheologic_data
        self.t_input = t_input
        self.m_input = m_input

    def __set_coordinate_data(self, initial_data):
        longitudes = initial_data[:, 0]
        latitudes = initial_data[:, 1]
        top_boundary = initial_data[:, 5]
        bottom_boundary = initial_data[:, 2]
        max_lon = self.__get_max(longitudes)
        min_lon = self.__get_min(longitudes)
        max_lat = self.__get_max(latitudes)
        min_lat = self.__get_min(latitudes)
        max_z = self.__get_max(top_boundary)
        min_z = self.__get_min(bottom_boundary)
        coordinate_data = {
            'longitudes': longitudes,
            'latitudes': latitudes,
            'max_lon': max_lon,
            'min_lon': min_lon,
            'max_lat': max_lat,
            'min_lat': min_lat,
            'max_z': max_z,
            'min_z': min_z
        }
        return DotMap(coordinate_data)

    def __set_geometric_data(self, initial_data):
        z_slab_lab = initial_data[:, 2]
        z_moho = initial_data[:, 3]
        z_icd = initial_data[:, 4]
        z_topo = initial_data[:, 5]
        geometric_data = {
            'z_slab_lab': z_slab_lab,
            'z_moho': z_moho,
            'z_icd': z_icd,
            'z_topo': z_topo
        }
        return DotMap(geometric_data)

    def get_cs_data(self):
        cs_data = {
            'coordinate_data': self.coordinate_data,
            # geometric data to get relevant_area
            'geometric_data': self.get_gm_data().geometric_data,
            'coast': self.coast
        }
        return DotMap(cs_data)

    def get_gm_data(self):
        gm_data = {
            'geometric_data': self.geometric_data,
            'slab_lab_areas': self.slab_lab_areas
        }
        return DotMap(gm_data)

    def get_tm_data(self):
        tm_data = {
            't_input': self.t_input,
            'trench_age': self.trench_age
        }
        return DotMap(tm_data)

    def get_mm_data(self):
        mm_data = {
            'm_input': self.m_input,
            'rheologic_data': self.rheologic_data
        }
        return DotMap(mm_data)

class CoordinateSystem(object):

    @staticmethod
    def round_to_step(x, step=1, prec=0):
        return (step * (np.array(x) / step).round()).round(prec)

    def __init__(self, cs_data, xy_step, z_step, z_precision=0):
        self.data = cs_data.coordinate_data
        self.geometries_not_rounded = cs_data.geometric_data
        self.coast = cs_data.coast
        self.xy_step = xy_step
        self.z_step = z_step
        self.z_precision = z_precision
        self.x_axis = self.__set_axis(self.data.max_lon, self.data.min_lon,
                                      self.xy_step)
        self.y_axis = self.__set_axis(self.data.max_lat, self.data.min_lat,
                                      self.xy_step, revert=True)
        self.z_axis = self.__set_axis(self.data.max_z, self.data.min_z,
                                      self.z_step, revert=True,
                                      precision=z_precision)
        self.grid_2D = self.__set_grid([self.x_axis, self.y_axis])
        self.grid_3D = self.__set_grid([self.x_axis, self.y_axis, self.z_axis])
        self.populated_area = self.__set_populated_area()
        self.sa_plate_area = self.__set_sa_plate_area()
        self.relevant_area = self.__set_relevant_area(
            [self.populated_area, self.sa_plate_area])
        self.relevant_volume = self.__set_relevant_volume(self.relevant_area)
        self.continent_area = self.__set_continent_area()
        self.relevant_area_eet = self.__set_relevant_area(
            [self.populated_area, self.sa_plate_area, self.continent_area])
        self.relevant_volume_eet = self.__set_relevant_volume(
            self.relevant_area_eet)

    def __set_axis(self, max_v, min_v, step, revert=False, precision=None):
        if precision is not None:
            max_v = self.round_to_step(max_v, step=step, prec=precision)
            min_v = self.round_to_step(min_v, step=step, prec=precision)
        else:
            max_v = np.ceil(max_v)
            min_v = np.floor(min_v)
        if revert is True:
            first_v, last_v = max_v, min_v
        else:
            first_v, last_v = min_v, max_v
        axis = np.linspace(first_v, last_v,
                           num=int(round((abs(last_v-first_v))/step+1)),
                           endpoint=True)
        if precision is not None:
            axis = self.round_to_step(axis, step=step, prec=precision)
        return axis

    def __set_grid(self, axes, mask=False):
        grid = np.meshgrid(*[n for n in axes], indexing='ij')
        return grid

    def __set_populated_area(self):
        slab_lab_valid = np.invert(np.isnan(self.reshape_data(
            self.geometries_not_rounded.z_slab_lab)))
        moho_valid = np.invert(np.isnan(self.reshape_data(
            self.geometries_not_rounded.z_moho)))
        icd_valid = np.invert(np.isnan(self.reshape_data(
            self.geometries_not_rounded.z_icd)))
        topo_valid = np.invert(np.isnan(self.reshape_data(
            self.geometries_not_rounded.z_topo)))
        populated_area = np.ones(self.get_2D_shape(), dtype=bool)
        populated_area[slab_lab_valid == 0] = 0
        populated_area[moho_valid == 0] = 0
        populated_area[icd_valid == 0] = 0
        populated_area[topo_valid == 0] = 0
        return populated_area

    def __set_sa_plate_area(self):
        """False (0) in Nazca Plate, True (1) in South American Plate"""
        slab_lab = self.reshape_data(self.geometries_not_rounded.z_slab_lab)
        g = np.gradient(slab_lab, axis=0)
        with np.errstate(invalid='ignore'):  # error_ignore
            high_g = np.absolute(g) > 1  # type: np.ndarray
        trench_start = np.argmax(high_g, axis=0)  # gets first true value
        i_idx = self.get_2D_indexes()[0]
        sa_plate_area = np.ones(self.get_2D_shape(), dtype=bool)
        sa_plate_area[i_idx < trench_start] = 0
        return sa_plate_area

    def __set_continent_area(self):
        lon = self.get_2D_grid()[0]
        continent_area = np.zeros(self.get_2D_shape())
        continent_area[lon >= self.coast[:, 0]] = 1
        return continent_area

    def __set_relevant_area(self, relevant_areas):
        relevant_area = np.ones(self.get_2D_shape(), dtype=bool)
        for area in relevant_areas:
            relevant_area[area == 0] = 0
        return relevant_area

    def __set_relevant_volume(self, relevant_area):
        relevant_volume = np.zeros(self.get_3D_shape(), dtype=bool)
        relevant_volume[:, :, :] = relevant_area[:, :, np.newaxis]
        return relevant_volume

    def round_depth_data(self, data):
        data = self.round_to_step(data, step=self.z_step, prec=self.z_precision)
        return data

    def reshape_data(self, data_column):
        return data_column.T.reshape(len(self.y_axis), len(self.x_axis)).T

    def get_xy_step(self):
        return self.xy_step

    def get_z_step(self):
        return self.z_step

    def get_axes(self):
        return [self.x_axis, self.y_axis, self.z_axis]

    def get_x_axis(self):
        return self.x_axis

    def get_y_axis(self):
        return self.y_axis

    def get_z_axis(self):
        return self.z_axis

    def get_xy_step(self):
        return self.xy_step

    def get_z_step(self):
        return self.z_step

    def get_2D_shape(self):
        return len(self.x_axis), len(self.y_axis)

    def get_3D_shape(self):
        return len(self.x_axis), len(self.y_axis), len(self.z_axis)

    def get_2D_grid(self):
        grid_2D = []
        for n in range(len(self.grid_2D)):
            grid_2D.append(SpatialArray2D(self.grid_2D[n],
                                          self).mask_irrelevant())
        return grid_2D

    #def get_3D_grid(self):
    #    grid_3D = []
    #    for n in range(len(self.grid_3D)):
    #        grid_3D.append(SpatialArray3D(self.grid_3D[n],
    #                                      self).mask_irrelevant())
    #    return grid_3D

    def get_3D_grid(self, masked=True):
        grid_3D = []
        for n in range(len(self.grid_3D)):
            n_grid = SpatialArray3D(self.grid_3D[n], self)
            if masked is True:
                n_grid = n_grid.mask_irrelevant()
            grid_3D.append(n_grid)
        return grid_3D

    #def get_2D_indexes(self):
    #    indexes_arrays = np.indices(self.get_2D_shape())
    #    indexes_2D = []
    #    for n in range(len(indexes_arrays)):
    #        indexes_2D = SpatialArray2D(indexes_arrays[n], self)
    #    return indexes_2D

    def get_2D_indexes(self):
        return np.indices(self.get_2D_shape())

    #def get_3D_indexes(self):
    #    indexes_arrays = np.indices(self.get_3D_shape())
    #    indexes_3D = []
    #    for n in range(len(indexes_arrays)):
    #        indexes_3D = SpatialArray3D(indexes_arrays[n], self)
    #    return indexes_3D

    def get_3D_indexes(self):
        return np.indices(self.get_3D_shape())

    def get_relevant_area(self):
        return SpatialArray2D(self.relevant_area, self)

    def get_relevant_volume(self):
        return SpatialArray3D(self.relevant_volume, self)

    def get_relevant_area_eet(self):
        return SpatialArray2D(self.relevant_area_eet, self)

    def get_relevant_volume_eet(self):
        return SpatialArray2D(self.relevant_volume_eet, self)


class SpatialArray(np.ndarray):

    def __new__(cls, input_array, coordinate_system=None):
        # print('new method of SpatialArray called')
        # print('cls is', cls)
        obj = np.asarray(input_array).view(cls)
        #Returns unmodified numpy array
        #obj.array = np.asarray(input_array).view(np.ndarray)
        obj.cs = coordinate_system
        return obj

    def __array_finalize__(self, obj):
        # print('__array_finalize__ method of SpatialArray called')
        # print('self type is', type(self))
        # print('obj type is', type(obj))
        if obj is None:
            return
        self.cs = getattr(obj, 'cs', None)
        #self.array = np.asarray(self).view(np.ndarray)
        #self.array = getattr(obj, 'array', None)

    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        reduce_tuple = super(SpatialArray, self).__reduce__()
        # Create our own state to pass to __setstate__
        new_state = reduce_tuple[2] + (self.cs,)
        # Return a tuple that replaces the parent's __setstate__ tuple
        return (reduce_tuple[0], reduce_tuple[1], new_state)

    def __setstate__(self, state):
        self.cs = state[-1]
        # Call the parent's __setstate__ with the other tuple elements
        super(SpatialArray, self).__setstate__(state[0:-1])

    @staticmethod
    def divide_array_by_areas(array, areas):
        if len(array.shape) == 2:
            array = array[:, :, np.newaxis]
        array_1 = array.copy()
        array_2 = array.copy()
        areas_3D = np.repeat(areas[:, :, np.newaxis], array.shape[2],
                             axis=2)
        array_1[np.invert(areas_3D)] = np.nan
        array_2[areas_3D] = np.nan
        if array.shape[2] == 1:
            array_1 = np.squeeze(array_1, axis=2)
            array_2 = np.squeeze(array_2, axis=2)
        return array_1, array_2

    @staticmethod
    def combine_arrays_by_areas(array_1, array_2, areas):
        if len(array_1.shape) == 2:
            array_1 = array_1[:, :, np.newaxis]
            array_2 = array_2[:, :, np.newaxis]
        combined_array = array_1.copy()
        areas_3D = np.repeat(areas[:, :, np.newaxis], array_1.shape[2],
                             axis=2)
        combined_array[np.invert(areas_3D)] = array_2[np.invert(areas_3D)]
        if combined_array.shape[2] == 1:
            combined_array = np.squeeze(combined_array, axis=2)
        return combined_array

    @staticmethod
    def get_values(array):
        if isinstance(array, (SpatialArray2D, SpatialArray3D)):
            array = array.get_array()
        return array


class SpatialArray2D(SpatialArray):

    def __new__(cls, input_array, coordinate_system=None):
        #print('__new__ method of SpatialArray2D called')
        #print('cls is', cls)
        obj = super(SpatialArray2D, cls).__new__(cls, input_array, coordinate_system)
        return obj

    def get_array(self):
        return self

    def mask_irrelevant(self, nan_fill=False):
        mask = np.invert(self.cs.get_relevant_area())
        # TODO: figure out if this .copy() is necessary
        masked_array = self.copy()
        masked_array[mask] = np.nan
        #masked_array = np.ma.array(self, mask=mask)
        #if nan_fill is True:
        #    masked_array = masked_array.filled(np.nan)
        return masked_array

    def mask_irrelevant_eet(self, nan_fill=False):
        mask = np.invert(self.cs.get_relevant_area_eet())
        masked_array = self.copy()
        masked_array[mask] = np.nan
        return masked_array

    def cross_section(self, latitude=None, longitude=None,
                      lat1=None, lon1=None, lat2=None, lon2=None):
        if self.ndim != 2:
            dims = self.ndim
            error = ("Array with 2 dimensions expected,",
                     " array with {} dimensions recieved").format(dims)
            raise ValueError(error)
        elif latitude:
            #index = list(self.cs.get_y_axis()).index(latitude)
            index = np.where(np.isclose(self.cs.get_y_axis(), latitude))[0][0]
            cross_section = self[:, index]
            return cross_section
        elif longitude:
            #index = list(self.cs.get_x_axis()).index(longitude)
            index = np.where(np.isclose(self.cs.get_x_axis(), longitude))[0][0]
            print("index:", index)
            cross_section = self[index, :]
        elif lat1 and lon1 and lat2 and lon2:
            cross_section = self
        return cross_section

    def extract_point(self, latitude=None, longitude=None):
        if self.ndim != 2:
            dims = self.ndim
            error = ("Array with 2 dimensions expected,",
                     " array with {} dimensions recieved").format(dims)
            raise ValueError(error)
        else:
            lat_index = np.where(np.isclose(self.cs.get_y_axis(), latitude))[0][0]
            lon_index = np.where(np.isclose(self.cs.get_x_axis(), longitude))[0][0]
            point = self[lon_index, lat_index]
        return point

    def divide_by_areas(self, areas):
        array_1, array_2 = super().divide_array_by_areas(self, areas)
        return array_1, array_2

    def replace_where_area_is_false(self, replacement_array, areas):
        combined_array = super().combine_arrays_by_areas(self,
                                                         replacement_array,
                                                         areas)
        return combined_array


class SpatialArray3D(SpatialArray):

    def __new__(cls, input_array, coordinate_system=None):
        #print('__new__ method of SpatialArray3D called')
        #print('cls is', cls)
        obj = super(SpatialArray3D, cls).__new__(cls, input_array, coordinate_system)
        return obj

    def mask_irrelevant(self, nan_fill=True):
        mask = np.invert(self.cs.get_relevant_volume())
        # TODO: figure out if this .copy() is necessary
        masked_array = self.copy()
        masked_array[mask] = np.nan
        #masked_array = np.ma.array(self, mask=mask)
        #if nan_fill is True:
        #    masked_array = masked_array.filled(np.nan)
        return masked_array

    def mask_irrelevant_eet(self, nan_fill=True):
        mask = np.invert(self.cs.get_relevant_volume_eet())
        masked_array = self.copy()
        masked_array[mask] = np.nan
        return masked_array

    def extract_surface(self, z_2D):
        z_2D = z_2D
        z_3D = self.cs.get_3D_grid()[2]
        # round to step to fix any floating precision errors
        a = z_3D != self.cs.round_depth_data(z_2D[:, :, np.newaxis])
        array_3D = self.copy()
        array_3D[a] = 0
        surface = np.sum(array_3D, axis=2)
        #surface2 = surface.copy()
        #surface2[surface2 == 0] = np.nan #Returns NaNs when input array has zeros
        # innecesary when using mask_irrelevant() down below
        surface = SpatialArray2D(surface, self.cs).mask_irrelevant()
        #np.testing.assert_equal(surface, surface2)
        return surface

    def crop(self, top_z='top', bottom_z='bottom'):
        z_3D = self.cs.get_3D_grid()[2]
        if isinstance(top_z, str):
            top_z = np.inf
        else:
            # round to step to fix any floating precision errors
            top_z = self.cs.round_depth_data(top_z[:, :, np.newaxis])
        if isinstance(bottom_z, str):
            bottom_z = -np.inf
        else:
            # round to step to fix any floating precision errors
            bottom_z = self.cs.round_depth_data(bottom_z[:, :, np.newaxis])
        with np.errstate(invalid='ignore'): # error_ignore
            a = top_z < z_3D
            b = bottom_z >= z_3D
        array_3D = self.copy()
        array_3D[a] = np.nan
        array_3D[b] = np.nan
        return array_3D

    def cross_section(self, latitude=None, longitude=None,
                      lat1=None, lon1=None, lat2=None, lon2=None):
        if self.ndim != 3:
            dims = self.ndim
            error = ("Array with 3 dimensions expected,",
                     " array with {} dimensions recieved").format(dims)
            raise ValueError(error)
        elif latitude:
            #index = list(self.cs.get_y_axis()).index(latitude)
            index = np.where(np.isclose(self.cs.get_y_axis(), latitude))[0][0]
            cross_section = self[:, index, :]
        elif longitude:
            #index = list(self.cs.get_x_axis()).index(longitude)
            index = np.where(np.isclose(self.cs.get_x_axis(), longitude))[0][0]
            cross_section = self[:, index, :]
        elif lat1 and lon1 and lat2 and lon2:
            cross_section = self
        return SpatialArray2D(cross_section, self.cs)

    def point_depth_profile(self, latitude=None, longitude=None):
        if self.ndim != 3:
            dims = self.ndim
            error = ("Array with 3 dimensions expected,",
                     " array with {} dimensions recieved").format(dims)
            raise ValueError(error)
        else:
            lat_index = np.where(np.isclose(self.cs.get_y_axis(), latitude))[0][0]
            lon_index = np.where(np.isclose(self.cs.get_x_axis(), longitude))[0][0]
            point_depth_profile = self[lon_index, lat_index, :]
        return point_depth_profile

    def divide_by_areas(self, areas):
        array_1, array_2 = super().divide_array_by_areas(self, areas)
        return array_1, array_2

    def replace_where_area_is_false(self, replacement_array, areas):
        combined_array = super().combine_arrays_by_areas(self,
                                                         replacement_array,
                                                         areas)
        return combined_array

    def get_array(self):
        return self


class GeometricModel(object):

    def __init__(self, gm_data, coordinate_system):
        self.data = gm_data.geometric_data
        self.slab_lab_areas = gm_data.slab_lab_areas
        self.cs = coordinate_system
        self.slab_lab = self.__set_boundary(self.data.z_slab_lab)
        self.moho = self.__set_boundary(self.data.z_moho)
        self.icd = self.__set_boundary(self.data.z_icd)
        self.topo = self.__set_boundary(self.data.z_topo)
        self.geo_model_3D = self.__set_3D_geometric_model()
        #TODO: Quizas mover metodos para obtener slab_lab_area a CoordinateSys.
        # Para generar un archivo de areas preliminares
        # self.slab_lab_int_index = self.__set_slab_lab_int_index()
        # self.slab_lab_int_area = self.__set_slab_lab_int_area(save=True)
        # Para utilizar un archivo de areas ya generado (areas.dat):
        self.slab_lab_int_index = self.__set_slab_lab_int_index_from_areas()
        self.slab_lab_int_area = self.__set_slab_lab_int_area()
        self.slab_lab_int_depth = self.__set_slab_lab_int_depth()
        self.slab_lab_int_topo = self.__set_slab_lab_int_topo()

    def __set_boundary(self, geometric_data):
        return SpatialArray2D(
            self.cs.reshape_data(self.cs.round_depth_data(geometric_data)),
            self.cs)

    def __set_3D_geometric_model(self):
        boundaries = self.get_boundaries()
        z = self.cs.get_3D_grid()[2]
        geo_model_3D = np.empty(z.shape)*np.nan
        for n in range(len(boundaries)):
            c = boundaries[n][:, :, np.newaxis]
            if n < (len(boundaries)-1):
                a = n + 1
            else:
                a = np.nan
            with np.errstate(invalid='ignore'): # error_ignore
                geo_model_3D[z <= c] = a
        return SpatialArray3D(geo_model_3D, self.cs).mask_irrelevant()

    def __set_layer_thickness(self):
        pass

    def __set_slab_lab_int_index(self):
        # Slab/Lab gradients array
        g = np.gradient(
            self.cs.reshape_data(self.cs.geometries_not_rounded.z_slab_lab),
            axis=0)
        # Min gradient Indexes
        min_g_idx = np.nanargmin(np.ma.masked_invalid(g), axis=0)
        # Positive gradient boolean array (True where gradient > 0 ...)
        with np.errstate(invalid='ignore'): # error_ignore
            pos_g = g > -4e-01 # (... or where gradient gets very close to 0)
        # Ignore (mask) the values to the left of the minimum gradient
        i_idx = self.cs.get_2D_indexes()[0]
        g_mask = i_idx < min_g_idx
        masked_pos_g = np.ma.array(pos_g, mask=g_mask)
        # Get index of first positive gradient (Slab/Lab intersection idx)
        sli_idx = np.argmax(masked_pos_g, axis=0)
        # Get visual boolean representation of sli_idx
        # Z = np.zeros(g.shape)
        # Z[sli_idx, self.cs.get_2D_indexes()[1][0]] = 1
        # return Z.T
        return sli_idx

    def __set_slab_lab_int_area(self, save=False):
        sli_idx = self.slab_lab_int_index
        i_idx = self.cs.get_2D_indexes()[0]
        sli_area = np.zeros(self.cs.get_2D_shape())
        # TODO: this moves sli_area 1 index to the right compared with areas.dat
        # is it necessary?
        sli_area[i_idx > sli_idx] = 1
        # TODO: Is sli_ara a typo? It is not returned.
        sli_ara = np.asarray(sli_area, dtype=bool)
        if save == True:
            # Agregar latitudes en la primera columna:
            # latitudes = self.cs.get_y_axis()
            # sli_ara = np.r_[latitudes[np.newaxis, :], sli_ara]
            # np.savetxt('data/preliminar_areas.dat', sli_ara.T, fmt='%2d')
            np.savetxt('data/preliminar_areas.dat', sli_ara.T, fmt='%.0f')
        return sli_area

    def __set_slab_lab_int_index_from_areas(self):
        sli_idx = np.argmax(self.slab_lab_areas, axis=1)
        # sli_idx = sli_idx - 1
        # TODO: Check what is this next line for?
        sli_idx[-1] = 0
        return sli_idx

    def __set_slab_lab_int_depth(self):
        sli_idx = self.slab_lab_int_index
        j_idx = self.cs.get_2D_indexes()[1][0]
        sli_depth = self.slab_lab[sli_idx, j_idx]
        return sli_depth

    def __set_slab_lab_int_topo(self):
        sli_idx = self.slab_lab_int_index
        j_idx = self.cs.get_2D_indexes()[1][0]
        sli_topo = self.topo[sli_idx, j_idx]
        return sli_topo

    def set_layer_property(self, cs, ci, ml):
        r = np.ones(self.geo_model_3D.shape) * np.nan
        if isinstance(cs, np.ndarray) and np.ndim(cs)==1:
            r = np.repeat(r[:,:,:,np.newaxis], len(cs), axis=3)
        r[self.geo_model_3D == 1] = cs
        r[self.geo_model_3D == 2] = ci
        r[self.geo_model_3D == 3] = ml
        return SpatialArray3D(r, self.cs)

    def get_boundaries(self):
        boundaries = [self.topo, self.icd, self.moho, self.slab_lab]
        return boundaries

    def get_topo(self):
        return self.topo

    def get_icd(self):
        return self.icd

    def get_moho(self):
        return self.moho

    def get_slab_lab(self):
        return self.slab_lab

    def get_3D_geometric_model(self):
        return self.geo_model_3D

    def get_slab_lab_int_index(self):
        return self.slab_lab_int_index

    def get_slab_lab_int_depth(self):
        return self.slab_lab_int_depth

    def get_slab_lab_int_topo(self):
        return self.slab_lab_int_topo

    def get_slab_lab_int_area(self):
        return self.slab_lab_int_area


class ThermalModel(object):

    @staticmethod
    def __calc_lab_temp(tp, g, depth):
        lab_temp = tp + g * abs(depth)*1e3
        return lab_temp

    @staticmethod
    def __calc_q_zero(k, tp, kappa, age):
        q_zero = (k * tp)/np.sqrt(np.pi * kappa * age)
        return q_zero

    @staticmethod
    def __calc_s(depth, topo, kappa, v, dip, b):
        s = 1. + (b * np.sqrt(((abs(depth-topo)*1.e3)
                              * v * abs(np.sin(dip)))/kappa))
        return s

    @staticmethod
    def __calc_slab_lab_int_sigma(sli_depth, sli_topo, sli_temp, sli_s, sli_k,
                                  sli_q_zero, v):
        sli_sigma = ((sli_temp * sli_s * sli_k)
                     / (v * abs(sli_depth-sli_topo)*1.e3)) - (sli_q_zero/v)
        return sli_sigma

    @staticmethod
    def __calc_slab_sigma(depth, topo, sli_depth, sli_topo, sli_sigma, d):
        mu = sli_sigma / (1. - np.exp(d))
        # print(abs(depth-topo)*1.e3 * d/ abs())
        slab_sigma = mu * (1. - np.exp(abs(depth-topo)*1.e3 * d
                                       / (abs(sli_depth-sli_topo)*1.e3)))
        return slab_sigma

    @staticmethod
    def __calc_slab_temp(depth, topo, q_zero, slab_sigma, v, k, s):
        slab_temp = ((q_zero + slab_sigma * v) * abs(depth-topo)*1.e3
                     / (k * s))
        return slab_temp

    @staticmethod
    def __calc_geotherm(h, delta, k, z, z_topo, z_sl, temp_sl):
        # TODO: optimize memory usage of this function
        #z = z*1.e3
        #z_topo = z_topo*1.e3
        #z_sl = z_sl*1.e3
        #base_temp = temp_sl-((h*delta**2)/k)*(np.exp(z_topo/delta)
        #                                      - np.exp(z_sl/delta))
        #rad_temp = ((h*delta**2)/k)*(np.exp(z_topo/delta)-np.exp(z/delta))
        #geotherm = rad_temp + (abs(z-z_topo)/abs(z_sl-z_topo))*base_temp

        geotherm = ((((h*delta**2)/k)*(np.exp((z_topo*1.e3)/delta)-np.exp((z*1.e3)/delta)))
                    +(abs((z*1.e3)-(z_topo*1.e3))/abs((z_sl*1.e3)-(z_topo*1.e3)))
                    * (temp_sl-((h*delta**2)/k)*(np.exp((z_topo*1.e3)/delta)
                                                 - np.exp((z_sl*1.e3)/delta))))
        return geotherm

    @staticmethod
    def __calc_surface_heat_flow(h, delta, k, z_topo, z_sl, temp_sl):
        z_topo = z_topo*1.e3
        z_sl = z_sl*1.e3
        base_temp = temp_sl-((h*delta**2)/k)*(np.exp(z_topo/delta)
                                              - np.exp(z_sl/delta))
        heat_flow = (-h*delta-k/(abs(z_sl-z_topo))*base_temp)
        #heat_flow2 = -(k*temp_sl)/abs(z_sl-z_topo) - (h*delta) + ((h*delta**2)/abs(z_sl-z_topo))*(np.exp(z_topo/delta)-np.exp(z_sl/delta))
        #if np.allclose(heat_flow, heat_flow2, equal_nan=True):
        #    print('Arrays are equal')
        return heat_flow

    @staticmethod
    def __calc_surface_heat_flow_tassa(h, delta, k, z_topo, z_sl, temp_sl):
        z_topo = z_topo*1.e3
        z_sl = z_sl*1.e3
        heat_flow = -k*(temp_sl/z_sl)- (h*delta)*(np.exp(z_topo/delta)-np.exp(z_sl/delta))

    @staticmethod
    def __calc_oceanic_plate_thickness(kappa, t):
        op_thickness = (2.32*np.sqrt(kappa*t))/1.e3
        return op_thickness

    def __init__(self, tm_data, geometric_model, coordinate_system):
        self.geo_model = geometric_model
        self.cs = coordinate_system
        self.vars = DotMap(self.__set_variables(tm_data.t_input,
                                                 tm_data.trench_age))
        self.slab_lab_temp = self.__set_slab_lab_temp()
        self.surface_heat_flow = self.__set_surface_heat_flow()
        self.geotherm = self.__set_geotherm()

    def __set_k_or_h(self, t_input, prop):
        boundaries = self.geo_model.get_boundaries()
        prop_v = None
        prop_fz = None
        if prop == 'k':
            prop_v = [t_input.k_cs, t_input.k_ci, t_input.k_ml]
            prop_fz = t_input.k_z
        elif prop == 'h':
            prop_v = [t_input.H_cs, t_input.H_ci, t_input.H_ml]
            prop_fz = t_input.H_z
        if prop_fz is True:
            r = self.geo_model.set_layer_property(prop_v[0], prop_v[1],
                                                  prop_v[2])
        else:
            rh = []
            h = []
            n = 0
            while n < len(boundaries)-1:
                h.append(boundaries[n] - boundaries[n+1])
                # Erase negative values due to data errors (e.g. z_moho < z_sl)
                h[n][h[n] < 0] = 0
                rh.append(prop_v[n]*h[n])
                n += 1
            r = (sum(rh) / sum(h))[:, :, np.newaxis]
            r = np.repeat(r, self.cs.get_3D_shape()[2], axis=2)
            geo_model_mask = np.isnan(self.geo_model.get_3D_geometric_model())
            #r = SpatialArray3D(r, self.cs).mask_irrelevant(irrelevant=geo_model_mask)
            r = np.ma.array(r, mask=geo_model_mask).filled(np.nan)
        return SpatialArray3D(r, self.cs)

    def __set_delta(self, t_input):
        boundaries = self.geo_model.get_boundaries()
        if t_input.delta_icd is True:
            delta = boundaries[0] - boundaries[1]
        else:
            delta = (t_input.delta)
        return delta

    def __set_trench_age(self, trench_age, t_input):
        if t_input.t_lat is True:
            trench_age = trench_age[:, 1]
        else:
            trench_age = t_input.t
        return trench_age

    def __set_variables(self, t_input, trench_age):
        t_input = DotMap(t_input)
        t_input_2 = t_input.copy()
        t_input_2['k_z'] = False
        t_input_2['H_z'] = False
        # TODO: find a way to save memory by not storing 3D arrays as k and h
        t_vars = {
            'k_cs': t_input['k_cs'],
            'k_ci': t_input['k_ci'],
            'k_ml': t_input['k_ml'],
            'k': self.__set_k_or_h(t_input, 'k'),
            'k_prom': self.__set_k_or_h(t_input_2, 'k'),
            'h_cs': t_input['H_cs'],
            'h_ci': t_input['H_ci'],
            'h_ml': t_input['H_ml'],
            'h': self.__set_k_or_h(t_input, 'h'),
            'h_prom': self.__set_k_or_h(t_input_2, 'h'),
            'delta': self.__set_delta(t_input)*1.e3,
            'tp': t_input['Tp'],
            'g': t_input['G'],
            'kappa': t_input['kappa'],
            'v': t_input['V']/(1.e6*365.*24.*60.*60.),
            'dip': t_input['dip'],
            'b': t_input['b'],
            't': self.__set_trench_age(trench_age, t_input)*(1.e6*365.*24.*60.*60.),
            'd': t_input['D']
        }
        return t_vars

    def __get_slab_lab_int_temperature(self):
        sli_depth = self.geo_model.get_slab_lab_int_depth()
        tp = self.vars.tp
        g = self.vars.g
        sli_temp = self.__calc_lab_temp(tp, g, sli_depth)
        return sli_temp

    def __get_slab_lab_int_sigma(self):
        sli_depth = self.geo_model.get_slab_lab_int_depth()
        sli_topo = self.geo_model.get_slab_lab_int_topo()
        sli_temp = self.__get_slab_lab_int_temperature()
        kappa = self.vars.kappa
        v = self.vars.v
        dip = self.vars.dip
        b = self.vars.b
        sli_s = self.__calc_s(sli_depth, sli_topo, kappa, v, dip, b)
        z_sl = self.geo_model.get_slab_lab()
        slab_k = self.vars.k.extract_surface(z_sl+self.cs.z_step)
        sli_idx = self.geo_model.get_slab_lab_int_index()
        j_idx = self.cs.get_2D_indexes()[1][0]
        sli_k = slab_k[sli_idx, j_idx]
        tp = self.vars.tp
        t = self.vars.t
        sli_q_zero = self.__calc_q_zero(sli_k, tp, kappa, t)
        sli_sigma = self.__calc_slab_lab_int_sigma(sli_depth, sli_topo,
                                                   sli_temp, sli_s, sli_k,
                                                   sli_q_zero, v)
        return sli_sigma

    def __get_slab_sigma(self):
        z_sl = self.geo_model.get_slab_lab()
        z_topo = self.geo_model.get_topo()
        sli_depth = self.geo_model.get_slab_lab_int_depth()
        sli_topo = self.geo_model.get_slab_lab_int_topo()
        sli_sigma = self.__get_slab_lab_int_sigma()
        d = self.vars.d
        slab_sigma = self.__calc_slab_sigma(z_sl, z_topo, sli_depth, sli_topo,
                                            sli_sigma, d)
        a = self.geo_model.get_slab_lab_int_area()
        b = a == 1
        slab_sigma = np.ma.array(slab_sigma, mask=b).filled(np.nan)
        return slab_sigma

    def __get_slab_temp(self):
        z_sl = self.geo_model.get_slab_lab()
        z_topo = self.geo_model.get_topo()
        slab_k = self.vars.k.extract_surface(z_sl + self.cs.z_step)
        tp = self.vars.tp
        kappa = self.vars.kappa
        v = self.vars.v
        t = self.vars.t
        q_zero = self.__calc_q_zero(slab_k, tp, kappa, t)
        dip = self.vars.dip
        b = self.vars.b
        s = self.__calc_s(z_sl, z_topo, kappa, v, dip, b)
        slab_sigma = self.__get_slab_sigma()
        slab_temp = self.__calc_slab_temp(z_sl, z_topo, q_zero, slab_sigma, v,
                                          slab_k, s)
        a = self.geo_model.get_slab_lab_int_area()
        b = a == 1
        slab_temp = np.ma.array(slab_temp, mask=b).filled(np.nan)
        return slab_temp

    def __get_lab_temp(self):
        z_sl = self.geo_model.get_slab_lab()
        # z_lab = self.geo_model.mask_slab(z_sl)
        tp = self.vars.tp
        g = self.vars.g
        lab_temp = self.__calc_lab_temp(tp, g, z_sl)
        a = self.geo_model.get_slab_lab_int_area()
        b = a == 0
        lab_temp = np.ma.array(lab_temp, mask=b).filled(np.nan)
        return lab_temp

    def __set_slab_lab_temp(self):
        a = self.geo_model.get_slab_lab_int_area()
        b = a == 0
        c = a == 1
        slab_temp = self.__get_slab_temp()  # type: np.ndarray
        lab_temp = self.__get_lab_temp()  # type: np.ndarray
        slab_lab_temp = np.zeros(self.cs.get_2D_shape())
        slab_lab_temp[b] = slab_temp[b]
        slab_lab_temp[c] = lab_temp[c]
        return slab_lab_temp

    def __set_surface_heat_flow(self):
        delta = self.vars.delta
        z_topo = self.geo_model.get_topo()
        z_sl = self.geo_model.get_slab_lab()
        #slab_lab_k = self.vars.k.extract_surface(z_sl+self.cs.z_step)
        #slab_lab_h = self.vars.h.extract_surface(z_sl+self.cs.z_step)
        topo_k = self.vars.k.extract_surface(z_topo-1)
        topo_h = self.vars.h.extract_surface(z_topo-1)
        temp_sl = self.slab_lab_temp
        heat_flow = self.__calc_surface_heat_flow(topo_h, delta, topo_k,
                                                  z_topo, z_sl, temp_sl)
        return heat_flow

    def __set_geotherm(self):
        k = self.vars.k
        h = self.vars.h
        if isinstance(self.vars.delta, float):
            delta = self.vars.delta
        else:
            delta = self.vars.delta[:, :, np.newaxis]
        z = self.cs.get_3D_grid()[2]
        z_topo = self.geo_model.get_topo()[:, :, np.newaxis]
        z_sl = self.geo_model.get_slab_lab()[:, :, np.newaxis]
        temp_sl = self.slab_lab_temp[:, :, np.newaxis]
        geotherm = self.__calc_geotherm(h, delta, k, z, z_topo, z_sl, temp_sl)
        return geotherm

    def get_surface_heat_flow(self, format='negative watts'):
        surface_heat_flow = SpatialArray2D(self.surface_heat_flow, self.cs)
        if format == 'positive milliwatts':
            surface_heat_flow = surface_heat_flow*-1*1.e3
        return surface_heat_flow

    def get_oceanic_plate_thickness(self):
        trench_age = self.vars.t
        kappa = self.vars.kappa
        op_thickness = self.__calc_oceanic_plate_thickness(kappa,trench_age)
        return op_thickness

    def get_geotherm(self):
        return SpatialArray3D(self.geotherm, self.cs)

    def get_geometric_model(self):
        return self.geo_model

    def get_coordinate_system(self):
        return self.cs

class MechanicModel(object):

    @staticmethod
    def __calc_brittle_yield_strength(bs, depth):
        bys = bs*(depth*1.e3)*1.e-6
        return bys

    @staticmethod
    def calc_ductile_yield_strength(es, n, a, h, r, temp):
        with np.errstate(over='ignore'):
            dys = (es/a)**(1/n)*np.exp(h/(n*r*(temp+273.15)))*1.e-6
        #result = np.empty(temp.shape)
        #stored_value = np.empty(temp.shape)
        #np.add(temp, 273.15, result)
        #np.multiply(r, result, result)
        #np.multiply(n, result, result)
        #np.divide(h, result, result)
        #np.exp(result, result)
        #np.multiply(1.e-6, result, stored_value)
        #np.divide(es, a, result)
        #with np.errstate(over='ignore'):
        #    np.power(result, (1/n), result)
        #np.multiply(result, stored_value,result)
        #return result
        return dys

    @staticmethod
    def __calc_depth_from_brittle_yield_strength(bs, bys):
        depth = bys/(bs*1.e3*1.e-6)
        return depth

    @staticmethod
    def calc_temperature_from_ductile_yield_strength(es, n, a, h, r, dys):
        #temp = 1/(np.log(dys/(((es/a)**(1/n))*1.e-6))*((n*r)/h))
        temp = h/(n*r*np.log(dys/(1.e-6*(es/a)**(1/n))))
        return temp - 273.15

    @staticmethod
    def __calc_eet_attached(attached_ths):
        attached_ths[np.isnan(attached_ths)] = 0
        attached_eet = (attached_ths[:, :, 0]
                        + attached_ths[:, :, 1]
                        + attached_ths[:, :, 2])
        return attached_eet

    @staticmethod
    def __calc_eet_detached(detached_ths):
        detached_ths[np.isnan(detached_ths)] = 0
        detached_eet = (detached_ths[:, :, 0]**3
                        + detached_ths[:, :, 1]**3
                        + detached_ths[:, :, 2]**3)**(1/3)
        return detached_eet

    def __init__(self, mm_data, geo_model, thermal_model, coordinate_system):
        self.geo_model = geo_model
        self.cs = coordinate_system
        self.thermal_model = thermal_model
        self.vars = DotMap(self.__set_variables(mm_data.m_input,
                                                 mm_data.rheologic_data))
        self.depth_from_topo = self.__get_depth_from_topo()
        self.bys_t, self.bys_c = self.__get_brittle_yield_strength()
        self.yse_t, self.yse_c = self.__set_yield_strength_envelope()
        #self.eet, self.eet_calc_data = self.__set_eet(self.yse_t)
        self.eet, self.eet_wrong, self.eet_calc_data = self.__set_eet_2()
        self.integrated_strength = self.__set_integrated_strength()

    def __get_rheologic_vars_from_model(self, rock_id):
        rock = RheologicModel.objects.get(name=rock_id)
        rock_dic = {
            'n': rock.n,
            'a': rock.A,
            'h': rock.H,
        }
        return DotMap(rock_dic)

    def __get_rheologic_vars(self, id_rh, rhe_data):
        rock = rhe_data[str(id_rh)]
        rock_dic = {
            'h': rock.H,
            'n': rock.n,
            'a': rock.A
        }
        return DotMap(rock_dic)

    def __get_polyphase_rheologic_vars(self, f1, id_rh_1, id_rh_2, rhe_data):
        # Based on Tullis et al., 1991
        f2 = 1 - f1
        rock_1 = rhe_data[str(id_rh_1)]
        n1 = rock_1.n
        a1 = rock_1.A
        h1 = rock_1.H
        rock_2 = rhe_data[str(id_rh_2)]
        n2 = rock_2.n
        a2 = rock_2.A
        h2 = rock_2.H
        n = 10**(f1*np.log10(n1)+f2*np.log10(n2))
        a = 10**((np.log10(a2)*(n-n1)-np.log10(a1)*(n-n2))/(n2-n1))
        h = (h2*(n-n1) - h1*(n-n2))/(n2-n1)
        rock_dic = {
            'h': h,
            'n': n,
            'a': a
        }
        #print(rock_dic)
        return DotMap(rock_dic)

    def __set_variables(self, m_input, rhe_data):
        m_input = DotMap(m_input)
        m_vars = {
            'slm': m_input['slm'],
            'bs_t': m_input['Bs_t'],
            'bs_c': m_input['Bs_c'],
            'e': m_input['e'],
            'r': m_input['R'],
            's_max': m_input['s_max'],
            'cs': self.__get_rheologic_vars(m_input['Cs'], rhe_data),
            'ci': self.__get_rheologic_vars(m_input['Ci'], rhe_data),
            'ml': self.__get_rheologic_vars(m_input['Ml'], rhe_data),
            'mla': self.__get_polyphase_rheologic_vars(
                m_input['spct'], m_input['Serp'], m_input['Ml'], rhe_data)
        }
        return m_vars

    def __get_depth_from_topo(self):
        depth_from_topo = -(self.cs.get_3D_grid()[2]
                            - self.geo_model.get_topo()[:, :, np.newaxis])
        return depth_from_topo

    def __get_brittle_yield_strength(self):
        bs_t = self.vars.bs_t
        bs_c = self.vars.bs_c
        #depth = self.cs.get_3D_indexes()[2]  # needs to be Z from topo
        depth_from_topo = self.__get_depth_from_topo()
        geo_model_mask = np.isnan(self.geo_model.get_3D_geometric_model())
        depth_from_topo = np.ma.array(depth_from_topo, mask=geo_model_mask).filled(np.nan)
        depth_from_topo = SpatialArray3D(depth_from_topo
                                         .astype(np.float64, copy=False),
                                                 self.cs).mask_irrelevant()
        #depth = self.cs.mask_3d_array(depth.astype(np.float64),
        #                                          nan_fill=True)
        bys_t = self.__calc_brittle_yield_strength(bs_t, depth_from_topo)
        bys_c = self.__calc_brittle_yield_strength(bs_c, depth_from_topo)
        return bys_t, bys_c

    def __get_ductile_yield_strength(self):
        e = self.vars.e
        r = self.vars.r
        temp = self.thermal_model.get_geotherm()
        cs = self.vars.cs
        ci = self.vars.ci
        ml = self.vars.ml
        mla = self.vars.mla
        # TODO: find a way to not store these next 3D arrays in memory
        """
        n = self.geo_model.set_layer_property(cs.n, ci.n, ml.n)
        a = self.geo_model.set_layer_property(cs.a, ci.a, ml.a)
        h = self.geo_model.set_layer_property(cs.h, ci.h, ml.h)
        #Antearco
        n = self.geo_model.set_layer_property(cs.n, ci.n, mla.n)
        a = self.geo_model.set_layer_property(cs.a, ci.a, mla.a)
        h = self.geo_model.set_layer_property(cs.h, ci.h, mla.h)
        """
        # 4-Dimensional Array
        rheo_array_backarc = self.geo_model.set_layer_property(
            np.array([cs.n, cs.a, cs.h]),
            np.array([ci.n, ci.a, ci.h]),
            np.array([ml.n, ml.a, ml.h]))
        if self.vars.slm is True:
            rheo_array_forearc = self.geo_model.set_layer_property(
                np.array([cs.n, cs.a, cs.h]),
                np.array([ci.n, ci.a, ci.h]),
                np.array([mla.n, mla.a, mla.h]))
            rheo_array = SpatialArray.combine_arrays_by_areas(
                rheo_array_backarc, rheo_array_forearc,
                self.geo_model.get_slab_lab_int_area().astype(bool))
        else:
            rheo_array = rheo_array_backarc
        n = rheo_array[:,:,:,0]
        a = rheo_array[:,:,:,1]
        h = rheo_array[:,:,:,2]
        dys = self.calc_ductile_yield_strength(e, n, a, h, r, temp)
        return dys

    def __set_yield_strength_envelope(self):
        bys_t, bys_c = self.__get_brittle_yield_strength()  # type: np.ndarray
        dys = self.__get_ductile_yield_strength()  # type: np.ndarray
        with np.errstate(invalid='ignore'):
            yse_t = np.where(bys_t < dys, bys_t, dys)
            yse_c = np.where(bys_c > -dys, bys_c, -dys)
        return yse_t, yse_c

    def __get_layer_elastic_tuple(self, elastic_z, layer):
        topo, icd, moho, slablab = self.geo_model.get_boundaries()
        if layer == 'uc':
            top_z = topo
            bottom_z = icd
        elif layer == 'lc':
            top_z = icd
            bottom_z = moho
        elif layer == 'lm':
            top_z = moho
            bottom_z = slablab
        #elastic_z = self.geo_model.data.cut_3d_array(elastic_z,
        #                                             top_z=top_z,
        #                                             bottom_z=bottom_z)
        elastic_z = elastic_z.crop(top_z=top_z, bottom_z=bottom_z)
        elastic_z = np.ma.array(elastic_z, mask=np.isnan(elastic_z))
        bottom_elastic_z = np.amax(elastic_z, axis=2).filled(np.nan)
        top_elastic_z = np.amin(elastic_z, axis=2).filled(np.nan)
        elastic_thickness = bottom_elastic_z - top_elastic_z
        layer_tuple = np.stack((top_elastic_z, bottom_elastic_z,
                                elastic_thickness), axis=2)
        return layer_tuple, elastic_z.filled(np.nan)

    def __set_eet(self, yse):
        #elastic_z = self.cs.get_3D_grid()[2].copy() # needs to be Z from topo
        elastic_z = self.__get_depth_from_topo().copy()
        with np.errstate(invalid='ignore'):
            elastic_z[yse < self.vars.s_max] = np.nan
        elastic_z[np.isnan(yse)] = np.nan
        uc_tuple, e_z_uc = self.__get_layer_elastic_tuple(elastic_z, 'uc')
        lc_tuple, e_z_lc = self.__get_layer_elastic_tuple(elastic_z, 'lc')
        lm_tuple, e_z_lm = self.__get_layer_elastic_tuple(elastic_z, 'lm')
        share_icd = uc_tuple[:, :, 1] + self.cs.z_step == lc_tuple[:, :, 0]
        share_moho = lc_tuple[:, :, 1] + self.cs.z_step == lm_tuple[:, :, 0]
        elastic_thickness = np.stack((uc_tuple[:, :, 2],
                                      lc_tuple[:, :, 2],
                                      lm_tuple[:, :, 2]),
                                     axis=2)
        attached_ths, detached_ths = SpatialArray.divide_array_by_areas(elastic_thickness, share_moho)
        detached_ths_2l, detached_ths_3l = SpatialArray.divide_array_by_areas(detached_ths, share_icd)
        detached_ths_2l[:, :, 0] = (detached_ths_2l[:, :, 0]
                                    + detached_ths_2l[:, :, 1])
        detached_ths_2l[:, :, 1] = np.nan
        detached_ths = SpatialArray.combine_arrays_by_areas(detached_ths_2l, detached_ths_3l, share_icd)
        attached_eet = self.__calc_eet_attached(attached_ths)
        detached_eet = self.__calc_eet_detached(detached_ths)
        eet = SpatialArray.combine_arrays_by_areas(attached_eet, detached_eet, share_moho)
        eet_calc_data = {
            'uc_tuple': uc_tuple,
            'lc_tuple': lc_tuple,
            'lm_tuple': lm_tuple,
            'e_z_uc': e_z_uc,
            'e_z_lc': e_z_lc,
            'e_z_lm': e_z_lm,
            'share_icd': SpatialArray2D(share_icd, self.cs),
            'share_moho': SpatialArray2D(share_moho, self.cs)
        }
        return eet, Box(eet_calc_data)

    def get_layer_elastic_tuple(self, bys_depth, dys_depth, layer):
        if layer == 'uc':
            top_boundary = self.geo_model.get_topo()
            bottom_boundary = self.geo_model.get_icd()
        elif layer == 'lc':
            top_boundary = self.geo_model.get_icd()
            bottom_boundary = self.geo_model.get_moho()
        elif layer == 'lm' or layer == 'flm':
            top_boundary = self.geo_model.get_moho()
            bottom_boundary = self.geo_model.get_slab_lab()
        top_elastic = np.minimum(top_boundary, bys_depth)
        bottom_elastic = np.maximum(bottom_boundary, dys_depth)
        thickness = top_elastic - bottom_elastic
        top_elastic[thickness < 0] = np.nan
        bottom_elastic[thickness < 0] = np.nan
        thickness[thickness < 0] = 0
        return np.stack((top_elastic, bottom_elastic, thickness), axis=2)

    def __set_eet_2(self):
        s_max = self.vars.s_max
        # Get depth at which bys reaches s_max
        bs_t = self.vars.bs_t
        bs_c = self.vars.bs_c
        bys_depth_t = self.__calc_depth_from_brittle_yield_strength(
            bs_t, s_max)
        bys_depth_c = self.__calc_depth_from_brittle_yield_strength(
            bs_c, -s_max)
        bys_depth_t = self.geo_model.get_topo() - bys_depth_t
        bys_depth_c = self.geo_model.get_topo() - bys_depth_c
        # Get temperature at which dys reaches s_max
        e = self.vars.e
        r = self.vars.r
        uc = self.vars.cs
        lc = self.vars.ci
        blm = self.vars.ml
        flm = self.vars.mla
        temp_uc = self.calc_temperature_from_ductile_yield_strength(
            e, uc.n, uc.a, uc.h, r, s_max)
        temp_lc = self.calc_temperature_from_ductile_yield_strength(
            e, lc.n, lc.a, lc.h, r, s_max)
        temp_blm = self.calc_temperature_from_ductile_yield_strength(
            e, blm.n, blm.a, blm.h, r, s_max)
        if self.vars.slm is True:
            temp_flm = self.calc_temperature_from_ductile_yield_strength(
                e, flm.n, flm.a, flm.h, r, s_max)
            dys_temp = np.array([temp_uc, temp_lc, temp_blm, temp_flm])
            dys_depth = np.empty((*self.cs.get_2D_shape(),4))
        else:
            dys_temp = np.array([temp_uc, temp_lc, temp_blm])
            dys_depth = np.empty((*self.cs.get_2D_shape(),3))
        # Interpolate depth to get exact depth of the ductile yield temperature
        # (Depth at wich dys reaches s_max)
        idxs = SpatialArray3D(self.cs.get_3D_indexes()[2], self.cs)
        idxs = idxs.astype(np.float).mask_irrelevant()
        top_idx = np.nan_to_num(
            idxs.extract_surface(self.geo_model.get_topo()))
        bottom_idx = np.nan_to_num(
            idxs.extract_surface(self.geo_model.get_slab_lab()))
        top_idx = top_idx.astype(int)
        bottom_idx = bottom_idx.astype(int)
        depth = self.cs.get_3D_grid()[2]
        geotherm = self.thermal_model.get_geotherm()
        dys_depth[:,:] = np.nan
        for i in np.arange(idxs.shape[0]):
            for j in np.arange(idxs.shape[1]):
                if (np.isnan(idxs[i,j,0]) == False
                        and (top_idx[i,j] < bottom_idx[i,j])):
                    dys_depth[i,j,:] = np.interp(dys_temp,
                        geotherm[i,j,top_idx[i,j]:bottom_idx[i,j]],
                        depth[i,j,top_idx[i,j]:bottom_idx[i,j]],
                        left=np.inf, right=-np.inf)
        # Get values of intersection between YSE with s_max,
        # and the elastic thickness for each litospheric layer
        uc_tuple = self.get_layer_elastic_tuple(
            bys_depth_c, dys_depth[:,:,0], 'uc')
        lc_tuple = self.get_layer_elastic_tuple(
            bys_depth_c, dys_depth[:,:,1], 'lc')
        blm_tuple = self.get_layer_elastic_tuple(
            bys_depth_c, dys_depth[:,:,2], 'lm')
        if self.vars.slm is True:
            flm_tuple = self.get_layer_elastic_tuple(
                bys_depth_c, dys_depth[:,:,3], 'lm')
            lm_tuple = SpatialArray.combine_arrays_by_areas(
                blm_tuple, flm_tuple,
                self.geo_model.get_slab_lab_int_area().astype(bool))
        else:
            lm_tuple = blm_tuple
        # Get the coupled and decoupled zones for moho and icd
        share_moho = lc_tuple[:, :, 1] == lm_tuple[:, :, 0]
        share_icd = uc_tuple[:, :, 1] == lc_tuple[:, :, 0]
        layers_thickness_i = np.stack(
            (uc_tuple[:,:,2], lc_tuple[:,:,2], lm_tuple[:,:,2]), axis=2)
        # If icd is coupled sum thickness of upper crust with lower crust
        coupled_icd_ths, decoupled_icd_ths = SpatialArray.divide_array_by_areas(
            layers_thickness_i, share_icd)
        coupled_icd_ths[:,:,1] = coupled_icd_ths[:,:,0] + coupled_icd_ths[:,:,1]
        coupled_icd_ths[:,:,0] = 0
        layers_thickness = SpatialArray.combine_arrays_by_areas(
            coupled_icd_ths, decoupled_icd_ths, share_icd)
        # If moho is coupled sum thickness of lower crust with litospheric mantle
        coupled_moho_ths, decoupled_moho_ths = SpatialArray.divide_array_by_areas(
            layers_thickness, share_moho)
        coupled_moho_ths[:,:,2] = coupled_moho_ths[:,:,1] + coupled_moho_ths[:,:,2]
        coupled_moho_ths[:,:,1] = 0
        layers_thickness = SpatialArray.combine_arrays_by_areas(
            coupled_moho_ths, decoupled_moho_ths, share_moho)
        ##### TEMPORAL #####
        # Metodo erroneo para calcular EET
        # If moho is decoupled separate the model in 3 layers or 2 layers (when icd coupled)
        coupled_ths_2, decoupled_ths_2 = SpatialArray.divide_array_by_areas(
            layers_thickness_i, share_moho)
        decoupled_ths_2l, decoupled_ths_3l = SpatialArray.divide_array_by_areas(
            decoupled_ths_2, share_icd)
        # If icd is coupled (2 layers) sum thickness of upper crust with lower crust
        decoupled_ths_2l[:,:,1] = decoupled_ths_2l[:,:,0] + decoupled_ths_2l[:,:,1]
        decoupled_ths_2l[:,:,0] = np.nan
        decoupled_ths_2 = SpatialArray.combine_arrays_by_areas(
            decoupled_ths_2l, decoupled_ths_3l, share_icd)
        # Calculate eet differently for the coupled and decoupled moho
        eet_coupled = self.__calc_eet_attached(coupled_ths_2)
        eet_decoupled = self.__calc_eet_detached(decoupled_ths_2)
        eet_wrong = SpatialArray.combine_arrays_by_areas(
            eet_coupled, eet_decoupled, share_moho)
        ##### FIN TEMPORAL #####
        # If topo is at less than 50 km. from slab, add oceanic plate thickness
        # to litospheric mantle thickness
        #diff = self.geo_model.get_topo() - self.geo_model.get_slab_lab()
        #op_points = diff <= 50
        #op_points_thickness = (
        #    op_points * self.thermal_model.get_oceanic_plate_thickness())
        #ths_plus_op, ths = SpatialArray.divide_array_by_areas(
        #    layers_thickness, op_points)
        #ths_plus_op[:,:,2] = np.nansum(
        #    np.dstack((ths_plus_op[:,:,2], op_points_thickness)), 2) 
        #layers_thickness = SpatialArray.combine_arrays_by_areas(
        #    ths_plus_op, ths, op_points)
        ## Calculate eet according to formula
        eet = self.__calc_eet_detached(layers_thickness)

        eet_calc_data = {
            'uc_tuple': uc_tuple,
            'lc_tuple': lc_tuple,
            'lm_tuple': lm_tuple,
            'share_icd': SpatialArray2D(share_icd, self.cs),
            'share_moho': SpatialArray2D(share_moho, self.cs),
            'bys_depth_t': bys_depth_t,
            'bys_depth_c': bys_depth_c,
            'dys_depth': dys_depth,
            'dys_temp': dys_temp,
            'layers_thickness': layers_thickness,
        }
        return eet, eet_wrong, Box(eet_calc_data)

    def __set_integrated_strength(self):
        yse_c = self.yse_c
        geo_model_mask = np.isnan(self.geo_model.get_3D_geometric_model())
        integrated_strength = np.trapz(
            np.ma.array(yse_c,mask=geo_model_mask), axis=2)
        return integrated_strength

    def get_yse(self):
        return (SpatialArray3D(self.yse_t, self.cs),
                SpatialArray3D(self.yse_c, self.cs))

    def get_eet(self):
        return SpatialArray2D(self.eet, self.cs).mask_irrelevant_eet()

    def get_eet_from_trench(self):
        return SpatialArray2D(self.eet, self.cs).mask_irrelevant()

    def get_eet_wrong(self):
        return SpatialArray2D(self.eet_wrong, self.cs).mask_irrelevant_eet()

    def get_integrated_strength(self):
        return SpatialArray2D(
            self.integrated_strength, self.cs).mask_irrelevant_eet()

    def get_geometric_model(self):
        return self.geo_model

    def get_coordinate_system(self):
        return self.cs


def compute(gm_data, slab_lab_areas, trench_age, rhe_data, coast,
        t_input, m_input):
    d = Data(gm_data, slab_lab_areas, trench_age, rhe_data, coast,
        t_input, m_input)
    cs_d = d.get_cs_data()
    gm_d = d.get_gm_data()
    tm_d = d.get_tm_data()
    mm_d = d.get_mm_data()
    cs = CoordinateSystem(cs_d, 0.2, 1)
    #cs = CoordinateSystem(cs_d, 0.2, 0.1, 1)
    gm = GeometricModel(gm_d, cs)
    tm = ThermalModel(tm_d, gm, cs)
    mm = MechanicModel(mm_d, gm, tm, cs)
    model = {
        'd': d,
        'cs': cs,
        'gm': gm,
        'tm': tm,
        'mm': mm
    }
    model = DotMap(model)
    return model
