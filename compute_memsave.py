import math
import numpy as np
from dotmap import DotMap
import resource

def mem():
    print('Memory usage         : % 2.2f MB' % round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1)
    )

class Data(object):

    @staticmethod
    def __get_max(values):
        max_v = math.ceil(np.nanmax(values))
        return max_v

    @staticmethod
    def __get_min(values):
        min_v = math.floor(np.nanmin(values))
        return min_v

    def __init__(self, initial_data, slab_lab_areas, trench_age, rheologic_data,
                 t_input, m_input):

        self.coordinate_data = self.__set_coordinate_data(initial_data)
        self.geometric_data = self.__set_geometric_data(initial_data)
        self.populated_points = np.invert(np.isnan(self
                                                   .geometric_data.z_slab_lab))
        self.slab_lab_areas = np.asarray(slab_lab_areas, dtype=bool)
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
            'populated_points': self.populated_points,
            'z_slab_lab': self.geometric_data.z_slab_lab  # to get plates_areas
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

    def __init__(self, cs_data, xy_step, z_step):
        self.data = cs_data.coordinate_data
        print("self.data")
        mem()
        self.xy_step = xy_step
        self.z_step = z_step
        self.x_axis = self.__set_axis(self.data.max_lon, self.data.min_lon,
                                      self.xy_step)
        print("self.x_axis")
        mem()
        self.y_axis = self.__set_axis(self.data.max_lat, self.data.min_lat,
                                      self.xy_step, revert=True)
        print("self.y_axis")
        mem()
        self.z_axis = self.__set_axis(self.data.max_z, self.data.min_z,
                                      self.z_step, revert=True)
        print("self.z_axis")
        mem()
        self.grid_2D = self.__set_grid([self.x_axis, self.y_axis])
        print("self.grid_2d")
        mem()
        self.grid_3D = self.__set_grid([self.x_axis, self.y_axis, self.z_axis])
        print("self.grid_3d")
        mem()
        self.populated_points = cs_data.populated_points
        print("self.populated_points")
        mem()
        self.plates_areas = self.__set_plates_areas(cs_data.z_slab_lab)
        print("self.paltes_areas")
        mem()
        self.relevant_area = self.__set_relevant_area(self.populated_points,
                                                      self.plates_areas)
        print("self.relevant_area")
        mem()
        self.relevant_volume = self.__set_relevant_volume(self.relevant_area)
        print("self.relevant_volume")
        mem()

    def __set_axis(self, max_v, min_v, step, revert=False):
        if revert is True:
            first_v, last_v = max_v, min_v
        else:
            first_v, last_v = min_v, max_v

        return np.linspace(first_v, last_v,
                           num=(abs(last_v-first_v))/step+1,
                           endpoint=True)

    def __set_grid(self, axes, mask=False):
        grid = np.meshgrid(*[n for n in axes], indexing='ij')
        return grid

    def __set_plates_areas(self, z_slab_lab):
        """False (0) in Nazca Plate, True (1) in South American Plate"""
        slab_lab = self.reshape_data(z_slab_lab)
        g = np.gradient(slab_lab, axis=0)
        with np.errstate(invalid='ignore'):  # error_ignore
            high_g = np.absolute(g) > 1  # type: np.ndarray
        trench_start = np.argmax(high_g, axis=0)  # gets first true value
        i_idx = self.get_2D_indexes()[0]
        plates_area = np.ones(self.get_2D_shape(), dtype=bool)
        plates_area[i_idx < trench_start] = 0
        return plates_area

    def __set_relevant_area(self, populated_points, plates_areas):
        populated_area = self.reshape_data(populated_points)
        relevant_area = np.ones(self.get_2D_shape(), dtype=bool)
        relevant_area[populated_area == 0] = 0
        relevant_area[plates_areas == 0] = 0
        return relevant_area

    def __set_relevant_volume(self, relevant_area):
        relevant_volume = np.zeros(self.get_3D_shape(), dtype=bool)
        relevant_volume[:, :, :] = relevant_area[:, :, np.newaxis]
        return relevant_volume

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

    def get_2D_shape(self):
        return len(self.x_axis), len(self.y_axis)

    def get_3D_shape(self):
        return len(self.x_axis), len(self.y_axis), len(self.z_axis)

    """
    def get_2D_grid(self):
        grid_2D = self.__set_grid([self.x_axis, self.y_axis])
        grid_2D_relevant = []
        for n in range(len(grid_2D)):
            grid_2D_relevant.append(SpatialArray2D(grid_2D[n],
                                                   self).mask_irrelevant())
        return grid_2D_relevant

    def get_3D_grid(self):
        grid_3D = self.__set_grid([self.x_axis, self.y_axis, self.z_axis])
        grid_3D_relevant = []
        for n in range(len(grid_3D)):
            grid_3D_relevant.append(SpatialArray3D(grid_3D[n],
                                                   self).mask_irrelevant())
        return grid_3D_relevant
    """

    def get_2D_grid(self):
        grid_2D = []
        for n in range(len(self.grid_2D)):
            grid_2D.append(SpatialArray2D(self.grid_2D[n],
                                          self).mask_irrelevant())
        return grid_2D

    def get_3D_grid(self):
        grid_3D = []
        for n in range(len(self.grid_3D)):
            grid_3D.append(SpatialArray3D(self.grid_3D[n],
                                          self).mask_irrelevant())
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


class SpatialArray(np.ndarray):

    def __new__(cls, input_array, coordinate_system=None):
        # print('new method of SpatialArray called')
        # print('cls is', cls)
        obj = np.asarray(input_array).view(cls)
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
        masked_array = self.copy()
        masked_array[mask] = np.nan
        #masked_array = np.ma.array(self, mask=mask)
        #if nan_fill is True:
        #    masked_array = masked_array.filled(np.nan)
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
            index = np.where(self.cs.get_y_axis() == latitude)[0][0]
            cross_section = self[:, index]
            return cross_section
        elif longitude:
            #index = list(self.cs.get_x_axis()).index(longitude)
            index = np.where(self.cs.get_x_axis() == longitude)[0][0]
            print("index:", index)
            cross_section = self[:, index]
        elif lat1 and lon1 and lat2 and lon2:
            cross_section = self
        return cross_section

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
        masked_array = self.copy()
        masked_array[mask] = np.nan
        #masked_array = np.ma.array(self, mask=mask)
        #if nan_fill is True:
        #    masked_array = masked_array.filled(np.nan)
        return masked_array

    def extract_surface(self, z_2D):
        z_2D = np.ceil(z_2D)
        z_3D = self.cs.get_3D_grid()[2]
        a = z_3D != z_2D[:, :, np.newaxis]
        array_3D = self.copy()
        array_3D[a] = 0
        surface = np.sum(array_3D, axis=2)
        surface[surface == 0] = np.nan
        return SpatialArray2D(surface, z_2D.cs)

    def crop(self, top_z='top', bottom_z='bottom'):
        z_3D = self.cs.get_3D_grid()[2]
        if isinstance(top_z, str):
            top_z = np.inf
        else:
            top_z = np.ceil(top_z)[:, :, np.newaxis]
        if isinstance(bottom_z, str):
            bottom_z = -np.inf
        else:
            bottom_z = np.ceil(bottom_z)[:, :, np.newaxis]
        with np.errstate(invalid='ignore'): # error_ignore
            a = top_z < z_3D
            b = bottom_z > z_3D
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
            index = np.where(self.cs.get_y_axis() == latitude)[0][0]
            cross_section = self[:, index, :]
            return cross_section
        elif longitude:
            #index = list(self.cs.get_x_axis()).index(longitude)
            index = np.where(self.cs.get_x_axis() == longitude)[0][0]
            cross_section = self[:, index, :]
        elif lat1 and lon1 and lat2 and lon2:
            cross_section = self
        return cross_section 

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
        print("self.data")
        mem()
        self.slab_lab_areas = gm_data.slab_lab_areas
        print("self.slab_lab_areas")
        mem()
        self.cs = coordinate_system
        print("self.cs")
        mem()
        self.slab_lab = self.__set_boundary(self.data.z_slab_lab)
        print("self.slab_lab")
        mem()
        self.moho = self.__set_boundary(self.data.z_moho)
        print("self.moho")
        mem()
        self.icd = self.__set_boundary(self.data.z_icd)
        print("self.icd")
        mem()
        self.topo = self.__set_boundary(self.data.z_topo)
        print("self.topo")
        mem()
        self.geo_model_3D = self.__set_3D_geometric_model()
        print("self.geo_model_3D")
        mem()
        # self.slab_lab_int_index = self.__set_slab_lab_int_index()
        self.slab_lab_int_index = self.__set_slab_lab_int_index_from_areas()
        print("self.slab_lab_int_index")
        mem()
        self.slab_lab_int_depth = self.__set_slab_lab_int_depth()
        print("self.slab_lab_int_depth")
        mem()
        self.slab_lab_int_topo = self.__set_slab_lab_int_topo()
        print("self.slab_lab_int_topo")
        mem()
        self.slab_lab_int_area = self.__set_slab_lab_int_area()
        print("self.slab_lab_int_area")
        mem()


    def __set_boundary(self, geometric_data):
        return SpatialArray2D(self.cs.reshape_data(geometric_data), self.cs)

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
                geo_model_3D[z < c] = a
        return SpatialArray3D(geo_model_3D, self.cs)#.mask_irrelevant()

    def sencos(self):
        A = self.cs.get_2D_grid()[0]
        B = self.cs.get_2D_grid()[1]
        C = A + B
        return C

    def sencos_sa(self):
        A = self.cs.get_2D_grid()[0]
        B = self.cs.get_2D_grid()[1]
        C = A + B
        return SpatialArray3D(C, self.cs)

    def __set_layer_thickness(self):
        pass

    def __set_slab_lab_int_index(self):
        # Slab/Lab gradients array
        g = np.gradient(self.slab_lab, axis=0)
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

    def __set_slab_lab_int_index_from_areas(self):
        sli_idx = np.argmax(self.slab_lab_areas, axis=1)
        # sli_idx = sli_idx - 1
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

    def __set_slab_lab_int_area(self):
        sli_idx = self.slab_lab_int_index
        i_idx = self.cs.get_2D_indexes()[0]
        sli_area = np.zeros(self.cs.get_2D_shape())
        sli_area[i_idx > sli_idx] = 1
        sli_ara = np.asarray(sli_area, dtype=bool)
        return sli_area

    def set_layer_property(self, cs, ci, ml):
        r = np.ones(self.geo_model_3D.shape) * np.nan
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
        lab_temp = tp + g * depth
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
                     / (v * abs(sli_depth-sli_topo)*1.e3)
                     - sli_q_zero/v)
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
        z = z*1.e3
        z_topo = z_topo*1.e3
        z_sl = z_sl*1.e3
        base_temp = temp_sl-((h*delta**2)/k)*(np.exp(z_topo/delta)
                                              - np.exp(z_sl/delta))
        rad_temp = ((h*delta**2)/k)*(np.exp(z_topo/delta)-np.exp(z/delta))
        geotherm = rad_temp + (abs(z-z_topo)/abs(z_sl-z_topo))*base_temp
        return geotherm

    @staticmethod
    def __calc_surface_heat_flow(h, delta, k, z_topo, z_sl, temp_sl):
        z_topo = z_topo*1.e3
        z_sl = z_sl*1.e3
        base_temp = temp_sl-((h*delta**2)/k)*(np.exp(z_topo/delta) -
                                              np.exp(z_sl/delta))
        heat_flow = (-h*delta-k/(z_sl-z_topo)*base_temp)
        return heat_flow

    def __init__(self, tm_data, geometric_model, coordinate_system):
        self.geo_model = geometric_model
        print("self.geo_model")
        mem()
        self.cs = coordinate_system
        print("self.cs")
        mem()
        self.vars = DotMap(self.__set_variables(tm_data.t_input,
                                                 tm_data.trench_age))
        print("self.vars")
        mem()
        self.slab_lab_temp = self.__set_slab_lab_temp()
        print("self.slab_lab_temp")
        mem()
        self.surface_heat_flow = self.__set_surface_heat_flow()
        print("self.surface_heat_flow")
        mem()
        self.geotherm = self.__set_geotherm()
        print("self.geotherm")
        mem()

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
                rh.append(prop_v[n]*h[n])
                n += 1
            r = (sum(rh) / sum(h))[:, :, np.newaxis]
            r = np.repeat(r, self.cs.get_3D_shape()[2], axis=2)
            geo_model_mask = np.isnan(self.geo_model.get_3D_geometric_model())
            r = np.ma.array(r, mask=geo_model_mask).filled(np.nan)
        return SpatialArray3D(r, self.cs)

    def __set_delta(self, t_input):
        boundaries = self.geo_model.get_boundaries()
        if t_input.delta_icd is True:
            delta = boundaries[0] - boundaries[1]
        else:
            delta = t_input.delta
        return delta

    def __set_trench_age(self, trench_age, t_input):
        if t_input.t_lat is True:
            trench_age = trench_age[:, 1]
        else:
            trench_age = t_input.t
        return trench_age

    def __set_variables(self, t_input, trench_age):
        t_input = DotMap(t_input)
        t_vars = {
            'k_cs': t_input['k_cs'],
            'k_ci': t_input['k_ci'],
            'k_ml': t_input['k_ml'],
            'k': self.__set_k_or_h(t_input, 'k'),
            'h_cs': t_input['H_cs'],
            'h_ci': t_input['H_ci'],
            'h_ml': t_input['H_ml'],
            'h': self.__set_k_or_h(t_input, 'h'),
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
        slab_k = self.vars.k.extract_surface(z_sl)
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
        slab_k = self.vars.k.extract_surface(z_sl)
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
        slab_temp = self.__get_slab_temp()
        lab_temp = self.__get_lab_temp()  # type: np.ndarray
        slab_lab_temp = np.zeros(self.cs.get_2D_shape())
        slab_lab_temp[b] = slab_temp[b]
        slab_lab_temp[c] = lab_temp[c]
        return slab_lab_temp

    def __set_surface_heat_flow(self):
        delta = self.vars.delta
        z_topo = self.geo_model.get_topo()
        z_sl = self.geo_model.get_slab_lab()
        slab_lab_k = self.vars.k.extract_surface(z_sl)
        slab_lab_h = self.vars.h.extract_surface(z_sl)
        temp_sl = self.slab_lab_temp
        heat_flow = self.__calc_surface_heat_flow(slab_lab_h, delta,
                                                  slab_lab_k, z_topo, z_sl,
                                                  temp_sl)
        return heat_flow

    def __set_geotherm(self):
        k = self.vars.k
        h = self.vars.h
        if self.vars.delta.shape:
            delta = self.vars.delta[:, :, np.newaxis]
        else:
            delta = self.vars.delta
        z = self.cs.get_3D_grid()[2]
        z_topo = self.geo_model.get_topo()[:, :, np.newaxis]
        z_sl = self.geo_model.get_slab_lab()[:, :, np.newaxis]
        temp_sl = self.slab_lab_temp[:, :, np.newaxis]
        geotherm = self.__calc_geotherm(h, delta, k, z, z_topo, z_sl, temp_sl)
        return geotherm

    def get_surface_heat_flow(self):
        return SpatialArray2D(self.surface_heat_flow, self.cs)

    def get_geotherm(self):
        return SpatialArray3D(self.geotherm, self.cs)


class MechanicModel(object):

    @staticmethod
    def __calc_brittle_yield_strength(bs, depth):
        bys = bs*(depth*1.e3)*1.e-6
        return bys

    @staticmethod
    def __calc_ductile_yield_strength(es, n, a, h, r, temp):
        with np.errstate(over='ignore'):
            dys = (es/a)**(1/n)*np.exp(h/(n*r*(temp+273.15)))*1.e-6
        return dys

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
        super().__init__()
        print("superinit")
        mem()
        self.geo_model = geo_model
        print("cs.geo_model")
        mem()
        self.cs = coordinate_system
        print("self.cs")
        mem()
        self.thermal_model = thermal_model
        print("self.thermal_model")
        mem()
        self.vars = DotMap(self.__set_variables(mm_data.m_input,
                                                 mm_data.rheologic_data))
        print("self.vars")
        mem()
        self.depth_from_topo = self.__set_depth_from_topo()
        print("self.depth_from_topo")
        mem()
        #self.bys_t, self.bys_c = self.__set_brittle_yield_strength()
        #print("self.bys_t, self.bys_c")
        #mem()
        #self.dys = self.__set_ductile_yield_strength()
        #print("self.dys")
        #mem()
        self.yse_t, self.yse_c = self.__set_yield_strength_envelope()
        print("self.yse_t, self.yse_c")
        mem()
        self.eet = self.__set_eet(self.yse_t)
        print("self.eet")
        mem()

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

    def __set_variables(self, m_input, rhe_data):
        m_input = DotMap(m_input)
        m_vars = {
            'bs_t': m_input['Bs_t'],
            'bs_c': m_input['Bs_c'],
            'e': m_input['e'],
            'r': m_input['R'],
            's_max': m_input['s_max'],
            'cs': self.__get_rheologic_vars(m_input['Cs'], rhe_data),
            'ci': self.__get_rheologic_vars(m_input['Ci'], rhe_data),
            'ml': self.__get_rheologic_vars(m_input['Ml'], rhe_data)
        }
        return m_vars

    def __set_depth_from_topo(self):
        depth_from_topo = -(self.cs.get_3D_grid()[2]
                            - self.geo_model.get_topo()[:, :, np.newaxis])
        return depth_from_topo

    def __get_brittle_yield_strength(self):
        bs_t = self.vars.bs_t
        bs_c = self.vars.bs_c
        #depth = self.cs.get_3D_indexes()[2]  # needs to be Z from topo
        depth_from_topo = self.depth_from_topo
        depth_from_topo = SpatialArray3D(depth_from_topo.astype(np.float64), self.cs).mask_irrelevant()
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
        n = self.geo_model.set_layer_property(cs.n, ci.n, ml.n)
        a = self.geo_model.set_layer_property(cs.a, ci.a, ml.a)
        h = self.geo_model.set_layer_property(cs.h, ci.h, ml.h)

        dys = self.__calc_ductile_yield_strength(e, n, a, h, r, temp)
        return dys

    def __set_yield_strength_envelope(self):
        bys_t, bys_c = self.__get_brittle_yield_strength()  # type: np.ndarray
        dys = self.__get_ductile_yield_strength()  # type: np.ndarray
        with np.errstate(invalid='ignore'):
            yse_t = np.where(bys_t < dys, bys_t, dys)
            yse_c = np.where(bys_c > dys, bys_c, dys)
        return SpatialArray3D(yse_t, self.cs), SpatialArray3D(yse_c, self.cs)

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
        top_elastic_z = np.amax(elastic_z, axis=2).filled(np.nan)
        bottom_elastic_z = np.amin(elastic_z, axis=2).filled(np.nan)
        elastic_thickness = top_elastic_z - bottom_elastic_z
        layer_tuple = np.stack((top_elastic_z, bottom_elastic_z,
                                elastic_thickness), axis=2)
        return layer_tuple, elastic_z.filled(np.nan)

    def __set_eet(self, yse):
        #elastic_z = self.cs.get_3D_grid()[2].copy() # needs to be Z from topo
        elastic_z = self.depth_from_topo.copy()
        with np.errstate(invalid='ignore'):
            elastic_z[yse < self.vars.s_max] = np.nan
        uc_tuple, e_z_uc = self.__get_layer_elastic_tuple(elastic_z, 'uc')
        lc_tuple, e_z_lc = self.__get_layer_elastic_tuple(elastic_z, 'lc')
        lm_tuple, e_z_lm = self.__get_layer_elastic_tuple(elastic_z, 'lm')
        share_icd = uc_tuple[:, :, 1] == lc_tuple[:, :, 0]
        share_moho = lc_tuple[:, :, 1] == lm_tuple[:, :, 0]
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

        return SpatialArray2D(eet, self.cs)

    def get_yse(self):
        return self.yse_t, self.yse_c

    def get_eet(self):
        return self.eet


def compute(gm_data, slab_lab_areas, trench_age, rhe_data, t_input, m_input):
    print("#####Initial M.S.:")
    mem()
    d = Data(gm_data, slab_lab_areas, trench_age, rhe_data, t_input, m_input)
    print("#####After D M.S:")
    mem()
    cs_d = d.get_cs_data()
    gm_d = d.get_gm_data()
    tm_d = d.get_tm_data()
    mm_d = d.get_mm_data()
    print("#####After D vars M.S:")
    mem()
    cs = CoordinateSystem(cs_d, 0.2, 1)
    print("#####After CS M.S:")
    mem()
    gm = GeometricModel(gm_d, cs)
    print("#####After GM M.S:")
    mem()
    tm = ThermalModel(tm_d, gm, cs)
    print("#####After TM M.S:")
    mem()
    mm = MechanicModel(mm_d, gm, tm, cs)
    print("#####After MM M.S:")
    mem()
    return d, cs, gm, tm, mm


