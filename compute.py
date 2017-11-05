""" Computations """

import math
import numpy as np
from utils import DotDict

class Data(object):

    @staticmethod
    def print_data(name, data, boolean=False):
        if boolean is True:
            np.savetxt(name, data, fmt="%1d", delimiter=" ")
        else:
            np.savetxt(name, data, fmt="%9.2e", delimiter="   ")

    def __init__(self, data, xy_step, z_step):
        self.data = data
        self.xy_step = xy_step
        self.z_step = z_step
        self.x_axis = self.__set_axis(values=self.data[:, 0])
        self.y_axis = self.__set_axis(values=self.data[:, 1], revert=True)
        self.z_axis = self.__set_axis(values=self.data[:, 5], revert=True,
                                      low_values=self.data[:, 2], z_axis=True)
        self.mask_2d = self.__set_2d_mask()
        self.mask_3d = self.__set_3d_mask()
        self.grid2d = self.__set_grid([self.x_axis, self.y_axis],
                                      mask=True)
        self.grid3d = self.__set_grid([self.x_axis, self.y_axis, self.z_axis],
                                      mask=True)
        pass

    def __set_axis(self, values, low_values=False, revert=False, z_axis=False):
        if low_values is False:
            low_values = values
        max_v = math.ceil(np.nanmax(values))
        min_v = math.floor(np.nanmin(low_values))
        if revert is True:
            first_v, last_v = max_v, min_v
        else:
            first_v, last_v = min_v, max_v
        if z_axis is True:
            step = self.z_step
        else:
            step = self.xy_step
        return np.linspace(first_v, last_v,
                           num=(abs(last_v-first_v))/step+1,
                           endpoint=True)

    def __set_2d_mask(self):
        lab_slab_data = self.reshape_2d_data(self.data[:, 2])
        g = np.gradient(lab_slab_data, axis=0)
        with np.errstate(invalid='ignore'):
            high_g = np.absolute(g) > 1  # type: np.ndarray
        trench_start = np.argmax(high_g, axis=0)
        i_idx = self.get_2d_indexes()[0]
        mask_2d = np.zeros(self.get_2d_shape())
        mask_2d[np.isnan(lab_slab_data)] = 1
        mask_2d[i_idx < trench_start] = 1
        return mask_2d

    def __set_3d_mask(self):
        # Use lab/slab data nan values as mask
        mask_2d = self.mask_2d
        mask_3d = np.zeros(self.get_3d_shape(), dtype=bool)
        mask_3d[:, :, :] = mask_2d[:, :, np.newaxis]
        return mask_3d

    def __set_grid(self, axes, mask=False):
        grid = np.meshgrid(*[n for n in axes], indexing='ij')
        masked_grid = []
        if mask is True:
            for n in grid:
                if len(n.shape) == 2:
                    masked_grid.append(self.mask_2d_array(n, nan_fill=True))
                else:
                    masked_grid.append(self.mask_3d_array(n, nan_fill=True))
            return masked_grid
        else:
            return grid

    def get_data(self):
        return self.data

    def get_axes(self):
        return [self.x_axis, self.y_axis, self.z_axis]

    def get_2d_shape(self):
        return len(self.x_axis), len(self.y_axis)

    def get_3d_shape(self):
        return len(self.x_axis), len(self.y_axis), len(self.z_axis)

    def get_2d_grid(self):
        return self.grid2d

    def get_3d_grid(self):
        return self.grid3d

    def get_2d_indexes(self):
        return np.indices(self.get_2d_shape())

    def get_3d_indexes(self):
        return np.indices(self.get_3d_shape())

    def reshape_2d_data(self, data2d):
        return data2d.T.reshape((len(self.y_axis), len(self.x_axis))).T

    def mask_2d_array(self, array_2d, nan_fill=False):
        mask_2d = np.ma.array(array_2d, mask=self.mask_2d)
        if nan_fill is True:
            mask_2d = mask_2d.filled(np.nan)
        return mask_2d

    def mask_3d_array(self, array_3d, nan_fill=False):
        mask_3d = np.ma.array(array_3d, mask=self.mask_3d)
        if nan_fill is True:
            mask_3d = mask_3d.filled(np.nan)
        return mask_3d

    def extract_2d_surface(self, z_2d, array_3d):
        z_2d = np.ceil(z_2d)
        z_3d = self.get_3d_grid()[2]
        a = z_3d != z_2d[:, :, np.newaxis]
        array_3d = np.copy(array_3d)
        array_3d[a] = 0
        surface_2d = np.sum(array_3d, axis=2)
        surface_2d[surface_2d == 0] = np.nan
        return surface_2d

    def cut_3d_array(self, array_3d, top_z='top', bottom_z='bottom'):
        z_3d = self.get_3d_grid()[2]
        if top_z == 'top':
            top_z = np.inf
        else:
            top_z = np.ceil(top_z)[:, :, np.newaxis]
        if bottom_z == 'bottom':
            bottom_z = -np.inf
        else:
            bottom_z = np.ceil(bottom_z)[:, :, np.newaxis]
        a = top_z < z_3d
        b = bottom_z > z_3d
        array_3d = np.copy(array_3d)
        array_3d[a] = np.nan
        array_3d[b] = np.nan
        return array_3d

    def delimit_area(self, idxs):
        rows = np.repeat(np.arange(len(idxs)), idxs)
        cols = np.ones(idxs.sum(), dtype=int)
        cols[np.cumsum(idxs)[:-1]] -= idxs[:-1]
        cols = np.cumsum(cols) - 1
        areas = np.ones((idxs.shape[0], self.grid2d[2].shape[1]))
        areas[rows, cols] = 0
        return areas


class Model3d(object):

    def divide_by_areas(self, array, areas):
        if len(array.shape) == 2:
            array = array[:, :, np.newaxis]
        array_1 = np.copy(array)
        array_2 = np.copy(array)
        areas_3d = np.repeat(areas[:, :, np.newaxis], array.shape[2],
                             axis=2)
        array_1[np.invert(areas_3d)] = np.nan
        array_2[areas_3d] = np.nan
        if array.shape[2] == 1:
            array_1 = np.squeeze(array_1, axis=2)
            array_2 = np.squeeze(array_2, axis=2)
        return array_1, array_2

    def combine_by_areas(self, array_1, array_2, areas):
        if len(array_1.shape) == 2:
            array_1 = array_1[:, :, np.newaxis]
            array_2 = array_2[:, :, np.newaxis]
        combined_array = np.copy(array_1)
        areas_3d = np.repeat(areas[:, :, np.newaxis], array_1.shape[2],
                             axis=2)
        combined_array[np.invert(areas_3d)] = array_2[np.invert(areas_3d)]
        if combined_array.shape[2] == 1:
            combined_array = np.squeeze(combined_array, axis=2)
        return combined_array

    def __init__(self):
        pass


class GeometricModel(Model3d):
    def __init__(self, data, areas):
        super().__init__()
        self.data = data
        self.z_sl = self.__set_boundary(self.data.get_data()[:, 2])
        self.z_moho = self.__set_boundary(self.data.get_data()[:, 3])
        self.z_icd = self.__set_boundary(self.data.get_data()[:, 4])
        self.z_topo = self.__set_boundary(self.data.get_data()[:, 5])
        self.geo_model_3d = self.__set_3d_geo_model()
        # self.slab_lab_int_index = self.__get_slab_lab_int_index()
        self.slab_lab_int_index = self.__set_slab_lab_int_index_from_areas(areas)
        self.slab_lab_int_depth = self.__set_slab_lab_int_depth()
        self.slab_lab_int_topo = self.__set_slab_lab_int_topo()
        self.slab_lab_int_area = self.__set_slab_lab_int_area()
        pass
        # self.GModel3d = self.__set3dGModel()

    def __set_boundary(self, depth_data):
        return self.data.reshape_2d_data(depth_data)

    def __set_3d_geo_model(self):
        boundaries = self.get_boundaries()
        z = self.data.get_3d_grid()[2]
        geo_model_3d = self.data.mask_3d_array(np.empty(z.shape)*np.nan,
                                               nan_fill=True)
        for n in range(len(boundaries)):
            c = boundaries[n][:, :, np.newaxis]
            if n < (len(boundaries)-1):
                a = n + 1
            else:
                a = np.nan
            with np.errstate(invalid='ignore'):  # ?? ###
                geo_model_3d[z < c] = a
        return geo_model_3d

    def __set_layer_thickness(self):
        pass

    def __get_slab_lab_int_index(self):
        # Gradients array
        g = np.gradient(self.z_sl, axis=0)
        # Min Gradient Indexes
        min_g_idx = np.nanargmin(np.ma.masked_invalid(g), axis=0)
        # Positive gradient boolean array (True where gradient > 0 ...)
        with np.errstate(invalid='ignore'):
            pos_g = g > -4e-01  # (... or where gradient gets very close to 0)
        # Ignore (mask) the values to the left of the minimum gradient
        i_idx = self.data.get_2d_indexes()[0]
        g_mask = i_idx < min_g_idx
        masked_pos_g = np.ma.array(pos_g, mask=g_mask)
        # Get index of first positive gradient (Slab/Lab idx)
        sli_idx = np.argmax(masked_pos_g, axis=0)
        # Get visual boolean representation of sli_idx
        # Z = np.zeros(g.shape)
        # Z[sli_idx, self.data.get_2d_indexes()[1][0]] = 1
        # return Z.T
        return sli_idx

    def __set_slab_lab_int_index_from_areas(self, areas):
        sli_idx = np.argmax(areas, axis=1)
        # sli_idx = sli_idx - 1
        sli_idx[-1] = 0
        return sli_idx

    def __set_slab_lab_int_depth(self):
        sli_idx = self.slab_lab_int_index
        j_idx = self.data.get_2d_indexes()[1][0]
        sli_depth = self.z_sl[sli_idx, j_idx]
        return sli_depth

    def __set_slab_lab_int_topo(self):
        sli_idx = self.slab_lab_int_index
        j_idx = self.data.get_2d_indexes()[1][0]
        sli_topo = self.z_topo[sli_idx, j_idx]
        return sli_topo

    def __set_slab_lab_int_area(self):
        sli_idx = self.slab_lab_int_index
        i_idx = self.data.get_2d_indexes()[0]
        sli_area = np.zeros(self.data.get_2d_shape())
        sli_area[i_idx > sli_idx] = 1
        return sli_area

    def get_3d_geo_model(self):
        return self.geo_model_3d

    def get_boundaries(self):
        return [self.z_topo, self.z_icd, self.z_moho, self.z_sl]

    def get_slab_lab_int_index(self):
        return self.slab_lab_int_index

    def get_slab_lab_int_depth(self):
        return self.slab_lab_int_depth

    def get_slab_lab_int_topo(self):
        return self.slab_lab_int_topo

    def get_slab_lab_int_area(self):
        return self.slab_lab_int_area

    def set_layer_property(self, cs, ci, ml):
        r = np.empty(self.geo_model_3d.shape) * np.nan
        r[self.geo_model_3d == 1] = cs
        r[self.geo_model_3d == 2] = ci
        r[self.geo_model_3d == 3] = ml
        return r


class ThermalModel(Model3d):

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

    def __init__(self, geo_model, t_input, ta_data):
        super().__init__()
        self.geo_model = geo_model
        self.vars = DotDict(self.__set_variables(t_input, ta_data))
        self.slab_lab_int_temp = self.__set_slab_lab_int_temperature()
        self.slab_lab_int_sigma = self.__set_slab_lab_int_sigma()
        self.slab_sigma = self.__set_slab_sigma()
        self.slab_temp = self.__set_slab_temp()
        self.lab_temp = self.__set_lab_temp()
        self.slab_lab_temp = self.__set_slab_lab_temp()
        self.geotherm = self.__set_geotherm()
        self.surface_heat_flow = self.__set_surface_heat_flow()

    def __set_k_or_h(self, t_input, prop):
        boundaries = self.geo_model.get_boundaries()
        prop_v = None
        prop_fz = None
        if prop == "k":
            prop_v = [t_input.k_cs, t_input.k_ci, t_input.k_ml]
            prop_fz = t_input.k_z
        elif prop == "h":
            prop_v = [t_input.H_cs, t_input.H_ci, t_input.H_ml]
            prop_fz = t_input.H_z
        if prop_fz is True:
            r = self.geo_model.set_layer_property(prop_v[0], prop_v[1],
                                                  prop_v[2])
        else:
            rh = []
            h = []  # h = thickness
            n = 0
            while n < len(boundaries)-1:
                h.append(boundaries[n] - boundaries[n+1])
                rh.append(prop_v[n]*h[n])
                n += 1
            r = (sum(rh) / sum(h))[:, :, np.newaxis]
            r = np.repeat(r, self.geo_model.data.get_3d_shape()[2], axis=2)
            geo_model_mask = np.isnan(self.geo_model.get_3d_geo_model())
            r = np.ma.array(r, mask=geo_model_mask).filled(np.nan)
        return r

    def __set_delta(self, t_input):
        # Work needs to be done here
        boundaries = self.geo_model.get_boundaries()
        if t_input.delta_icd is True:
            delta = boundaries[0] - boundaries[1]
        else:
            delta = t_input.delta
        return delta

    def __set_trench_age(self, t_input, ta_data):
        if t_input.t_lat is True:
            trench_age = ta_data[:, 1]
        else:
            trench_age = t_input.t
        return trench_age

    def __set_variables(self, t_input, ta_data):
        t_input = DotDict(t_input)
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
            't': self.__set_trench_age(t_input, ta_data)*(1.e6*365.*24.*60.*60.),
            'd': t_input['D']
        }
        return t_vars

    def __set_slab_lab_int_temperature(self):
        sli_depth = self.geo_model.get_slab_lab_int_depth()
        tp = self.vars.tp
        g = self.vars.g
        sli_temp = self.__calc_lab_temp(tp, g, sli_depth)
        return sli_temp

    def __set_slab_lab_int_sigma(self):
        sli_depth = self.geo_model.get_slab_lab_int_depth()
        sli_topo = self.geo_model.get_slab_lab_int_topo()
        sli_temp = self.slab_lab_int_temp
        kappa = self.vars.kappa
        v = self.vars.v
        dip = self.vars.dip
        b = self.vars.b
        sli_s = self.__calc_s(sli_depth, sli_topo, kappa, v, dip, b)
        z_sl = self.geo_model.get_boundaries()[3]
        slab_k = self.geo_model.data.extract_2d_surface(z_sl, self.vars.k)
        sli_idx = self.geo_model.get_slab_lab_int_index()
        j_idx = self.geo_model.data.get_2d_indexes()[1][0]
        sli_k = slab_k[sli_idx, j_idx]
        tp = self.vars.tp
        t = self.vars.t
        sli_q_zero = self.__calc_q_zero(sli_k, tp, kappa, t)
        sli_sigma = self.__calc_slab_lab_int_sigma(sli_depth, sli_topo,
                                                   sli_temp, sli_s, sli_k,
                                                   sli_q_zero, v)
        return sli_sigma

    def __set_slab_sigma(self):
        z_sl = self.geo_model.get_boundaries()[3]
        z_topo = self.geo_model.get_boundaries()[0]
        sli_depth = self.geo_model.get_slab_lab_int_depth()
        sli_topo = self.geo_model.get_slab_lab_int_topo()
        sli_sigma = self.slab_lab_int_sigma
        d = self.vars.d
        slab_sigma = self.__calc_slab_sigma(z_sl, z_topo, sli_depth, sli_topo,
                                            sli_sigma, d)
        a = self.geo_model.get_slab_lab_int_area()
        b = a == 1
        slab_sigma = np.ma.array(slab_sigma, mask=b).filled(np.nan)
        return slab_sigma

    def __set_slab_temp(self):
        z_sl = self.geo_model.get_boundaries()[3]
        z_topo = self.geo_model.get_boundaries()[0]
        slab_k = self.geo_model.data.extract_2d_surface(z_sl, self.vars.k)
        tp = self.vars.tp
        kappa = self.vars.kappa
        v = self.vars.v
        t = self.vars.t
        q_zero = self.__calc_q_zero(slab_k, tp, kappa, t)
        dip = self.vars.dip
        b = self.vars.b
        s = self.__calc_s(z_sl, z_topo, kappa, v, dip, b)
        slab_sigma = self.slab_sigma
        slab_temp = self.__calc_slab_temp(z_sl, z_topo, q_zero, slab_sigma, v,
                                          slab_k, s)
        a = self.geo_model.get_slab_lab_int_area()
        b = a == 1
        slab_temp = np.ma.array(slab_temp, mask=b).filled(np.nan)
        return slab_temp

    def __set_lab_temp(self):
        z_sl = self.geo_model.get_boundaries()[3]
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
        slab_temp = self.slab_temp  # type: np.ndarray
        lab_temp = self.lab_temp  # type: np.ndarray
        slab_lab_temp = np.zeros(self.geo_model.data.get_2d_shape())
        slab_lab_temp[b] = slab_temp[b]
        slab_lab_temp[c] = lab_temp[c]
        return slab_lab_temp

    def __set_surface_heat_flow(self):
        delta = self.vars.delta
        z_topo = self.geo_model.get_boundaries()[0]
        z_sl = self.geo_model.get_boundaries()[3]
        slab_lab_k = self.geo_model.data.extract_2d_surface(z_sl,
                                                            self.vars.k)
        slab_lab_h = self.geo_model.data.extract_2d_surface(z_sl,
                                                            self.vars.h)
        temp_sl = self.slab_lab_temp
        heat_flow = self.__calc_surface_heat_flow(slab_lab_h, delta,
                                                  slab_lab_k, z_topo, z_sl,
                                                  temp_sl)
        return heat_flow

    def __set_geotherm(self):
        k = self.vars.k
        h = self.vars.h
        print('delta:', self.vars.delta)
        if self.vars.delta.shape:
            delta = self.vars.delta[:, :, np.newaxis]
        else:
            delta = self.vars.delta
        z = self.geo_model.data.get_3d_grid()[2]
        z_topo = self.geo_model.get_boundaries()[0][:, :, np.newaxis]
        z_sl = self.geo_model.get_boundaries()[3][:, :, np.newaxis]
        temp_sl = self.slab_lab_temp[:, :, np.newaxis]
        geotherm = self.__calc_geotherm(h, delta, k, z, z_topo, z_sl, temp_sl)
        return geotherm

    def get_geotherm(self):
        return self.geotherm


class MechanicModel(Model3d):

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


    def __init__(self, geo_model, thermal_model, m_input, rhe_data):
        super().__init__()
        self.geo_model = geo_model
        self.thermal_model = thermal_model
        self.vars = DotDict(self.__set_variables(m_input, rhe_data))
        self.bys_t, self.bys_c = self.__set_brittle_yield_strength()
        self.dys = self.__set_ductile_yield_strength()
        self.yse_t, self.yse_c = self.__set_yield_strength_envelope()
        self.uc,self.lc,self.lm,self.eet = self.__set_eet(self.yse_t)
        pass

    def __get_rheologic_vars_from_model(self, rock_id):
        rock = RheologicModel.objects.get(name=rock_id)
        rock_dic = {
            'n': rock.n,
            'a': rock.A,
            'h': rock.H,
        }
        return DotDict(rock_dic)

    def __get_rheologic_vars(self, id_rh, rhe_data):
        rock = rhe_data[str(id_rh)]
        rock_dic = {
            'h': rock[1],
            'n': rock[2],
            'a': rock[3]
        }
        return DotDict(rock_dic)

    def __set_variables(self, m_input, rhe_data):
        m_input = DotDict(m_input)
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

    def __set_brittle_yield_strength(self):
        bs_t = self.vars.bs_t
        bs_c = self.vars.bs_c
        depth = self.geo_model.data.get_3d_indexes()[2]  # depth from surface
        depth = self.geo_model.data.mask_3d_array(depth.astype(np.float64),
                                                  nan_fill=True)
        bys_t = self.__calc_brittle_yield_strength(bs_t, depth)
        bys_c = self.__calc_brittle_yield_strength(bs_c, depth)
        return bys_t, bys_c

    def __set_ductile_yield_strength(self):
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
        bys_t = self.bys_t  # type: np.ndarray
        bys_c = self.bys_c  # type: np.ndarray
        dys = self.dys  # type: np.ndarray
        with np.errstate(invalid='ignore'):
            yse_t = np.where(bys_t < dys, bys_t, dys)
            yse_c = np.where(bys_c > dys, bys_c, dys)
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
        elastic_z = self.geo_model.data.cut_3d_array(elastic_z,
                                                     top_z=top_z,
                                                     bottom_z=bottom_z)
        elastic_z = np.ma.array(elastic_z, mask=np.isnan(elastic_z))
        top_elastic_z = np.amax(elastic_z, axis=2).filled(np.nan)
        bottom_elastic_z = np.amin(elastic_z, axis=2).filled(np.nan)
        elastic_thickness = top_elastic_z - bottom_elastic_z
        layer_tuple = np.stack((top_elastic_z, bottom_elastic_z,
                                elastic_thickness), axis=2)
        return layer_tuple, elastic_z.filled(np.nan)

    def __set_eet(self, yse):
        elastic_z = np.copy(self.geo_model.data.get_3d_grid()[2])
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
        attached_ths, detached_ths = self.divide_by_areas(elastic_thickness,
                                                          share_moho)
        detached_ths_2l, detached_ths_3l = self.divide_by_areas(detached_ths,
                                                                share_icd)
        detached_ths_2l[:, :, 0] = (detached_ths_2l[:, :, 0]
                                    + detached_ths_2l[:, :, 1])
        detached_ths_2l[:, :, 1] = np.nan

        detached_ths = self.combine_by_areas(detached_ths_2l, detached_ths_3l,
                                             share_icd)

        attached_eet = self.__calc_eet_attached(attached_ths)
        detached_eet = self.__calc_eet_detached(detached_ths)

        eet = self.combine_by_areas(attached_eet, detached_eet, share_moho)

        return uc_tuple, lc_tuple, lm_tuple, eet


    def get_yse(self):
        return self.yse_t, self.yse_c
        # ##WORk IN PROGRESS###

def compute(gm_data, ta_data, rhe_data, areas, t_input, m_input):
    D = Data(gm_data, 0.2, 0.1)
    GM = GeometricModel(D, areas)
    TM = ThermalModel(GM, t_input, ta_data)
    # D = C.get_geotherm()[45, 35, 50]
    MM = MechanicModel(GM, TM, m_input, rhe_data)
    return D, GM, TM, MM

def get_old_output(t_input, m_input, areas):
    A = Data(gm_data, 0.2, 1)
    B = GeometricModel(A, areas)
    C = ThermalModel(B, t_input)
    D = MechanicModel(B, C, m_input)
    os.chdir('../misc/prints')
    if not os.path.exists('perfs-new'):
        os.makedirs('perfs-new/Input', exist_ok=True)
        os.makedirs('perfs-new/Output/0_Termal/Perfiles', exist_ok=True)
        os.makedirs('perfs-new/Output/0_Termal/0_Mecanico/Perfiles',
                    exist_ok=True)
    os.chdir('perfs-new')
    for n in range(A.get_2d_shape()[1]):
        x_axis = A.get_axes()[0][:, np.newaxis]
        y_axis = np.repeat((-10-(n)*0.2), len(x_axis))[:, np.newaxis]
        print(y_axis)
        sl = B.get_boundaries()[3][:, n][:, np.newaxis]
        moho = B.get_boundaries()[2][:, n][:, np.newaxis]
        icd = B.get_boundaries()[1][:, n][:, np.newaxis]
        topo = B.get_boundaries()[0][:, n][:, np.newaxis]
        areas = A.mask_2d_array(B.get_slab_lab_int_area(), nan_fill=True)
        area = areas[:, n][:, np.newaxis]
        area[area == 1] = 2
        area[area == 0] = 1
        perfil = np.append(x_axis, y_axis, axis=1)
        perfil = np.append(perfil, sl, axis=1)
        perfil = np.append(perfil, moho, axis=1)
        perfil = np.append(perfil, icd, axis=1)
        perfil = np.append(perfil, topo, axis=1)
        perfil = np.append(perfil, area, axis=1)
        q = np.repeat(np.nan, len(x_axis))[:, np.newaxis]
        s = np.repeat(np.nan, len(x_axis))[:, np.newaxis]
        geotherm = C.get_geotherm()[:, n, :]
        perfilT = np.append(perfil, q, axis=1)
        perfilT = np.append(perfilT, s, axis=1)
        perfilT = np.append(perfilT, geotherm, axis=1)
        yse_T = D.get_yse()[0][:, n, :]
        yse_C = D.get_yse()[1][:, n, :]
        perfilM_T = np.append(perfil, yse_T, axis=1)
        perfilM_C = np.append(perfil, yse_C, axis=1)
        num = 1000 + n * 20
        os.chdir('Input')
        np.savetxt('perf-' + str(num), perfil, fmt='%11.4e', delimiter="   ")
        os.chdir('../Output/0_Termal/Perfiles')
        copy2('../../../../../t_data.pkl', '../')
        np.savetxt('perf-' + str(num), perfilT, fmt='%11.4e', delimiter="   ")
        os.chdir('../../../Output/0_Termal/0_Mecanico/Perfiles')
        copy2('../../../../../../m_data.pkl', '../')
        np.savetxt('perfT-' + str(num), perfilM_T, fmt='%11.4e',
                   delimiter="   ")
        np.savetxt('perfC-' + str(num), perfilM_C, fmt='%11.4e',
                   delimiter="   ")
        os.chdir('../../../../')
        if n == 174:
            break


# get_old_output(t_data, m_data, areas)
#B = compute(t_data, m_data)
#geotherm = B.thermal_model.get_geotherm()
#yse_t, yse_c = B.get_yse()
#imageToVTK("./thermomecanic", spacing=(10.0, 10.0, 1.0),
#           pointData={"temperature": geotherm, "yield_strength": yse_t})
