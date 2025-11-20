# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 15:53:00 2025

@author: Zirdi
"""

import numpy as np
from one_ring.base.pint_units import ureg
from one_ring.base.pint_units import default_unit
from scipy import constants
import one_ring.base.h5py_tools as h5t
import one_ring.base.fitting_tools as fitt
import one_ring.base.pandas_tools as pdt
import pickle

#constants useful for polariation measurements
bohr_magneton = constants.physical_constants['Bohr magneton'][0] * ureg.J/ureg.T
hbar = constants.hbar * ureg.J * ureg.s / ureg.rad
r_e = constants.physical_constants['classical electron radius'][0] * ureg.m
c = constants.c * ureg.m / ureg.s
f_K = 0.338
f_Rb = 0.004
g_factor_K = 0.5006
gyromagnetic_ratio_K = (bohr_magneton / hbar * g_factor_K).to(ureg.kHz/ureg.uT*ureg.rad)
#potassium absorption cross section 

#temperature and frequency calibration for probe laser at 55mA
probe_calib = (-35.579346 * ureg.delta_degC / ureg.THz) #
potassium_f0_T = (16.048 * ureg.delta_degC) 
HWHM = 0.170 * ureg.delta_degC / probe_calib
rot_calib = 1.1274997 * ureg.V/ureg.deg #at 55mA
density_K = 3.7e11 * ureg.cm**-3

#wulfgang data
wulf_HWHM = 0.00427070309215070 * ureg.THz #from fit_absorption_spectrum measurement by junyi a while ago

class polarization_data():
    
    def __init__(self,file_name, group_name, style, data_mask):
        """
        Initializes a class representing a polarization measurement. After initializing with
        the appropraite h5 data location, the class has various methods and variable calculations
        useful for the experiment.

        Parameters:
        -----------
        file_name : str
            A string with the the file location.
        group_name : str
            A string of the group that the data is within the h5 hierarchy
        style : str
            Because data has been taken either with pdt.save_dataframe or pdt.save_dataset, data is either
            stores as a dataframe os a np.ndarray. Write style = 'numpy' or 'pandas' in order to specify
            which it is for the selected dataset
        data_mask : tuple, (#, #)
            tuple consisting of two numbers s.t. it only considers datapoints taken
            at pump powers above data_mask[0] and below data_mask[1]
        """
        self._file_name = file_name
        self._group_name = group_name
        self._style = style
        self._data_mask = data_mask
        self._fit_result = None
        self._polarization = None
        
        if (style == 'pandas'):
            self._total_dataframe = pdt.load_dataframe(file_name, group_name)
            self._metadata = pickle.loads(h5t.get_node_metadata(file_name, group_name)['metadata'])
            
            if (data_mask != None):
                lower_bound = data_mask[0]
                upper_bound = data_mask[1]
                mask = (self._total_dataframe.index > lower_bound) & (self._total_dataframe.index < upper_bound)
                self._data = self._total_dataframe[mask]
                self._excluded_data  = self._total_dataframe[~mask]       
            else:
                self._data = self._total_dataframe
                
            self._R = self._data['R [V]'].to_numpy() * ureg.V
            self._Power = self._data.index.to_numpy() * ureg.mW

            self._excludedR = self._excluded_data['R [V]'].to_numpy() * ureg.V
            self._excludedPower = self._excluded_data.index.to_numpy() * ureg.mW

        
        elif (style == 'numpy'):
            self._data = h5t.load_group(file_name, group_name)
            self._metadata = h5t.get_node_metadata(file_name, group_name)

            if (data_mask != None):
                lower_bound = data_mask[0]
                upper_bound = data_mask[1]
            else: 
                lower_bound = -9999
                upper_bound = 9999
            
            self._mask = (self._data['pump_pwr'].m > lower_bound) & (self._data['pump_pwr'].m < upper_bound)
            self._Power = self._data['pump_pwr'][self._mask]
            self._X = self._data['x'][self._mask]
            self._Y = self._data['y'][self._mask]
            self._R = np.sqrt(self._X**2 + self._Y**2)

            self._excludedR = np.sqrt(self._data['x'][~self._mask]**2 + self._data['y'][~self._mask]**2)
            self._excludedPower = self._data['pump_pwr'][~self._mask]
            

        #defining properties if it can be pulled from the metadata,
        #otherwise defining them as none
        self._probe_temp = default_unit(self._metadata['probe_temp'].to(ureg.degC).m, ureg.degC) if ('probe_temp' in self._metadata) else None
        self._Bx_field_AC_amp = default_unit( self._metadata['Bx_field_AC_amp'].to(ureg.uT).m, ureg.uT) if ('Bx_field_AC_amp' in self._metadata) else None
        self._By_field = default_unit( self._metadata['By_field'].to(ureg.uT).m, ureg.uT) if ('By_field' in self._metadata) else None
        self._Bz_field = default_unit( self._metadata['Bz_field'].to(ureg.uT).m, ureg.uT) if ('Bz_field' in self._metadata) else None
        self._field_tilt = np.arctan( self._Bx_field_AC_amp / self._By_field )
        self._detune = detune(self.probe_temp.m * ureg.delta_degC) if (self._probe_temp is not None) else None
        self._cell_temp = self._metadata['cell_temperature'] if ('cell_temperature' in self._metadata) else None
 
    @property
    def probe_temp(self):
        return self._probe_temp

    @property
    def field_tilt(self):
        return self._field_tilt
    
    @property
    def R(self):
        return self._R
    
    @property
    def Power(self):
        return self._Power

    @property 
    def metadata(self):
        return self._metadata

    @property
    def file_location(self):
        return self._file_name, self._group_name

    @property
    def cell_temperature(self):
        return self._cell_temp

    def plot_fit(self, ax=None, excluded_data=False, dlabel='data', flabel='fit', guess=np.array([0.5,80,1])):
        """
        Function that will help plot the raw data and the best fit.

        Parameters:
        -----------
        ax : str
            The ax from matplotlib that will be used to plot. Defaults to plt
        excluded_data : Boolean
            If true, will plot excluded data that has been masked by self._mask
        dlabel : str
            data label
        flabel : str
            fit label
        """
        if ax is None:
            ax = plt
        
        #fit data
        pump_fitter = fitt.PolarizationPumpFitter()
        fit_result = pump_fitter.fit( self._Power.m, self._R.m, guess=guess)
        x_interp_arr = np.arange(np.min(self._Power.m),np.max(self._Power.m),0.5) * ureg.mW
        fitt_data = pump_fitter.fit_function(x_interp_arr.m, **fit_result.fit_parameters)
        
        polarization_estimate = fit_result.polarization(x_interp_arr.m)
        #plot data
        ax.plot(self._Power.m, self._R.m, label=dlabel,marker='o')
        ax.plot(x_interp_arr.m, fitt_data,label=flabel)
        if(excluded_data):
            ax.errorbar(self._excludedPower.m, self._excludedR.m,
                        marker='^', ls='none', label='excluded data')
        ax.set_xlabel('Pump Power [mW]')
        ax.set_ylabel('Optical Rotation [' + str(self._R.u) + ']')
        ax.set_title('Raw Measurements and Fit')
        ax.legend()

        self._fit_result = fit_result
    
    def plot_potassium_polarization(self, ax=None, label=None ,guess=np.array([0.5,2,-1])):
        """
        Function that will help plot the polarization based on the best fit.

        Parameters:
        -----------
        ax : str
            The ax from matplotlib that will be used to plot. Defaults to plt
        label : str
            data label
        guess : str
            guessing parameters for the fit for A, R and R0
        """
        if ax is None:
            ax = plt
        
        if self._fit_result is None:
            fit_result = pump_fitter.fit( self._Power, self._R, guess=guess)
            self._fit_result = fit_result

        polarization = self._fit_result.polarization(self._Power.m)
        self._polarization = polarization
        ax.plot(self._Power.m, polarization,label=label)
        ax.set_xlabel('Pump Power [mW]')
        ax.set_ylabel('K-Polarization')
        ax.set_title('K-Polarization')
        
        if label is not None:
            ax.legend()

def detune(nu):   
    nu = default_unit(nu, ureg.delta_degC).to(ureg.delta_degC) / probe_calib
    delta_nu = nu - potassium_f0_T / probe_calib #temperature of f0
    return np.abs(delta_nu / (HWHM**2 + delta_nu**2))

def polarization_K_calib(data, temp, tilt, rot_calib, density, Lx):
    d = detune(temp)
    optical_rotation = np.sqrt(2) * data / rot_calib
    denom = density * r_e * c * f_K * Lx * tilt * d
    return (2 * optical_rotation / denom).to_base_units()

def estimate_cell_length(data, temp, tilt, rot_calib, density, Pz, noMean = False):
    """
    Given some parameters, estimate the length of the cell
    Parameters
    ----
    data : np.ndarray, or ureg.Quantity
        Raw optical rotation data from lock in, units of Volts
    temp : float, or ureg.Quantity
        probe laser temperature in units of degC
    tilt : float, or ureg.Quantity
        angle of field tipping in radians
    rot_calib : float, or ureg.Quantity
        calibration that converts lock in voltage to rotation angle, in units of rad/V
    density : float, or ureg.Quantity
        density of atoms in units of 1/cm^3
    Pz : float, or ureg.Quantity
        fitted polarization from method 1

    Returns
    --------
    Lx : ureg.Quantity
        Effective length of cell in units of cm
        
    """
    rot_calib = default_unit(rot_calib, ureg.rad/ureg.V)
    density = default_unit(density, ureg.cm**-3)
    
    d = detune(temp)
    optical_rotation = np.sqrt(2) * data * rot_calib
    denom = density * r_e * c * f_K * Pz * tilt * d
    if (noMean):
        return (2 * optical_rotation / denom).to(ureg.cm)
    return (2 * optical_rotation / denom).to(ureg.cm).mean()


    