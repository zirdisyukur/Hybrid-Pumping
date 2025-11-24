# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 15:53:00 2025

@author: Zirdi
"""

import numpy as np
from one_ring.base.pint_units import ureg
from one_ring.base.pint_units import default_unit
import one_ring.base.h5py_tools as h5t
import one_ring.base.fitting_tools as fitt
from scipy import constants

r_e = constants.physical_constants['classical electron radius'][0] * ureg.m
c = constants.c * ureg.m / ureg.s
f_K = 0.338
f_Rb = 0.004

class absorption_data():
    def __init__(self,file_name, group_name, f_osc, style='numpy', data_mask=None):
        """
        Initializes a class representing a absorption measurement. After initializing with
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
        data_mask : ndarray
            ndarrauy of the same shape as the set. It consists of booleans in every entry where True is values that
            is to be included.
        """

        self._file_name = file_name
        self._group_name = group_name
        self._style = style
        self._data_mask = data_mask
        self._fit_result = None
        self._fit_params = None
        self._probe_calib_intercept = 389.74896028618315 * ureg.THz
        self._reference_scan = False
        self._normalized_data = False
        self._fit_xunit = None
        self._f_osc = f_osc
        
        if (style == 'pandas'):
            return 'dataframes not usually used for this type of measurment'

        elif (style == 'numpy'):
            self._rawdata = h5t.load_group(file_name, group_name)
            self._metadata = h5t.get_node_metadata(file_name, group_name)
            self._mask = data_mask
            
            if (data_mask != None):
                self._data = self._rawdata[self._mask]
                self._excluded_data = self._rawdata[~self._mask]
                self._excludedTemperature = self._excluded_data['temperature']
                self._excludedPower = self._excluded_data['power']
            else: 
                self._data = self._rawdata
                self._excludedTemperature = None
                self._excludedPower = None
            
            
            self._Temperature = self._data['temperature']
            self._Power = self._data['power']

            

        #defining properties if it can be pulled from the metadata,
        #otherwise defining them as none
        self._probe_temp = default_unit(self._metadata['probe_temp'].to(ureg.degC).m, ureg.degC) if ('probe_temp' in self._metadata) else None
        self._detune = detune(self.probe_temp.m * ureg.delta_degC) if (self._probe_temp is not None) else None
        self._cell_temp = self._metadata['cell_temperature'] if ('cell_temperature' in self._metadata) else None
        self._probe_laser_current = self._metadata['probe_laser_current'] if ('probe_laser_current' in self._metadata) else None
        if ('probe_calib' in self._metadata):
            self._probe_calib = self._matadata['probe_calib']
        elif self._probe_laser_current == 55*ureg.mA:   
            self._probe_calib = (-0.028106193969950993 * ureg.THz / ureg.delta_degC) #from old measurement with WS7 i think
        else:
            self._probe_calib = None
        if (self._cell_temp.to(ureg.degC).m < 30):
            self._reference_scan = True
        
    @property 
    def metadata(self):
        return self._metadata

    @property
    def file_location(self):
        return self._file_name, self._group_name

    @property
    def cell_temperature(self):
        return self._cell_temp

    @property 
    def Y(self):
        return self._Power
        
    def X(self, unit='temp'):
        """
        write unit
        """
        if (unit == 'temp'):
            return self._Temperature.to(ureg.degC)
        elif (unit == 'delta_temp'):
            return self._Temperature.to(ureg.degC).m * ureg.delta_degC
        elif (unit == 'freq'):
            return (self._Temperature.to(ureg.degC).m * ureg.delta_degC) * self._probe_calib + self._probe_calib_intercept

    def fit_absorption(self, ref, xunit='freq', ax=None, excluded_data=False, dlabel='data', flabel='fit', tlabel='Raw Measurements and Fit', guess=None):
        """
        Function that will help plot the raw data and the best fit.

        Parameters:
        ----------------
        ref : absorption_data object
            absorption data object that has reference scan data loaded
        ax : str
            The ax from matplotlib that will be used to plot. Defaults to plt
        xunit : str
            options include 'delta_temp' or 'freq' for the units of the x-axis 
        excluded_data : Boolean
            If true, will plot excluded data that has been masked by self._mask
        dlabel : str
            data label
        flabel : str
            fit label
        """
        if (self._reference_scan):
            return "not hot cell data"
            
        if ax is None:
            ax = plt
            
        if (guess == None):
            if (xunit=='temp'):
                guess = [0.008,16,0.121,0.0]
            if (xunit=='freq'):
                guess = [0.008,389.3,0.1,0.0]
        else:
            guess = guess

        normalized_data = self.Y / ref.Y
        self._normalized_data = normalized_data
        absorption_data = -np.log(normalized_data)
            
        #fit data
        abs_fitter = fitt.AbsorptiveLorentzianFitter()
        fit_result = abs_fitter.fit( self.X(xunit).m, (absorption_data).m, guess=guess)
        self._fit_xunit = xunit
        xmax = np.max(self.X(xunit))
        xmin = np.min(self.X(xunit))
        x_step = (xmax - xmin) / 200
        x_interp_arr = np.arange(xmin.m, xmax.m, x_step.m ) * self.X(xunit).u

        fit_params = fit_result.fit_parameters
        fit_params['C'] = fit_params.pop('Offset')
        #plot data
        ax.plot(self.X(xunit).m, (absorption_data.m), label=dlabel,marker='o')
        ax.plot(self.X(xunit).m, fit_result._fit_function(self.X(xunit).m,**fit_params), label=flabel)

        #if(excluded_data):
        #    ax.errorbar(self._excludedPower.m, self._excludedR.m,
       #                 marker='^', ls='none', label='excluded data')
        if (xunit == 'temp'):
            ax.set_xlabel('Probe Temperature [degC]')
        if (xunit == 'freq'):
            ax.set_xlabel('Probe frequency [{}]'.format(self.X(xunit).u))
        ax.set_ylabel('Absorption')
        ax.set_title(tlabel)
        ax.legend()

        self._fit_params = fit_params
        self._fit_result = fit_result
        
    def absorption(self):
        """
        given a absorption spectrum fit result, use the fit parameters to estimate the absorption of this signal,
        a nominal improvement over just taking the absorption = log(data.min/data.max).
    
        returns absorption, unitless
        
        """
        A = self._fit_params['A']
        f0 = self._fit_params['f0']
        HWHM = self._fit_params['HWHM']
        Offset = self._fit_params['C']
        absorption = (A / HWHM) - A * HWHM / (f0**2 + HWHM**2)
    
        return absorption

    def cross_section(self):
        """
        Given a absorption spectrum's HWHM, and a particular transition with a known oscillator strength,
        this method computes the absorption cross section.
        
        Parameters
        ----------
        HWHM : float or ureg.Quantity, optional
            half width half max of absorption signal, assumed to be in THz
        transition: string
            options include 'K/D1' or 'Rb/4S_4P0.5'
    
        Returns
        -------
        cs : ureg.Quantity
            absorption cross section in units of cm**2
            
        """
        if (self._fit_xunit == 'freq'):
            HWHM = np.abs(default_unit(self._fit_result.fit_parameters['HWHM'], ureg.THz).to(ureg.THz))
        if (self._fit_xunit == 'temp'):
            HWHM_temp = default_unit(self._fit_result.fit_parameters['HWHM'], ureg.delta_degC).to(ureg.delta_degC)
            HWHM =  np.abs((HWHM_temp) * self._probe_calib)
        
        cs = (c * r_e * self._f_osc / HWHM).to(ureg.cm**2)
        return cs

    def density(self, L, order_notation=False):
        """
        Given a the density of the absorption measurement
        
        Parameters
        ----------
        L : float or ureg.Quantity
            length of optical path, defaults to ureg.cm units
            
        Returns
        -------
        density : ureg.Quantity
            number density of atoms in units of 1/cm^3
            
        """
        L = default_unit(L, ureg.cm)
        density_K = self.absorption() / (L * self.cross_section())
        if(order_notation):
            return print('{:.3e}'.format(density_K))
        return density_K.to(1/ureg.cm**3)

def potassium_density(T, order_notation=False):
    """
    Returns the density of potassium at a given temperature in units of kelvin.
    Assumes that cell only has potassium.
    
    Parameters
    ----------
    T : float or ureg.Quantity, optional
        Temperature, assumed to be in celsius

    Returns
    -------
    n : ureg.Quantity
        number density of potassium in units 1/cm^3
        
    """

    T = default_unit(T, ureg.degC).to(ureg.K)
    
    logp = 7.4077 - ureg.Quantity(4453,ureg.K) / T #nist, reports in mbar
    p = ureg.Quantity(10**logp, ureg.mbar)
    k_b = constants.physical_constants['Boltzmann constant'][0] * ureg.J / ureg.K
    n = p / k_b / T
    if(order_notation):
            return print('{:.3e}'.format(n.to(ureg.cm**-3)))
    return n.to(ureg.cm**-3)