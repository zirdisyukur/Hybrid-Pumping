# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 15:53:00 2025

@author: Zirdi
"""

import numpy as np
from one_ring.base.pint_units import ureg
from one_ring.base.pint_units import default_unit
from scipy import constants
def absorption(fit_result):
    """
    given a absorption spectrum fit result, use the fit parameters to estimate the absorption of this signal,
    a nominal improvement over just taking the absorption = log(data.min/data.max).

    returns absorption, unitless
    
    """
    A = fit_result.fit_parameters['A']
    f0 = fit_result.fit_parameters['f0']
    HWHM = fit_result.fit_parameters['HWHM']
    Offset = fit_result.fit_parameters['Offset']
    absorption = (A / HWHM) - A * HWHM / (f0**2 + HWHM**2)

    return absorption

def cross_section(HWHM, transition):
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
    HWHM = default_unit(HWHM, ureg.THz).to(ureg.THz)
    if transition == 'K/D1':
        f_osc = 0.338
    if transition == 'Rb/4S_4P0.5':
        f_osc = 0.004
    cs = (c * r_e * f_osc / HWHM).to(ureg.cm**2)
    return cs

def potassium_density(T):
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
    return n.to(ureg.cm**-3)