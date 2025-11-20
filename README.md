# Hybrid Pumping Project
this repository is a log of the progress of the hybrid pumping project. Developed by Zirdi.

***

# Experimental Overview
Hybrid pumping is a technique where high polarization of a atomic species in a vapor cell is achieved by optically pumping a technically easier "donor" species, then relying on spin exchange to transfer that polarization to a "receiver" species. When the system is in spin-temperature equilibrium, the polarization of the two species should be equal.

The goal of the project is to experimentally demonstrate how the relative density of the two species prevents spin-temperature equilibrium in an atomic gas. Our donor species of choice is Potassium and receiver species is Rubidium. We intend to do this by manufacturing atomic cells filled with different density ratios $D = \frac{\rho_{Rb}}{\rho_{K}}$ and measuring the polarization of the two species. 

Where $D$ is small, we expect the polarization of the two species to be to equal, $P_K = P_{Rb}$. 

Where $D$ is large, we expect the polarization to diverge, and the $P_{K} > P_{Rb}$

## Cells Overview

|Cell Name|RbCl [g]|KCl [g] |$P_{N2}$[Torr]|$P_{He}$[Torr]| $D_{est}$ [T=190] |
|-------- |------- |------  |------|------|------|
| Elsa    |6.1     |0.55    |150   |500   |43.31 |
| Anna    |6.1     |0.55    |150   |500   |43.31 | 
| Olaf    |5.6042  |0.104   |-     |-     |210.45|

$D$ was estimated using the following formula [citation needed], Then $\rho_K$, and $\rho_{Rb}$ were calculated using Raoult's law. See HybridPumping_DataAnalysis.ipynb for more details
$$
D(T) = 10^{\frac{413 K}{T} - 0.09} * \frac{m_{rb}}{m_k}
$$

At a given temperature, the density of the two species is fixed, irrespective of mole fractions. These are also empirically calculated from the same source of the above equation
$$
\rho_{Rb}(T=190\degree C) = 6.162\cdot10^{14} { cm}^{-3}\\
\rho_K(T=190\degree C)= 9.728\cdot10^{13} { cm}^{-3}
$$

# Density Measurements
Density measurements, broadly speaking is done through absorption spectroscopy. By scanning the frequency of a known transition, while monitoring the power transmission across a hot vapor, we can plot an absorption profile, fit it to extract its half-width-half-max (HWHM), and Absorbance (A), andthen estimate the density of the species using Beer-Lambert's law.

$$
\rho_{Rb} = \frac{A}{\sigma L} = \frac{log(P/P_0)}{cr_ef_{osc}} \cdot \text{HWHM}
$$

Where $\sigma$ is the cros section given by the speed of light $c$, $r_e$ is the classical radius of the electron, $f_{osc}$ is the oscillator strength of a given transition. 

## Potassium Absorption

We excite the Potassium's D1, which is $4S \rightarrow 4P_{1/2}$ transition which has an oscillator strength $f_{osc-rb}$ = 0.338.

The transition should lie in 390 THz, or 770 nm in vacuum.

I measured the potassium density with Olaf. This measurement was done with an average calibration curves consisting of 18 scans done when the cell was cold, measured at different points over the course of 6 days. This was necessary because there is a discrepancy in response between the reference powermeter and the transmission powermeter due to their different models. This helped minimize how much power-fluctuations from the drift and scanning of the laser would affect our transmission measurements. The transmission measurements were 8 different measurements across 2 days that all used the same average calibration curve. All measurements done at a probe current of 55mA.

The fitt and estimated parameters are as follows.

![alt text](<k_olaf_abs.png>)

### Measurement results:
$
\overline{\rho_K} = (3.7 \pm 0.2)10^{11} \:\text{ cm}^{-3}\\\\
\overline{{HWHM}_K} = (4.8 \pm 0.2) \:\text{GHz}\\\\
\overline{f_{K}} = (389.28428 \pm 0.00007) \:\text{THz}
$

As you can see, the potassium density is nearly $10^3$ orders of magnitude smaller than expected. This may be attributed to a poor mixing procedure when making the cells, preventing a sufficient amount of K entering the cell. This is one of the current issues of this project. If we consider this at face value, then the $D$ of this cell is 55,000!

We have progressed by assuming this is true, and trying to observe the optical rotation of this cell, and comparing the polarizations.

## Rubidium Absorption

The Rubidium is dense and too optically thick on the D1 transition, so instead we use the $5S_{1/2} \rightarrow 6P_{1/2}$ transition which has an oscillator strength $f_{osc-rb}$ = 0.004.

The transition should lie in 710 THz, or 421 nm in vacuum.

I measured the rubidium absorption using a beam splitter configuration and two different powermeters. A transmission and reference scan when the cell was cold and then another transmission and refernece scan when the cell was hot. This allowed us to gain match the signals to unity absorption off-resonance, and calibration for mismatched power-sensitivity from the two different powermeters. Because the Rb laser had rather trong modehopes periodically across the scan, we had to omit those regions from the data analysis, which explains some of the residual inconsistencies off resonance.

![alt text](<rb_olaf_abs.png>)

### Measurement results:
$
\overline{\rho_{Rb}} = (2.64 \pm ??)10^{14} \:\text{ cm}^{-3}\\\\
\overline{{HWHM}_{Rb}} = (15.34 \pm 0.08) \:\text{GHz}\\\\
\overline{f_{Rb}} = (710.96165 \pm 0.00005) \:\text{THz}
$

# Polarization measurements
## Method 1
This method can be generally characterized by measuring the optical rotation as a function of pumping power. I can then fit this curve to a fitt that estimates P ~= P_up / ( P_up + P_dwn ), where P_up is proportional to the pumping power and Pd_dwn is proportional to the relaxation rate. Notice in absence of any relaxation, polarization equals to one. This method works well for non-hybrid cells as the pumping power is the pumping of the one speices and the relaxation power is the same species. In hybrid cells this is not the case as only the donor's pumping power is input to the fit, but the estimated relaxation is that of the receiver species. 

To do this measurement follow these steps
1. Ensure that you have the probe and pump optics aligned with the cell. Use the iris to align the cell to the oven
2. Ensure that the lasers also have a good beam shape, i.e. comparable to the size of the cell
3. Ensure that the balanced photodetector has maximal amount of light into both of its photodiodes (optimize by monitoring voltage of each diode). 
4. Optics have been finalized and for the purpose of this measurement, they should not change.
5. Perform optical rotation calibration
    1. Calibration is dependent on probe power (standardized to power just before it hits the cell). If you are using a laser that drifts in power on a timescale on the order of the time it takes to do this measurement, it is wise to measure the optical calibration at different powers and make a fit out of that. This calibration will hold for any power up to any changes in optical geometry (step 4). If not though, it is sufficient to standardize to one power level, whatever power level minimizes probe realaxation
    1. Make sure to turn off any pumping and field compensation, no need to heat up the cell.
    2. To calculate the calibration you must use the vernier rotation mount for a polarizer , adjust a half wave plate in front of it in order to normalize power across the measurement. 
    3. Measure the differential voltage of the photodiode at that particular angle. (try to get the full -10V to 10V range for the calibration). 3-4 points should be sufficient for each probe power.
6. Perform absorption measurement
    1. Heat up the cell, dont' turn on any pump beam or field compensation.
    1. If not done absorption measurement already, it is important to know the temperature of the probe that corresponds to the transition (especially if it has drifted away), and to also know the HWHM of the signal, this will become usefull in calculating the polarization later. This measurement can be done just once as it usually is irrespective of optical geometry.
    2. First heat up the cell to the desired The absorption measurement process can vary depending on technicalities but overall consider doing: 
        Submethod 1: a single photodiode measuring a reference beam in front of oven and a transimssion beam after oven and normalize. Good for cells with strong densities
        Submethod 2: a reference beam into a photodiode aftet the oven, done with cold atoms, Then a transmission beam into the same photodiode after the oven, but done with hot atoms. Good eliminating any scattering from the oven
        Submethod 3: use two photodiodes, one with a reference beamsample before oven and one after the oven. You may potentially need to calibrate for the two photodiode responses if they havef diffeent responsitivities. This is done by making many many cold cell calibration measurements. 
7. Perform Polarization method 1 and 2
    1. Set probe power to that which was done in optical rotation calibration and take note of power
    2. Turn on cell heaters, pump beam, field compensation (both cancellation, storng quantization axis field, and weak adiabatic oscillation in the transverse plane). Ensure that pump circular polarization is good, and that the strength of the magnetic fields are good (usually a 50uT DC and 1 uT AC field). You'll want to ensure that the pump is capable of going to around 40mW of power, and when you start to measure R (optical rotation) from the lock in, you can also maximize R by varying the pump temperature. Note down all of these quantities that are constant throughout the experiment.
    3. Insert the *right* polarizer in pump optics so that you can attenuate the pump power. 
    4. Ensure the probe temperature is maybe around 5HWHM away from the transition, you may have to scan across the cell in a quick absorption measurement to find the peak. Also note down the length of cell (estimate can be comparing to the estimate given from method 2, otherwise it should be just some value slightly smaller than 1). You may also need the probe laser temperature to frequency calibration for the detuning parameter. After finding the transition temperature, just use the calibration to find where 5HWHM away from it is. 
    5. Polarization method 1: Vary the pump power from 1mW to 40mW and record the R value. Use one_ring and fitt the curve in order to get polarization (only for potassium polarization measurement)
    6. Polarization method 2: Take raw data from method one, and calibrate using cell length, optical calibration, detuning parameter, the field strengths, and density. (this can be compared to a good reference done with method one to determine things like the effective cell length)


## Authors and acknowledgment
M. Zirdi Syukur, Lee Junyi, and the consultations of Md. Danial Afiq Bin Abdullah and Ran Qiandong

## License
--

## Citations
[1] = Chen et. al. Spin-exchagne optical pumping of 3He with Rb-K mixtures and pure K (2007)
## Project status
This project is still ongoing.
