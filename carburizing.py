#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 17:38:11 2023

This script attempts to simulate carbon diffusion depth in austenite 
with the pydiffusion library by Zhangqi Chen. The original library is 
very handy and useful, but a carburizing step needs a source term at the 
surface (here: left side of interval). 
This can be accomplished with a simple modification of the original sphSim 
function (see function sphSim_source below).
However, the script does not take into account carbon flux from the
atmosphere into the surface, which is dependend on more parameters and
very likely requires experiments. 

This script can easily adapted to more steps, also at different temperatures.

Based on diffusion.py in pydiffsion examples, but with sphSim and modified
sphSim_source

Carbon diffusion coefficients D(T) in austenite from:
Thibaux, P., Métenier, A. & Xhoffer, C. Carbon Diffusion Measurement in 
Austenite in the Temperature Range 500 °C to 900 °C. Metall Mater Trans 
A 38, 1169–1176 (2007). https://doi.org/10.1007/s11661-007-9150-5

@author: emefff
"""

import numpy as np
import matplotlib.pyplot as plt
from pydiffusion import DiffSystem
from pydiffusion import step, mesh
from pydiffusion.simulation import sphSim
from pydiffusion import profileplot #, DCplot
# from pydiffusion.io import save_csv

# for spSim_source_at_surf
from pydiffusion.core import DiffProfile
from scipy.interpolate import splev

###############################################################################
def sphSim_source(profile, c_source_at_surf, diffsys, time, \
                  output=True, name=''):
    """
    Single-Phase Diffusion Simulation with source at left side of interval.
    Modified version of sphSim in pydiffusion library.

    Parameters
    ----------
    profile : DiffProfile
        Initial diffusion profile before simulation.
    c_source_at_surf: float
        concentration at surface during diffusion. It acts as a source for
        the diffusing species.
    diffsys : DiffSystem
        Diffusion coefficients.
    time : float
        time in seconds.
    output : boolean, optional
        Print simulation progress, default = True.
    name : str, optional
        Name the output DiffProfile.

    Returns
    -------
    profile : DiffProfile
        Simulated diffusion profile

    """
    if name == '':
        name = diffsys.name+'_%.1fh' % (time/3600)

    dis, Xs = profile.dis.copy()/1e6, profile.X.copy()
    fD = diffsys.Dfunc
    d = dis[1:]-dis[:-1]
    dj = 0.5*(d[1:]+d[:-1])
    t, m = 0.0, 0
    while t < time:
        Xm = 0.5*(Xs[1:]+Xs[:-1])
        DCs = np.exp(splev(Xm, fD[0]))
        dt = min(d**2/DCs/2)
        dt = time-t if t+dt > time else dt*0.95
        t += dt
        m += 1
        Jf = -DCs*(Xs[1:]-Xs[:-1])/d
        Xs[1:-1] = Xs[1:-1]-dt*(Jf[1:]-Jf[:-1])/dj
        # Xs[0] -= Jf[0]/d[0]*dt*2 # ORIGINAL LINE
        #####################################################################
        Xs[0] = c_source_at_surf  # NEW LINE, C is constant the whole time
        #####################################################################
        Xs[-1] += Jf[-1]/d[-1]*dt*2
        if output and np.mod(m, 3e4) == 0:
            print('%.3f/%.3f hrs simulated' % (t/3600, time/3600))
    if output:
        print('Simulation Complete')
    return DiffProfile(dis*1e6, Xs, name=name)

##############################################################################

def main():
    """
    Generates 4 diffusion profiles from 4 separate steps with the pydiffusion
    library:
        
        Carburizing 1
        Diffusion 1
        Carburizing 2
        Diffusion 2
    
    Each step takes the resulting profile from the predecessing step as input. 
    In the two carburizing steps we use our own modified version of the 
    sphSim in the pydisffusion library. It leaves the concentration at the
    surface at a constant value (can be altered to a time-dependent function
    too), just like in real-life carburizing. 
    
    Carbon diffusion coefficients D(T) in austenite from:
    Thibaux, P., Métenier, A. & Xhoffer, C. Carbon Diffusion Measurement in 
    Austenite in the Temperature Range 500 °C to 900 °C. Metall Mater Trans 
    A 38, 1169–1176 (2007). https://doi.org/10.1007/s11661-007-9150-5

    Returns
    -------
    Saves the two diagrams as figures in your working directory. 
    Can also save .csvs of the profiles. Uncomment at your convenience.

    """
    # calculate diffusion coefficient from temperature with below expression
    # (see paper)
    temp_diffusion = 920  # °C
    D = 1.23E-6 * np.exp(-15050 / (temp_diffusion + 273.15))  # m²/s
    print(
        f"\nDiffusion coeff. of C in austenite at {temp_diffusion}°C is {D}m²/s\n")

    # Create diffusion system with constant DC
    # D can be concentration dependend D(c)
    diff_sys = DiffSystem(Xr=[0., 1], X=[0., 1], DC=[D, D], \
                          name=f'Constant D({temp_diffusion}°C)')

    # Create initial step profile,
    # it is somewhat useless when sphSim_source is used, but sphSim needs
    # a profile_init.
    # Xlim is the initial profile at distance
    # length_interval is the total depth
    # We should use mol/m³ as unit for concentration. But in this
    # case, we can use weight-% directly.
    # c_at_surface1_mol_m3 = 9728.45 # 1.5% in mol/m³
    # c_at_surface1_mol_m3 = 3266.15 # 0.5%
    # c_in_core_mol_m3 = 1178.5058 # 0.18%
    num_points = 501
    c_at_surface1 = 1.5 # C-concentration @ first surface 'layer' during carburizing step 1
    c_at_surface2 = 0.5 # C-concentration @ first surface 'layer' during carburizing step 2

    c_in_core = 0.18 # C-concentration of carburizing steel, 18CrNiMo7-6 --> 0.18
    distance = 0  # µm
    length_interval = 3000  # µm, interval to look at, diffusion depth should be
    # lower
    dist_mesh = mesh(0, length_interval, num_points)
    profile_init = step(dist_mesh, distance, diff_sys, Xlim=[
                        c_at_surface1, c_in_core], name='Intitial step profile')
    print("")
    
    # we leave the plot of the initial profile commented for now, because
    # we use sphSim_source for genertating an initial profile
    # fig = plt.figure(figsize=(16, 6))
    # ax1, ax2 = fig.add_subplot(121), fig.add_subplot(122)
    # ax1.set_title('Diffusion Coefficients', fontsize=15)
    # ax2.set_title('Initial Step Profile', fontsize=15)
    # DCplot(diff_sys, ax1)
    # profileplot(profile_init, ax2)
    # ax2.set_ylabel('C percentage / %') # we overwrite the default labels, because
    # ax2.set_xlabel('distance / µm') # they are not correct and/or not pretty.

    # Diffusion simulation using the setups with source at surface
    # Carburizing step, also called boost step.
    time_1 = 4 * 3600  # time in seconds
    profile_intermediate1 = sphSim_source(profile_init, c_at_surface1,
                                          diff_sys, time_1)
    print("Carburizing step 1 finished..... \n")
    # Diffusion simulation without source
    # Diffusion step
    time_2 = 4 * 3600  # time in seconds
    profile_intermediate2 = sphSim(profile_intermediate1, diff_sys, time_2)
    print("Diffusion step 1 finished..... \n")

    # Diffusion simulation using the setups with source at surface
    # Carburizing step 2
    time_3 = 2. * 3600  # time in seconds
    profile_intermediate3 = sphSim_source(profile_intermediate2, c_at_surface2,
                                          diff_sys, time_3)
    print("Carburizing step 2 finished..... \n")

    # Diffusion simulation without source
    # Diffusion step 2
    time_4 = 1.2 * 3600  # time in seconds
    profile_final = sphSim(profile_intermediate3, diff_sys, time_4)
    print("Diffusion step 2 finished..... \n")

    total_duration = time_1 + time_2 + time_3 + time_4
    print(f"Total process duration is {total_duration/3600:.1f} hours.")

    # Plotting the profiles
    fig = plt.figure(figsize=(12, 6))
    ax1, ax2 = fig.add_subplot(121), fig.add_subplot(122)
    ax1.set_title('Carburizing step 1', fontsize=12)
    ax2.set_title('Diffusion step 1', fontsize=12)
    profileplot(profile_intermediate1, ax1, c='g', ls='solid')
    profileplot(profile_intermediate2, ax2, c='r', ls='solid')
    ax1.set_ylabel('C percentage / %')
    ax2.set_ylabel('C percentage / %')
    ax1.set_xlabel('distance / µm')
    ax2.set_xlabel('distance / µm')
    fig.savefig('Carburizing+Diffusion_steps1.png')

    fig2 = plt.figure(figsize=(12, 6))
    ax3, ax4 = fig2.add_subplot(121), fig2.add_subplot(122)
    ax3.set_title('Carburizing step 2', fontsize=12)
    ax4.set_title('Diffusion step 2', fontsize=12)
    profileplot(profile_intermediate3, ax3, c='g', ls='solid')
    profileplot(profile_final, ax4, c='r', ls='solid')
    ax3.set_ylabel('C percentage / %')
    ax4.set_ylabel('C percentage / %')
    ax3.set_xlabel('distance / µm')
    ax4.set_xlabel('distance / µm')
    fig2.savefig('Carburizing+Diffusion_steps2.png')
    ###########################################################################

    # =========================================================================
    # # Save results
    # save_csv(f"Carburizing_step_1_{time_1/3600:.0f}hours.csv",
    #          profile=profile_intermediate1, diffsys=diffsys)
    # save_csv(f"Diffusion_step_1_{time_2/3600:.0f}hours.csv",
    #          profile=profile_intermediate2, diffsys=diffsys)
    # save_csv(f"Carburizing_step_2_{time_3/3600:.0f}hours.csv",
    #          profile=profile_intermediate3, diffsys=diffsys)
    # save_csv(f"Diffusion_step_2_{time_4/3600:.0f}hours.csv",
    #          profile=profile_final, diffsys=diffsys)
    # =========================================================================

if __name__ == "__main__":
    main()
    
