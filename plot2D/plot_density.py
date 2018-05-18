#***********************************************************************!
#                       The Zeltron code project.                       !
#***********************************************************************!
# Copyright (C) 2012-2015. Authors: Benoit Cerutti & Greg Werner        !
#                                                                       !
# This program is free software: you can redistribute it and/or modify  !
# it under the terms of the GNU General Public License as published by  !
# the Free Software Foundation, either version 3 of the License, or     !
# (at your option) any later version.                                   !
#                                                                       !
# This program is distributed in the hope that it will be useful,       !
# but WITHOUT ANY WARRANTY; without even the implied warranty of        !
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
# GNU General Public License for more details.                          !
#                                                                       !
# You should have received a copy of the GNU General Public License     !
# along with this program. If not, see <http://www.gnu.org/licenses/>.  !
#***********************************************************************!
#
# This python script will draw the electron or ion density and the 
# in-plane magnetic field lines as function of x and y, at a given time 
# step it.
#
# To execute, type for instance:$ python plot_density.py 0 ions &
# This command will draw the ions density at time step 0.
#
#***********************************************************************!

import numpy
import math
import matplotlib.pyplot as plt
import sys

def plot_density(it,spec):

    if it=='':
       it='0'

    if spec=='':
       spec='electrons'

    #===============================================================================
    # Parameters of the simulation

    params1=numpy.loadtxt(".././Zeltron2D/data/phys_params.dat",skiprows=1)
    params2=numpy.loadtxt(".././Zeltron2D/data/input_params.dat",skiprows=1)

    # Nominal Larmor radius
    rho=2.9979246e+10/params1[3]

    #===============================================================================
    # The grid
    x=numpy.loadtxt(".././Zeltron2D/data/x.dat")
    y=numpy.loadtxt(".././Zeltron2D/data/y.dat")

    # The density
    mapxy0=numpy.loadtxt(".././Zeltron2D/data/densities/mapxy_"+spec+"_drift0.dat")
    map_bg=numpy.loadtxt(".././Zeltron2D/data/densities/mapxy_"+spec+"_bg"+it+".dat")
    map_drift=numpy.loadtxt(".././Zeltron2D/data/densities/mapxy_"+spec+"_drift"+it+".dat")
    
    mapxy=(map_bg+map_drift)/numpy.max(mapxy0)

    plt.pcolormesh(x/rho,y/rho,mapxy)
    plt.xlabel(r'$x/\rho$',fontsize=18)
    plt.ylabel(r'$y/\rho$',fontsize=18)
    
    plt.colorbar().ax.set_ylabel(r'$n/n_0$',fontsize=18)
    
    plt.title("Time step="+it+", Species="+spec,fontsize=18)
    
    #===============================================================================
    # Plotting magnetic field lines
    
    xf=numpy.loadtxt(".././Zeltron2D/data/xfield.dat")
    yf=numpy.loadtxt(".././Zeltron2D/data/yfield.dat")
    
    # The fields
    Bx=numpy.loadtxt(".././Zeltron2D/data/fields/Bx"+it+".dat")
    By=numpy.loadtxt(".././Zeltron2D/data/fields/By"+it+".dat")
    
    # The magnetic flux function
    psi=numpy.empty((len(yf),len(xf)))

    dx=xf[1]-xf[0]
    dy=yf[1]-yf[0]
    
    for iy in range(0,len(yf)-1):
        psi[iy+1,0]=psi[iy,0]+Bx[iy,0]*dy
        
        for ix in range(0,len(xf)-1):
            psi[iy,ix+1]=psi[iy,ix]-By[iy,ix]*dx
            
    # Periodic boundary conditions
    psi[:,len(xf)-1]=psi[:,0]
    psi[len(yf)-1,:]=psi[0,:]
    
    # To set the negative contour as solid lines (dashed is default)
    plt.rcParams['contour.negative_linestyle']='solid'
    plt.contour(xf/rho,yf/rho,psi,25,lw=1.0,colors='white')
    
    #===============================================================================

    plt.show()
    
    #===============================================================================

plot_density(sys.argv[1],sys.argv[2])
