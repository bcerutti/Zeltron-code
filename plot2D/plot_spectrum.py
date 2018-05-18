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
# This python script will draw the total electron or ion energy spectrum:
# dN/du(u), where u is the 4 velocity, at a given time step it.
#
# To execute, type for instance:$ python plot_spectrum.py 0 ions &
# This command will draw the ions spectrum at time step 0.
#
#***********************************************************************!

import matplotlib
import numpy
import math
import matplotlib.pyplot as plt
import sys

def plot_spectrum(it,spec):

    if it=='':
       it='0'

    if spec=='':
       spec='electrons'

    #===============================================================================
    # 4-velocity
    u=numpy.loadtxt(".././Zeltron2D/data/u.dat")
    u=u[0:len(u)-1]
    
    # Spectrum
    dNdud=numpy.loadtxt(".././Zeltron2D/data/spectrum_"+spec+"_drift"+it+".dat")
    dNdub=numpy.loadtxt(".././Zeltron2D/data/spectrum_"+spec+"_bg"+it+".dat")
    
    # Total spectrum
    dNdu=dNdud+dNdub
    
    plt.plot(u,dNdu,color='blue',lw=2)
    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'$\gamma\beta$',fontsize=20)
    plt.ylabel(r'$\frac{dN}{d(\gamma\beta)}$',fontsize=20)
    plt.ylim([1e-3*numpy.max(dNdu),3.0*numpy.max(dNdu)])
    
    plt.title("Time step="+it+", Species="+spec, fontsize=18)
    
    #===============================================================================

    plt.show()
    
    #===============================================================================

plot_spectrum(sys.argv[1],sys.argv[2])
