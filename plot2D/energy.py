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
# This python script will draw the total magnetic, electric, particle, 
# and radiative energies as a function of time.
#
# To execute, type for instance:$ python energy.py &
#
#***********************************************************************!

import numpy
import math
import matplotlib.pyplot as plt
import pylab

def energy():

    #===============================================================================
    # Magnetic and electric energies

    data=numpy.loadtxt(".././Zeltron2D/data/Eem.dat")
    emag=data[:,0]
    eelc=data[:,1]

    plt.plot(emag,color='blue',lw=2)
    plt.plot(eelc,color='red',lw=2)

    plt.xlabel('Time step',fontsize=18)
    plt.ylabel('Energy',fontsize=18)

    #===============================================================================
    # Particle kinetic energies

    # Electrons
    ekineb=numpy.loadtxt(".././Zeltron2D/data/Ekin_electrons_bg.dat")
    ekined=numpy.loadtxt(".././Zeltron2D/data/Ekin_electrons_drift.dat")

    # Ions
    ekinpb=numpy.loadtxt(".././Zeltron2D/data/Ekin_ions_bg.dat")
    ekinpd=numpy.loadtxt(".././Zeltron2D/data/Ekin_ions_drift.dat")

    plt.plot(ekineb+ekined+ekinpb+ekinpd,ls='--',lw=2,color='green')

    #===============================================================================
    # Radiative energies

    # Synchrotron
    # Electrons
    esyne=numpy.loadtxt(".././Zeltron2D/data/Esyn_electrons.dat")
    # Ions
    esynp=numpy.loadtxt(".././Zeltron2D/data/Esyn_electrons.dat")

    plt.plot(esyne+esynp,ls='--',lw=2,color='magenta')

    # Inverse Compton
    # Electrons
    eicse=numpy.loadtxt(".././Zeltron2D/data/Eics_electrons.dat")
    # Ions
    eicsp=numpy.loadtxt(".././Zeltron2D/data/Eics_electrons.dat")

    plt.plot(eicse+eicsp,ls=':',lw=2,color='magenta')

    #===============================================================================
    # Total energy

    etot=emag+eelc+ekineb+ekined+ekinpb+ekinpd+esyne+esynp+eicse+eicsp
    plt.plot(etot,ls='-',lw=3,color='black')

    error=(etot[len(etot)-1]-etot[0])/etot[0]*100.0

    # Relative error
    print ""
    print "********************************"
    print "Total energy relative error:"
    print error,"%"
    print "********************************"

    # Plot legend
    pylab.legend(('Magnetic','Electric','Particles','Synchrotron','Inv. Compton','Total'),shadow=False,loc=(0.65,0.27))
    ltext = pylab.gca().get_legend().get_texts()
    pylab.setp(ltext[0],fontsize=15)
    pylab.setp(ltext[1],fontsize=15)
    pylab.setp(ltext[2],fontsize=15)
    pylab.setp(ltext[3],fontsize=15)
    pylab.setp(ltext[4],fontsize=15)
    pylab.setp(ltext[5],fontsize=15)
    
    plt.title("Error= %+2.3f " % error + "%", fontsize=20)

    #===============================================================================

    plt.show()
    
    #===============================================================================
    
energy()
