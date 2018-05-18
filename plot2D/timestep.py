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
# This python script will draw the duration of each time step as a
# function of time.
#
# To execute, type for instance:$ python timestep.py &
#
#***********************************************************************!

import numpy
import math
import matplotlib.pyplot as plt

def timestep():

    #===============================================================================
    time=numpy.loadtxt(".././Zeltron2D/data/timestep.dat")
    
    plt.plot(time,color='blue',lw=2)
    plt.yscale('log')

    plt.xlabel('Time step',fontsize=18)
    plt.ylabel('Duration time step [s]',fontsize=18)
    
    average=numpy.mean(time)

    plt.plot([0,len(time)],[average,average],lw=2,ls='--',color='black')
    
    print ""
    print "********************************"
    print "Average: ", average," seconds"
    print "********************************"

    plt.title("Average time step = %2.3f " % average + "seconds", fontsize=20)

    #===============================================================================

    plt.show()
    
    #===============================================================================
    
timestep()
