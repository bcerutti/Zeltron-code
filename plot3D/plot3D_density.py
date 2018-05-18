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
# This python script will draw the 3D electron or ion density isosurfaces
# as function of x, y, z, at a given time step it.
#
# To execute, type for instance:$ python plot3D_density.py 0 ions &
# This command will draw the 3D ions density at time step 0.
#
# Notice: You need the Mayavi library to run this script
#***********************************************************************!

import numpy
import math
from mayavi.mlab import *
from mayavi import mlab
import sys

def plot3D_density(it,spec):

    if it=='':
       it='0'

    if spec=='':
       spec='electrons'

    #===============================================================================
    # Parameters of the simulation

    params1=numpy.loadtxt(".././Zeltron3D/data/phys_params.dat",skiprows=1)
    params2=numpy.loadtxt(".././Zeltron3D/data/input_params.dat",skiprows=1)

    # Nominal cyclotron frequency
    rho=2.9979246e+10/params1[3]

    #===============================================================================
    # The grid
    x=numpy.loadtxt(".././Zeltron3D/data/x.dat")
    y=numpy.loadtxt(".././Zeltron3D/data/y.dat")
    z=numpy.loadtxt(".././Zeltron3D/data/z.dat")
   
    # The density
    map0=numpy.loadtxt(".././Zeltron3D/data/densities/mapxyz_"+spec+"_drift0.dat")
    mapd=numpy.loadtxt(".././Zeltron3D/data/densities/mapxyz_"+spec+"_drift"+it+".dat")
    mapb=numpy.loadtxt(".././Zeltron3D/data/densities/mapxyz_"+spec+"_bg"+it+".dat")

    nx=len(x)-1
    ny=len(y)-1
    nz=len(z)-1

    dx,ptx=x[1]-x[0],nx*1j
    dy,pty=y[1]-y[0],ny*1j
    dz,ptz=z[1]-z[0],nz*1j

    x,y,z=numpy.ogrid[-dx:dx:ptx,-dy:dy:pty,-dz:dz:ptz]

    mapxyz=numpy.empty((nx,ny,nz))

    for ix in range(0,nx):
        for iy in range(0,ny):
            for iz in range(0,nz):   
                mapxyz[ix,iy,iz]=mapd[iy+iz*ny,ix]+mapb[iy+iz*ny,ix]
      
    mapxyz=mapxyz/numpy.max(map0)
    
    contour3d(mapxyz,contours=10,transparent=True,opacity=0.9)
    mlab.outline()
    mlab.colorbar(title='n/n0',orientation='vertical',nb_labels=10)
    mlab.axes()
    mlab.title("Time="+it+", "+spec)
    
    #===============================================================================

    mlab.show()
    
    #===============================================================================

plot3D_density(sys.argv[1],sys.argv[2])
