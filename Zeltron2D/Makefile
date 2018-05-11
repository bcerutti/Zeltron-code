#***********************************************************************!
#                       The Zeltron code project.                       !
#***********************************************************************!
# Copyright (C) 2012-2015. Authors: Benoît Cerutti & Greg Werner        !
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

# ---------------------Macro-Defs-----------------------------
## Common options for F90 are: ftn, h5pfc, mpif90
F90=mpif90

## for debugging, OPT=-g -fbacktrace -fbounds-check -O0
OPT=-g -O3

##-----------------I/O---------------
## IO should be "h5" for hdf5 output (hdf5 modules need to be loaded,
##   or include and library files available), or "txt" for text output
##   files, in which case no hdf5 capability need be present.
IO?=h5

ifeq ($(IO),h5)
	modio=mod_io_hdf5
	## The following compiler/linker options should be set for each machine,
	# except for supercomputers where they are set through module loading
	# and the correct compiler command (e.g., ftn or h5pfc).
	LPATH=-L /scr_verus/wernerg/contrib2/hdf5-1.8.11-par/lib
	LIB=-lhdf5_fortran -lhdf5_f90cstub -lhdf5
	INCLUDE=-I /scr_verus/wernerg/contrib2/hdf5-1.8.11-par/include/fortran
	## On supercomputers with hdf5 modules loaded 
	#  (or systems with h5pfc defined), uncomment the following 3 lines,
	#  and comment the above 3 lines.  
	#LPATH=
	#LIB=
	#INCLUDE=
else
	modio=mod_io_txt
	LPATH=
	LIB=
	INCLUDE=
endif

# -------------------End-macro-Defs---------------------------

all:main testExtcmd

# Here is the link step 
main:mod_input.o mod_extcmd.o mod_initial.o mod_mpi.o mod_log.o mod_interp.o $(modio).o mod_fields.o mod_motion.o mod_rhoj.o mod_sync.o mod_track.o mod_analysis.o main.o
	 $(F90) $(OPT) -o zeltron.exe $(LPATH) mod_input.o mod_extcmd.o mod_initial.o mod_mpi.o mod_log.o mod_interp.o $(modio).o mod_fields.o mod_motion.o mod_rhoj.o mod_sync.o mod_track.o mod_analysis.o main.o $(LIB) 

testExtcmd:extcmd.o mod_extcmd.o
	$(F90) $(OPT) mod_extcmd.o extcmd.o -o testExtcmd

# Here are the compile steps

mod_extcmd.o:mod_extcmd.f90
	$(F90) $(OPT) -c mod_extcmd.f90

extcmd.o:extcmd.f90
	$(F90) $(OPT) -c extcmd.f90
 
mod_input.o:mod_input.f90  
	 $(F90) $(OPT) -c mod_input.f90 

mod_initial.o:mod_initial.f90  
	 $(F90) $(OPT) -c mod_initial.f90 

mod_mpi.o:mod_mpi.f90  
	 $(F90) $(OPT) -c mod_mpi.f90 

mod_log.o:mod_log.f90  
	 $(F90) $(OPT) -c mod_log.f90 

mod_interp.o:mod_interp.f90  
	 $(F90) $(OPT) -c mod_interp.f90 

$(modio).o:$(modio).f90  
	 $(F90) $(OPT) $(INCLUDE) -c $(modio).f90 

mod_fields.o:mod_fields.f90  
	 $(F90) $(OPT) -c mod_fields.f90 

mod_motion.o:mod_motion.f90  
	 $(F90) $(OPT) -c mod_motion.f90

mod_rhoj.o:mod_rhoj.f90  
	 $(F90) $(OPT) -c mod_rhoj.f90
 
mod_sync.o:mod_sync.f90  
	 $(F90) $(OPT) -c mod_sync.f90

mod_track.o:mod_track.f90  
	 $(F90) $(OPT) -c mod_track.f90 

mod_analysis.o:mod_analysis.f90  
	 $(F90) $(OPT) -c mod_analysis.f90
 
main.o:main.f90 mod_input.o mod_initial.o mod_mpi.o mod_log.o mod_interp.o $(modio).o mod_fields.o mod_motion.o mod_rhoj.o mod_sync.o mod_track.o mod_analysis.o 
	 $(F90) $(OPT) -c main.f90

