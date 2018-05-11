!***********************************************************************!
!                       The Zeltron code project.                       !
!***********************************************************************!
! Copyright (C) 2012-2015. Authors: Benoît Cerutti & Greg Werner        !
!                                                                       !
! This program is free software: you can redistribute it and/or modify  !
! it under the terms of the GNU General Public License as published by  !
! the Free Software Foundation, either version 3 of the License, or     !
! (at your option) any later version.                                   !
!                                                                       !
! This program is distributed in the hope that it will be useful,       !
! but WITHOUT ANY WARRANTY; without even the implied warranty of        !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
! GNU General Public License for more details.                          !
!                                                                       !
! You should have received a copy of the GNU General Public License     !
! along with this program. If not, see <http://www.gnu.org/licenses/>.  !
!***********************************************************************!

MODULE MOD_INPUT

IMPLICIT NONE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ CONSTANTS +++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Speed of light [cm/s]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: c=299792458d2
! Fundamental charge [esu]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: e=4.8032068d-10
! Boltzmann constant [erg/K]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: k=1.380658d-16
! Planck constant [erg.s]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: h=6.6260755d-27
! Mass of the electron [g]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: me=9.1093897d-28
! Mass of the proton [g]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: mp=1.6726231d-24
! 1 eV in erg [erg]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: evtoerg=1.602177d-12
! pi
DOUBLE PRECISION, PARAMETER, PUBLIC           :: pi=dacos(-1d0)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ DIMENSION ++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! spatial dimension
INTEGER, PARAMETER, PUBLIC :: NDIM=3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ SAVE/RESTORE ++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Save particle and field data for future restoration of the simulation
LOGICAL, PARAMETER, PUBLIC :: CHECKPOINT=.FALSE.

! Restore a simulation where it stopped
LOGICAL, PARAMETER, PUBLIC :: RESTORE=.FALSE.

! Saving frequency, in seconds (the simulation will try to save its 
! current state *before* elapsed time n*FSAVE, for integer n).
! N.B. It is highly recommended to have more than 1 checkpoint during
! the simulation, unless it is certain that the simulation will end
! naturally (as opposed to being killed before it is finished).
DOUBLE PRECISION, PARAMETER, PUBLIC :: FSAVE=0.

! Give the time step from which the simulation should restart
INTEGER, PARAMETER, PUBLIC :: time_ref=0

! Whether to test restore capability: this restores data from step 
!   time_ref, and then immediately write data (as if dumping
!   for a checkpoint) under the next time step; the simulation then halts.
LOGICAL, PARAMETER, PUBLIC :: TEST_RESTORE = .FALSE.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ BOUNDARY CONDITIONS ++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Specify the boundary conditions for the fields:
! 1. "PERIODIC": Periodic boundary conditions
! 2. "METAL": Perfect metal with infinite conductivity

CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_XMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_XMAX="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_YMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_YMAX="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_ZMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_ZMAX="PERIODIC"

! Specify the boundary conditions for the particles:
! 1. "PERIODIC": Periodic boundary conditions
! 2. "REFLECT": Particles are elastically reflected at the wall
! 3. "ABSORB": Particles are absorbed at the wall

CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_XMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_XMAX="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_YMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_YMAX="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_ZMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_ZMAX="PERIODIC"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ INPUT PARAMETERS ++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Number of cells in X
INTEGER*8, PARAMETER, PUBLIC :: NCX=32

! Number of cells in Y
INTEGER*8, PARAMETER, PUBLIC :: NCY=32

! Number of cells in Z
INTEGER*8, PARAMETER, PUBLIC :: NCZ=16

! Number of particles per cell per species
INTEGER*8, PARAMETER, PUBLIC :: PPC=1

! Number of process (domain decomposition in the X- Y- and Z-directions)
INTEGER, PARAMETER, PUBLIC :: NPX=1
INTEGER, PARAMETER, PUBLIC :: NPY=1
INTEGER, PARAMETER, PUBLIC :: NPZ=1

! Mass ratio IONS/ELECTRONS
DOUBLE PRECISION, PARAMETER, PUBLIC :: mass_ratio=1d0

! Spatial boundaries in the X-direction
DOUBLE PRECISION, PARAMETER, PUBLIC :: xmin=0d0,xmax=32d0

! Spatial boundaries in the Y-direction
DOUBLE PRECISION, PARAMETER, PUBLIC :: ymin=0d0,ymax=32d0

! Spatial boundaries in the Z-direction
DOUBLE PRECISION, PARAMETER, PUBLIC :: zmin=0d0,zmax=16d0

! Dump data frequency in terms of timesteps
INTEGER, PARAMETER, PUBLIC :: FDUMP=200

! Number of data dumps
INTEGER, PARAMETER, PUBLIC :: NDUMP=5

! Number of time steps
INTEGER, PARAMETER, PUBLIC :: NT=FDUMP*NDUMP

! Poisson solver calling frequency in terms of timesteps
INTEGER, PARAMETER, PUBLIC :: FREQ_POISSON=25

! Number of iterations to solve Poisson's equation
INTEGER, PARAMETER, PUBLIC :: NIT=500

! Switch ON/OFF the radiation reaction force (synchrotron radiation and
! inverse Compton scattering in the Thomson regime only)
LOGICAL, PARAMETER, PUBLIC :: RAD_FORCE=.FALSE.

! Switch ON/OFF spatial filtering the electric field
LOGICAL, PARAMETER, PUBLIC :: FILTER=.FALSE.

! Smoothing parameter
DOUBLE PRECISION, PARAMETER, PUBLIC :: alpha=1d0/64d0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Number of grid cells in X
INTEGER, PARAMETER, PUBLIC :: NX=NCX+1

! Number of grid cells in Y
INTEGER, PARAMETER, PUBLIC :: NY=NCY+1

! Number of grid cells in Z
INTEGER, PARAMETER, PUBLIC :: NZ=NCZ+1

! Number of cell per process in the X-direction
INTEGER, PARAMETER, PUBLIC :: NCXP=NCX/NPX

! Number of cell per process in the Y-direction
INTEGER, PARAMETER, PUBLIC :: NCYP=NCY/NPY

! Number of cell per process in the Z-direction
INTEGER, PARAMETER, PUBLIC :: NCZP=NCZ/NPZ

! Number of grid cells per process in X
INTEGER, PARAMETER, PUBLIC :: NXP=NCXP+1

! Number of grid cells per process in Y
INTEGER, PARAMETER, PUBLIC :: NYP=NCYP+1

! Number of grid cells per process in Z
INTEGER, PARAMETER, PUBLIC :: NZP=NCZP+1

! Total number of particles per species NP=NCX*NCY*NCZ*PPC
INTEGER*8, PARAMETER, PUBLIC :: NP=NCX*NCY*NCZ*PPC

! Spatial step
DOUBLE PRECISION, PARAMETER, PUBLIC :: dx=(xmax-xmin)/NCX
DOUBLE PRECISION, PARAMETER, PUBLIC :: dy=(ymax-ymin)/NCY
DOUBLE PRECISION, PARAMETER, PUBLIC :: dz=(zmax-zmin)/NCZ

! Time step (Courant-Friedrichs-Lewy timestep)
DOUBLE PRECISION, PARAMETER, PUBLIC :: dt=0.99*1d0/sqrt(1d0/dx**2d0+&
                                          1d0/dy**2d0+1d0/dz**2d0)/c

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Initial (tearing-mode) perturbation amplitude
DOUBLE PRECISION, PARAMETER, PUBLIC :: perturb_amp=0.0

! Ratio of Bz/B0 (Guide field / reconnecting field)
DOUBLE PRECISION, PARAMETER, PUBLIC :: guide_field=0.0
                                          
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Temperature in unit of me*c^2/k of the DRIFTING ELECTRONS in the co-moving frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: thde=1d0

! Temperature in unit of mi*c^2/k of the DRIFTING IONS in the co-moving frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: thdi=1d0

! Temperature in unit of me*c^2/k of the BACKGROUND ELECTRONS in the lab frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: thbe=1d0

! Temperature in unit of mi*c^2/k of the BACKGROUND IONS in the lab frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: thbi=1d0

! Minimum Larmor radius of the electrons
DOUBLE PRECISION, PARAMETER, PUBLIC :: rhoc=2d0

! Drift velocity betad=vdrift/c
DOUBLE PRECISION, PARAMETER, PUBLIC :: betad=0.6d0

! Density ratio between initial DRIFTING and BACKGROUND particles
DOUBLE PRECISION, PARAMETER, PUBLIC :: density_ratio=0.1d0

! Energy density ratio between external radiation field and the magnetic field
! udens_ratio=Uph/Ub, where Ub=B0^2/8*pi
DOUBLE PRECISION, PARAMETER, PUBLIC :: udens_ratio=0d0

! 4-velocity boundaries for the particles
DOUBLE PRECISION, PARAMETER, PUBLIC :: umin=1d-2,umax=1d2

! 4-velocity boundaries for the drifting electrons
DOUBLE PRECISION, PARAMETER, PUBLIC :: udemin=-1d2,udemax=1d2
! 4-velocity boundaries for the drifting ions
DOUBLE PRECISION, PARAMETER, PUBLIC :: udpmin=-1d2,udpmax=1d2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ DEFINE PARTICLE ARRAYS +++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! pcl_species(1,:)=x      : x position of the particle
! pcl_species(2,:)=y      : y position of the particle
! pcl_species(3,:)=z      : z position of the particle
! pcl_species(4,:)=ux     : x-component of the reduced particle 4-velocity
! pcl_species(5,:)=uy     : y-component of the reduced particle 4-velocity
! pcl_species(6,:)=uz     : z-component of the reduced particle 4-velocity
! pcl_species(7,:)=weight : weight of the particle
!
! Defition of the particle distribution data function: it gives extra
! information regarding the particles, e.g., electric/magnetic field/force
! at the location of the particle.
!
! pcl_data_species(1,:)=Ell : Electric field parallel to the particle velocity
! pcl_data_species(2,:)=Bpp : Magnetic field perp. to the particle velocity
! pcl_data_species(3,:)=Fr  : Synchrotron radiation reaction force
! pcl_data_species(4,:)=Fe  : Electric force
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! DRIFTING ELECTRONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_ed(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_ed(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: taged(:)

! DRIFTING IONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_pd(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_pd(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagpd(:)

! BACKGROUND ELECTRONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_eb(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_eb(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tageb(:)

! BACKGROUND IONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_pb(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_pb(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagpb(:)

! Temporary particle distribution function
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_f(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_f(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagf(:)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ ANALYSIS ++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Data formatting (for data in text files) for integers and doubles
CHARACTER(LEN=6), PARAMETER, PUBLIC  :: FMT_INT="I26"
CHARACTER(LEN=11), PARAMETER, PUBLIC :: FMT_DBL="E26.16E3"

! For text output of fields/density/pressure/current: (irrelevant for hdf5 output)
! If .TRUE., each MPI process will write to a separate file.
! This can be faster (for large simulations), but may create too many files.
LOGICAL, PARAMETER, PUBLIC :: DUMP_FIELDS_PARALLEL=.FALSE.

! Switch ON/OFF write fields to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_FIELDS=.TRUE.

! Switch ON/OFF write the charge and current densities to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_RHOJ=.TRUE.

! For text output of particle data: (irrelevant for hdf5 output)
! If .TRUE., each MPI process will write to its own file.
! This can be faster (for large simulations), but may create too many files.
LOGICAL, PARAMETER, PUBLIC :: DUMP_PARTICLES_PARALLEL=.FALSE.

! Switch ON/OFF write raw particle data to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_PARTICLES=.FALSE.

! Switch ON/OFF particle tracker
LOGICAL, PARAMETER, PUBLIC :: TRACK_PARTICLES=.FALSE.

! Switch ON/OFF analysis of particle energy and angular distributions
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_DISTRIBUTIONS=.TRUE.

! Switch ON/OFF analysis the plasma density
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_DENSITIES=.TRUE.

! Switch ON/OFF analysis the relativistic macroscopic fluid quantities
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_FLUID=.FALSE.

! Switch ON/OFF analysis the radiation spectra and anisotropies
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_RADIATION=.FALSE.

! Total number of tracked particles
INTEGER, PARAMETER, PUBLIC :: NSAMPLE=100

! Give the time step from which the tracker should restart
INTEGER, PARAMETER, PUBLIC :: time_track=1

! The tracking frequency (timesteps between recording of tracking data)
! E.g., 1 = record every step; 2 = record every other step; ....
INTEGER, PARAMETER, PUBLIC :: FTRACK=1

! Particle's Lorentz factor grid for the spectrum
INTEGER, PARAMETER, PUBLIC :: NU=100

! Spatial grid size for dumping EM-related fields (E,B,J)
! Spatial grid size for dumping EM-related fields
! If NFIELDX = NCX, etc., Zeltron dumps exactly the (nodal) field it uses.
INTEGER, PARAMETER, PUBLIC :: NFIELDX=NCX,NFIELDY=NCY,NFIELDZ=NCZ

! Spatial grid size for dumping particle-related fields (densities, etc.)
INTEGER, PARAMETER, PUBLIC :: NDX=NCX,NDY=NCY,NDZ=NCZ

! Angular distribution grid
INTEGER, PARAMETER, PUBLIC :: NPHI=50,NLBA=90

! Synchrotron radiation frequency
INTEGER, PARAMETER, PUBLIC :: NNU=100

! Frequency boundaries for the synchrotron radiation
DOUBLE PRECISION, PARAMETER, PUBLIC :: nmin=1d9,nmax=1d14

! Latitude boundaries for the angular distribution (phi)
DOUBLE PRECISION, PARAMETER, PUBLIC :: pmin=-89.99d0,pmax=89.99d0

! Longitude boundaries for the angular distribution (lambda)
DOUBLE PRECISION, PARAMETER, PUBLIC :: lmin=-179.99d0,lmax=179.99d0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ TUNING ++++++++++++++++++++++++++++++++++++
! You might want to adjust these parameters, but they should mostly  
! be left as they are. 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! (Inverse) frequency, in steps, with which all ranks communicate
! general info, such as whether to checkpoint, or the contents of the
! command file.
INTEGER, PARAMETER, PUBLIC :: FGENCOMM=20

! Fraction of FSAVE by which a checkpoint time can be shifted (earlier)
! so that it occurs as the same time as a dump.
DOUBLE PRECISION, PARAMETER, PUBLIC :: FSAVE_SHIFT_FRAC=0.25

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ DEBUGGING +++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! For normal runs, RANDOMIZE = true ensures a different seed each run
! For testing new code, can set to false for reproducibility
LOGICAL, PARAMETER, PUBLIC :: RANDOMIZE=.TRUE.

! whether to write one diagnostic file per rank
LOGICAL, PARAMETER, PUBLIC :: writePerRankFiles=.FALSE.

! the name of the per-rank file
CHARACTER(len=18), PUBLIC  :: perRankFile

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE MOD_INPUT
