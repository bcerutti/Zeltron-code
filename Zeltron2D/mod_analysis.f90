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

MODULE MOD_ANALYSIS

USE MOD_INPUT
USE MOD_MPI
USE MOD_LOG
USE MOD_INTERP
USE MOD_IO
USE MOD_SYNC
USE MOD_INITIAL
USE MOD_TRACK

IMPLICIT NONE

! =======================================================
! keep track of coordinates that need be dumped only once
! u (gamma-velocity)
LOGICAL :: uDumped = .FALSE.
! phi (latitude for angular velocity)
LOGICAL :: phiDumped = .FALSE.
! lambda (longitude for angular velocity)
LOGICAL :: lambdaDumped = .FALSE.
! nu (photon spectrum)
LOGICAL :: nuDumped = .FALSE.
! =======================================================

PRIVATE

PUBLIC :: START_LOGS ! set up log files
PUBLIC :: WRITE_LOGS ! force log data to disk
PUBLIC :: END_LOGS   ! write final log data do disk
PUBLIC :: SPECTRUM_ANGULAR ! Particles' spectrum and angular distribution
PUBLIC :: MAP_XY ! Particles' spatial distribution
PUBLIC :: MAP_FLUID ! Computes the macroscopic relativistic fluid quantities
PUBLIC :: DUMP_COORDS ! dump coordinate data to disk (once per simulation)
PUBLIC :: DUMP_FIELD_NOUGC ! DUMP_FIELD, without global upper guard cells
PUBLIC :: REDUCE_GRID ! Compute a reduced grid, to save space
PUBLIC :: REDUCE_GRID_1D ! Compute a reduced grid in 1D
PUBLIC :: REDUCE_FIELDS ! Sample fields on smaller (actually arbitrary) grid
PUBLIC :: DUMP_REDUCED_FIELDS ! Dump reduced E, B
PUBLIC :: DUMP_REDUCED_RHOJ ! Writes the charge and current arrays to disk
PUBLIC :: DUMP_GRID ! Dump grid points
PUBLIC :: EM_ENERGY ! Electric and magnetic energies
PUBLIC :: KIN_ENERGY ! Particles' kinetic energy
PUBLIC :: RAD_ENERGY ! Total radiative energy losses
PUBLIC :: ANALYSIS_SYNC ! Synchrotron radiation spectrum and angular distrib.

TYPE(DataLogObj) :: EemLog, Eics_eLog, Eics_pLog, Esyn_eLog, Esyn_pLog
TYPE(DataLogObj) :: Ekin_ebLog, Ekin_edLog, Ekin_pbLog, Ekin_pdLog

 CONTAINS

!***********************************************************************
! Subroutine START_LOGS
! Create logs: write after 100 entries or 120 seconds, whichever comes first
! This must be called by mpi rank 0; may be called by other ranks as well.
! INPUT:
! - id: the mpi rank
!
!***********************************************************************
SUBROUTINE START_LOGS(id)

IMPLICIT NONE

! mpi rank
INTEGER :: id

INTEGER          :: maxEntries = 100
DOUBLE PRECISION :: maxTime = 120.
!***********************************************************************

IF (id == 0) THEN
CALL START_LOG(EemLog, "./data/Eem.dat", 2, maxEntries, maxTime)
CALL START_LOG(Eics_eLog, "./data/Eics_electrons.dat", 1, maxEntries, maxTime)
CALL START_LOG(Eics_pLog, "./data/Eics_ions.dat", 1, maxEntries, maxTime)
CALL START_LOG(Esyn_eLog, "./data/Esyn_electrons.dat", 1, maxEntries, maxTime)
CALL START_LOG(Esyn_pLog, "./data/Esyn_ions.dat", 1, maxEntries, maxTime)
CALL START_LOG(Ekin_ebLog, "./data/Ekin_electrons_bg.dat", 1, maxEntries, maxTime)
CALL START_LOG(Ekin_edLog, "./data/Ekin_electrons_drift.dat", 1, maxEntries, maxTime)
CALL START_LOG(Ekin_pbLog, "./data/Ekin_ions_bg.dat", 1, maxEntries, maxTime)
CALL START_LOG(Ekin_pdLog, "./data/Ekin_ions_drift.dat", 1, maxEntries, maxTime)
ENDIF

END SUBROUTINE START_LOGS

!***********************************************************************
! Subroutine WRITE_LOGS
! Force all log data to be written to disk
! This must be called by mpi rank 0; may be called by other ranks as well.
! INPUT:
! - id: the mpi rank
!
!***********************************************************************
SUBROUTINE WRITE_LOGS(id)

IMPLICIT NONE

! mpi rank
INTEGER :: id
!***********************************************************************

IF (id == 0) THEN
CALL WRITE_LOG_TO_DISK(EemLog)
CALL WRITE_LOG_TO_DISK(Eics_eLog)
CALL WRITE_LOG_TO_DISK(Eics_pLog)
CALL WRITE_LOG_TO_DISK(Esyn_eLog)
CALL WRITE_LOG_TO_DISK(Esyn_pLog)
CALL WRITE_LOG_TO_DISK(Ekin_ebLog)
CALL WRITE_LOG_TO_DISK(Ekin_edLog)
CALL WRITE_LOG_TO_DISK(Ekin_pbLog)
CALL WRITE_LOG_TO_DISK(Ekin_pdLog)
ENDIF

END SUBROUTINE WRITE_LOGS

!***********************************************************************
! Subroutine END_LOGS
! Finish writing logs, forever
! This must be called by mpi rank 0; may be called by other ranks as well.
! INPUT:
! - id: the mpi rank
!
!***********************************************************************
SUBROUTINE END_LOGS(id)

IMPLICIT NONE

! mpi rank
INTEGER :: id

IF (id == 0) THEN
CALL END_LOG(EemLog)
CALL END_LOG(Eics_eLog)
CALL END_LOG(Eics_pLog)
CALL END_LOG(Esyn_eLog)
CALL END_LOG(Esyn_pLog)
CALL END_LOG(Ekin_ebLog)
CALL END_LOG(Ekin_edLog)
CALL END_LOG(Ekin_pbLog)
CALL END_LOG(Ekin_pdLog)
ENDIF

END SUBROUTINE END_LOGS

!***********************************************************************
! Subroutine SPECTRUM_ANGULAR
! This subroutine computes the particles' spectrum and angular distribution.
!
! INPUT:
! 
! - pcl: Particle distribution function
! - it:  Timestep number
! - spec: Particle species
! - sym: Name of the particles
!
! OUTPUT: u,phi,lambda,dN/dOmega/du,dNdu
!***********************************************************************

SUBROUTINE SPECTRUM_ANGULAR(pcl,NPP,it,spec,sym,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                               :: id,COMM,ierr
INTEGER*8                                             :: ip,NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)                :: pcl
DOUBLE PRECISION                                      :: xp,yp,zp,ux,uy,uz,wp
DOUBLE PRECISION, DIMENSION(1:NU)                     :: uref
DOUBLE PRECISION, DIMENSION(1:NU-1)                   :: du
DOUBLE PRECISION, DIMENSION(1:NU-1)                   :: spectrap,spectra
DOUBLE PRECISION, DIMENSION(1:NPHI)                   :: pref
DOUBLE PRECISION, DIMENSION(1:NPHI-1)                 :: dSinPhi
DOUBLE PRECISION, DIMENSION(1:NLBA)                   :: lref
DOUBLE PRECISION, DIMENSION(1:NLBA-1)                 :: dLrad
DOUBLE PRECISION, DIMENSION(1:NU-1,1:NPHI-1,1:NLBA-1) :: angular,angularp
DOUBLE PRECISION                                      :: yhalf
DOUBLE PRECISION                                      :: utot,clba,slba,lambda
DOUBLE PRECISION                                      :: cphi,sphi,phi
INTEGER                                               :: iu,iph,il
INTEGER                                               :: iu2,iph2,il2,it
CHARACTER(len=10)                                     :: cit,spec,sym
!***********************************************************************

! Convert the integer it into a string cit
 write(cit,'(i10)') it

! This step left-justify the string
 cit=adjustl(cit)
 spec=adjustl(spec)
 sym=adjustl(sym)
 
! Define the y separating the 2 layers
yhalf=(ymax-ymin)/2.0+ymin

DO iu=1,NU
uref(iu)=1d1**((iu-1)*1d0/(NU-1d0)*(log10(umax)-&
         log10(umin))+log10(umin))
ENDDO
du = uref(2:NU) - uref(1:NU-1)

! Latitude PHI in degrees
DO iph=1,NPHI
pref(iph)=(iph-1)*1d0/(NPHI-1d0)*(pmax-pmin)+pmin
ENDDO

!DO iph=1,NPHI
!dSinPhi(iph)= sin(pref(iph+1)*pi/180.0)-sin(pref(iph)*pi/180.0)
!ENDDO

dSinPhi = sin(pref(2:NPHI)*pi/180.0) - sin(pref(1:NPHI-1)*pi/180.0)

! Longitude LAMBDA in degrees
DO il=1,NLBA
lref(il)=(il-1)*1d0/(NLBA-1d0)*(lmax-lmin)+lmin
ENDDO
dLrad = pi/180. * (lref(2:NLBA) - lref(1:NLBA-1))

! Angular distribution per process
angularp=0.0

! Angular distribution for all process
angular=0.0

! Particle energy spectrum per process
spectrap=0.0

! Particle energy spectrum for all process
spectra=0.0

!===================================================================
DO ip=1,NPP

  xp=pcl(1,ip)
  yp=pcl(2,ip)
  zp=pcl(3,ip)
  ux=pcl(4,ip)
  uy=pcl(5,ip)
  uz=pcl(6,ip)
  wp=pcl(7,ip)
  
  !IF (yp.LT.yhalf) THEN

    ! 4-velocity
    utot=sqrt(ux*ux+uy*uy+uz*uz)

    clba=uz/sqrt(ux*ux+uz*uz)
    slba=ux/sqrt(ux*ux+uz*uz)

      IF (slba.GT.0d0) THEN
      lambda=acos(clba)
      ELSE
      lambda=-acos(clba)
      END IF

    lambda=lambda*180.0/pi

    cphi=sqrt(ux*ux+uz*uz)/sqrt(ux*ux+uy*uy+uz*uz)
    sphi=uy/sqrt(ux*ux+uy*uy+uz*uz)

    phi=asin(sphi)
    phi=phi*180.0/pi
      
    iu2=FLOOR((LOG10(utot)-LOG10(umin))/((LOG10(umax)-LOG10(umin))/(NU-1)))+1
    iph2=FLOOR((phi-pmin)/((pmax-pmin)/(NPHI-1)))+1
    il2=FLOOR((lambda-lmin)/((lmax-lmin)/(NLBA-1)))+1

    IF ((iph2<1) .OR. (iph2>NPHI-1) .OR. (il2<1) .OR. (il2>NLBA-1) .OR. &
      (iu2<1) .OR. (iu2>NU-1)) THEN
      ! do nothing (should we warn?)
    ELSE
      ! dN/cos(phi)dphi*dlambda*du
      angularp(iu2,iph2,il2)=angularp(iu2,iph2,il2)+wp/( &
        (dSinPhi(iph2) * dLrad(il2)) * du(iu2))
    
      ! dN/du
      spectrap(iu2)=spectrap(iu2)+wp/du(iu2)
    ENDIF

ENDDO

! Total angular map
CALL MPI_REDUCE(angularp,angular,(NU-1)*(NPHI-1)*(NLBA-1),MPI_DOUBLE_PRECISION,&
                MPI_SUM,0,COMM,ierr)

! Total spectrum
CALL MPI_REDUCE(spectrap,spectra,NU-1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,ierr)
  
! Transpose array
CALL DUMP_GLOBAL_FIELD_3D("./data/angular_" // trim(spec) // "_" // &
  trim(sym), "field", INT8( (/ NU-1, NPHI-1, NLBA-1 /) ),angular, 0, &
  it, .TRUE., COMM, uref, "|u|", pref, "phi", lref, "lambda")

CALL DUMP_GLOBAL_FIELD_1D("./data/spectrum_" // trim(spec) // "_" // &
  trim(sym), "field", (/ INT8(NU-1) /), spectra/(4.0*pi), 0, &
  it, COMM, uref, "|u|")

CALL DUMP_COORDS("u", uref, uDumped, COMM)

CALL DUMP_COORDS("phi", pref, phiDumped, COMM)

CALL DUMP_COORDS("lambda", lref, lambdaDumped, COMM)

END SUBROUTINE SPECTRUM_ANGULAR

!***********************************************************************
! Subroutine MAP_XY
! This subroutine computes the total particles' spatial distribution
!
! INPUT:
! 
! - pcl: Particle distribution function
! - it: Timestep number
! - xr,yr: (global) reduced grid points
! - nrp: the shape of (local) reduced grid points on this domain 
! - irp: the index of the first red. grid point on this domain (in each dir)
! - spec: Particle species
! - sym: Name of the particles
!
! OUTPUT: writes dNxy to disk.
!***********************************************************************

SUBROUTINE MAP_XY(pcl,NPP, xr,yr,nrp, irp, it, spec, sym, COMM)

IMPLICIT NONE

! INPUT
INTEGER*8                                :: NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
DOUBLE PRECISION                         :: xr(:)
DOUBLE PRECISION                         :: yr(:)
INTEGER, DIMENSION(NDIM)                 :: nrp, irp
CHARACTER(len=*)                         :: spec, sym
INTEGER                                  :: it, COMM

! LOCAL
INTEGER*8                                 :: i
INTEGER                                   :: ix, iy
!INTEGER, DIMENSION(:), ALLOCATABLE        :: ix, iy
DOUBLE PRECISION                          :: idxr, idyr
DOUBLE PRECISION, ALLOCATABLE             :: dNxyp(:,:),dNxy(:,:)
!***********************************************************************

! Just in case, allocate 1 lower grid cell, and 2 upper grid cells
ALLOCATE(dNxyp(nrp(1)+3, nrp(2)+3)) !{
idxr = 1./(xr(2)-xr(1))
idyr = 1./(yr(2)-yr(1))

! a particle in the first local cell (e.g., x = xr(irp(1)) + tiny)
! gets mapped to the second bin of dNxyp

!ix = FLOOR((pcl(1,:) - xmin)*idxr) + (1 - irp(1) + 2)
!iy = FLOOR((pcl(2,:) - ymin)*idyr) + (1 - irp(2) + 2)
  
dNxyp = 0.
DO i=1,NPP
  ix = FLOOR((pcl(1,i) - xmin)*idxr) + (1 - irp(1) + 2)
  iy = FLOOR((pcl(2,i) - ymin)*idyr) + (1 - irp(2) + 2)
  dNxyp(ix,iy)=dNxyp(ix,iy)+pcl(7,i)*idxr*idyr
ENDDO

CALL DUMP_FIELD_NOUGC('./data/densities/mapxy_' // trim(spec) // "_" // trim(sym),&
  'field', dNxyp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")

DEALLOCATE(dNxyp) !}

END SUBROUTINE MAP_XY

!***********************************************************************
! Subroutine MAP_FLUID
! This subroutine computes the relativistic fluid quantities: the plasma flux,
! the pressure tensor, the energy density and the momentum density.
!
! INPUT:
! 
! - mass: particle mass
! - pcl: Particle distribution function
! - it: Timestep number
! - xr,yr: (global) reduced grid points
! - nrp: the shape of (local) reduced grid points on this domain 
! - irp: the index of the first red. grid point on this domain (in each dir)
! - spec: Particle species
! - sym: Name of the particles
!
! OUTPUT: writes Fx-y-z, Ue, Upx-y-z, and Pij to disk.
!***********************************************************************

SUBROUTINE MAP_FLUID(mass,pcl,NPP, xr,yr,nrp, irp, it, spec, sym, COMM)

IMPLICIT NONE

! INPUT
DOUBLE PRECISION                         :: mass
INTEGER*8                                :: NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
DOUBLE PRECISION                         :: xr(:)
DOUBLE PRECISION                         :: yr(:)
INTEGER, DIMENSION(NDIM)                 :: nrp, irp
CHARACTER(len=*)                         :: spec, sym
INTEGER                                  :: it, COMM

! LOCAL
INTEGER*8                                 :: i
INTEGER                                   :: ix, iy
!INTEGER, DIMENSION(:), ALLOCATABLE        :: ix, iy
DOUBLE PRECISION                          :: idxr, idyr, iArea
DOUBLE PRECISION                          :: xp, yp, zp, ux, uy, uz, wt
DOUBLE PRECISION                          :: gam
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Fxp, Fyp, Fzp 
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Uep, Upxp, Upyp, Upzp 
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Pxxp,Pyyp,Pzzp
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Pxyp,Pxzp,Pyzp

!***********************************************************************

! Just in case, allocate 1 lower grid cell, and 2 upper grid cells
ALLOCATE(Fxp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Fyp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Fzp(nrp(1)+3,nrp(2)+3)) !{

ALLOCATE(Uep(nrp(1)+3,nrp(2)+3)) !{

ALLOCATE(Upxp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Upyp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Upzp(nrp(1)+3,nrp(2)+3)) !{

ALLOCATE(Pxxp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Pyyp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Pzzp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Pxyp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Pxzp(nrp(1)+3,nrp(2)+3)) !{
ALLOCATE(Pyzp(nrp(1)+3,nrp(2)+3)) !{

Fxp=0.0
Fyp=0.0
Fzp=0.0
Uep=0.0
Upxp=0.0
Upyp=0.0
Upzp=0.0
Pxxp=0.0
Pyyp=0.0
Pzzp=0.0
Pxyp=0.0
Pxzp=0.0
Pyzp=0.0

idxr = 1./(xr(2)-xr(1))
idyr = 1./(yr(2)-yr(1))
iArea = idxr*idyr

! a particle in the first local cell (e.g., x = xr(irp(1)) + tiny)
! gets mapped to the second bin of dNxyp

!ix = FLOOR((pcl(1,:) - xmin)*idxr) + (1 - irp(1) + 2)
!iy = FLOOR((pcl(2,:) - ymin)*idyr) + (1 - irp(2) + 2)
  
!===================================================================
DO i=1,NPP

  xp=pcl(1,i)
  yp=pcl(2,i)
  ux=pcl(4,i)
  uy=pcl(5,i)
  uz=pcl(6,i)
  wt=pcl(7,i)
  
  ix = FLOOR((xp - xmin)*idxr) + (1 - irp(1) + 2)
  iy = FLOOR((yp - ymin)*idyr) + (1 - irp(2) + 2)

  gam=sqrt(1d0+ux*ux+uy*uy+uz*uz)
  
  ! Particle flux
  Fxp(ix,iy)=Fxp(ix,iy)+(wt*c*ux/gam)
  Fyp(ix,iy)=Fyp(ix,iy)+(wt*c*uy/gam)
  Fzp(ix,iy)=Fzp(ix,iy)+(wt*c*uz/gam)
  
  ! Plasma energy density Ue=<gamma>*mass*c^2*n
  Uep(ix,iy)=Uep(ix,iy)+(wt*gam*mass*c*c)
  
  ! Plasma momentum density Up=<u>*mass*c*n
  Upxp(ix,iy)=Upxp(ix,iy)+(wt*ux*mass*c)
  Upyp(ix,iy)=Upyp(ix,iy)+(wt*uy*mass*c)
  Upzp(ix,iy)=Upzp(ix,iy)+(wt*uz*mass*c)
  
  ! Plasma pressure tensor components
  Pxxp(ix,iy)=Pxxp(ix,iy)+(wt*mass*c*c/gam)*ux*ux
  Pyyp(ix,iy)=Pyyp(ix,iy)+(wt*mass*c*c/gam)*uy*uy
  Pzzp(ix,iy)=Pzzp(ix,iy)+(wt*mass*c*c/gam)*uz*uz
  Pxyp(ix,iy)=Pxyp(ix,iy)+(wt*mass*c*c/gam)*ux*uy
  Pxzp(ix,iy)=Pxzp(ix,iy)+(wt*mass*c*c/gam)*ux*uz
  Pyzp(ix,iy)=Pyzp(ix,iy)+(wt*mass*c*c/gam)*uy*uz
  
ENDDO

Fxp = Fxp * iArea
Fyp = Fyp * iArea
Fzp = Fzp * iArea
Uep = Uep * iArea
Upxp = Upxp * iArea
Upyp = Upyp * iArea
Upzp = Upzp * iArea
Pxxp = Pxxp * iArea
Pyyp = Pyyp * iArea
Pzzp = Pzzp * iArea
Pxyp = Pxyp * iArea
Pxzp = Pxzp * iArea
Pyzp = Pyzp * iArea

CALL DUMP_FIELD_NOUGC('./data/densities/Fx_' // trim(spec) // "_" // trim(sym),&
  'field', Fxp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Fy_' // trim(spec) // "_" // trim(sym),&
  'field', Fyp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Fz_' // trim(spec) // "_" // trim(sym),&
  'field', Fzp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Ue_' // trim(spec) // "_" // trim(sym),&
  'field', Uep(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Upx_' // trim(spec) // "_" // trim(sym),&
  'field', Upxp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Upy_' // trim(spec) // "_" // trim(sym),&
  'field', Upyp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Upz_' // trim(spec) // "_" // trim(sym),&
  'field', Upzp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Pxx_' // trim(spec) // "_" // trim(sym),&
  'field', Pxxp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Pyy_' // trim(spec) // "_" // trim(sym),&
  'field', Pyyp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Pzz_' // trim(spec) // "_" // trim(sym),&
  'field', Pzzp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Pxy_' // trim(spec) // "_" // trim(sym),&
  'field', Pxyp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Pxz_' // trim(spec) // "_" // trim(sym),&
  'field', Pxzp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD_NOUGC('./data/densities/Pyz_' // trim(spec) // "_" // trim(sym),&
  'field', Pyzp(2:nrp(1)+1, 2:nrp(2)+1), it, COMM,xr,"x", yr,"y")

DEALLOCATE(Fxp) !}
DEALLOCATE(Fyp) !}
DEALLOCATE(Fzp) !}
DEALLOCATE(Uep) !}
DEALLOCATE(Upxp) !}
DEALLOCATE(Upyp) !}
DEALLOCATE(Upzp) !}
DEALLOCATE(Pxxp) !}
DEALLOCATE(Pyyp) !}
DEALLOCATE(Pzzp) !}
DEALLOCATE(Pxyp) !}
DEALLOCATE(Pxzp) !}
DEALLOCATE(Pyzp) !}

END SUBROUTINE MAP_FLUID

!***********************************************************************
! Subroutine DUMP_COORDS
! This subroutine writes a 1D array of coordinates to text file,
!   for files dumped once per simulation.
! This should be called by all ranks (currently only rank 0 does anything).
!
! INPUT:
! 
! - fileBaseName the name of the file, to which
!     "./data/" will be prepended, and ".dat" appended
! - x: array of (global) coordinate values
! - wasDumped: whether these coordinates have been dumped
!        (if .TRUE., nothing will be done, if .FALSE., x will be dumped)
! - COMM: mpiComm
!
! OUTPUT:
! - wasDumped: .TRUE. upon return
!
!***********************************************************************

SUBROUTINE DUMP_COORDS(fileBaseName, x, wasDumped, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
CHARACTER(len=*)                         :: fileBaseName
DOUBLE PRECISION, DIMENSION(:)           :: x
LOGICAL, INTENT(INOUT)                   :: wasDumped
INTEGER                                  :: COMM
! LOCAL
INTEGER :: i, rank, mpiErr
!***********************************************************************

IF (.NOT. wasDumped) THEN

CALL MPI_COMM_RANK(COMM, rank, mpiErr)

IF ((rank == 0)) THEN
  OPEN(9,FILE="./data/" // fileBaseName // ".dat")
  DO i=1,SIZE(x)
  WRITE(9,"(1" // FMT_DBL // ")") (x(i))
  ENDDO
  CLOSE(9)
ENDIF

ENDIF

wasDumped = .TRUE.

END SUBROUTINE DUMP_COORDS

!***********************************************************************
! Subroutine DUMP_FIELD_NOUGC
! This subroutine writes a 2D array, decomposed in 2D across all ranks,
! to disk -- without the (global) uppermost guard cells 
!
! INPUT:
! 
! - fileBaseName the base name of the file, 
!     to which "_DUMPNUM.h5" will be appended
! - datasetname the name of the field dataset
! - F: array of (local domain) field values
! - xgp,ygp: coordinates of (local domain) field values
! - it: Timestep number
! - COMM: mpi comm
!
! OUTPUT:
! - ierr: error code
!
!***********************************************************************

SUBROUTINE DUMP_FIELD_NOUGC(fileBaseName, datasetname, F, it,COMM,&
  xgp, xLabel, ygp, yLabel, zgp, zLabel)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER                       :: FDIM = NDIM

! INPUT
CHARACTER(len=*)                         :: fileBaseName
CHARACTER(len=*)                         :: datasetname
! following must have rank=NDIM
DOUBLE PRECISION, DIMENSION(:,:)         :: F
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xgp, ygp, zgp
CHARACTER(len=*), OPTIONAL               :: xLabel, yLabel, zLabel
INTEGER                                  :: it
INTEGER                                  :: COMM

! LOCAL
INTEGER, DIMENSION(FDIM) :: ubs
INTEGER                  :: d, mpiErr
INTEGER, DIMENSION(FDIM) :: domainDecomp, domainIndex
!***********************************************************************

CALL FILL_DOMAIN_DECOMP(COMM, mpiErr, domainDecomp, domainIndex)

ubs = SHAPE(F)
DO d = 1, FDIM
  IF (domainIndex(d)+1==domainDecomp(d)) ubs(d) = ubs(d) - 1
ENDDO

CALL DUMP_FIELD(fileBaseName, datasetname, F(1:ubs(1),1:ubs(2)), it,COMM,&
  xgp, xLabel, ygp, yLabel, zgp, zLabel)

END SUBROUTINE DUMP_FIELD_NOUGC

!***********************************************************************
! Subroutine REDUCE_GRID_1D
! This subroutine calculates a reduced grid.
!
! INPUT:
! 
! - dir: the direction (1=x, 2=y, 3=z)
! - xg: global grid points (x represents x,y, or z, depending on dir)
! - xminp, xmaxp: the edges of this domain
!   xmaxp must equal the xminp of the next-higher domain (in x)
!   exactly.
! - method
!   0: Include (global) upper guard cell; exclude local upper guard cell.
!      Spaces the reduced points evenly (with dxr) within the range of xg, 
!      with xr(1) = xpg(1), and xr(last) = xg(last)
!      Local reduced points do not include upper guard cell, except
!      for the global upper guard cell in included in the uppermost domain.
!   1: Include (global and local) upper guard cell.
!      Spaces the reduced points evenly (with dxr) within the range of xg, 
!      with xr(1) = xpg(1), and xr(last) = xg(last)
!      Local reduced points include upper guard cell.
!   2: Exclude (global and local) upper guard cell.
!      Spaces the reduced points evenly (with dxr) within the range of xg, 
!      with xr(1) = xpg(1), and xr(last) = xg(last) - dxr
!      Local reduced points do not include upper guard cell.
!   4, 5: Decrease the number of reduced grid points until it equals a
!      multiple of the number of mpi ranks (in direction dir) plus 1
!      (for the global upper guard cell).  
!      Then use method 0 (if method=4) or 1 (method=5).
!   6: Decrease the number of reduced grid points until it equals a
!      multiple of the number of mpi ranks (in direction dir).
!      Then use method 2.
! - COMM: mpi comm
! - xr: reduced grid points (xr(n) is undefined for n > nr)
! - nr: the number of reduced grid points actually given
! - nrp: the number of reduced grid points on this mpi domain 
! - ixp: xr(ixp) is the first grid point on this mpi domain
!
!***********************************************************************

SUBROUTINE REDUCE_GRID_1D(dir, xg, xminp, xmaxp, method, COMM, &
  nr, nrp, ixp, xr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: dir
DOUBLE PRECISION                         :: xg(:)
INTEGER, INTENT(IN)                      :: method
INTEGER                                  :: COMM
DOUBLE PRECISION, INTENT(OUT)            :: xr(:)
INTEGER, INTENT(OUT)                     :: nr, nrp, ixp

! LOCAL
INTEGER                  :: methodCopy
INTEGER                  :: mpiErr, i
INTEGER                  :: nx, nrpOwned
DOUBLE PRECISION         :: L, dxr, xminp, xmaxp
INTEGER, DIMENSION(NDIM) :: domainDecomp, domainIndex
INTEGER                  :: np
INTEGER                  :: ixpRa(1), nrpRa(1)
!***********************************************************************

CALL FILL_DOMAIN_DECOMP(COMM, mpiErr, domainDecomp, domainIndex)

np = domainDecomp(dir)
nx = SIZE(xg)
nr = SIZE(xr)

methodCopy = method
IF ((methodCopy==4).OR.(methodCopy==5)) THEN
  nr = ((nr-1)/np) * np + 1
  nrpOwned = (nr-1)/np
  methodCopy = methodCopy - 4
ELSEIF (methodCopy==6) THEN
  nr = (nr/np) * np
  nrpOwned = nr
  methodCopy = methodCopy - 4
ENDIF

L = xg(nx) - xg(1)

IF (method == 2) THEN
  dxr = L / nr 
ELSE
  dxr = L / (nr-1)
ENDIF

DO i = 1, nr
  xr(i) = xg(1) + (i-1) * dxr
ENDDO
! Assure that xr(nr) = xg(nx), in spite of finite-precision arithmetic
IF (methodCopy /= 2) THEN
  xr(nr) = xg(nx)
ENDIF

! figure out which reduced points belong to this domain
IF (method >= 4) THEN
  ! this is easy
  ! nrpOwned is already set
  ixp = domainIndex(dir) * nrpOwned + 1
  ! I have seen a case for 6 domains, 240 cells (40 cells/domain)
  ! interpolated to a reduced grid of 351, reduced to 348 cells/domain,
  ! for L=100.
  ! where the first point (on global and reduced grid) on domain 5
  ! should be 83.3333333333333, but the reduced point was slightly 
  ! smaller:   83.333333333333329 (red) vs. 83.333333333333343 (global)
  xr(ixp) = MAX(xr(ixp), xminp)
  ! :TODO: worry about upper point
ELSE
  ! I worry about finite precision; if a reduced point should be
  ! exactly on a domain boundary, then finite precision can result in
  ! it being a little on either side.
  ixpRa = MINLOC( xr - xminp, MASK = (xr - xminp >= 0) )
  ixp = ixpRa(1)
  IF (domainIndex(dir) + 1 == np) THEN
    IF (methodCopy == 2) THEN
      nrpOwned = nr + 1 - ixp ! include uppermost
    ELSE 
      nrpOwned = nr - ixp ! exclude uppermost
    ENDIF
  ELSE
    nrpRa = MINLOC( xr - xmaxp, MASK = (xr - xmaxp >= 0) ) - ixp
    nrpOwned = nrpRa(1)
  ENDIF

ENDIF
 
! does not include upper guard cell
nrp = nrpOwned
IF (methodCopy == 0) THEN
  ! include upper guard cell only for uppermost domain
  IF (domainIndex(dir) + 1 == np) nrp = nrp + 1
ELSEIF (methodCopy == 1) THEN
  ! include upper guard for every domain
  nrp = nrp + 1
ENDIF

IF (.FALSE. .AND. writePerRankFiles) THEN
  OPEN(9, FILE=trim(perRankFile), POSITION='APPEND')
  WRITE(9, *) "GridRed1D, dir=", dir
  WRITE(9, *) 'SIZE(xg) =', SIZE(xg)
  WRITE(9, *) 'xg =', xg
  WRITE(9, *) 'xr =', xr
  WRITE(9, *) "ixp, nrp, nr:", ixp, nrp, nr
ENDIF

END SUBROUTINE REDUCE_GRID_1D

!***********************************************************************
! Subroutine REDUCE_GRID
! This subroutine calculates a reduced grid.
!
! INPUT:
! - xg, yg: global grid points 
! - minp, maxp: the edges of this domain (in each dir)
! - COMM: mpi comm
! OUTPUT:
! - xr, yr: global reduced grid points (xr(n) is undefined for n > nr)
! - nr: the number of reduced grid points actually given
! - nrp: the number of reduced grid points on this mpi domain 
! - ip: xr(ixp(1)), yr(ipy(2)) is the first grid point on this mpi domain
!
!***********************************************************************

SUBROUTINE REDUCE_GRID(xg, yg, minp, maxp, COMM, nr, nrp, ip, xr, yr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: dir
DOUBLE PRECISION, DIMENSION(:)           :: xg, yg
DOUBLE PRECISION, DIMENSION(NDIM)        :: minp, maxp
INTEGER                                  :: COMM
DOUBLE PRECISION, INTENT(OUT)            :: xr(:), yr(:)
INTEGER, DIMENSION(NDIM), INTENT(OUT)    :: nr, nrp, ip

INTEGER                                  :: method = 0
!***********************************************************************

CALL REDUCE_GRID_1D(1, xg, minp(1), maxp(1), method, COMM, &
  nr(1), nrp(1), ip(1), xr)
CALL REDUCE_GRID_1D(2, yg, minp(2), maxp(2), method, COMM, &
  nr(2), nrp(2), ip(2), yr)

END SUBROUTINE REDUCE_GRID

!***********************************************************************
! Subroutine REDUCE_FIELDS
! This subroutine reduces a field (or fields) to a desired size,
! to use less disk space.  Reducing multiple fields at once has the
! advantage of calculating interpolation weights once for all fields.
!
! INPUT:
! 
! - xgp,ygp: grid positions (on local domain) where field values are located
! - xgpNew,ygpNew: grid values where the reduced field values are located
!   N.B. each xgpNew must be between xgp(1) and xpg(NXP), same for ygpNew.
!   Otherwise simulation may crash, or worse.
! - F1, F2...: (local) field values at xgp,ygp 
! - F1r, F2r...: field values at xgpNew,ygpNew
!   - these must have shape (SIZE(xgpNew), SIZE(ygpNew))
!   - If F2 is given, then F2r must be given, etc.
!
!***********************************************************************

SUBROUTINE REDUCE_FIELDS(xgp,ygp, xgpNew, ygpNew, F1, F1r, &
  F2, F2r, F3, F3r, F4, F4r, F5, F5r, F6, F6r)

IMPLICIT NONE

INCLUDE 'mpif.h'

DOUBLE PRECISION, DIMENSION(1:NXP)                 :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)                 :: ygp
DOUBLE PRECISION, DIMENSION(:)                     :: xgpNew, ygpNew
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP)           :: F1
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP), OPTIONAL :: F2, F3, F4, F5, F6

DOUBLE PRECISION, DIMENSION(:,:)           :: F1r
DOUBLE PRECISION, DIMENSION(:,:), OPTIONAL :: F2r, F3r, F4r, F5r, F6r

! Raw fields per domain

DOUBLE PRECISION, DIMENSION(1:6)         :: Fields
INTEGER                                  :: i,j, numFields

!***********************************************************************

!print *, "rank ", globalMpiRank, "interpolating from (local) y ", ygp, " to ", ygpNew 

numFields = 1
IF (PRESENT(F2r)) numFields = 2
IF (PRESENT(F3r)) numFields = 3
IF (PRESENT(F4r)) numFields = 4
IF (PRESENT(F5r)) numFields = 5
IF (PRESENT(F6r)) numFields = 6

! Note: it would be interesting to know if the compiler will optimize
! out "IF PRESENT()" statements.  I suspect not...also I'm not sure how
! f90 handles passing optional arguments to other functions taking
! optional arguments.  So playing it safe for now, and writing loop 6 times.

SELECT CASE(numFields)
  CASE(1)
    DO j=1,SIZE(ygpNew)
      DO i=1,SIZE(xgpNew)
      CALL BILINEAR_FIELDS(xgp,ygp,xgpNew(i),ygpNew(j), Fields, F1)
      F1r(i,j)=Fields(1)
      ENDDO
    ENDDO
  CASE(2)
    DO j=1,SIZE(ygpNew)
      DO i=1,SIZE(xgpNew)
      CALL BILINEAR_FIELDS(xgp,ygp,xgpNew(i),ygpNew(j), Fields, F1, F2)
      F1r(i,j)=Fields(1)
      F2r(i,j)=Fields(2)
      ENDDO
    ENDDO
  CASE(3)
    DO j=1,SIZE(ygpNew)
      DO i=1,SIZE(xgpNew)
      CALL BILINEAR_FIELDS(xgp,ygp,xgpNew(i),ygpNew(j), Fields, F1, F2, F3)
      F1r(i,j)=Fields(1)
      F2r(i,j)=Fields(2)
      F3r(i,j)=Fields(3)
      ENDDO
    ENDDO
  CASE(4)
    DO j=1,SIZE(ygpNew)
      DO i=1,SIZE(xgpNew)
      CALL BILINEAR_FIELDS(xgp,ygp,xgpNew(i),ygpNew(j), Fields, F1, F2, F3, &
        F4)
      F1r(i,j)=Fields(1)
      F2r(i,j)=Fields(2)
      F3r(i,j)=Fields(3)
      F4r(i,j)=Fields(4)
      ENDDO
    ENDDO
  CASE(5)
    DO j=1,SIZE(ygpNew)
      DO i=1,SIZE(xgpNew)
      CALL BILINEAR_FIELDS(xgp,ygp,xgpNew(i),ygpNew(j), Fields, F1, F2, F3, &
        F4, F5)
      F1r(i,j)=Fields(1)
      F2r(i,j)=Fields(2)
      F3r(i,j)=Fields(3)
      F4r(i,j)=Fields(4)
      F5r(i,j)=Fields(5)
      ENDDO
    ENDDO
  CASE(6)
    DO j=1,SIZE(ygpNew)
      DO i=1,SIZE(xgpNew)
      CALL BILINEAR_FIELDS(xgp,ygp,xgpNew(i),ygpNew(j), Fields, F1, F2, F3, &
        F4, F5, F6)
      F1r(i,j)=Fields(1)
      F2r(i,j)=Fields(2)
      F3r(i,j)=Fields(3)
      F4r(i,j)=Fields(4)
      F5r(i,j)=Fields(5)
      F6r(i,j)=Fields(6)
      ENDDO
    ENDDO
ENDSELECT

END SUBROUTINE REDUCE_FIELDS

!***********************************************************************
! Subroutine DUMP_REDUCED_FIELDS
! This subroutine interpolates the E,B fields to a reduced grid and
! saves it to disk (using a reduced grid saves disk space).
!
! INPUT:
! 
! - xgp,ygp: Grid of domain
! - xr,yr: (global) reduced grid points
! - nrp: the shape of (local) reduced grid points on this domain 
! - irp: the index of the first red. grid point on this domain (in each dir)
! - Eg,Bg: Fields at nodes
! - it: Timestep number
! - coords: coordinates of each process in the cartesian topology
!
!***********************************************************************

SUBROUTINE DUMP_REDUCED_FIELDS(Exg,Eyg,Ezg,Bxg,Byg,Bzg,xgp,ygp, &
  xr,yr,nrp, irp, it, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: COMM
INTEGER                                  :: it

! Raw grid per domain
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION                         :: xr(:)
DOUBLE PRECISION                         :: yr(:)
INTEGER, DIMENSION(NDIM)                 :: nrp, irp

! Raw fields per domain
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bxg,Byg,Bzg,Exg,Eyg,Ezg

! Reduced fields per domain
DOUBLE PRECISION, ALLOCATABLE            :: Bxrp(:,:),Byrp(:,:),Bzrp(:,:)
DOUBLE PRECISION, ALLOCATABLE            :: Exrp(:,:),Eyrp(:,:),Ezrp(:,:)

! Reduced fields per domain in the global reduced grid
DOUBLE PRECISION, ALLOCATABLE            :: xrp(:), yrp(:)

INTEGER                                  :: mpiErr
INTEGER                                  :: i,j
INTEGER                                  :: NXRP,NYRP
INTEGER                                  :: NXRPM,NYRPM
!***********************************************************************

ALLOCATE(xrp(nrp(1)))
ALLOCATE(yrp(nrp(2)))
xrp = xr(irp(1):irp(1)+nrp(1)-1)
yrp = yr(irp(2):irp(2)+nrp(2)-1)

NXRPM = SIZE(xrp)
NYRPM = SIZE(yrp)
ALLOCATE(Bxrp(NXRPM, NYRPM))
ALLOCATE(Byrp(NXRPM, NYRPM))
ALLOCATE(Bzrp(NXRPM, NYRPM))
ALLOCATE(Exrp(NXRPM, NYRPM))
ALLOCATE(Eyrp(NXRPM, NYRPM))
ALLOCATE(Ezrp(NXRPM, NYRPM))

CALL REDUCE_FIELDS(xgp,ygp, xrp,yrp, Bxg,Bxrp, Byg,Byrp, Bzg,Bzrp, &
  Exg, Exrp, Eyg, Eyrp, Ezg, Ezrp)

! :TODO: psi used to be calculated and dumped when the fields were
! gathered on one rank.  Do we need to do this?

! For backward compatibility, we need the upper guard field values;
! but with the recent refactoring, we don't have those.  Either get them,
! or decide we don't need them.
! Don't write the upper guard cell
!NXRP = NXRPM
!NYRP = NYRPM
!IF (DUMP_FIELDS_PARALLEL .AND. .NOT. USE_HDF5) THEN
!! Except, for backward compatibility, when dumping in text format
!! with one file per mpi rank, include the upper guard cell in each file
!  NXRP = NXRPM+1
!  NYRP = NYRPM+1
!ENDIF

CALL DUMP_FIELD('./data/fields/Bx', 'field', Bxrp,it, COMM,xr,"x", yr,"y")
DEALLOCATE(Bxrp)
CALL DUMP_FIELD('./data/fields/By', 'field', Byrp,it, COMM,xr,"x", yr,"y")
DEALLOCATE(Byrp)
CALL DUMP_FIELD('./data/fields/Bz', 'field', Bzrp,it, COMM,xr,"x", yr,"y")
DEALLOCATE(Bzrp)
CALL DUMP_FIELD('./data/fields/Ex', 'field', Exrp,it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD('./data/fields/Ey', 'field', Eyrp,it, COMM,xr,"x", yr,"y")
CALL DUMP_FIELD('./data/fields/Ez', 'field', Ezrp,it, COMM,xr,"x", yr,"y")
DEALLOCATE(Exrp,Eyrp,Ezrp)

DEALLOCATE(xrp,yrp)

END SUBROUTINE DUMP_REDUCED_FIELDS

!***********************************************************************
! Subroutine DUMP_REDUCED_RHOJ
! This subroutine writes the E,B fields array to disk. Each array has a size fixed
! by NFIELDX and NFIELDY, in order to avoid memory issue for large simulations.
!
! INPUT:
! 
! - xgp,ygp: Grid of domain
! - xr,yr: (global) reduced grid points
! - nrp: the shape of (local) reduced grid points on this domain 
! - irp: the index of the first red. grid point on this domain (in each dir)
! - rho,J: Fields at nodes
! - it: Timestep number
! - coords: coordinates of each process in the cartesian topology
!
!***********************************************************************

SUBROUTINE DUMP_REDUCED_RHOJ(rho,Jx,Jy,Jz,xgp,ygp, &
  xr,yr,nrp, irp, spec, it, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: COMM
CHARACTER(len=*)                         :: spec
INTEGER                                  :: it

! Raw grid per domain
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION                         :: xr(:)
DOUBLE PRECISION                         :: yr(:)
INTEGER, DIMENSION(NDIM)                 :: nrp, irp

! Raw fields per domain
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: rho,Jx,Jy,Jz

! Reduced fields per domain
DOUBLE PRECISION, ALLOCATABLE            :: rhorp(:,:)
DOUBLE PRECISION, ALLOCATABLE            :: Jxrp(:,:),Jyrp(:,:),Jzrp(:,:)

! Reduced fields per domain in the global reduced grid
DOUBLE PRECISION, ALLOCATABLE            :: xrp(:), yrp(:)

INTEGER                                  :: i,j
INTEGER                                  :: NXRP,NYRP
INTEGER                                  :: NXRPM,NYRPM

!***********************************************************************

spec = adjustl(spec)

ALLOCATE(xrp(nrp(1)))
ALLOCATE(yrp(nrp(2)))
xrp = xr(irp(1):irp(1)+nrp(1)-1)
yrp = yr(irp(2):irp(2)+nrp(2)-1)

NXRPM = SIZE(xrp)
NYRPM = SIZE(yrp)
ALLOCATE(rhorp(NXRPM, NYRPM))
ALLOCATE(Jxrp(NXRPM, NYRPM))
ALLOCATE(Jyrp(NXRPM, NYRPM))
ALLOCATE(Jzrp(NXRPM, NYRPM))

CALL REDUCE_FIELDS(xgp,ygp, xrp,yrp, rho,rhorp, Jx,Jxrp, Jy,Jyrp, Jz,Jzrp)

! For backward compatibility, we need the upper guard field values;
! but with the recent refactoring, we don't have those.  Either get them,
! or decide we don't need them.
! Don't write the upper guard cell
!NXRP = NXRPM
!NYRP = NYRPM
!IF (DUMP_FIELDS_PARALLEL .AND. .NOT. USE_HDF5) THEN
!! Except, for backward compatibility, when dumping in text format
!! with one file per mpi rank, include the upper guard cell in each file
!  NXRP = NXRPM+1
!  NYRP = NYRPM+1
!ENDIF

CALL DUMP_FIELD('./data/currents/rho_' // trim(spec), 'field', &
  rhorp,it, COMM,xr,"x", yr,"y")
DEALLOCATE(rhorp)
CALL DUMP_FIELD('./data/currents/Jx_' // trim(spec), 'field', &
  Jxrp,it, COMM,xr,"x", yr,"y")
DEALLOCATE(Jxrp)
CALL DUMP_FIELD('./data/currents/Jy_' // trim(spec), 'field', &
  Jyrp,it, COMM,xr,"x", yr,"y")
DEALLOCATE(Jyrp)
CALL DUMP_FIELD('./data/currents/Jz_' // trim(spec), 'field', &
  Jzrp,it, COMM,xr,"x", yr,"y")
DEALLOCATE(Jzrp)

DEALLOCATE(xrp,yrp)

END SUBROUTINE DUMP_REDUCED_RHOJ

!***********************************************************************
! Subroutine DUMP_GRID
!  Writes grid data to disk.
!
! INPUT:
! 
! - gridName: x-data will be saved to "x" + name + ".dat", etc.
! - xs, ys: Global grid points
! - COMM mpi comm
!
!***********************************************************************

SUBROUTINE DUMP_GRID(gridName, xs, ys, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

CHARACTER(len=*)                         :: gridName
DOUBLE PRECISION, DIMENSION(:)           :: xs, ys
INTEGER                                  :: COMM

INTEGER                                  :: i,j, rank, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)

IF (rank.EQ.0) THEN

OPEN(9,FILE="./data/x" // gridName // ".dat")
DO i=1,SIZE(xs)
WRITE(9,"(1" // FMT_DBL // ")") xs(i)
ENDDO
CLOSE(9)

OPEN(9,FILE="./data/y" // gridName // ".dat")
DO j=1,SIZE(ys)
WRITE(9,"(1" // FMT_DBL // ")") ys(j)
ENDDO
CLOSE(9)

END IF

END SUBROUTINE DUMP_GRID

!***********************************************************************
! Subroutine EM_ENERGY
! This subroutine computes the total electric and magnetic energies in 
! each domain.
!
! INPUT:
! 
! - Bx: x-component of the magnetic field at nodes
! - By: y-component of the magnetic field at nodes
! - Bz: z-component of the magnetic field at nodes
! - Ex: x-component of the magnetic field at nodes
! - Ey: y-component of the magnetic field at nodes
! - Ez: z-component of the magnetic field at nodes
!
! OUTPUT: Emag,Eelc in erg
!***********************************************************************

SUBROUTINE EM_ENERGY(Bx,By,Bz,Ex,Ey,Ez,xgp,ygp,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: id,COMM,ierr
INTEGER                                  :: ix,iy
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bx,By,Bz,Ex,Ey,Ez
DOUBLE PRECISION                         :: Emag,Eelc,Emagp,Eelcp
!***********************************************************************

! Total energies
Emag=0.0
Eelc=0.0

! Energies per process
Emagp=0.0
Eelcp=0.0

DO ix=1,NXP-1
  DO iy=1,NYP-1

  ! Magnetic energy per domain
  Emagp=Emagp+(Bx(ix,iy)*Bx(ix,iy)+By(ix,iy)*By(ix,iy)+&
  Bz(ix,iy)*Bz(ix,iy))/(8.0*pi)*(xgp(ix+1)-xgp(ix))*(ygp(iy+1)-ygp(iy))

  !Electric energy per domain
  Eelcp=Eelcp+(Ex(ix,iy)*Ex(ix,iy)+Ey(ix,iy)*Ey(ix,iy)+&
  Ez(ix,iy)*Ez(ix,iy))/(8.0*pi)*(xgp(ix+1)-xgp(ix))*(ygp(iy+1)-ygp(iy))

  ENDDO
ENDDO

! Total magnetic energy
CALL MPI_REDUCE(Emagp,Emag,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,ierr)

! Total electric energy
CALL MPI_REDUCE(Eelcp,Eelc,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,ierr)

! Writing the data to disk. 
! The option 'APPEND' enables to append data to the existing file
IF (id.EQ.0) THEN
  CALL ADD_LOG_DATA(EemLog, (/ Emag, Eelc /) )
END IF

END SUBROUTINE EM_ENERGY

!***********************************************************************
! Subroutine KIN_ENERGY
! This subroutine computes the total particles' kinetic energy per species
!
! INPUT:
! 
! - mass: Mass of the particle
! - pcl:  Particle distribution function
! - spec: Species of the particles (i.e., electrons or ions)
! - sym:  Type of particles (i.e., bg or drift)
!
! OUTPUT: Ekin_spec
!***********************************************************************

SUBROUTINE KIN_ENERGY(mass,pcl,NPP,spec,sym,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER*8                              :: ip,NPP
INTEGER                                :: id,COMM,ierr
DOUBLE PRECISION, DIMENSION(1:7,1:NPP) :: pcl
DOUBLE PRECISION                       :: Ekin,Ekinp,mass,gam
CHARACTER(LEN=10)                      :: spec,sym
!***********************************************************************

spec=adjustl(spec)
sym=adjustl(sym)

! Total kinetic energy
Ekin=0.0

! Kinetic energy per process
Ekinp=0.0

DO ip=1,NPP
gam=sqrt(1.0+pcl(4,ip)*pcl(4,ip)+pcl(5,ip)*pcl(5,ip)+pcl(6,ip)*pcl(6,ip))
Ekinp=Ekinp+pcl(7,ip)*gam*mass*c*c
ENDDO

! Total kinetic energy background particles
CALL MPI_REDUCE(Ekinp,Ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,ierr)

! Writing the data to disk. 
IF (id.EQ.0) THEN
  IF (trim(spec) == "electrons" .AND. trim(sym) == "bg") THEN
    CALL ADD_LOG_DATA(Ekin_ebLog, (/ Ekin /) )
  ELSEIF (trim(spec) == "electrons" .AND. trim(sym) == "drift") THEN
    CALL ADD_LOG_DATA(Ekin_edLog, (/ Ekin /) )
  ELSEIF (trim(spec) == "ions" .AND. trim(sym) == "bg") THEN
    CALL ADD_LOG_DATA(Ekin_pbLog, (/ Ekin /) )
  ELSEIF (trim(spec) == "ions" .AND. trim(sym) == "drift") THEN
    CALL ADD_LOG_DATA(Ekin_pdLog, (/ Ekin /) )
  ENDIF
ENDIF

END SUBROUTINE KIN_ENERGY

!***********************************************************************
! Subroutine RAD_ENERGY
! This subroutine computes the total radiative energy losses per processes
!
! INPUT:
! 
! - Erad_bg:    Total radiative energy lost by the bg particles per domain
! - Erad_drift: Total radiative energy lost by the drift particles per domain
! - Erad: Total radiative energy lost by all the particles
! - spec: Particle species
! - rad_process: Name of the radiative process considered
!
! OUTPUT: Erad at time t+dt
!
!***********************************************************************

SUBROUTINE RAD_ENERGY(Erad_bg,Erad_drift,Erad,spec,rad_process,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                             :: id,COMM,ierr
DOUBLE PRECISION                    :: Erad_bg,Erad_drift,Eradd,Erad,Etemp
CHARACTER(LEN=10)                   :: spec,rad_process
!***********************************************************************

spec=adjustl(spec)
rad_process=adjustl(rad_process)

! Total radiative energy lost between t and t+dt
Etemp=0.0

! Total radiative energy per domain lost between t and t+dt
Eradd=Erad_bg+Erad_drift

! Reduction
CALL MPI_REDUCE(Eradd,Etemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,ierr)

! Total energy lost at time t+dt
Erad=Erad+Etemp

! Writing the data to disk. 
IF (id.EQ.0) THEN
  IF (trim(rad_process) == "syn" .AND. trim(spec) == "electrons") THEN
    CALL ADD_LOG_DATA(Esyn_eLog, (/ Erad /) )
  ELSEIF (trim(rad_process) == "syn" .AND. trim(spec) == "ions") THEN
    CALL ADD_LOG_DATA(Esyn_pLog, (/ Erad /) )
  ELSEIF (trim(rad_process) == "ics" .AND. trim(spec) == "electrons") THEN
    CALL ADD_LOG_DATA(Eics_eLog, (/ Erad /) )
  ELSEIF (trim(rad_process) == "ics" .AND. trim(spec) == "ions") THEN
    CALL ADD_LOG_DATA(Eics_pLog, (/ Erad /) )
  ENDIF
END IF

END SUBROUTINE RAD_ENERGY

!***********************************************************************
! Subroutine ANALYSIS_SYNC
! This subroutine computes the total synchrotron radiation spectrum, and
! the radiation angular distribution.
!
! INPUT:
! 
! - mass: Mass of the particle
! - pcl: Particle distribution function
! - Bxg: x-component of B at the nodes at t
! - Byg: y-component of B at the nodes at t
! - Bzg: z-component of B at the nodes at t
! - it: Timestep number
! - spec: Particle species
! - sym: Name of the particles
!
! OUTPUT: nu(Hz),nuFnu(erg/s)
!***********************************************************************

SUBROUTINE ANALYSIS_SYNC(mass,pcl,Bxg,Byg,Bzg,xgp,ygp,it,spec,sym,NPP,&
                         id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                          :: id,COMM,ierr
INTEGER*8                                        :: ip,NPP
DOUBLE PRECISION, DIMENSION(1:NNU)               :: nu
DOUBLE PRECISION, DIMENSION(1:NPHI)              :: pref
DOUBLE PRECISION, DIMENSION(1:NLBA)              :: lref
DOUBLE PRECISION, DIMENSION(1:NNU)               :: nuFnup_iso,nuFnu_iso
DOUBLE PRECISION, DIMENSION(1:NNU)               :: dNdnucp_iso,dNdnuc_iso
DOUBLE PRECISION, DIMENSION(1:NNU,1:NPHI,1:NLBA) :: nuFnup,nuFnu
DOUBLE PRECISION, DIMENSION(1:NNU,1:NPHI,1:NLBA) :: dNdnucp,dNdnuc
DOUBLE PRECISION                                 :: mass,yhalf
DOUBLE PRECISION                                 :: gam,clba,slba,lambda
DOUBLE PRECISION                                 :: cphi,sphi,phi
DOUBLE PRECISION                                 :: kern
DOUBLE PRECISION                                 :: alpha,calpha,Btot,nuc
DOUBLE PRECISION, DIMENSION(1:3)                 :: BField
DOUBLE PRECISION                                 :: Bxi,Byi,Bzi
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)           :: pcl
DOUBLE PRECISION                                 :: x,y,ux,uy,uz,wt
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP)         :: Bxg,Byg,Bzg
DOUBLE PRECISION, DIMENSION(1:NXP)               :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)               :: ygp
DOUBLE PRECISION, DIMENSION(1)                   :: minnu
DOUBLE PRECISION, DIMENSION(1)                   :: minp
DOUBLE PRECISION, DIMENSION(1)                   :: minl
CHARACTER(len=10)                                :: cit,sym,spec,CNLBA
INTEGER                                          :: in,iph,il
INTEGER                                          :: in2,iph2,il2,it
!***********************************************************************

! Convert the integer it into a string cit
 WRITE(cit,'(i10)') it
 WRITE(CNLBA,'(i10)') NLBA

! This step left-justify the string
 cit=adjustl(cit)
 spec=adjustl(spec)
 sym=adjustl(sym)
 CNLBA=adjustl(CNLBA)

! Define the y separating the 2 layers
yhalf=(ymax-ymin)/2.0+ymin

! nu
DO in=1,NNU
nu(in)=1d1**((in-1)*1d0/((NNU-1)*1d0)*(log10(nmax)-log10(nmin))+log10(nmin))
ENDDO

! Latitude PHI in degrees
DO iph=1,NPHI
pref(iph)=(iph-1)*1d0/((NPHI-1)*1d0)*(pmax-pmin)+pmin
ENDDO

! Longitude LAMBDA in degrees
DO il=1,NLBA
lref(il)=(il-1)*1d0/((NLBA-1)*1d0)*(lmax-lmin)+lmin
ENDDO

! Anisotropic photon Flux nuFnu per process
nuFnup=0.0

! Anisotropic photon Flux nuFnu for all processes
nuFnu=0.0

! Isotropic photon Flux nuFnu per process
nuFnup_iso=0.0

! Isotropic photon Flux nuFnu for all processes
nuFnu_iso=0.0

! Anisotropic critical synchrotron frequency spectrum per process
dNdnucp=0.0

! Anisotropic critical synchrotron frequency spectrum for all process
dNdnuc=0.0

! Critical synchrotron frequency spectrum per process
dNdnucp_iso=0.0

! Critical synchrotron frequency spectrum for all process
dNdnuc_iso=0.0

DO ip=1,NPP

  x=pcl(1,ip)
  y=pcl(2,ip)
  ux=pcl(4,ip)
  uy=pcl(5,ip)
  uz=pcl(6,ip)
  wt=pcl(7,ip)

  IF (y.LT.yhalf) THEN

    clba=uz/sqrt(ux*ux+uz*uz)
    slba=ux/sqrt(ux*ux+uz*uz)

      IF (slba.GT.0d0) THEN
      lambda=acos(clba)
      ELSE
      lambda=-acos(clba)
      END IF

    lambda=lambda*180.0/pi

    cphi=sqrt(ux*ux+uz*uz)/sqrt(ux*ux+uy*uy+uz*uz)
    sphi=uy/sqrt(ux*ux+uy*uy+uz*uz)

    phi=asin(sphi)
    phi=phi*180.0/pi

    IF (phi.LT.pmin.OR.phi.GT.pmax.OR.lambda.LT.lmin.OR.lambda.GT.lmax) THEN
      nuFnup=nuFnup
      nuFnup_iso=nuFnup_iso
      dNdnucp=dNdnucp
      dNdnucp_iso=dNdnucp_iso
    ELSE

      minp=minloc(abs(pref-phi))
      minl=minloc(abs(lref-lambda))

      iph2=minp(1)
      il2=minl(1)

      IF (pref(iph2).GT.phi) THEN
        iph2=iph2-1
      ENDIF

      IF (lref(il2).GT.lambda) THEN
        il2=il2-1
      ENDIF

      CALL BILINEAR_FIELDS(xgp,ygp,x,y, BField, Bxg,Byg,Bzg)
      
      Bxi=BField(1)
      Byi=BField(2)
      Bzi=BField(3)

      calpha=(ux*Bxi+uy*Byi+uz*Bzi)
      calpha=calpha/sqrt(ux*ux+uy*uy+uz*uz)
      calpha=calpha/sqrt(Bxi*Bxi+Byi*Byi+Bzi*Bzi)

      alpha=acos(calpha)
      Btot=sqrt(Bxi*Bxi+Byi*Byi+Bzi*Bzi)
      gam=sqrt(1.0+ux*ux+uy*uy+uz*uz)
      nuc=3.0*e*Btot*gam*gam*sin(alpha)/(4.0*pi*mass*c)

      ! Spectrum nuFnu
      DO in=1,NNU

      CALL SYNC(mass,kern,h*nu(in)/evtoerg,Btot,gam,alpha)

      nuFnup(in,iph2,il2)=nuFnup(in,iph2,il2)+wt*(h*nu(in))**2.0*kern/&
      ((sin(pref(iph2+1)*pi/180.0)-sin(pref(iph2)*pi/180.0))*&
      (lref(il2+1)-lref(il2))*pi/180.0)

      nuFnup_iso(in)=nuFnup_iso(in)+wt*(h*nu(in))**2.0*kern

      ENDDO

      ! Distribution of critical synchrotron frequency nuc
      IF (nuc.LT.nmin.OR.nuc.GT.nmax) THEN
      dNdnucp=dNdnucp
      dNdnucp_iso=dNdnucp_iso
      ELSE

      minnu=minloc(abs(nu-nuc))
      in2=minnu(1)

      IF (nu(in2).GT.nuc) THEN
        in2=in2-1
      ENDIF

      dNdnucp(in2,iph2,il2)=dNdnucp(in2,iph2,il2)+wt/&
      ((sin(pref(iph2+1)*pi/180.0)-sin(pref(iph2)*pi/180.0))*&
      (lref(il2+1)-lref(il2))*pi/180.0)/(nu(in2+1)-nu(in2))

      dNdnucp_iso(in2)=dNdnucp_iso(in2)+wt/(nu(in2+1)-nu(in2))

      END IF

    END IF

  ELSE
  nuFnup=nuFnup
  nuFnup_iso=nuFnup_iso
  dNdnucp=dNdnucp
  dNdnucp_iso=dNdnucp_iso
  END IF

ENDDO

!===================================================================

! Copy of the number of particles at lambda=-180deg to lambda=180deg
nuFnup(:,:,NLBA)=nuFnup(:,:,1)
dNdnucp(:,:,NLBA)=dNdnucp(:,:,1)

! Special case: phi=+-90deg
nuFnup(:,NPHI,:)=nuFnup(:,NPHI-1,:)
dNdnucp(:,NPHI,:)=dNdnucp(:,NPHI-1,:)

! Total isotropic spectra
CALL MPI_REDUCE(nuFnup_iso,nuFnu_iso,NNU,MPI_DOUBLE_PRECISION,MPI_SUM,&
                0,COMM,ierr)

! Total anisotropic spectra
CALL MPI_REDUCE(nuFnup,nuFnu,NNU*NPHI*NLBA,MPI_DOUBLE_PRECISION,MPI_SUM,&
                0,COMM,ierr)

! Total nuc isotropic spectrum
CALL MPI_REDUCE(dNdnucp_iso,dNdnuc_iso,NNU,MPI_DOUBLE_PRECISION,MPI_SUM,&
                0,COMM,ierr)

! Total nuc anisotropic spectrum
CALL MPI_REDUCE(dNdnucp,dNdnuc,NNU*NPHI*NLBA,MPI_DOUBLE_PRECISION,MPI_SUM,&
                0,COMM,ierr)

IF (id.EQ.0) THEN

  OPEN(9,FILE="./data/synchrotron_" // trim(spec) // "_" // trim(sym) &
              // trim(cit) // ".dat")
  DO in=1,NNU
    DO iph=1,NPHI
    WRITE(9,"(" // trim(CNLBA) // FMT_DBL // ")") (nuFnu(in,iph,il),il=1,NLBA)
    ENDDO
  ENDDO
  CLOSE(9)

  ! Export the result in a file
  OPEN(9,FILE="./data/synchrotron_iso_" // trim(spec) // "_" // trim(sym) &
              // trim(cit) // ".dat")
  DO in=1,NNU
  WRITE(9,"(1" // FMT_DBL // ")") (nuFnu_iso(in)/(4.0*pi))
  ENDDO
  CLOSE(9)

  OPEN(9,FILE="./data/ecut_" // trim(spec) // "_" // trim(sym) &
              // trim(cit) // ".dat")
  DO in=1,NNU
    DO iph=1,NPHI
    WRITE(9,"(" // trim(CNLBA) // FMT_DBL // ")") (dNdnuc(in,iph,il),il=1,NLBA)
    ENDDO
  ENDDO
  CLOSE(9)

  ! Export the result in a file
  OPEN(9,FILE="./data/ecut_iso_" // trim(spec) // "_" // trim(sym) &
              // trim(cit) // ".dat")
  DO in=1,NNU
  WRITE(9,"(1" // FMT_DBL // ")") (dNdnuc_iso(in)/(4.0*pi))
  ENDDO
  CLOSE(9)

  CALL DUMP_COORDS("nu", nu, nuDumped, COMM)
  CALL DUMP_COORDS("phi", pref, phiDumped, COMM)
  CALL DUMP_COORDS("lambda", lref, lambdaDumped, COMM)

ENDIF

END SUBROUTINE ANALYSIS_SYNC

!***********************************************************************

END MODULE MOD_ANALYSIS
