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

MODULE MOD_INITIAL

USE MOD_INPUT

IMPLICIT NONE

PRIVATE

PUBLIC :: COM_TOPOLOGY ! Initialize the virtual topology of the processes
PUBLIC :: SET_FIELDS ! Generation of the initial fields E and B
PUBLIC :: WEIGHT ! This subroutine weight each super-particles
PUBLIC :: SET_TAG ! This subroutine assigns a tag to each super-particles
PUBLIC :: SET_MAXWELLIAN ! Generation of a relativistic Maxwellian distribution function
PUBLIC :: DIST_MAXWELL ! Integration of the Maxwell distribution function
PUBLIC :: GEN_U ! Generates the particle 4-velocity
PUBLIC :: SET_DRIFT_MAXWELLIAN ! Generation of a drifting Maxwellian distribution function
PUBLIC :: INIT_DRIFT_MAXWELLIAN ! Computes the cumulative distribution functions
PUBLIC :: GEN_UP ! Generates particle parallel momentum
PUBLIC :: GEN_PS ! Generates particle perpendicular momentum
PUBLIC :: INIT_RANDOM_SEED ! Avoid repeating series of random numbers

 CONTAINS

!***********************************************************************
! This subroutine initializes a virtual cartesian topology of the processes.
!
! INPUT:
!
! - ngh: neighbor array
! - coords: coordinates of the process id in the cartesian topology
! - NPROC: Total number of processes
! - id: process rank
! - COMM: communicator
! - ierr: error code
!
! OUTPUT: ngh,coords,COMM
!
!    NW----N-----NE
!    |     |     | 
!    W-----id----E
!    |     |     |
!    SW----S-----SE
!
!***********************************************************************

SUBROUTINE COM_TOPOLOGY(ngh,coords,NPROC,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
INTEGER                             :: NPROC,COMM,direction,step,ierr,id
LOGICAL                             :: reorder
INTEGER, PARAMETER                  :: North=1,East=2,South=3,West=4
INTEGER, PARAMETER                  :: NEast=5,SEast=6,SWest=7,NWest=8
INTEGER, DIMENSION(2)               :: dims,coords
LOGICAL, DIMENSION(2)               :: periods
INTEGER, DIMENSION(8)               :: ngh

!***********************************************************************

! Initialization of the cartesian topology
periods(1)=.TRUE.
periods(2)=.TRUE.
reorder=.FALSE.
dims(1)=NPX
dims(2)=NPY

! Creation of the dimension in each direction
CALL MPI_DIMS_CREATE(NPROC,2,dims,ierr)

! Creation of the topology
CALL MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,reorder,COMM,ierr)

! To obtain the ID number of each process
CALL MPI_COMM_RANK(COMM,id,ierr)

! To obtain the coordinates of the process
CALL MPI_CART_COORDS(COMM,id,2,coords,ierr)

!PRINT*,'*******************************************************'
!PRINT*,'Processus: ',id,', Coordinates: ',coords(1),coords(2)
!PRINT*,'*******************************************************'

! Search neighbors

! West-East
direction=0
step=1

CALL MPI_CART_SHIFT(COMM,direction,step,ngh(West),ngh(East),ierr)

! North-South
direction=1
step=1

CALL MPI_CART_SHIFT(COMM,direction,step,ngh(South),ngh(North),ierr)

! Diagonal neighbors

! South-West
CALL MPI_SENDRECV(ngh(South),1,MPI_INTEGER,ngh(East),1,ngh(SWest),1,&
MPI_INTEGER,ngh(West),1,COMM,stat,ierr)

! North-West
CALL MPI_SENDRECV(ngh(North),1,MPI_INTEGER,ngh(East),2,ngh(NWest),1,&
MPI_INTEGER,ngh(West),2,COMM,stat,ierr)

! South-East
CALL MPI_SENDRECV(ngh(South),1,MPI_INTEGER,ngh(West),3,ngh(SEast),1,&
MPI_INTEGER,ngh(East),3,COMM,stat,ierr)

! North-East
CALL MPI_SENDRECV(ngh(North),1,MPI_INTEGER,ngh(West),4,ngh(NEast),1,&
MPI_INTEGER,ngh(East),4,COMM,stat,ierr)

!PRINT*,'*******************************************************'
!PRINT*,'Processus: ',id,', my neighbor W is  : ',ngh(West)
!PRINT*,'Processus: ',id,', my neighbor E is  : ',ngh(East)
!PRINT*,'Processus: ',id,', my neighbor S is  : ',ngh(South)
!PRINT*,'Processus: ',id,', my neighbor N is  : ',ngh(North)
!PRINT*,'Processus: ',id,', my neighbor NE is : ',ngh(NEast)
!PRINT*,'Processus: ',id,', my neighbor SE is : ',ngh(SEast)
!PRINT*,'Processus: ',id,', my neighbor SW is : ',ngh(SWest)
!PRINT*,'Processus: ',id,', my neighbor NW is : ',ngh(NWest)
!PRINT*,'*******************************************************'

END SUBROUTINE COM_TOPOLOGY

!***********************************************************************
! Subroutine SET_FIELDS
! This subroutine generates the initial electromagnetic fields E and B and 
! initial current density J in the Yee lattice.
!
! INPUT:
! 
! - Bx: x-component of B
! - By: y-component of B
! - Bz: z-component of B
! - Ex: x-component of E
! - Ey: y-component of E
! - Ez: z-component of E
! - Jx: x-component of J
! - Jy: y-component of J
! - Jz: z-component of J
!
! OUTPUT: Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz
!***********************************************************************

SUBROUTINE SET_FIELDS(nd0,delta,de,B0,grad,sigma,Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,&
                      xgp,ygp,xyeep,yyeep)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bx,By,Bz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Jx,Jy,Jz

! Nodal grid in each domain
DOUBLE PRECISION, DIMENSION(1:NXP) :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP) :: ygp

! Yee grid in each domain
DOUBLE PRECISION, DIMENSION(1:NXP) :: xyeep
DOUBLE PRECISION, DIMENSION(1:NYP) :: yyeep

! Other plasma paramters
DOUBLE PRECISION :: delta,de,B0,J0,nd0,gamd,Az,sigma,grad
DOUBLE PRECISION :: y12,y14,y34,x14,x34,xm,ym

! Loop indexes
INTEGER :: ix,iy,id

!***********************************************************************
! Initial plasma parameters according to the Harris equilibrium

! Bulk Lorentz factor drifting particles
gamd=1.0/sqrt(1.0-betad*betad)

! Upstream initial reconnecting magnetic field
B0=thde*me*c*c/(e*rhoc)

! Initial layer thickness
delta=2.0*thde*me*c*c/(betad*gamd*e*B0)

! Electron density of the drifting particles in the layer.
nd0=(thde*me*c*c)/(4.0*pi*e*e*betad*betad*gamd*delta*delta)

! Magnetization parameter
sigma=B0*B0/(4.0*pi*me*c*c*2.0*nd0*density_ratio)

! Electron skin depth.
de=sqrt(thde*me*c*c/(4.0*pi*nd0*e*e))

! Initial current in the +/-z-direction
J0=(c*thde*me*c*c)/(2.0*pi*e*betad*gamd*delta*delta)

! Electron maximum radiation reaction limit energy
grad=sqrt(3.0*e/(2.0*e**4.0/(me**2.0*c**4.0)*B0))

DO ix=1,NXP
DO iy=1,NYP

! Initial magnetic field components

y12=(ymax-ymin)/2.0+ymin
y14=(ymax-ymin)/4.0+ymin
y34=3.0*(ymax-ymin)/4.0+ymin
x14=(xmax-xmin)/4.0+xmin
x34=3.0*(xmax-xmin)/4.0+xmin
xm=xmax-xmin
ym=ymax-ymin

IF (ygp(iy).LT.y12) THEN

Bx(ix,iy)=-B0*tanh((yyeep(iy)-y14)/delta)*(1.0+perturb_amp*cos(2.0*pi*(xgp(ix)-&
x14)/xm)*(cos(2.0*pi*(yyeep(iy)-y14)/ym))**2.0)-2.0*(B0*delta*&
log(cosh(y14/delta))-delta*B0*log(cosh((yyeep(iy)-y14)/delta)))*&
perturb_amp*cos(2.0*pi*(xgp(ix)-x14)/xm)*2.0*pi/ym*cos(2.0*pi*(yyeep(iy)-&
y14)/ym)*sin(2.0*pi*(yyeep(iy)-y14)/ym)

By(ix,iy)=(B0*delta*log(cosh(y14/delta))-delta*B0*log(cosh((ygp(iy)-&
y14)/delta)))*perturb_amp*2.0*pi/xm*sin(2.0*pi*(xyeep(ix)-&
x14)/xm)*(cos(2.0*pi*(ygp(iy)-y14)/ym))**2.0

ELSE

Bx(ix,iy)=B0*tanh((yyeep(iy)-y34)/delta)*(1.0+perturb_amp*cos(2.0*pi*(xgp(ix)-&
x34)/xm)*(cos(2.0*pi*(yyeep(iy)-y34)/ym))**2.0)+2.0*(B0*delta*&
log(cosh(y14/delta))-delta*B0*log(cosh((yyeep(iy)-y34)/delta)))*&
perturb_amp*cos(2.0*pi*(xgp(ix)-x34)/xm)*2.0*pi/ym*cos(2.0*pi*(yyeep(iy)-&
y34)/ym)*sin(2.0*pi*(yyeep(iy)-y34)/ym)

By(ix,iy)=-(B0*delta*log(cosh(y14/delta))-delta*B0*log(cosh((ygp(iy)-&
y34)/delta)))*perturb_amp*2.0*pi/xm*sin(2.0*pi*(xyeep(ix)-&
x34)/xm)*(cos(2.0*pi*(ygp(iy)-y34)/ym))**2.0

END IF

Bz(ix,iy)=B0*guide_field

! Initial electric field components
Ex(ix,iy)=0.0
Ey(ix,iy)=0.0
Ez(ix,iy)=0.0

! Initial current density components
Jx(ix,iy)=0.0
Jy(ix,iy)=0.0

IF (ygp(iy).LT.y12) THEN
! z-component potential vecteur
Az=(thde*me*c*c)/(e*betad*gamd)*log((cosh((ygp(iy)-y14)/delta))**(-2.0))
Jz(ix,iy)=J0*exp(e*betad*gamd*Az/(thde*me*c*c))
ELSE
! z-component potential vecteur
Az=(thde*me*c*c)/(e*betad*gamd)*log((cosh((ygp(iy)-y34)/delta))**(-2.0))
Jz(ix,iy)=-1.0*J0*exp(e*betad*gamd*Az/(thde*me*c*c))
END IF

ENDDO
ENDDO

END SUBROUTINE SET_FIELDS

!***********************************************************************
! Subroutine WEIGHT
! This subroutine weight each super-particles
!
! INPUT:
! 
! - pcl0: Initial particle distribution function
! - delta: Initial layer thickness
! - NPP: Number of initial particles per process
!
! OUTPUT: wght=pcl0(7,:)
!***********************************************************************

SUBROUTINE WEIGHT(pcl0,delta,NPP)

IMPLICIT NONE

INTEGER*8                              :: ip,NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP) :: pcl0
DOUBLE PRECISION                       :: delta,y,y12,y14,y34

!***********************************************************************

y12=(ymax-ymin)/2.0+ymin
y14=(ymax-ymin)/4.0+ymin
y34=3.0*(ymax-ymin)/4.0+ymin

! From the Harris equilibirum

DO ip=1,NPP

y=pcl0(2,ip)

IF (y.LT.y12) THEN
pcl0(7,ip)=1.0*(cosh((y-y14)/delta))**(-2.0)
ELSE
pcl0(7,ip)=1.0*(cosh((y-y34)/delta))**(-2.0)
END IF

ENDDO

END SUBROUTINE WEIGHT

!***********************************************************************
! Subroutine SET_TAG
! This subroutine assigns a tag to each super-particles
!
! INPUT:
! 
! - tag: tag of the super-particle
! - id: process rank
! - NPP: Number of initial particles per process
!
! OUTPUT: tag
!***********************************************************************

SUBROUTINE SET_TAG(tag,id,NPP)

IMPLICIT NONE

INTEGER*8                   :: ip,NPP
INTEGER*8, DIMENSION(1:NPP) :: tag
INTEGER                     :: id

!***********************************************************************

tag=0

DO ip=1,NPP
tag(ip)=id*NPP+ip ! tag must be > 0 so -tag means a tracked particle
ENDDO

END SUBROUTINE SET_TAG

!***********************************************************************
! Subroutine SET_MAXWELLIAN
! This subroutine generates the initial distribution function of the
! particles with a Maxwellian distribution in energy.
!
! INPUT: 
! - pcl0: Initial distribution function of the particles
! - theta: Temperature of plasma in mc^2/k units
! - xminp,yminp: Lower spatial boundary for each domain
! - NPP: Number of initial particles per process
!
! OUTPUT: Particle distribution function
!***********************************************************************

SUBROUTINE SET_MAXWELLIAN(pcl0,theta,xminp,yminp,NPP)

IMPLICIT NONE

! Temperature in mc^2/k units
DOUBLE PRECISION :: theta,umint,umaxt

! Cosine of the polar angle between the z-axis and the velocity
DOUBLE PRECISION, PARAMETER           :: cmin=-1.0,cmax=1.0

! Azimuth between x-axis and the projected velocity in the xy-plane
DOUBLE PRECISION, PARAMETER           :: pmin=0.0,pmax=2.0*pi

! Number of elements for integration distribution function
INTEGER, PARAMETER :: ND=500

! Number of initial particles per process
INTEGER*8 :: NPP

! Lower spatial boundary for each domain
DOUBLE PRECISION :: xminp,yminp

DOUBLE PRECISION, DIMENSION(1:ND)           :: ud,gFu
DOUBLE PRECISION, DIMENSION(1:PPC)          :: x0c,y0c,cth0,phi0
DOUBLE PRECISION, DIMENSION(1:PPC)          :: u0,ux0c,uy0c,uz0c
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)      :: pcl0

! Loop indexes
INTEGER   :: ic,id,ix,iy
INTEGER*8 :: ip

!***********************************************************************

CALL INIT_RANDOM_SEED()

ip=1

pcl0=0.0
ud=0.0

! Cumulative distribution function
gFu=0.0

umint=1d-2*sqrt(theta)
umaxt=1d2*sqrt(theta)

DO id=1,ND
ud(id)=1d1**((id-1)*1d0/((ND-1)*1d0)*(log10(umaxt)-log10(umint))+log10(umint))
ENDDO

! Numerical calculation of the cumulative distribution function
CALL DIST_MAXWELL(ud,gFu,theta,ND)

DO ix=1,NCXP
DO iy=1,NCYP

!***********************************************************************
! Definition of the initial position x0 as a uniform random value

x0c=0.0

CALL RANDOM_NUMBER(x0c)

x0c=xminp+(ix-1)*dx+x0c*dx

! Definition of the initial position y0 as a uniform random value

y0c=0.0

CALL RANDOM_NUMBER(y0c)

y0c=yminp+(iy-1)*dy+y0c*dy

!***********************************************************************
! Definition of the initial particle 4-velocity u0

u0=0.0

CALL GEN_U(u0,ud,gFu,ND)

!***********************************************************************
! Definition of the initial phi0 as a uniform random value

phi0=0.0

CALL RANDOM_NUMBER(phi0)

 phi0=phi0*(pmax-pmin)+pmin

! Definition of the initial cth0=cos(theta0) as a uniform random value

 cth0=0.0

CALL RANDOM_NUMBER(cth0)

 cth0=cth0*(cmax-cmin)+cmin

! Initial 4-velocity components
ux0c=u0*sqrt(1.0-cth0*cth0)*cos(phi0)
uy0c=u0*sqrt(1.0-cth0*cth0)*sin(phi0)
uz0c=u0*cth0

DO ic=1,PPC
pcl0(1,ip)=x0c(ic)
pcl0(2,ip)=y0c(ic)
pcl0(4,ip)=ux0c(ic)
pcl0(5,ip)=uy0c(ic)
pcl0(6,ip)=uz0c(ic)
ip=ip+1
ENDDO

ENDDO
ENDDO

END SUBROUTINE SET_MAXWELLIAN

!***********************************************************************
! Subroutine DIST_MAXWELL
! This subroutine generates a random distribution of particle Lorentz 
! factor gamma following a relativistic Maxwellian distribution.
!
! INPUT: 
! - ud: reference 4-velocity array
! - gFu: pre-calculated integrated distribution function
! - theta: temperature of plasma in mc^2/k units
! - ND: number of elements in the array ud and gFu
!
! OUTPUT: Integrated distribution function gFu
!***********************************************************************

SUBROUTINE DIST_MAXWELL(ud,gFu,theta,ND)

IMPLICIT NONE

! Input parameter
DOUBLE PRECISION :: theta
INTEGER          :: ND
DOUBLE PRECISION, DIMENSION(1:ND) :: ud,gFu,fu,udp,fup
DOUBLE PRECISION :: norm,hu,upmax,hup

! Loop indexes
INTEGER :: id1,id2

!***********************************************************************

! Relativistic Maxwellian distribution
fu=ud*ud*exp(-(sqrt(1.0+ud*ud)-1.0)/theta)*ud

norm=0.0

DO id1=1,ND-1
hu=(log(ud(id1+1))-log(ud(id1)))
norm=norm+hu/2.0*(fu(id1+1)+fu(id1))
ENDDO

DO id1=1,ND

upmax=ud(id1)

  DO id2=1,ND
  udp(id2)=1d1**((id2-1)*1d0/((ND-1)*1d0)*(log10(upmax)-&
           log10(minval(ud)))+log10(minval(ud)))
  ENDDO

  ! Relativistic Maxwellian distribution
  fup=udp*udp*exp(-(sqrt(1.0+udp*udp)-1.0)/theta)*udp

  DO id2=1,ND-1
  hup=(log(udp(id2+1))-log(udp(id2)))
  gFu(id1)=gFu(id1)+hup/2.0*(fup(id2+1)+fup(id2))
  ENDDO

ENDDO

gFu=gFu/norm

END SUBROUTINE DIST_MAXWELL

!***********************************************************************
! Subroutine GEN_U
! This subroutine generates a random distribution of particle 4-velocity
!
! INPUT: 
! - u0: 4-velocity array
! - ud: reference 4-velocity array
! - uFg: pre-calculated integrated distribution function
! - theta: temperature of plasma in mc^2/k units
! - ND: number of elements in the array gd and gFg
!
! OUTPUT: distribution of u0
!***********************************************************************

SUBROUTINE GEN_U(u0,ud,gFu,ND)

IMPLICIT NONE

INTEGER :: ND
DOUBLE PRECISION, DIMENSION(1:ND)  :: ud,gFu
DOUBLE PRECISION, DIMENSION(1:PPC) :: u0,Ru
DOUBLE PRECISION                   :: gFu1,gFu2,u1,u2
INTEGER, DIMENSION(1)              :: minu

! Loop indexes
INTEGER :: ic,iu

!***********************************************************************

Ru=0.0

CALL RANDOM_NUMBER(Ru)
Ru=Ru*0.9999

DO ic=1,PPC

minu=minloc(abs(gFu-Ru(ic)))
iu=minu(1)

IF (iu.EQ.ND) THEN
iu=iu-1
END IF

gFu1=gFu(iu)
gFu2=gFu(iu+1)

u1=ud(iu)
u2=ud(iu+1)

u0(ic)=exp((Ru(ic)*(log(u2)-log(u1))-log(u2)*gFu1+log(u1)*gFu2)/(gFu2-gFu1))

ENDDO

END SUBROUTINE GEN_U

!***********************************************************************
! Subroutine SET_DRIFT_MAXWELLIAN
! This subroutine generates the initial distribution function of the
! particles with a drifting Maxwellian distribution in energy.
!
! INPUT: 
! - q: electric charge sign
! - pcl0: Initial distribution function of the particles
! - theta: Temperature of plasma in mc^2/k units
! - up,gFp,ps,gFs,ND: Particle distribution generator quantities
! - xminp,yminp: Lower spatial boundary for each domain
! - NPP: Number of initial particles per process
!
! OUTPUT: Particle distribution function
!***********************************************************************

SUBROUTINE SET_DRIFT_MAXWELLIAN(q,pcl0,theta,up,gFp,ps,gFs,xminp,yminp,ND,NPP)

IMPLICIT NONE

! Temperature in mc^2/k units
DOUBLE PRECISION :: theta,q

! Angle defined in the plane perpendicular to the direction of the bulk motion
DOUBLE PRECISION, PARAMETER           :: pmin=0.0,pmax=2.0*pi

! Number of elements in the pre-calculated cumulative distribution functions
INTEGER :: ND

! Number of initial particles per process
INTEGER*8 :: NPP

! Lower spatial boundary for each domain
DOUBLE PRECISION :: xminp,yminp,gamd,y12

DOUBLE PRECISION, DIMENSION(1:ND)           :: up,gFp,ps
DOUBLE PRECISION, DIMENSION(1:ND,1:ND)      :: gFs
DOUBLE PRECISION, DIMENSION(1:PPC)          :: x0c,y0c,phi0
DOUBLE PRECISION, DIMENSION(1:PPC)          :: up0,ps0,uperp0,ux0c,uy0c,uz0c
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)      :: pcl0

! Loop indexes
INTEGER   :: ic,ix,iy
INTEGER*8 :: ip

!***********************************************************************

! Drift Lorentz factor
gamd=1.0/sqrt(1.0-betad*betad)

! Mid-box
y12=(ymax-ymin)/2.0+ymin

! Initialization of the seed for the generation of random variables
CALL INIT_RANDOM_SEED()

ip=1

pcl0=0.0

DO ix=1,NCXP
DO iy=1,NCYP

!***********************************************************************
! Definition of the initial position x0 as a uniform random value

x0c=0.0

CALL RANDOM_NUMBER(x0c)

x0c=xminp+(ix-1)*dx+x0c*dx

! Definition of the initial position y0 as a uniform random value

y0c=0.0

CALL RANDOM_NUMBER(y0c)

y0c=yminp+(iy-1)*dy+y0c*dy

!***********************************************************************
! Definition of the initial momentum parallel to the drift velocity

up0=0.0

CALL GEN_UP(up0,up,gFp,ND)

!***********************************************************************
! Definition of the initial perpendicular momentum

ps0=0.0

CALL GEN_PS(ps0,ps,up0,up,gFs,ND)

uperp0=ps0*sqrt(1.0+up0*up0)
  
!***********************************************************************
! Definition random angle in the xy-plane

phi0=0.0

CALL RANDOM_NUMBER(phi0)

 phi0=phi0*(pmax-pmin)+pmin

! Initial 4-velocity
ux0c=uperp0*cos(phi0)
uy0c=uperp0*sin(phi0)
uz0c=up0

DO ic=1,PPC
pcl0(1,ip)=x0c(ic)
pcl0(2,ip)=y0c(ic)
pcl0(4,ip)=ux0c(ic)
pcl0(5,ip)=uy0c(ic)

  IF (y0c(ic).LT.y12) THEN
  pcl0(6,ip)=q*uz0c(ic)
  ELSE
  pcl0(6,ip)=-q*uz0c(ic)
  END IF

ip=ip+1
ENDDO

ENDDO
ENDDO

END SUBROUTINE SET_DRIFT_MAXWELLIAN

!***********************************************************************
! Subroutine INIT_DRIFT_MAXWELLIAN
! This subroutine computes the cumulative distribution functions for
! generating a drifting relativistic maxwellian distribution of particles.
! The distributions are taken from Swisdak (2013), PoP.
!
! INPUT:
! - up: parallel momentum array
! - gFp: cumulative distribution function for the parallel momentum
! - ps=uperp/sqrt(1+up*up)
! - gFs: conditional distribution function for ps knowing up
! - theta: temperature of plasma in mc^2/k units
! - ND: number of elements for all arrays
! - COMM: mpi comm
!
! OUTPUT: Integrated distribution function gFp and gFs
!***********************************************************************

SUBROUTINE INIT_DRIFT_MAXWELLIAN(up,gFp,ps,gFs,theta,ND, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! Input parameter
DOUBLE PRECISION :: theta
INTEGER          :: ND, COMM
DOUBLE PRECISION, DIMENSION(1:ND)      :: up,gFp,fg1,up2,fg2,g1,g2
DOUBLE PRECISION, DIMENSION(1:ND)      :: ps,gs1,gs2,ps2,fs1,fs2
DOUBLE PRECISION, DIMENSION(1:ND,1:ND) :: gFs, gFs2
DOUBLE PRECISION :: norm1,norm2,hp,hp2,gamd,pu,umint,hs,hs2

INTEGER :: rank, numRanks, mpiErr
! Loop indexes
INTEGER :: i1,i2,i3

!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)

! Drift Lorentz factor
gamd=1.0/sqrt(1.0-betad*betad)

! Drift 4-velocity
pu=gamd*betad

g1=sqrt(1.0+up*up)

! f(p||) see Eq.(12) Swisdak (2013)
fg1=(1.0+gamd*g1/theta)*exp(-(up-pu)**2.0/(g1*gamd+up*pu+1d0)/theta)

norm1=0.0

DO i1=1,ND-1
hp=up(i1+1)-up(i1)
norm1=norm1+hp/2.0*(fg1(i1+1)+fg1(i1))
ENDDO

umint=minval(up)

!***********************************************************************
! CALCULATION OF THE PARALLEL MOMENTUM CUMULATIVE DISTRIBUTION gFp
!***********************************************************************

gFp=0.0

DO i1=1,ND

  DO i2=1,ND
  up2(i2)=(i2-1)*1d0/((ND-1)*1d0)*(up(i1)-umint)+umint
  ENDDO

  g2=sqrt(1.0+up2*up2)
  
  ! f(p||) see Eq.(12) Swisdak (2013)
  fg2=(1.0+gamd*g2/theta)*exp(-(up2-pu)**2.0/(g2*gamd+up2*pu+1d0)/theta)
  
  DO i2=1,ND-1
  hp2=up2(i2+1)-up2(i2)
  gFp(i1)=gFp(i1)+hp2/2.0*(fg2(i2+1)+fg2(i2))
  ENDDO

ENDDO

gFp=gFp/norm1

!***********************************************************************
! CALCULATION OF THE CONDITIONAL PERPENDICULAR MOMENTUM CUMULATIVE
! DISTRIBUTION gFs
!***********************************************************************

gs1=sqrt(1.0+ps*ps)
gFs2=0.0

! GRW: this triple loop is time-consuming; it decreases parallel 
! efficiency when run on many processors, but I believe each outer-loop
! is independent, so farm it out.
! N.B. For ND=1000, on verus, this subroutine took about 97 s, in serial.
!   After parallelizing, it took about 3 s on 48 procs.

DO i1=1,ND

  IF (MOD(i1-1, numRanks) == rank) THEN !{
  
  ! f(ps|up) see Eq.(15) Swisdak (2013)
  fs1=ps*exp(-((up(i1)-pu)**2.0+g1(i1)*g1(i1)*gamd*gamd*ps*ps)/&
              (g1(i1)*gamd*gs1+up(i1)*pu+1.0)/theta)
  
  norm2=0.0
  
  DO i2=1,ND-1
    hs=ps(i2+1)-ps(i2)
    norm2=norm2+hs/2.0*(fs1(i2+1)+fs1(i2))
  ENDDO
  
  !***********************************************************************

  DO i2=1,ND

     DO i3=1,ND
     ps2(i3)=(i3-1)*1d0/((ND-1)*1d0)*(ps(i2)-0.0)+0.0
     ENDDO
     
     gs2=sqrt(1.0+ps2*ps2)
     
     ! f(ps|p||) see Eq.(15) Swisdak (2013)
     fs2=ps2*exp(-((up(i1)-pu)**2.0+g1(i1)*g1(i1)*gamd*gamd*ps2*ps2)/&
                  (g1(i1)*gamd*gs2+up(i1)*pu+1.0)/theta)
     
     DO i3=1,ND-1
     hs2=ps2(i3+1)-ps2(i3)
     gFs2(i2,i1)=gFs2(i2,i1)+hs2/2.0*(fs2(i3+1)+fs2(i3))
     ENDDO

     gFs2(i2,i1)=gFs2(i2,i1)/norm2     

  ENDDO

  ENDIF !}
     
ENDDO

IF (numRanks > 1) THEN
  CALL MPI_ALLREDUCE(gFs2, gFs, SIZE(gFs), MPI_DOUBLE_PRECISION, MPI_SUM,&
    COMM, mpiErr)
ELSE
  gFs = gFs2
ENDIF

END SUBROUTINE INIT_DRIFT_MAXWELLIAN

!***********************************************************************
! Subroutine GEN_UP
! This subroutine generates a random distribution of particle momenta 
! parallel to the drift velocity.
!
! INPUT: 
! - up0: Particle momentum parallel to the drift velocity
! - up: reference momentum array
! - gFp: pre-calculated integrated distribution function
! - ND: number of elements in the array ug and gFg
!
! OUTPUT: distribution of up0
!***********************************************************************

SUBROUTINE GEN_UP(up0,up,gFp,ND)

IMPLICIT NONE

INTEGER :: ND
DOUBLE PRECISION, DIMENSION(1:ND)  :: up,gFp
DOUBLE PRECISION, DIMENSION(1:PPC) :: up0,Rp
DOUBLE PRECISION                   :: gFp1,gFp2,u1,u2
INTEGER, DIMENSION(1)              :: minu

! Loop indexes
INTEGER :: ic,iu

!***********************************************************************

Rp=0.0

CALL RANDOM_NUMBER(Rp)
Rp=Rp*0.999

DO ic=1,PPC

minu=minloc(abs(gFp-Rp(ic)))
iu=minu(1)

IF (iu.EQ.ND) THEN
iu=iu-1
END IF

gFp1=gFp(iu)
gFp2=gFp(iu+1)

u1=up(iu)
u2=up(iu+1)

up0(ic)=(Rp(ic)*(u2-u1)-u2*gFp1+u1*gFp2)/(gFp2-gFp1)

ENDDO

END SUBROUTINE GEN_UP

!***********************************************************************
! Subroutine GEN_PS
! This subroutine generates a random distribution of particle momenta ps.
!
! INPUT: 
! - ps0: Particle momentum ps=uperp/sqrt(1+up*up)
! - ps: reference momentum array
! - up0: Particle momentum parallel to the drift velocity
! - up: reference momentum array
! - gFs: pre-calculated integrated distribution function
! - ND: number of elements in the array ug and gFg
!
! OUTPUT: distribution of ps0
!***********************************************************************

SUBROUTINE GEN_PS(ps0,ps,up0,up,gFs,ND)

IMPLICIT NONE

INTEGER :: ND
DOUBLE PRECISION, DIMENSION(1:ND)      :: ps,up
DOUBLE PRECISION, DIMENSION(1:ND,1:ND) :: gFs
DOUBLE PRECISION, DIMENSION(1:PPC)     :: ps0,up0,Rp
DOUBLE PRECISION                       :: F11,F12,F21,F22,fp,fq
INTEGER, DIMENSION(1)                  :: minu,minf

! Loop indexes
INTEGER :: ic,i1,i2

!***********************************************************************

Rp=0.0

CALL RANDOM_NUMBER(Rp)
Rp=Rp*0.999

DO ic=1,PPC

minu=minloc(abs(up-up0(ic)))
i1=minu(1)

minf=minloc(abs(gFs(:,i1)-Rp(ic)))
i2=minf(1)

IF (i1.EQ.ND) THEN
i1=i1-1
END IF

IF (i2.EQ.ND) THEN
i2=i2-1
END IF

F11=gFs(i2,i1)
F12=gFs(i2,i1+1)
F21=gFs(i2+1,i1)
F22=gFs(i2+1,i1+1)

fq=(up0(ic)-up(i1))/(up(i1+1)-up(i1))
fp=(Rp(ic)-(1.0-fq)*F11-fq*F12)/((1.0-fq)*(F21-F11)+fq*(F22-F12))

ps0(ic)=fp*(ps(i2+1)-ps(i2))+ps(i2)

ENDDO

END SUBROUTINE GEN_PS

!***********************************************************************
! init_random_seed() subroutine enables to avoid the repeating series of
! number given by random_number.
! Reference: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!***********************************************************************

SUBROUTINE INIT_RANDOM_SEED()

INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
!***********************************************************************

IF (RANDOMIZE.EQV..TRUE.) THEN
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
        
CALL SYSTEM_CLOCK(COUNT=clock)
         
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
          
DEALLOCATE(seed)
ENDIF

END SUBROUTINE INIT_RANDOM_SEED

!***********************************************************************

END MODULE MOD_INITIAL
