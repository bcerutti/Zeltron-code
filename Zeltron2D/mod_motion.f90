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

MODULE MOD_MOTION

USE MOD_INPUT
USE MOD_INTERP

IMPLICIT NONE

PRIVATE

PUBLIC :: INITIAL_PUSH ! Computes the 4-veocity at t=0-dt/2
PUBLIC :: BORIS_PUSH ! Boris push on the particle 4-velocity
PUBLIC :: PUSH_PARTICLES ! Update the positions of the particles
PUBLIC :: COUNT_ESCAPE ! Counts the number of particles escaping the domain
PUBLIC :: BOUNDARIES_PARTICLES ! Applies boundary conditions to particles
PUBLIC :: COM_PARTICLES ! Exchange of part. at the boundaries btw processes

 CONTAINS

!***********************************************************************
! Subroutine INITIAL_PUSH
! This subroutine computes the initial half push of the particles' 
! 4-velocity of the particles, using the Boris algorithm. 
! INCLUDES THE SYNCHROTRON AND THE INVERSE COMPTION RADIATION
! REACTION FORCES, following the procedure by Tamburini et al., New Journal of
! Physics, 12 (2010) 123005.
!
! INPUT: 
! - q: sign electric charge
! - mass: Mass of the particle
! - pcl: Particle distribution function at t=0
! - Bxg: x-component of B at the nodes at t=0
! - Byg: y-component of B at the nodes at t=0
! - Bzg: z-component of B at the nodes at t=0
! - Exg: x-component of E at the nodes at t=0
! - Eyg: y-component of E at the nodes at t=0
! - Ezg: z-component of E at the nodes at t=0
! - Uph: External photon energy density
! - xgp,ygp: spatial grid
! - NPP: Number of particles per process
!
! OUTPUT: ux,uy,uz at t=0-dt/2
!***********************************************************************

SUBROUTINE INITIAL_PUSH(q,mass,pcl,Bxg,Byg,Bzg,Exg,Eyg,Ezg,Uph,xgp,ygp,NPP)

IMPLICIT NONE

INTEGER*8 :: ip,NPP
DOUBLE PRECISION :: q,mass,gam,Bxi,Byi,Bzi,Exi,Eyi,Ezi,Psyn,Pics,Uph
DOUBLE PRECISION :: uxp,uyp,uzp,ux0,uy0,uz0,uxL,uyL,uzL
DOUBLE PRECISION :: sx,sy,sz,tx,ty,tz
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
DOUBLE PRECISION                         :: x,y,z,ux,uy,uz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bxg,Byg,Bzg
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Exg,Eyg,Ezg
DOUBLE PRECISION, DIMENSION(1:6)         :: Fields

!***********************************************************************

! Switch ON/OFF of the radiation reaction force:

!***********************************************************************
! WITH radiation reaction force
!***********************************************************************

IF (RAD_FORCE.EQV..TRUE.) THEN

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)
z=pcl(3,ip)
ux=pcl(4,ip)
uy=pcl(5,ip)
uz=pcl(6,ip)

!***********************************************************************
! Interpolation of E and B fields from the nodes to the particles at time t
CALL BILINEAR_FIELDS(xgp,ygp,x,y, Fields, Bxg,Byg,Bzg,Exg,Eyg,Ezg)

Bxi=Fields(1)
Byi=Fields(2)
Bzi=Fields(3)

Exi=Fields(4)
Eyi=Fields(5)
Ezi=Fields(6)

!***********************************************************************
! Computation of u(-1/2)

! gamma at t=0
gam=sqrt(1.0+ux*ux+uy*uy+uz*uz)

! Calculation of t-vector
tx=q*e*Bxi*dt/(4.0*gam*mass*c)
ty=q*e*Byi*dt/(4.0*gam*mass*c)
tz=q*e*Bzi*dt/(4.0*gam*mass*c)

ux0=-ux+uy*tz-uz*ty
uy0=-uy+uz*tx-ux*tz
uz0=-uz+ux*ty-uy*tx

! Calculation of u+
sx=2.0*tx/(1.0+tx*tx+ty*ty+tz*tz)
sy=2.0*ty/(1.0+tx*tx+ty*ty+tz*tz)
sz=2.0*tz/(1.0+tx*tx+ty*ty+tz*tz)

uxp=ux+uy0*sz-uz0*sy
uyp=uy+uz0*sx-ux0*sz
uzp=uz+ux0*sy-uy0*sx

! uL at t=-dt/2
uxL=uxp-q*e*Exi*dt/(2.0*mass*c)
uyL=uyp-q*e*Eyi*dt/(2.0*mass*c)
uzL=uzp-q*e*Ezi*dt/(2.0*mass*c)

! Total synchrotron radiative power losses at t=0
Psyn=(2.0/3.0)*e**4.0/(mass**2.0*c**4.0)*c*((gam*Exi+uy*Bzi-uz*Byi)**2.0+&
     (gam*Eyi+uz*Bxi-ux*Bzi)**2.0+(gam*Ezi+ux*Byi-uy*Bxi)**2.0-&
     (ux*Exi+uy*Eyi+uz*Ezi)**2.0)

! Radiative inverse Compton energy losses (Thomson regime) at t=0
Pics=(32.0/9.0)*pi*e**4.0/(mass**2.0*c**4.0)*c*Uph*(ux*ux+uy*uy+uz*uz)

! u at t=-dt/2
ux=uxL+dt*(Psyn+Pics)/(2.0*mass*c*c)*ux/gam
uy=uyL+dt*(Psyn+Pics)/(2.0*mass*c*c)*uy/gam
uz=uzL+dt*(Psyn+Pics)/(2.0*mass*c*c)*uz/gam

!***********************************************************************

pcl(4,ip)=ux
pcl(5,ip)=uy
pcl(6,ip)=uz

ENDDO

!***********************************************************************
! WITHOUT radiation reaction force
!***********************************************************************

ELSE

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)
z=pcl(3,ip)
ux=pcl(4,ip)
uy=pcl(5,ip)
uz=pcl(6,ip)

!***********************************************************************
! Interpolation of E and B fields from the nodes to the particles at time t
CALL BILINEAR_FIELDS(xgp,ygp,x,y, Fields, Bxg,Byg,Bzg,Exg,Eyg,Ezg)

Bxi=Fields(1)
Byi=Fields(2)
Bzi=Fields(3)

Exi=Fields(4)
Eyi=Fields(5)
Ezi=Fields(6)

!***********************************************************************
! Computation of u(-1/2)

! gamma at t=0
gam=sqrt(1.0+ux*ux+uy*uy+uz*uz)

! Calculation of t-vector
tx=q*e*Bxi*dt/(4.0*gam*mass*c)
ty=q*e*Byi*dt/(4.0*gam*mass*c)
tz=q*e*Bzi*dt/(4.0*gam*mass*c)

ux0=-ux+uy*tz-uz*ty
uy0=-uy+uz*tx-ux*tz
uz0=-uz+ux*ty-uy*tx

! Calculation of u+
sx=2.0*tx/(1.0+tx*tx+ty*ty+tz*tz)
sy=2.0*ty/(1.0+tx*tx+ty*ty+tz*tz)
sz=2.0*tz/(1.0+tx*tx+ty*ty+tz*tz)

uxp=ux+uy0*sz-uz0*sy
uyp=uy+uz0*sx-ux0*sz
uzp=uz+ux0*sy-uy0*sx

! u at t=-dt/2
ux=uxp-q*e*Exi*dt/(2.0*mass*c)
uy=uyp-q*e*Eyi*dt/(2.0*mass*c)
uz=uzp-q*e*Ezi*dt/(2.0*mass*c)

!***********************************************************************

pcl(4,ip)=ux
pcl(5,ip)=uy
pcl(6,ip)=uz

ENDDO

ENDIF

END SUBROUTINE INITIAL_PUSH

!***********************************************************************
! Subroutine BORIS_PUSH
! This subroutine computes the 4-velocity of the particles, using the
! Boris algorithm. INCLUDES THE SYNCHROTRON AND THE INVERSE COMPTION RADIATION
! REACTION FORCES, following the procedure by Tamburini et al., New Journal of
! Physics, 12 (2010) 123005.
!
! INPUT: 
! - q: sign electric charge
! - mass: Mass of the particle
! - pcl: Particle distribution function
! - pcl_data: Extra data about the particle distribution function
! - Bxg: x-component of B at the nodes at t
! - Byg: y-component of B at the nodes at t
! - Bzg: z-component of B at the nodes at t
! - Exg: x-component of E at the nodes at t
! - Eyg: y-component of E at the nodes at t
! - Ezg: z-component of E at the nodes at t
! - Uph: External photon energy density
! - xgp,ygp: spatial grid
! - Esyn: Total synchrotron energy losses between t and t+dt
! - Eics: Total inverse Compton energy losses between t and t+dt
! - NPP: Number of particles per process
!
! OUTPUT: ux,uy,uz at t=t+dt/2
!***********************************************************************

SUBROUTINE BORIS_PUSH(q,mass,pcl,pcl_data,Bxg,Byg,Bzg,Exg,Eyg,Ezg,Uph,&
                      xgp,ygp,Esyn,Eics,NPP)

IMPLICIT NONE

INTEGER*8 :: ip,NPP
DOUBLE PRECISION :: q,mass,gam,Bxi,Byi,Bzi,Exi,Eyi,Ezi,Psyn,Pics,Esyn,Eics,Uph
DOUBLE PRECISION :: uxm,uym,uzm,uxp,uyp,uzp,ux0,uy0,uz0,uxL,uyL,uzL,uxt,uyt,uzt
DOUBLE PRECISION :: sx,sy,sz,tx,ty,tz
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
DOUBLE PRECISION, DIMENSION(1:4,1:NPP)   :: pcl_data
DOUBLE PRECISION                         :: x,y,z,ux,uy,uz,wgt,El,Bp,Fr,Fe
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bxg,Byg,Bzg
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Exg,Eyg,Ezg
DOUBLE PRECISION, DIMENSION(1:6)         :: Fields

!***********************************************************************

Esyn=0.0
Eics=0.0

! Switch ON/OFF of the radiation reaction force:

!***********************************************************************
! WITH radiation reaction force
!***********************************************************************

IF (RAD_FORCE.EQV..TRUE.) THEN

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)
z=pcl(3,ip)
ux=pcl(4,ip)
uy=pcl(5,ip)
uz=pcl(6,ip)
wgt=pcl(7,ip)

!***********************************************************************
! Interpolation of E and B fields from the nodes to the particles at time t
CALL BILINEAR_FIELDS(xgp,ygp,x,y, Fields, Bxg,Byg,Bzg,Exg,Eyg,Ezg)

Bxi=Fields(1)
Byi=Fields(2)
Bzi=Fields(3)

Exi=Fields(4)
Eyi=Fields(5)
Ezi=Fields(6)

! Parallel electric field
pcl_data(1,ip)=(Exi*ux+Eyi*uy+Ezi*uz)/sqrt(ux*ux+uy*uy+uz*uz)

! Perpendicular magnetic field
pcl_data(2,ip)=sqrt((Byi*uz-Bzi*uy)**2.0+(Bzi*ux-Bxi*uz)**2.0+&
               (Bxi*uy-Byi*ux)**2.0)/sqrt(ux*ux+uy*uy+uz*uz)

! Electric force
pcl_data(4,ip)=e*sqrt(Exi*Exi+Eyi*Eyi+Ezi*Ezi)

!***********************************************************************
! Computation of uL(t+dt/2) using the Boris algorithm

! Calculation of u-
uxm=ux+q*e*Exi*dt/(2.0*mass*c)
uym=uy+q*e*Eyi*dt/(2.0*mass*c)
uzm=uz+q*e*Ezi*dt/(2.0*mass*c)

gam=sqrt(1.0+uxm*uxm+uym*uym+uzm*uzm)

! Calculation of u'
tx=q*e*Bxi*dt/(2.0*gam*mass*c)
ty=q*e*Byi*dt/(2.0*gam*mass*c)
tz=q*e*Bzi*dt/(2.0*gam*mass*c)

ux0=uxm+uym*tz-uzm*ty
uy0=uym+uzm*tx-uxm*tz
uz0=uzm+uxm*ty-uym*tx

! Calculation of u+
sx=2.0*tx/(1.0+tx*tx+ty*ty+tz*tz)
sy=2.0*ty/(1.0+tx*tx+ty*ty+tz*tz)
sz=2.0*tz/(1.0+tx*tx+ty*ty+tz*tz)

uxp=uxm+uy0*sz-uz0*sy
uyp=uym+uz0*sx-ux0*sz
uzp=uzm+ux0*sy-uy0*sx

! uL at t+dt/2
uxL=uxp+q*e*Exi*dt/(2.0*mass*c)
uyL=uyp+q*e*Eyi*dt/(2.0*mass*c)
uzL=uzp+q*e*Ezi*dt/(2.0*mass*c)

! uL at t
uxt=(uxL+ux)/2.0
uyt=(uyL+uy)/2.0
uzt=(uzL+uz)/2.0

! gamma at time t
gam=sqrt(1.0+uxt*uxt+uyt*uyt+uzt*uzt)

! Total synchrotron radiative power losses at t
Psyn=(2.0/3.0)*e**4.0/(mass**2.0*c**4.0)*c*((gam*Exi+uyt*Bzi-uzt*Byi)**2.0+&
     (gam*Eyi+uzt*Bxi-uxt*Bzi)**2.0+(gam*Ezi+uxt*Byi-uyt*Bxi)**2.0-&
     (uxt*Exi+uyt*Eyi+uzt*Ezi)**2.0)

! Synchrotron radiation reaction force
pcl_data(3,ip)=Psyn/(gam*c)*sqrt(uxt*uxt+uyt*uyt+uzt*uzt)

! Radiative inverse Compton energy losses (Thomson regime)
Pics=(32.0/9.0)*pi*e**4.0/(mass**2.0*c**4.0)*c*Uph*(uxt*uxt+uyt*uyt+uzt*uzt)

! u at t+dt/2
ux=uxL-dt*(Psyn+Pics)/(mass*c*c)*uxt/gam
uy=uyL-dt*(Psyn+Pics)/(mass*c*c)*uyt/gam
uz=uzL-dt*(Psyn+Pics)/(mass*c*c)*uzt/gam

! Total energy lost via synchrotron radiation between dt
Esyn=Esyn+wgt*Psyn*dt

! Total energy lost via inverse Compton scattering between dt
Eics=Eics+wgt*Pics*dt

!***********************************************************************

pcl(4,ip)=ux
pcl(5,ip)=uy
pcl(6,ip)=uz

ENDDO

!***********************************************************************
! WITHOUT radiation reaction force
!***********************************************************************

ELSE

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)
z=pcl(3,ip)
ux=pcl(4,ip)
uy=pcl(5,ip)
uz=pcl(6,ip)
wgt=pcl(7,ip)

!***********************************************************************
! Interpolation of E and B fields from the nodes to the particles at time t

CALL BILINEAR_FIELDS(xgp,ygp,x,y, Fields, Bxg,Byg,Bzg,Exg,Eyg,Ezg)

Bxi=Fields(1)
Byi=Fields(2)
Bzi=Fields(3)

Exi=Fields(4)
Eyi=Fields(5)
Ezi=Fields(6)

! Parallel electric field
pcl_data(1,ip)=(Exi*ux+Eyi*uy+Ezi*uz)/sqrt(ux*ux+uy*uy+uz*uz)

! Perpendicular magnetic field
pcl_data(2,ip)=sqrt((Byi*uz-Bzi*uy)**2.0+(Bzi*ux-Bxi*uz)**2.0+&
               (Bxi*uy-Byi*ux)**2.0)/sqrt(ux*ux+uy*uy+uz*uz)

! Electric force
pcl_data(4,ip)=e*sqrt(Exi*Exi+Eyi*Eyi+Ezi*Ezi)

! Synchrotron radiation reaction force
pcl_data(3,ip)=0.0

!***********************************************************************
! Computation of u(t+dt/2) using the Boris algorithm

! Calculation of u-
uxm=ux+q*e*Exi*dt/(2.0*mass*c)
uym=uy+q*e*Eyi*dt/(2.0*mass*c)
uzm=uz+q*e*Ezi*dt/(2.0*mass*c)

gam=sqrt(1.0+uxm*uxm+uym*uym+uzm*uzm)

! Calculation of u'
tx=q*e*Bxi*dt/(2.0*gam*mass*c)
ty=q*e*Byi*dt/(2.0*gam*mass*c)
tz=q*e*Bzi*dt/(2.0*gam*mass*c)

ux0=uxm+uym*tz-uzm*ty
uy0=uym+uzm*tx-uxm*tz
uz0=uzm+uxm*ty-uym*tx

! Calculation of u+
sx=2.0*tx/(1.0+tx*tx+ty*ty+tz*tz)
sy=2.0*ty/(1.0+tx*tx+ty*ty+tz*tz)
sz=2.0*tz/(1.0+tx*tx+ty*ty+tz*tz)

uxp=uxm+uy0*sz-uz0*sy
uyp=uym+uz0*sx-ux0*sz
uzp=uzm+ux0*sy-uy0*sx

! t-dt/2 --> t+dt/2
ux=uxp+q*e*Exi*dt/(2.0*mass*c)
uy=uyp+q*e*Eyi*dt/(2.0*mass*c)
uz=uzp+q*e*Ezi*dt/(2.0*mass*c)

!***********************************************************************

pcl(4,ip)=ux
pcl(5,ip)=uy
pcl(6,ip)=uz

ENDDO

ENDIF

END SUBROUTINE BORIS_PUSH

!***********************************************************************
! Subroutine PUSH_PARTICLES
! This subroutine updates the positions of particles from the 4-velocity
! u known at t+dt/2
!
! INPUT: 
! - pcl: Particle distribution function
! - NPP: Number of particles per process
!
! OUTPUT: x,y,z at t=t+dt
!***********************************************************************

SUBROUTINE PUSH_PARTICLES(pcl,NPP)

IMPLICIT NONE

INTEGER*8 :: ip,NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
DOUBLE PRECISION                         :: x,y,z,ux,uy,uz,gam

!***********************************************************************

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)
z=pcl(3,ip)
ux=pcl(4,ip)
uy=pcl(5,ip)
uz=pcl(6,ip)

! Particle's Lorentz factor
gam=sqrt(1d0+ux*ux+uy*uy+uz*uz)

!***********************************************************************
! Computation of x(t+dt)
!***********************************************************************

! t --> t+dt
x=x+(c*dt/gam)*ux
y=y+(c*dt/gam)*uy
z=z+(c*dt/gam)*uz

pcl(1,ip)=x
pcl(2,ip)=y
pcl(3,ip)=z
pcl(4,ip)=ux
pcl(5,ip)=uy
pcl(6,ip)=uz

ENDDO

END SUBROUTINE PUSH_PARTICLES

!***********************************************************************
! Subroutine COUNT_ESCAPE
! This subroutine counts the number of particles escaping the domain 
! towards the neighbors.

! INPUT: 
! - pcl: particle distribution function
! - xminp,xmaxp,yminp,ymaxp: spatial boundaries of the domain
! - NPP: Number of macroparticles per domain
! - NESC: Number of particles escaping the domain towards neighbors
!
! OUTPUT: NESC at time t+dt
!***********************************************************************

SUBROUTINE COUNT_ESCAPE(pcl,xminp,xmaxp,yminp,ymaxp,NPP,NESC)

IMPLICIT NONE

INTEGER*8                              :: ip,NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP) :: pcl
DOUBLE PRECISION                       :: x,y,xminp,xmaxp,yminp,ymaxp
INTEGER, DIMENSION(8)                  :: NESC

!***********************************************************************

NESC=0

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)

! Case 1: the particle escapes north
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.GT.ymaxp) THEN
NESC(1)=NESC(1)+1
END IF

! Case 2: the particle escapes East
IF (x.GT.xmaxp.AND.y.LE.ymaxp.AND.y.GE.yminp) THEN
NESC(2)=NESC(2)+1
END IF

! Case 3: the particle escapes South
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LT.yminp) THEN
NESC(3)=NESC(3)+1
END IF

! Case 4: the particle escapes West
IF (x.LT.xminp.AND.y.LE.ymaxp.AND.y.GE.yminp) THEN
NESC(4)=NESC(4)+1
END IF

! Case 5: the particle escapes North-East
IF (x.GT.xmaxp.AND.y.GT.ymaxp) THEN
NESC(5)=NESC(5)+1
END IF

! Case 6: the particle escapes South-East
IF (x.GT.xmaxp.AND.y.LT.yminp) THEN
NESC(6)=NESC(6)+1
END IF

! Case 7: the particle escapes South-West
IF (x.LT.xminp.AND.y.LT.yminp) THEN
NESC(7)=NESC(7)+1
END IF

! Case 8: the particle escapes North-West
IF (x.LT.xminp.AND.y.GT.ymaxp) THEN
NESC(8)=NESC(8)+1
END IF

ENDDO

END SUBROUTINE COUNT_ESCAPE

!***********************************************************************
! Subroutine BOUNDARIES_PARTICLES
! This subroutine applies the boundary conditions chosen in the input 
! file to all the particles leaving the boundaries of the box.
!
! INPUT: 
! - pcl: particle distribution function
! - pcl: particle data distribution function
! - tag: particle tags
! - NPP: Number of particle per domain
!
! OUTPUT: Updated particle distribution function at time t+dt
!***********************************************************************

SUBROUTINE BOUNDARIES_PARTICLES(pcl,pcl_data,tag,NPP)

IMPLICIT NONE

INTEGER*8 :: ip,is,NPP,NTEMP
DOUBLE PRECISION, ALLOCATABLE :: pcl(:,:)
DOUBLE PRECISION, ALLOCATABLE :: pcl_data(:,:)
INTEGER*8, ALLOCATABLE        :: tag(:)
DOUBLE PRECISION              :: x,y,ux,uy,uz,wt

!***********************************************************************

NTEMP=0

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)
ux=pcl(4,ip)
uy=pcl(5,ip)
uz=pcl(6,ip)
wt=pcl(7,ip)

!***********************************************************************
! Case 1: x>xmax
!***********************************************************************
IF (x.GT.xmax) THEN

   ! Elastic reflection
   IF (BOUND_PART_XMAX.EQ."REFLECT") THEN
   x=2.0*xmax-x
   ux=-ux
   END IF
   
   ! Absorption
   IF (BOUND_PART_XMAX.EQ."ABSORB") THEN
   wt=0d0
   END IF
   
END IF

!***********************************************************************
! Case 2: x<xmin
!***********************************************************************
IF (x.LT.xmin) THEN

   ! Elastic reflection
   IF (BOUND_PART_XMIN.EQ."REFLECT") THEN
   x=2.0*xmin-x
   ux=-ux
   END IF
   
   ! Absorption
   IF (BOUND_PART_XMIN.EQ."ABSORB") THEN
   wt=0d0
   END IF
   
END IF

!***********************************************************************
! Case 3: y>ymax
!***********************************************************************
IF (y.GT.ymax) THEN

   ! Elastic reflection
   IF (BOUND_PART_YMAX.EQ."REFLECT") THEN
   y=2.0*ymax-y
   uy=-uy
   END IF
   
   ! Absorption
   IF (BOUND_PART_YMAX.EQ."ABSORB") THEN
   wt=0d0
   END IF

END IF

!***********************************************************************
! Case 4: y<ymin
!***********************************************************************
IF (y.LT.ymin) THEN

   ! Elastic reflection
   IF (BOUND_PART_YMIN.EQ."REFLECT") THEN
   y=2.0*ymin-y
   uy=-uy
   END IF
   
   ! Absorption
   IF (BOUND_PART_YMIN.EQ."ABSORB") THEN
   wt=0d0
   END IF

END IF

pcl(1,ip)=x
pcl(2,ip)=y
pcl(4,ip)=ux
pcl(5,ip)=uy
pcl(6,ip)=uz
pcl(7,ip)=wt

!***********************************************************************
! Counting particles remaining in the box
!***********************************************************************
IF (wt.NE.0d0) THEN
NTEMP=NTEMP+1
END IF

ENDDO

!***********************************************************************

ALLOCATE(pcl_f(1:7,1:NTEMP),pcl_data_f(1:4,1:NTEMP))
ALLOCATE(tagf(1:NTEMP))

!***********************************************************************
! Removing the leaking particles
!***********************************************************************

IF (NTEMP.NE.NPP) THEN

is=0

DO ip=1,NPP

   IF (pcl(7,ip).NE.0d0) THEN
   
   is=is+1
   
   pcl_f(:,is)=pcl(:,ip)
   pcl_data_f(:,is)=pcl_data(:,ip)
   tagf(is)=tag(ip)
   
   ENDIF
   
ENDDO

ELSE

pcl_f(:,:)=pcl(:,:)
pcl_data_f(:,:)=pcl_data(:,:)
tagf(:)=tag(:)

END IF

DEALLOCATE(pcl,pcl_data)
DEALLOCATE(tag)

!***********************************************************************
! Transfer of memory and content FROM pcl_f TO pcl

CALL MOVE_ALLOC(pcl_f,pcl)
CALL MOVE_ALLOC(pcl_data_f,pcl_data)
CALL MOVE_ALLOC(tagf,tag)

NPP=NTEMP

!***********************************************************************

END SUBROUTINE BOUNDARIES_PARTICLES

!***********************************************************************
! This subroutine exchanges particles at the boundaries between processes
! in the virtual topology.
!
! INPUT:
! - pcl: Particle distribution function at time t+dt
! - pcl_data: Extra data of the particle distribution function at time t+dt
! - tagi: Particle tags at time t+dt
! - xminp,xmaxp,yminp,ymaxp: Spatial boundaries of the domain
! - NPI: number of particles per process at time t (input) and t+dt (output)
! - NESC: Number of particles escaping the domain towards neighbors
! - id: process rank
! - ngh: neighbor array
! - COMM: communicator
! - ierr: error code
!
! OUTPUT: Particles in the domain at time t+dt
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE COM_PARTICLES(pcl,pcl_data,tagi,xminp,xmaxp,yminp,ymaxp,NPI,NESC,&
                         id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)    :: stat
DOUBLE PRECISION, ALLOCATABLE          :: pcl(:,:)
DOUBLE PRECISION, ALLOCATABLE          :: pcl_data(:,:)
INTEGER*8, ALLOCATABLE                 :: tagi(:)
DOUBLE PRECISION                       :: x,y
DOUBLE PRECISION                       :: xminp,xmaxp,yminp,ymaxp
INTEGER, DIMENSION(8)                  :: ngh,NESC,NINC
INTEGER, PARAMETER                     :: tag1=11,tag2=21,tag3=31,tag4=41
INTEGER, PARAMETER                     :: tag5=51,tag6=61,tag7=71,tag8=81
INTEGER                                :: id,COMM,ierr

DOUBLE PRECISION, ALLOCATABLE :: bufSN(:,:),bufSE(:,:),bufSS(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSW(:,:),bufSNE(:,:),bufSSE(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSSW(:,:),bufSNW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRN(:,:),bufRE(:,:),bufRS(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRW(:,:),bufRNE(:,:),bufRSE(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRSW(:,:),bufRNW(:,:)

INTEGER*8, ALLOCATABLE          :: buf_tagSN(:),buf_tagSE(:)
INTEGER*8, ALLOCATABLE          :: buf_tagSS(:),buf_tagSW(:)
INTEGER*8, ALLOCATABLE          :: buf_tagSNE(:),buf_tagSSE(:)
INTEGER*8, ALLOCATABLE          :: buf_tagSSW(:),buf_tagSNW(:)

INTEGER*8, ALLOCATABLE          :: buf_tagRN(:),buf_tagRE(:)
INTEGER*8, ALLOCATABLE          :: buf_tagRS(:),buf_tagRW(:)
INTEGER*8, ALLOCATABLE          :: buf_tagRNE(:),buf_tagRSE(:)
INTEGER*8, ALLOCATABLE          :: buf_tagRSW(:),buf_tagRNW(:)

INTEGER*8 :: NPI,NPO,NPT,NPT2,Ntemp
INTEGER*8 :: ip,ipo
INTEGER*8 :: i1,i2,i3,i4,i5,i6,i7,i8
!***********************************************************************

i1=0
i2=0
i3=0
i4=0
i5=0
i6=0
i7=0
i8=0

NPT2=NPI

ALLOCATE(bufSN(NESC(1),11),bufSE(NESC(2),11),bufSS(NESC(3),11),bufSW(NESC(4),11))
ALLOCATE(bufSNE(NESC(5),11),bufSSE(NESC(6),11),bufSSW(NESC(7),11))
ALLOCATE(bufSNW(NESC(8),11))

ALLOCATE(buf_tagSN(NESC(1)),buf_tagSE(NESC(2)),buf_tagSS(NESC(3)))
ALLOCATE(buf_tagSW(NESC(4)),buf_tagSNE(NESC(5)),buf_tagSSE(NESC(6)))
ALLOCATE(buf_tagSSW(NESC(7)),buf_tagSNW(NESC(8)))

!***********************************************************************
! ESCAPING PARTICLE ARRAY TO BE SEND
!***********************************************************************

DO ip=1,NPI

x=pcl(1,ip)
y=pcl(2,ip)

! Case 1: the particle escapes north
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.GT.ymaxp) THEN
i1=i1+1
bufSN(i1,1)=pcl(1,ip)
bufSN(i1,2)=pcl(2,ip)
bufSN(i1,3)=pcl(3,ip)
bufSN(i1,4)=pcl(4,ip)
bufSN(i1,5)=pcl(5,ip)
bufSN(i1,6)=pcl(6,ip)
bufSN(i1,7)=pcl(7,ip)
bufSN(i1,8)=pcl_data(1,ip)
bufSN(i1,9)=pcl_data(2,ip)
bufSN(i1,10)=pcl_data(3,ip)
bufSN(i1,11)=pcl_data(4,ip)
buf_tagSN(i1)=tagi(ip)

  IF (bufSN(i1,2).GT.ymax) THEN
  bufSN(i1,2)=bufSN(i1,2)-(ymax-ymin)
  END IF

END IF

! Case 2: the particle escapes East
IF (x.GT.xmaxp.AND.y.LE.ymaxp.AND.y.GE.yminp) THEN
i2=i2+1
bufSE(i2,1)=pcl(1,ip)
bufSE(i2,2)=pcl(2,ip)
bufSE(i2,3)=pcl(3,ip)
bufSE(i2,4)=pcl(4,ip)
bufSE(i2,5)=pcl(5,ip)
bufSE(i2,6)=pcl(6,ip)
bufSE(i2,7)=pcl(7,ip)
bufSE(i2,8)=pcl_data(1,ip)
bufSE(i2,9)=pcl_data(2,ip)
bufSE(i2,10)=pcl_data(3,ip)
bufSE(i2,11)=pcl_data(4,ip)
buf_tagSE(i2)=tagi(ip)

  IF (bufSE(i2,1).GT.xmax) THEN
  bufSE(i2,1)=bufSE(i2,1)-(xmax-xmin)
  END IF

END IF

! Case 3: the particle escapes South
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LT.yminp) THEN
i3=i3+1
bufSS(i3,1)=pcl(1,ip)
bufSS(i3,2)=pcl(2,ip)
bufSS(i3,3)=pcl(3,ip)
bufSS(i3,4)=pcl(4,ip)
bufSS(i3,5)=pcl(5,ip)
bufSS(i3,6)=pcl(6,ip)
bufSS(i3,7)=pcl(7,ip)
bufSS(i3,8)=pcl_data(1,ip)
bufSS(i3,9)=pcl_data(2,ip)
bufSS(i3,10)=pcl_data(3,ip)
bufSS(i3,11)=pcl_data(4,ip)
buf_tagSS(i3)=tagi(ip)

  IF (bufSS(i3,2).LT.ymin) THEN
  bufSS(i3,2)=bufSS(i3,2)+(ymax-ymin)
  END IF

END IF

! Case 4: the particle escapes West
IF (x.LT.xminp.AND.y.LE.ymaxp.AND.y.GE.yminp) THEN
i4=i4+1
bufSW(i4,1)=pcl(1,ip)
bufSW(i4,2)=pcl(2,ip)
bufSW(i4,3)=pcl(3,ip)
bufSW(i4,4)=pcl(4,ip)
bufSW(i4,5)=pcl(5,ip)
bufSW(i4,6)=pcl(6,ip)
bufSW(i4,7)=pcl(7,ip)
bufSW(i4,8)=pcl_data(1,ip)
bufSW(i4,9)=pcl_data(2,ip)
bufSW(i4,10)=pcl_data(3,ip)
bufSW(i4,11)=pcl_data(4,ip)
buf_tagSW(i4)=tagi(ip)

  IF (bufSW(i4,1).LT.xmin) THEN
  bufSW(i4,1)=bufSW(i4,1)+(xmax-xmin)
  END IF

END IF

! Case 5: the particle escapes North-East
IF (x.GT.xmaxp.AND.y.GT.ymaxp) THEN
i5=i5+1
bufSNE(i5,1)=pcl(1,ip)
bufSNE(i5,2)=pcl(2,ip)
bufSNE(i5,3)=pcl(3,ip)
bufSNE(i5,4)=pcl(4,ip)
bufSNE(i5,5)=pcl(5,ip)
bufSNE(i5,6)=pcl(6,ip)
bufSNE(i5,7)=pcl(7,ip)
bufSNE(i5,8)=pcl_data(1,ip)
bufSNE(i5,9)=pcl_data(2,ip)
bufSNE(i5,10)=pcl_data(3,ip)
bufSNE(i5,11)=pcl_data(4,ip)
buf_tagSNE(i5)=tagi(ip)

  IF (bufSNE(i5,1).GT.xmax) THEN
  bufSNE(i5,1)=bufSNE(i5,1)-(xmax-xmin)
  END IF

  IF (bufSNE(i5,2).GT.ymax) THEN
  bufSNE(i5,2)=bufSNE(i5,2)-(ymax-ymin)
  END IF

END IF

! Case 6: the particle escapes South-East
IF (x.GT.xmaxp.AND.y.LT.yminp) THEN
i6=i6+1
bufSSE(i6,1)=pcl(1,ip)
bufSSE(i6,2)=pcl(2,ip)
bufSSE(i6,3)=pcl(3,ip)
bufSSE(i6,4)=pcl(4,ip)
bufSSE(i6,5)=pcl(5,ip)
bufSSE(i6,6)=pcl(6,ip)
bufSSE(i6,7)=pcl(7,ip)
bufSSE(i6,8)=pcl_data(1,ip)
bufSSE(i6,9)=pcl_data(2,ip)
bufSSE(i6,10)=pcl_data(3,ip)
bufSSE(i6,11)=pcl_data(4,ip)
buf_tagSSE(i6)=tagi(ip)

  IF (bufSSE(i6,1).GT.xmax) THEN
  bufSSE(i6,1)=bufSSE(i6,1)-(xmax-xmin)
  END IF

  IF (bufSSE(i6,2).LT.ymin) THEN
  bufSSE(i6,2)=bufSSE(i6,2)+(ymax-ymin)
  END IF

END IF

! Case 7: the particle escapes South-West
IF (x.LT.xminp.AND.y.LT.yminp) THEN
i7=i7+1
bufSSW(i7,1)=pcl(1,ip)
bufSSW(i7,2)=pcl(2,ip)
bufSSW(i7,3)=pcl(3,ip)
bufSSW(i7,4)=pcl(4,ip)
bufSSW(i7,5)=pcl(5,ip)
bufSSW(i7,6)=pcl(6,ip)
bufSSW(i7,7)=pcl(7,ip)
bufSSW(i7,8)=pcl_data(1,ip)
bufSSW(i7,9)=pcl_data(2,ip)
bufSSW(i7,10)=pcl_data(3,ip)
bufSSW(i7,11)=pcl_data(4,ip)
buf_tagSSW(i7)=tagi(ip)

  IF (bufSSW(i7,1).LT.xmin) THEN
  bufSSW(i7,1)=bufSSW(i7,1)+(xmax-xmin)
  END IF

  IF (bufSSW(i7,2).LT.ymin) THEN
  bufSSW(i7,2)=bufSSW(i7,2)+(ymax-ymin)
  END IF

END IF

! Case 8: the particle escapes North-West
IF (x.LT.xminp.AND.y.GT.ymaxp) THEN
i8=i8+1
bufSNW(i8,1)=pcl(1,ip)
bufSNW(i8,2)=pcl(2,ip)
bufSNW(i8,3)=pcl(3,ip)
bufSNW(i8,4)=pcl(4,ip)
bufSNW(i8,5)=pcl(5,ip)
bufSNW(i8,6)=pcl(6,ip)
bufSNW(i8,7)=pcl(7,ip)
bufSNW(i8,8)=pcl_data(1,ip)
bufSNW(i8,9)=pcl_data(2,ip)
bufSNW(i8,10)=pcl_data(3,ip)
bufSNW(i8,11)=pcl_data(4,ip)
buf_tagSNW(i8)=tagi(ip)

  IF (bufSNW(i8,1).LT.xmin) THEN
  bufSNW(i8,1)=bufSNW(i8,1)+(xmax-xmin)
  END IF

  IF (bufSNW(i8,2).GT.ymax) THEN
  bufSNW(i8,2)=bufSNW(i8,2)-(ymax-ymin)
  END IF

END IF

ENDDO

!***********************************************************************
! EXCHANGE OF PARTICLES BETWEEN PROCESSES
!***********************************************************************

! Exchange sizes of the OUTGOING buffered arrays.

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(NESC(1),1,MPI_INTEGER,ngh(1),tag1,&
                  NINC(3),1,MPI_INTEGER,ngh(3),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(2),1,MPI_INTEGER,ngh(2),tag2,&
                  NINC(4),1,MPI_INTEGER,ngh(4),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(3),1,MPI_INTEGER,ngh(3),tag3,&
                  NINC(1),1,MPI_INTEGER,ngh(1),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(4),1,MPI_INTEGER,ngh(4),tag4,&
                  NINC(2),1,MPI_INTEGER,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(5),1,MPI_INTEGER,ngh(5),tag5,&
                  NINC(7),1,MPI_INTEGER,ngh(7),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(6),1,MPI_INTEGER,ngh(6),tag6,&
                  NINC(8),1,MPI_INTEGER,ngh(8),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(7),1,MPI_INTEGER,ngh(7),tag7,&
                  NINC(5),1,MPI_INTEGER,ngh(5),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(8),1,MPI_INTEGER,ngh(8),tag8,&
                  NINC(6),1,MPI_INTEGER,ngh(6),tag8,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(NESC(1),1,MPI_INTEGER,ngh(1),tag1,&
                  NINC(3),1,MPI_INTEGER,ngh(3),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(2),1,MPI_INTEGER,ngh(2),tag2,&
                  NINC(4),1,MPI_INTEGER,ngh(4),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(3),1,MPI_INTEGER,ngh(3),tag3,&
                  NINC(1),1,MPI_INTEGER,ngh(1),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(4),1,MPI_INTEGER,ngh(4),tag4,&
                  NINC(2),1,MPI_INTEGER,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(5),1,MPI_INTEGER,ngh(5),tag5,&
                  NINC(7),1,MPI_INTEGER,ngh(7),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(6),1,MPI_INTEGER,ngh(6),tag6,&
                  NINC(8),1,MPI_INTEGER,ngh(8),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(7),1,MPI_INTEGER,ngh(7),tag7,&
                  NINC(5),1,MPI_INTEGER,ngh(5),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(8),1,MPI_INTEGER,ngh(8),tag8,&
                  NINC(6),1,MPI_INTEGER,ngh(6),tag8,COMM,stat,ierr)

ENDIF

!***********************************************************************
! Total number of unchanged particles
NPT=NPI-(NESC(1)+NESC(2)+NESC(3)+NESC(4)+NESC(5)+NESC(6)+NESC(7)+NESC(8))

! Total number of particles (+incoming-outcoming)
NPO=NPT+(NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8))
!***********************************************************************

! Exchange of the OUTCOMING and INCOMING buffered arrays.

ALLOCATE(bufRN(NINC(1),11),bufRE(NINC(2),11),bufRS(NINC(3),11),bufRW(NINC(4),11))
ALLOCATE(bufRNE(NINC(5),11),bufRSE(NINC(6),11),bufRSW(NINC(7),11))
ALLOCATE(bufRNW(NINC(8),11))

ALLOCATE(buf_tagRN(NINC(1)),buf_tagRE(NINC(2)),buf_tagRS(NINC(3)))
ALLOCATE(buf_tagRW(NINC(4)),buf_tagRNE(NINC(5)),buf_tagRSE(NINC(6)))
ALLOCATE(buf_tagRSW(NINC(7)),buf_tagRNW(NINC(8)))

! Exchange particle positions, velocity and weight

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSN,NESC(1)*11,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
     bufRS,NINC(3)*11,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSE,NESC(2)*11,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
     bufRW,NINC(4)*11,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS,NESC(3)*11,MPI_DOUBLE_PRECISION,ngh(3),tag3,&
     bufRN,NINC(1)*11,MPI_DOUBLE_PRECISION,ngh(1),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NESC(4)*11,MPI_DOUBLE_PRECISION,ngh(4),tag4,&
     bufRE,NINC(2)*11,MPI_DOUBLE_PRECISION,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNE,NESC(5)*11,MPI_DOUBLE_PRECISION,ngh(5),tag5,&
     bufRSW,NINC(7)*11,MPI_DOUBLE_PRECISION,ngh(7),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSE,NESC(6)*11,MPI_DOUBLE_PRECISION,ngh(6),tag6,&
     bufRNW,NINC(8)*11,MPI_DOUBLE_PRECISION,ngh(8),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSW,NESC(7)*11,MPI_DOUBLE_PRECISION,ngh(7),tag7,&
     bufRNE,NINC(5)*11,MPI_DOUBLE_PRECISION,ngh(5),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNW,NESC(8)*11,MPI_DOUBLE_PRECISION,ngh(8),tag8,&
     bufRSE,NINC(6)*11,MPI_DOUBLE_PRECISION,ngh(6),tag8,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(bufSN,NESC(1)*11,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
     bufRS,NINC(3)*11,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSE,NESC(2)*11,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
     bufRW,NINC(4)*11,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS,NESC(3)*11,MPI_DOUBLE_PRECISION,ngh(3),tag3,&
     bufRN,NINC(1)*11,MPI_DOUBLE_PRECISION,ngh(1),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NESC(4)*11,MPI_DOUBLE_PRECISION,ngh(4),tag4,&
     bufRE,NINC(2)*11,MPI_DOUBLE_PRECISION,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNE,NESC(5)*11,MPI_DOUBLE_PRECISION,ngh(5),tag5,&
     bufRSW,NINC(7)*11,MPI_DOUBLE_PRECISION,ngh(7),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSE,NESC(6)*11,MPI_DOUBLE_PRECISION,ngh(6),tag6,&
     bufRNW,NINC(8)*11,MPI_DOUBLE_PRECISION,ngh(8),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSW,NESC(7)*11,MPI_DOUBLE_PRECISION,ngh(7),tag7,&
     bufRNE,NINC(5)*11,MPI_DOUBLE_PRECISION,ngh(5),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNW,NESC(8)*11,MPI_DOUBLE_PRECISION,ngh(8),tag8,&
     bufRSE,NINC(6)*11,MPI_DOUBLE_PRECISION,ngh(6),tag8,COMM,stat,ierr)

ENDIF

! Exchange particle tags

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(buf_tagSN,NESC(1),MPI_INTEGER8,ngh(1),tag1,&
                  buf_tagRS,NINC(3),MPI_INTEGER8,ngh(3),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSE,NESC(2),MPI_INTEGER8,ngh(2),tag2,&
                  buf_tagRW,NINC(4),MPI_INTEGER8,ngh(4),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSS,NESC(3),MPI_INTEGER8,ngh(3),tag3,&
                  buf_tagRN,NINC(1),MPI_INTEGER8,ngh(1),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSW,NESC(4),MPI_INTEGER8,ngh(4),tag4,&
                  buf_tagRE,NINC(2),MPI_INTEGER8,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSNE,NESC(5),MPI_INTEGER8,ngh(5),tag5,&
                  buf_tagRSW,NINC(7),MPI_INTEGER8,ngh(7),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSSE,NESC(6),MPI_INTEGER8,ngh(6),tag6,&
                  buf_tagRNW,NINC(8),MPI_INTEGER8,ngh(8),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSSW,NESC(7),MPI_INTEGER8,ngh(7),tag7,&
                  buf_tagRNE,NINC(5),MPI_INTEGER8,ngh(5),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSNW,NESC(8),MPI_INTEGER8,ngh(8),tag8,&
                  buf_tagRSE,NINC(6),MPI_INTEGER8,ngh(6),tag8,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(buf_tagSN,NESC(1),MPI_INTEGER8,ngh(1),tag1,&
                  buf_tagRS,NINC(3),MPI_INTEGER8,ngh(3),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSE,NESC(2),MPI_INTEGER8,ngh(2),tag2,&
                  buf_tagRW,NINC(4),MPI_INTEGER8,ngh(4),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSS,NESC(3),MPI_INTEGER8,ngh(3),tag3,&
                  buf_tagRN,NINC(1),MPI_INTEGER8,ngh(1),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSW,NESC(4),MPI_INTEGER8,ngh(4),tag4,&
                  buf_tagRE,NINC(2),MPI_INTEGER8,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSNE,NESC(5),MPI_INTEGER8,ngh(5),tag5,&
                  buf_tagRSW,NINC(7),MPI_INTEGER8,ngh(7),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSSE,NESC(6),MPI_INTEGER8,ngh(6),tag6,&
                  buf_tagRNW,NINC(8),MPI_INTEGER8,ngh(8),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSSW,NESC(7),MPI_INTEGER8,ngh(7),tag7,&
                  buf_tagRNE,NINC(5),MPI_INTEGER8,ngh(5),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSNW,NESC(8),MPI_INTEGER8,ngh(8),tag8,&
                  buf_tagRSE,NINC(6),MPI_INTEGER8,ngh(6),tag8,COMM,stat,ierr)

END IF

!***********************************************************************
! New particle arrays
!***********************************************************************

ALLOCATE(pcl_f(1:7,1:NPO))
ALLOCATE(pcl_data_f(1:4,1:NPO))
ALLOCATE(tagf(NPO))

pcl_f=0.0
pcl_data_f=0.0
tagf=0

ipo=0

DO ip=1,NPT2

x=pcl(1,ip)
y=pcl(2,ip)

IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LE.ymaxp.AND.y.GE.yminp) THEN
ipo=ipo+1
pcl_f(1,ipo)=pcl(1,ip)
pcl_f(2,ipo)=pcl(2,ip)
pcl_f(3,ipo)=pcl(3,ip)
pcl_f(4,ipo)=pcl(4,ip)
pcl_f(5,ipo)=pcl(5,ip)
pcl_f(6,ipo)=pcl(6,ip)
pcl_f(7,ipo)=pcl(7,ip)
pcl_data_f(1,ipo)=pcl_data(1,ip)
pcl_data_f(2,ipo)=pcl_data(2,ip)
pcl_data_f(3,ipo)=pcl_data(3,ip)
pcl_data_f(4,ipo)=pcl_data(4,ip)
tagf(ipo)=tagi(ip)
END IF
END DO

Ntemp=NPT

DO ip=1,NINC(1)
pcl_f(1,ip+Ntemp)=bufRN(ip,1)
pcl_f(2,ip+Ntemp)=bufRN(ip,2)
pcl_f(3,ip+Ntemp)=bufRN(ip,3)
pcl_f(4,ip+Ntemp)=bufRN(ip,4)
pcl_f(5,ip+Ntemp)=bufRN(ip,5)
pcl_f(6,ip+Ntemp)=bufRN(ip,6)
pcl_f(7,ip+Ntemp)=bufRN(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRN(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRN(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRN(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRN(ip,11)
tagf(ip+Ntemp)=buf_tagRN(ip)
ENDDO

Ntemp=NPT+NINC(1)

DO ip=1,NINC(2)
pcl_f(1,ip+Ntemp)=bufRE(ip,1)
pcl_f(2,ip+Ntemp)=bufRE(ip,2)
pcl_f(3,ip+Ntemp)=bufRE(ip,3)
pcl_f(4,ip+Ntemp)=bufRE(ip,4)
pcl_f(5,ip+Ntemp)=bufRE(ip,5)
pcl_f(6,ip+Ntemp)=bufRE(ip,6)
pcl_f(7,ip+Ntemp)=bufRE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRE(ip,11)
tagf(ip+Ntemp)=buf_tagRE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)

DO ip=1,NINC(3)
pcl_f(1,ip+Ntemp)=bufRS(ip,1)
pcl_f(2,ip+Ntemp)=bufRS(ip,2)
pcl_f(3,ip+Ntemp)=bufRS(ip,3)
pcl_f(4,ip+Ntemp)=bufRS(ip,4)
pcl_f(5,ip+Ntemp)=bufRS(ip,5)
pcl_f(6,ip+Ntemp)=bufRS(ip,6)
pcl_f(7,ip+Ntemp)=bufRS(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRS(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRS(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRS(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRS(ip,11)
tagf(ip+Ntemp)=buf_tagRS(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)

DO ip=1,NINC(4)
pcl_f(1,ip+Ntemp)=bufRW(ip,1)
pcl_f(2,ip+Ntemp)=bufRW(ip,2)
pcl_f(3,ip+Ntemp)=bufRW(ip,3)
pcl_f(4,ip+Ntemp)=bufRW(ip,4)
pcl_f(5,ip+Ntemp)=bufRW(ip,5)
pcl_f(6,ip+Ntemp)=bufRW(ip,6)
pcl_f(7,ip+Ntemp)=bufRW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRW(ip,11)
tagf(ip+Ntemp)=buf_tagRW(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)

DO ip=1,NINC(5)
pcl_f(1,ip+Ntemp)=bufRNE(ip,1)
pcl_f(2,ip+Ntemp)=bufRNE(ip,2)
pcl_f(3,ip+Ntemp)=bufRNE(ip,3)
pcl_f(4,ip+Ntemp)=bufRNE(ip,4)
pcl_f(5,ip+Ntemp)=bufRNE(ip,5)
pcl_f(6,ip+Ntemp)=bufRNE(ip,6)
pcl_f(7,ip+Ntemp)=bufRNE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRNE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRNE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRNE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRNE(ip,11)
tagf(ip+Ntemp)=buf_tagRNE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)

DO ip=1,NINC(6)
pcl_f(1,ip+Ntemp)=bufRSE(ip,1)
pcl_f(2,ip+Ntemp)=bufRSE(ip,2)
pcl_f(3,ip+Ntemp)=bufRSE(ip,3)
pcl_f(4,ip+Ntemp)=bufRSE(ip,4)
pcl_f(5,ip+Ntemp)=bufRSE(ip,5)
pcl_f(6,ip+Ntemp)=bufRSE(ip,6)
pcl_f(7,ip+Ntemp)=bufRSE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRSE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRSE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRSE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRSE(ip,11)
tagf(ip+Ntemp)=buf_tagRSE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)

DO ip=1,NINC(7)
pcl_f(1,ip+Ntemp)=bufRSW(ip,1)
pcl_f(2,ip+Ntemp)=bufRSW(ip,2)
pcl_f(3,ip+Ntemp)=bufRSW(ip,3)
pcl_f(4,ip+Ntemp)=bufRSW(ip,4)
pcl_f(5,ip+Ntemp)=bufRSW(ip,5)
pcl_f(6,ip+Ntemp)=bufRSW(ip,6)
pcl_f(7,ip+Ntemp)=bufRSW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRSW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRSW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRSW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRSW(ip,11)
tagf(ip+Ntemp)=buf_tagRSW(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)

DO ip=1,NINC(8)
pcl_f(1,ip+Ntemp)=bufRNW(ip,1)
pcl_f(2,ip+Ntemp)=bufRNW(ip,2)
pcl_f(3,ip+Ntemp)=bufRNW(ip,3)
pcl_f(4,ip+Ntemp)=bufRNW(ip,4)
pcl_f(5,ip+Ntemp)=bufRNW(ip,5)
pcl_f(6,ip+Ntemp)=bufRNW(ip,6)
pcl_f(7,ip+Ntemp)=bufRNW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRNW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRNW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRNW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRNW(ip,11)
tagf(ip+Ntemp)=buf_tagRNW(ip)
ENDDO

DEALLOCATE(bufSN,bufSE,bufSS,bufSW,bufSNE,bufSSE,bufSSW,bufSNW)
DEALLOCATE(bufRN,bufRE,bufRS,bufRW,bufRNE,bufRSE,bufRSW,bufRNW)

DEALLOCATE(buf_tagSN,buf_tagSE,buf_tagSS,buf_tagSW)
DEALLOCATE(buf_tagSNE,buf_tagSSE,buf_tagSSW,buf_tagSNW)

DEALLOCATE(buf_tagRN,buf_tagRE,buf_tagRS,buf_tagRW)
DEALLOCATE(buf_tagRNE,buf_tagRSE,buf_tagRSW,buf_tagRNW)

DEALLOCATE(pcl,pcl_data)
DEALLOCATE(tagi)

! Transfer of memory and content FROM pcl_f TO pcl
CALL MOVE_ALLOC(pcl_f,pcl)
CALL MOVE_ALLOC(pcl_data_f,pcl_data)
CALL MOVE_ALLOC(tagf,tagi)

! Update the number of particles per process at t+dt
NPI=NPO

END SUBROUTINE COM_PARTICLES

!***********************************************************************

END MODULE MOD_MOTION
