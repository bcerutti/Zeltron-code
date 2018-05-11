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
! - pcl: Particle distribution function
! - Bxg: x-component of B at the nodes at t
! - Byg: y-component of B at the nodes at t
! - Bzg: z-component of B at the nodes at t
! - Exg: x-component of E at the nodes at t
! - Eyg: y-component of E at the nodes at t
! - Ezg: z-component of E at the nodes at t
! - NPP: Number of particles per process
! - Uph: External photon energy density
!
! OUTPUT: ux,uy,uz at t=-dt/2
!***********************************************************************

SUBROUTINE INITIAL_PUSH(q,mass,pcl,Bxg,Byg,Bzg,Exg,Eyg,Ezg,Uph,xgp,ygp,zgp,NPP)

IMPLICIT NONE

INTEGER*8 :: ip,NPP
DOUBLE PRECISION :: q,mass,gam,Bxi,Byi,Bzi,Exi,Eyi,Ezi,Psyn,Pics,Uph
DOUBLE PRECISION :: uxp,uyp,uzp,ux0,uy0,uz0,uxL,uyL,uzL
DOUBLE PRECISION :: sx,sy,sz,tx,ty,tz
DOUBLE PRECISION, DIMENSION(1:6)               :: Fields
DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)         :: pcl
DOUBLE PRECISION                               :: x,y,z,ux,uy,uz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bxg,Byg,Bzg
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Exg,Eyg,Ezg

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

! Interpolation of the fields E,B at the location of the particle at t=0
CALL TRILINEAR_FIELDS(xgp,ygp,zgp,x,y,z, Fields, Bxg,Byg,Bzg,Exg,Eyg,Ezg)

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

! Interpolation of the fields E,B at the location of the particle at t=0
CALL TRILINEAR_FIELDS(xgp,ygp,zgp,x,y,z, Fields, Bxg,Byg,Bzg,Exg,Eyg,Ezg)

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

END IF

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
! - NPP: Number of particles per process
! - Esyn: Total synchrotron energy losses between t and t+dt
! - Eics: Total inverse Compton energy losses between t and t+dt
! - Uph: External photon energy density
!
! OUTPUT: ux,uy,uz at t=t+dt/2
!***********************************************************************

SUBROUTINE BORIS_PUSH(q,mass,pcl,pcl_data,Bxg,Byg,Bzg,Exg,Eyg,Ezg,Uph,&
                      xgp,ygp,zgp,Esyn,Eics,NPP)

IMPLICIT NONE

INTEGER*8 :: ip,NPP
DOUBLE PRECISION :: q,mass,gam,Bxi,Byi,Bzi,Exi,Eyi,Ezi,Psyn,Pics,Esyn,Eics,Uph
DOUBLE PRECISION :: uxm,uym,uzm,uxp,uyp,uzp,ux0,uy0,uz0,uxL,uyL,uzL,uxt,uyt,uzt
DOUBLE PRECISION :: sx,sy,sz,tx,ty,tz
DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)         :: pcl
DOUBLE PRECISION, DIMENSION(1:4,1:NPP)         :: pcl_data
DOUBLE PRECISION                               :: x,y,z,ux,uy,uz,wgt,El,Bp,Fr,Fe
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bxg,Byg,Bzg
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Exg,Eyg,Ezg
DOUBLE PRECISION, DIMENSION(1:6)               :: Fields

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
CALL TRILINEAR_FIELDS(xgp,ygp,zgp,x,y,z, Fields, Bxg,Byg,Bzg,Exg,Eyg,Ezg)

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

! Total radiative power losses at t
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

CALL TRILINEAR_FIELDS(xgp,ygp,zgp,x,y,z, Fields, Bxg,Byg,Bzg,Exg,Eyg,Ezg)

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
! - xminp,xmaxp,yminp,ymaxp,zminp,zmaxp: spatial boundaries of the domain
! - NESC: Number of particles escaping the domain towards neighbors
!
! OUTPUT: NESC at time t+dt
!
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE COUNT_ESCAPE(pcl,xminp,xmaxp,yminp,ymaxp,zminp,zmaxp,NPP,NESC)

IMPLICIT NONE

INTEGER*8                              :: ip,NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP) :: pcl
DOUBLE PRECISION       :: x,y,z,xminp,xmaxp,yminp,ymaxp,zminp,zmaxp
INTEGER, DIMENSION(26) :: NESC

!***********************************************************************

NESC=0

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)
z=pcl(3,ip)

! Case 1: the particle escapes Back-North
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.GT.ymaxp.AND.z.LT.zminp) THEN
NESC(1)=NESC(1)+1
END IF

! Case 2: the particle escapes Back-East
IF (y.LE.ymaxp.AND.y.GE.yminp.AND.x.GT.xmaxp.AND.z.LT.zminp) THEN
NESC(2)=NESC(2)+1
END IF

! Case 3: the particle escapes Back-South
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LT.yminp.AND.z.LT.zminp) THEN
NESC(3)=NESC(3)+1
END IF

! Case 4: the particle escapes Back-West
IF (y.LE.ymaxp.AND.y.GE.yminp.AND.x.LT.xminp.AND.z.LT.zminp) THEN
NESC(4)=NESC(4)+1
END IF

! Case 5: the particle escapes Back-North-East
IF (x.GT.xmaxp.AND.y.GT.ymaxp.AND.z.LT.zminp) THEN
NESC(5)=NESC(5)+1
END IF

! Case 6: the particle escapes Back-South-East
IF (x.GT.xmaxp.AND.y.LT.yminp.AND.z.LT.zminp) THEN
NESC(6)=NESC(6)+1
END IF

! Case 7: the particle escapes Back-South-West
IF (x.LT.xminp.AND.y.LT.yminp.AND.z.LT.zminp) THEN
NESC(7)=NESC(7)+1
END IF

! Case 8: the particle escapes Back-North-West
IF (x.LT.xminp.AND.y.GT.ymaxp.AND.z.LT.zminp) THEN
NESC(8)=NESC(8)+1
END IF

! Case 9: the particle escapes Back Center
IF (x.GE.xminp.AND.x.LE.xmaxp.AND.y.GE.yminp.AND.y.LE.ymaxp.AND.z.LT.zminp) THEN
NESC(9)=NESC(9)+1
END IF

! Case 10: the particle escapes north
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.GT.ymaxp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
NESC(10)=NESC(10)+1
END IF

! Case 11: the particle escapes East
IF (x.GT.xmaxp.AND.y.LE.ymaxp.AND.y.GE.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
NESC(11)=NESC(11)+1
END IF

! Case 12: the particle escapes South
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LT.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
NESC(12)=NESC(12)+1
END IF

! Case 13: the particle escapes West
IF (x.LT.xminp.AND.y.LE.ymaxp.AND.y.GE.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
NESC(13)=NESC(13)+1
END IF

! Case 14: the particle escapes North-East
IF (x.GT.xmaxp.AND.y.GT.ymaxp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
NESC(14)=NESC(14)+1
END IF

! Case 15: the particle escapes South-East
IF (x.GT.xmaxp.AND.y.LT.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
NESC(15)=NESC(15)+1
END IF

! Case 16: the particle escapes South-West
IF (x.LT.xminp.AND.y.LT.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
NESC(16)=NESC(16)+1
END IF

! Case 17: the particle escapes North-West
IF (x.LT.xminp.AND.y.GT.ymaxp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
NESC(17)=NESC(17)+1
END IF

! Case 18: the particle escapes Front-North
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.GT.ymaxp.AND.z.GT.zmaxp) THEN
NESC(18)=NESC(18)+1
END IF

! Case 19: the particle escapes Front-East
IF (y.LE.ymaxp.AND.y.GE.yminp.AND.x.GT.xmaxp.AND.z.GT.zmaxp) THEN
NESC(19)=NESC(19)+1
END IF

! Case 20: the particle escapes Front-South
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LT.yminp.AND.z.GT.zmaxp) THEN
NESC(20)=NESC(20)+1
END IF

! Case 21: the particle escapes Front-West
IF (y.LE.ymaxp.AND.y.GE.yminp.AND.x.LT.xminp.AND.z.GT.zmaxp) THEN
NESC(21)=NESC(21)+1
END IF

! Case 22: the particle escapes Front-North-East
IF (x.GT.xmaxp.AND.y.GT.ymaxp.AND.z.GT.zmaxp) THEN
NESC(22)=NESC(22)+1
END IF

! Case 23: the particle escapes Front-South-East
IF (x.GT.xmaxp.AND.y.LT.yminp.AND.z.GT.zmaxp) THEN
NESC(23)=NESC(23)+1
END IF

! Case 24: the particle escapes Front-South-West
IF (x.LT.xminp.AND.y.LT.yminp.AND.z.GT.zmaxp) THEN
NESC(24)=NESC(24)+1
END IF

! Case 25: the particle escapes Front-North-West
IF (x.LT.xminp.AND.y.GT.ymaxp.AND.z.GT.zmaxp) THEN
NESC(25)=NESC(25)+1
END IF

! Case 26: the particle escapes Front Center
IF (x.GE.xminp.AND.x.LE.xmaxp.AND.y.GE.yminp.AND.y.LE.ymaxp.AND.z.GT.zmaxp) THEN
NESC(26)=NESC(26)+1
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
DOUBLE PRECISION              :: x,y,z,ux,uy,uz,wt

!***********************************************************************

NTEMP=0

DO ip=1,NPP

x=pcl(1,ip)
y=pcl(2,ip)
z=pcl(3,ip)
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

!***********************************************************************
! Case 5: z>zmax
!***********************************************************************
IF (z.GT.zmax) THEN

   ! Elastic reflection
   IF (BOUND_PART_ZMAX.EQ."REFLECT") THEN
   z=2.0*zmax-z
   uz=-uz
   END IF
   
   ! Absorption
   IF (BOUND_PART_ZMAX.EQ."ABSORB") THEN
   wt=0d0
   END IF

END IF

!***********************************************************************
! Case 6: z<zmin
!***********************************************************************
IF (z.LT.zmin) THEN

   ! Elastic reflection
   IF (BOUND_PART_ZMIN.EQ."REFLECT") THEN
   z=2.0*zmin-z
   uz=-uz
   END IF
   
   ! Absorption
   IF (BOUND_PART_ZMIN.EQ."ABSORB") THEN
   wt=0d0
   END IF

END IF

pcl(1,ip)=x
pcl(2,ip)=y
pcl(3,ip)=z
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
! - NPI: number of particles per process at time t
! - NESC: Number of particles escaping the domain towards neighbors
! - id: process rank
! - ngh: neighbor array
! - COMM: communicator
! - ierr: error code
!
! OUTPUT: Particles in the domain at time t+dt
!
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE COM_PARTICLES(pcl,pcl_data,tagi,xminp,xmaxp,yminp,ymaxp,zminp,zmaxp,NPI,&
                         NESC,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)    :: stat
DOUBLE PRECISION, ALLOCATABLE          :: pcl(:,:)
DOUBLE PRECISION, ALLOCATABLE          :: pcl_data(:,:)
INTEGER*8, ALLOCATABLE                 :: tagi(:)
DOUBLE PRECISION                       :: x,y,z
DOUBLE PRECISION                       :: xminp,xmaxp,yminp,ymaxp,zminp,zmaxp
INTEGER, DIMENSION(26)                 :: ngh,NESC,NINC
INTEGER                                :: id,COMM,ierr

INTEGER, PARAMETER                     :: tag1=1,tag2=2,tag3=3,tag4=4
INTEGER, PARAMETER                     :: tag5=5,tag6=6,tag7=7,tag8=8
INTEGER, PARAMETER                     :: tag9=9,tag10=10,tag11=11,tag12=12
INTEGER, PARAMETER                     :: tag13=13,tag14=14,tag15=15,tag16=16
INTEGER, PARAMETER                     :: tag17=17,tag18=18,tag19=19,tag20=20
INTEGER, PARAMETER                     :: tag21=21,tag22=22,tag23=23,tag24=24
INTEGER, PARAMETER                     :: tag25=25,tag26=26

DOUBLE PRECISION, ALLOCATABLE :: bufSBN(:,:),bufSBE(:,:),bufSBS(:,:),bufSBW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSBNE(:,:),bufSBSE(:,:),bufSBSW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSBNW(:,:),bufSB(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSN(:,:),bufSE(:,:),bufSS(:,:),bufSW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSNE(:,:),bufSSE(:,:),bufSSW(:,:),bufSNW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSFN(:,:),bufSFE(:,:),bufSFS(:,:),bufSFW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSFNE(:,:),bufSFSE(:,:),bufSFSW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufSFNW(:,:),bufSF(:,:)

DOUBLE PRECISION, ALLOCATABLE :: bufRBN(:,:),bufRBE(:,:),bufRBS(:,:),bufRBW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRBNE(:,:),bufRBSE(:,:),bufRBSW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRBNW(:,:),bufRB(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRN(:,:),bufRE(:,:),bufRS(:,:),bufRW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRNE(:,:),bufRSE(:,:),bufRSW(:,:),bufRNW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRFN(:,:),bufRFE(:,:),bufRFS(:,:),bufRFW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRFNE(:,:),bufRFSE(:,:),bufRFSW(:,:)
DOUBLE PRECISION, ALLOCATABLE :: bufRFNW(:,:),bufRF(:,:)

INTEGER*8, ALLOCATABLE :: buf_tagSBN(:),buf_tagSBE(:),buf_tagSBS(:),buf_tagSBW(:)
INTEGER*8, ALLOCATABLE :: buf_tagSBNE(:),buf_tagSBSE(:),buf_tagSBSW(:)
INTEGER*8, ALLOCATABLE :: buf_tagSBNW(:),buf_tagSB(:)
INTEGER*8, ALLOCATABLE :: buf_tagSN(:),buf_tagSE(:),buf_tagSS(:),buf_tagSW(:)
INTEGER*8, ALLOCATABLE :: buf_tagSNE(:),buf_tagSSE(:),buf_tagSSW(:),buf_tagSNW(:)
INTEGER*8, ALLOCATABLE :: buf_tagSFN(:),buf_tagSFE(:),buf_tagSFS(:),buf_tagSFW(:)
INTEGER*8, ALLOCATABLE :: buf_tagSFNE(:),buf_tagSFSE(:),buf_tagSFSW(:)
INTEGER*8, ALLOCATABLE :: buf_tagSFNW(:),buf_tagSF(:)

INTEGER*8, ALLOCATABLE :: buf_tagRBN(:),buf_tagRBE(:),buf_tagRBS(:),buf_tagRBW(:)
INTEGER*8, ALLOCATABLE :: buf_tagRBNE(:),buf_tagRBSE(:),buf_tagRBSW(:)
INTEGER*8, ALLOCATABLE :: buf_tagRBNW(:),buf_tagRB(:)
INTEGER*8, ALLOCATABLE :: buf_tagRN(:),buf_tagRE(:),buf_tagRS(:),buf_tagRW(:)
INTEGER*8, ALLOCATABLE :: buf_tagRNE(:),buf_tagRSE(:),buf_tagRSW(:),buf_tagRNW(:)
INTEGER*8, ALLOCATABLE :: buf_tagRFN(:),buf_tagRFE(:),buf_tagRFS(:),buf_tagRFW(:)
INTEGER*8, ALLOCATABLE :: buf_tagRFNE(:),buf_tagRFSE(:),buf_tagRFSW(:)
INTEGER*8, ALLOCATABLE :: buf_tagRFNW(:),buf_tagRF(:)

INTEGER*8 :: NPI,NPO,NPT,NPT2,Ntemp
INTEGER*8 :: ip,ipo
INTEGER*8 :: i1,i2,i3,i4,i5,i6,i7,i8
INTEGER*8 :: i9,i10,i11,i12,i13,i14,i15,i16
INTEGER*8 :: i17,i18,i19,i20,i21,i22,i23,i24
INTEGER*8 :: i25,i26

!***********************************************************************

i1=0
i2=0
i3=0
i4=0
i5=0
i6=0
i7=0
i8=0
i9=0
i10=0
i11=0
i12=0
i13=0
i14=0
i15=0
i16=0
i17=0
i18=0
i19=0
i20=0
i21=0
i22=0
i23=0
i24=0
i25=0
i26=0

NPT2=NPI

ALLOCATE(bufSBN(NESC(1),11),bufSBE(NESC(2),11),bufSBS(NESC(3),11))
ALLOCATE(bufSBW(NESC(4),11))
ALLOCATE(bufSBNE(NESC(5),11),bufSBSE(NESC(6),11),bufSBSW(NESC(7),11))
ALLOCATE(bufSBNW(NESC(8),11),bufSB(NESC(9),11),bufSN(NESC(10),11))
ALLOCATE(bufSE(NESC(11),11))
ALLOCATE(bufSS(NESC(12),11),bufSW(NESC(13),11),bufSNE(NESC(14),11))
ALLOCATE(bufSSE(NESC(15),11))
ALLOCATE(bufSSW(NESC(16),11),bufSNW(NESC(17),11),bufSFN(NESC(18),11))
ALLOCATE(bufSFE(NESC(19),11),bufSFS(NESC(20),11),bufSFW(NESC(21),11))
ALLOCATE(bufSFNE(NESC(22),11),bufSFSE(NESC(23),11),bufSFSW(NESC(24),11))
ALLOCATE(bufSFNW(NESC(25),11),bufSF(NESC(26),11))

ALLOCATE(buf_tagSBN(NESC(1)),buf_tagSBE(NESC(2)),buf_tagSBS(NESC(3)))
ALLOCATE(buf_tagSBW(NESC(4)),buf_tagSBNE(NESC(5)),buf_tagSBSE(NESC(6)))
ALLOCATE(buf_tagSBSW(NESC(7)),buf_tagSBNW(NESC(8)),buf_tagSB(NESC(9)))
ALLOCATE(buf_tagSN(NESC(10)),buf_tagSE(NESC(11)),buf_tagSS(NESC(12)))
ALLOCATE(buf_tagSW(NESC(13)),buf_tagSNE(NESC(14)),buf_tagSSE(NESC(15)))
ALLOCATE(buf_tagSSW(NESC(16)),buf_tagSNW(NESC(17)),buf_tagSFN(NESC(18)))
ALLOCATE(buf_tagSFE(NESC(19)),buf_tagSFS(NESC(20)),buf_tagSFW(NESC(21)))
ALLOCATE(buf_tagSFNE(NESC(22)),buf_tagSFSE(NESC(23)),buf_tagSFSW(NESC(24)))
ALLOCATE(buf_tagSFNW(NESC(25)),buf_tagSF(NESC(26)))

!***********************************************************************
! ESCAPING PARTICLE ARRAY TO BE SEND
!***********************************************************************

DO ip=1,NPI

x=pcl(1,ip)
y=pcl(2,ip)
z=pcl(3,ip)

! Case 1: the particle escapes Back-North
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.GT.ymaxp.AND.z.LT.zminp) THEN
i1=i1+1

bufSBN(i1,1)=pcl(1,ip)
bufSBN(i1,2)=pcl(2,ip)
bufSBN(i1,3)=pcl(3,ip)
bufSBN(i1,4)=pcl(4,ip)
bufSBN(i1,5)=pcl(5,ip)
bufSBN(i1,6)=pcl(6,ip)
bufSBN(i1,7)=pcl(7,ip)
bufSBN(i1,8)=pcl_data(1,ip)
bufSBN(i1,9)=pcl_data(2,ip)
bufSBN(i1,10)=pcl_data(3,ip)
bufSBN(i1,11)=pcl_data(4,ip)
buf_tagSBN(i1)=tagi(ip)

  IF (bufSBN(i1,2).GT.ymax) THEN
  bufSBN(i1,2)=bufSBN(i1,2)-(ymax-ymin)
  END IF

  IF (bufSBN(i1,3).LT.zmin) THEN
  bufSBN(i1,3)=bufSBN(i1,3)+(zmax-zmin)
  END IF

END IF

! Case 2: the particle escapes Back-East
IF (y.LE.ymaxp.AND.y.GE.yminp.AND.x.GT.xmaxp.AND.z.LT.zminp) THEN
i2=i2+1

bufSBE(i2,1)=pcl(1,ip)
bufSBE(i2,2)=pcl(2,ip)
bufSBE(i2,3)=pcl(3,ip)
bufSBE(i2,4)=pcl(4,ip)
bufSBE(i2,5)=pcl(5,ip)
bufSBE(i2,6)=pcl(6,ip)
bufSBE(i2,7)=pcl(7,ip)
bufSBE(i2,8)=pcl_data(1,ip)
bufSBE(i2,9)=pcl_data(2,ip)
bufSBE(i2,10)=pcl_data(3,ip)
bufSBE(i2,11)=pcl_data(4,ip)
buf_tagSBE(i2)=tagi(ip)

  IF (bufSBE(i2,1).GT.xmax) THEN
  bufSBE(i2,1)=bufSBE(i2,1)-(xmax-xmin)
  END IF

  IF (bufSBE(i2,3).LT.zmin) THEN
  bufSBE(i2,3)=bufSBE(i2,3)+(zmax-zmin)
  END IF

END IF

! Case 3: the particle escapes Back-South
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LT.yminp.AND.z.LT.zminp) THEN
i3=i3+1

bufSBS(i3,1)=pcl(1,ip)
bufSBS(i3,2)=pcl(2,ip)
bufSBS(i3,3)=pcl(3,ip)
bufSBS(i3,4)=pcl(4,ip)
bufSBS(i3,5)=pcl(5,ip)
bufSBS(i3,6)=pcl(6,ip)
bufSBS(i3,7)=pcl(7,ip)
bufSBS(i3,8)=pcl_data(1,ip)
bufSBS(i3,9)=pcl_data(2,ip)
bufSBS(i3,10)=pcl_data(3,ip)
bufSBS(i3,11)=pcl_data(4,ip)
buf_tagSBS(i3)=tagi(ip)

  IF (bufSBS(i3,2).LT.ymin) THEN
  bufSBS(i3,2)=bufSBS(i3,2)+(ymax-ymin)
  END IF

  IF (bufSBS(i3,3).LT.zmin) THEN
  bufSBS(i3,3)=bufSBS(i3,3)+(zmax-zmin)
  END IF

END IF

! Case 4: the particle escapes Back-West
IF (y.LE.ymaxp.AND.y.GE.yminp.AND.x.LT.xminp.AND.z.LT.zminp) THEN
i4=i4+1

bufSBW(i4,1)=pcl(1,ip)
bufSBW(i4,2)=pcl(2,ip)
bufSBW(i4,3)=pcl(3,ip)
bufSBW(i4,4)=pcl(4,ip)
bufSBW(i4,5)=pcl(5,ip)
bufSBW(i4,6)=pcl(6,ip)
bufSBW(i4,7)=pcl(7,ip)
bufSBW(i4,8)=pcl_data(1,ip)
bufSBW(i4,9)=pcl_data(2,ip)
bufSBW(i4,10)=pcl_data(3,ip)
bufSBW(i4,11)=pcl_data(4,ip)
buf_tagSBW(i4)=tagi(ip)

  IF (bufSBW(i4,1).LT.xmin) THEN
  bufSBW(i4,1)=bufSBW(i4,1)+(xmax-xmin)
  END IF

  IF (bufSBW(i4,3).LT.zmin) THEN
  bufSBW(i4,3)=bufSBW(i4,3)+(zmax-zmin)
  END IF

END IF

! Case 5: the particle escapes Back-North-East
IF (x.GT.xmaxp.AND.y.GT.ymaxp.AND.z.LT.zminp) THEN
i5=i5+1

bufSBNE(i5,1)=pcl(1,ip)
bufSBNE(i5,2)=pcl(2,ip)
bufSBNE(i5,3)=pcl(3,ip)
bufSBNE(i5,4)=pcl(4,ip)
bufSBNE(i5,5)=pcl(5,ip)
bufSBNE(i5,6)=pcl(6,ip)
bufSBNE(i5,7)=pcl(7,ip)
bufSBNE(i5,8)=pcl_data(1,ip)
bufSBNE(i5,9)=pcl_data(2,ip)
bufSBNE(i5,10)=pcl_data(3,ip)
bufSBNE(i5,11)=pcl_data(4,ip)
buf_tagSBNE(i5)=tagi(ip)

  IF (bufSBNE(i5,1).GT.xmax) THEN
  bufSBNE(i5,1)=bufSBNE(i5,1)-(xmax-xmin)
  END IF

  IF (bufSBNE(i5,2).GT.ymax) THEN
  bufSBNE(i5,2)=bufSBNE(i5,2)-(ymax-ymin)
  END IF

  IF (bufSBNE(i5,3).LT.zmin) THEN
  bufSBNE(i5,3)=bufSBNE(i5,3)+(zmax-zmin)
  END IF

END IF

! Case 6: the particle escapes Back-South-East
IF (x.GT.xmaxp.AND.y.LT.yminp.AND.z.LT.zminp) THEN
i6=i6+1

bufSBSE(i6,1)=pcl(1,ip)
bufSBSE(i6,2)=pcl(2,ip)
bufSBSE(i6,3)=pcl(3,ip)
bufSBSE(i6,4)=pcl(4,ip)
bufSBSE(i6,5)=pcl(5,ip)
bufSBSE(i6,6)=pcl(6,ip)
bufSBSE(i6,7)=pcl(7,ip)
bufSBSE(i6,8)=pcl_data(1,ip)
bufSBSE(i6,9)=pcl_data(2,ip)
bufSBSE(i6,10)=pcl_data(3,ip)
bufSBSE(i6,11)=pcl_data(4,ip)
buf_tagSBSE(i6)=tagi(ip)

  IF (bufSBSE(i6,1).GT.xmax) THEN
  bufSBSE(i6,1)=bufSBSE(i6,1)-(xmax-xmin)
  END IF

  IF (bufSBSE(i6,2).LT.ymin) THEN
  bufSBSE(i6,2)=bufSBSE(i6,2)+(ymax-ymin)
  END IF

  IF (bufSBSE(i6,3).LT.zmin) THEN
  bufSBSE(i6,3)=bufSBSE(i6,3)+(zmax-zmin)
  END IF

END IF

! Case 7: the particle escapes Back-South-West
IF (x.LT.xminp.AND.y.LT.yminp.AND.z.LT.zminp) THEN
i7=i7+1

bufSBSW(i7,1)=pcl(1,ip)
bufSBSW(i7,2)=pcl(2,ip)
bufSBSW(i7,3)=pcl(3,ip)
bufSBSW(i7,4)=pcl(4,ip)
bufSBSW(i7,5)=pcl(5,ip)
bufSBSW(i7,6)=pcl(6,ip)
bufSBSW(i7,7)=pcl(7,ip)
bufSBSW(i7,8)=pcl_data(1,ip)
bufSBSW(i7,9)=pcl_data(2,ip)
bufSBSW(i7,10)=pcl_data(3,ip)
bufSBSW(i7,11)=pcl_data(4,ip)
buf_tagSBSW(i7)=tagi(ip)

  IF (bufSBSW(i7,1).LT.xmin) THEN
  bufSBSW(i7,1)=bufSBSW(i7,1)+(xmax-xmin)
  END IF

  IF (bufSBSW(i7,2).LT.ymin) THEN
  bufSBSW(i7,2)=bufSBSW(i7,2)+(ymax-ymin)
  END IF

  IF (bufSBSW(i7,3).LT.zmin) THEN
  bufSBSW(i7,3)=bufSBSW(i7,3)+(zmax-zmin)
  END IF

END IF

! Case 8: the particle escapes Back-North-West
IF (x.LT.xminp.AND.y.GT.ymaxp.AND.z.LT.zminp) THEN
i8=i8+1

bufSBNW(i8,1)=pcl(1,ip)
bufSBNW(i8,2)=pcl(2,ip)
bufSBNW(i8,3)=pcl(3,ip)
bufSBNW(i8,4)=pcl(4,ip)
bufSBNW(i8,5)=pcl(5,ip)
bufSBNW(i8,6)=pcl(6,ip)
bufSBNW(i8,7)=pcl(7,ip)
bufSBNW(i8,8)=pcl_data(1,ip)
bufSBNW(i8,9)=pcl_data(2,ip)
bufSBNW(i8,10)=pcl_data(3,ip)
bufSBNW(i8,11)=pcl_data(4,ip)
buf_tagSBNW(i8)=tagi(ip)

  IF (bufSBNW(i8,1).LT.xmin) THEN
  bufSBNW(i8,1)=bufSBNW(i8,1)+(xmax-xmin)
  END IF

  IF (bufSBNW(i8,2).GT.ymax) THEN
  bufSBNW(i8,2)=bufSBNW(i8,2)-(ymax-ymin)
  END IF

  IF (bufSBNW(i8,3).LT.zmin) THEN
  bufSBNW(i8,3)=bufSBNW(i8,3)+(zmax-zmin)
  END IF

END IF

! Case 9: the particle escapes Back Center
IF (x.GE.xminp.AND.x.LE.xmaxp.AND.y.GE.yminp.AND.y.LE.ymaxp.AND.z.LT.zminp) THEN
i9=i9+1

bufSB(i9,1)=pcl(1,ip)
bufSB(i9,2)=pcl(2,ip)
bufSB(i9,3)=pcl(3,ip)
bufSB(i9,4)=pcl(4,ip)
bufSB(i9,5)=pcl(5,ip)
bufSB(i9,6)=pcl(6,ip)
bufSB(i9,7)=pcl(7,ip)
bufSB(i9,8)=pcl_data(1,ip)
bufSB(i9,9)=pcl_data(2,ip)
bufSB(i9,10)=pcl_data(3,ip)
bufSB(i9,11)=pcl_data(4,ip)
buf_tagSB(i9)=tagi(ip)

  IF (bufSB(i9,3).LT.zmin) THEN
  bufSB(i9,3)=bufSB(i9,3)+(zmax-zmin)
  END IF

END IF

! Case 10: the particle escapes north
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.GT.ymaxp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
i10=i10+1
bufSN(i10,1)=pcl(1,ip)
bufSN(i10,2)=pcl(2,ip)
bufSN(i10,3)=pcl(3,ip)
bufSN(i10,4)=pcl(4,ip)
bufSN(i10,5)=pcl(5,ip)
bufSN(i10,6)=pcl(6,ip)
bufSN(i10,7)=pcl(7,ip)
bufSN(i10,8)=pcl_data(1,ip)
bufSN(i10,9)=pcl_data(2,ip)
bufSN(i10,10)=pcl_data(3,ip)
bufSN(i10,11)=pcl_data(4,ip)
buf_tagSN(i10)=tagi(ip)

  IF (bufSN(i10,2).GT.ymax) THEN
  bufSN(i10,2)=bufSN(i10,2)-(ymax-ymin)
  END IF

END IF

! Case 11: the particle escapes East
IF (x.GT.xmaxp.AND.y.LE.ymaxp.AND.y.GE.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
i11=i11+1
bufSE(i11,1)=pcl(1,ip)
bufSE(i11,2)=pcl(2,ip)
bufSE(i11,3)=pcl(3,ip)
bufSE(i11,4)=pcl(4,ip)
bufSE(i11,5)=pcl(5,ip)
bufSE(i11,6)=pcl(6,ip)
bufSE(i11,7)=pcl(7,ip)
bufSE(i11,8)=pcl_data(1,ip)
bufSE(i11,9)=pcl_data(2,ip)
bufSE(i11,10)=pcl_data(3,ip)
bufSE(i11,11)=pcl_data(4,ip)
buf_tagSE(i11)=tagi(ip)

  IF (bufSE(i11,1).GT.xmax) THEN
  bufSE(i11,1)=bufSE(i11,1)-(xmax-xmin)
  END IF

END IF

! Case 12: the particle escapes South
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LT.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
i12=i12+1
bufSS(i12,1)=pcl(1,ip)
bufSS(i12,2)=pcl(2,ip)
bufSS(i12,3)=pcl(3,ip)
bufSS(i12,4)=pcl(4,ip)
bufSS(i12,5)=pcl(5,ip)
bufSS(i12,6)=pcl(6,ip)
bufSS(i12,7)=pcl(7,ip)
bufSS(i12,8)=pcl_data(1,ip)
bufSS(i12,9)=pcl_data(2,ip)
bufSS(i12,10)=pcl_data(3,ip)
bufSS(i12,11)=pcl_data(4,ip)
buf_tagSS(i12)=tagi(ip)

  IF (bufSS(i12,2).LT.ymin) THEN
  bufSS(i12,2)=bufSS(i12,2)+(ymax-ymin)
  END IF

END IF

! Case 13: the particle escapes West
IF (x.LT.xminp.AND.y.LE.ymaxp.AND.y.GE.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
i13=i13+1
bufSW(i13,1)=pcl(1,ip)
bufSW(i13,2)=pcl(2,ip)
bufSW(i13,3)=pcl(3,ip)
bufSW(i13,4)=pcl(4,ip)
bufSW(i13,5)=pcl(5,ip)
bufSW(i13,6)=pcl(6,ip)
bufSW(i13,7)=pcl(7,ip)
bufSW(i13,8)=pcl_data(1,ip)
bufSW(i13,9)=pcl_data(2,ip)
bufSW(i13,10)=pcl_data(3,ip)
bufSW(i13,11)=pcl_data(4,ip)
buf_tagSW(i13)=tagi(ip)

  IF (bufSW(i13,1).LT.xmin) THEN
  bufSW(i13,1)=bufSW(i13,1)+(xmax-xmin)
  END IF

END IF

! Case 14: the particle escapes North-East
IF (x.GT.xmaxp.AND.y.GT.ymaxp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
i14=i14+1
bufSNE(i14,1)=pcl(1,ip)
bufSNE(i14,2)=pcl(2,ip)
bufSNE(i14,3)=pcl(3,ip)
bufSNE(i14,4)=pcl(4,ip)
bufSNE(i14,5)=pcl(5,ip)
bufSNE(i14,6)=pcl(6,ip)
bufSNE(i14,7)=pcl(7,ip)
bufSNE(i14,8)=pcl_data(1,ip)
bufSNE(i14,9)=pcl_data(2,ip)
bufSNE(i14,10)=pcl_data(3,ip)
bufSNE(i14,11)=pcl_data(4,ip)
buf_tagSNE(i14)=tagi(ip)

  IF (bufSNE(i14,1).GT.xmax) THEN
  bufSNE(i14,1)=bufSNE(i14,1)-(xmax-xmin)
  END IF

  IF (bufSNE(i14,2).GT.ymax) THEN
  bufSNE(i14,2)=bufSNE(i14,2)-(ymax-ymin)
  END IF

END IF

! Case 15: the particle escapes South-East
IF (x.GT.xmaxp.AND.y.LT.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
i15=i15+1
bufSSE(i15,1)=pcl(1,ip)
bufSSE(i15,2)=pcl(2,ip)
bufSSE(i15,3)=pcl(3,ip)
bufSSE(i15,4)=pcl(4,ip)
bufSSE(i15,5)=pcl(5,ip)
bufSSE(i15,6)=pcl(6,ip)
bufSSE(i15,7)=pcl(7,ip)
bufSSE(i15,8)=pcl_data(1,ip)
bufSSE(i15,9)=pcl_data(2,ip)
bufSSE(i15,10)=pcl_data(3,ip)
bufSSE(i15,11)=pcl_data(4,ip)
buf_tagSSE(i15)=tagi(ip)

  IF (bufSSE(i15,1).GT.xmax) THEN
  bufSSE(i15,1)=bufSSE(i15,1)-(xmax-xmin)
  END IF

  IF (bufSSE(i15,2).LT.ymin) THEN
  bufSSE(i15,2)=bufSSE(i15,2)+(ymax-ymin)
  END IF

END IF

! Case 16: the particle escapes South-West
IF (x.LT.xminp.AND.y.LT.yminp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
i16=i16+1
bufSSW(i16,1)=pcl(1,ip)
bufSSW(i16,2)=pcl(2,ip)
bufSSW(i16,3)=pcl(3,ip)
bufSSW(i16,4)=pcl(4,ip)
bufSSW(i16,5)=pcl(5,ip)
bufSSW(i16,6)=pcl(6,ip)
bufSSW(i16,7)=pcl(7,ip)
bufSSW(i16,8)=pcl_data(1,ip)
bufSSW(i16,9)=pcl_data(2,ip)
bufSSW(i16,10)=pcl_data(3,ip)
bufSSW(i16,11)=pcl_data(4,ip)
buf_tagSSW(i16)=tagi(ip)

  IF (bufSSW(i16,1).LT.xmin) THEN
  bufSSW(i16,1)=bufSSW(i16,1)+(xmax-xmin)
  END IF

  IF (bufSSW(i16,2).LT.ymin) THEN
  bufSSW(i16,2)=bufSSW(i16,2)+(ymax-ymin)
  END IF

END IF

! Case 17: the particle escapes North-West
IF (x.LT.xminp.AND.y.GT.ymaxp.AND.z.GE.zminp.AND.z.LE.zmaxp) THEN
i17=i17+1
bufSNW(i17,1)=pcl(1,ip)
bufSNW(i17,2)=pcl(2,ip)
bufSNW(i17,3)=pcl(3,ip)
bufSNW(i17,4)=pcl(4,ip)
bufSNW(i17,5)=pcl(5,ip)
bufSNW(i17,6)=pcl(6,ip)
bufSNW(i17,7)=pcl(7,ip)
bufSNW(i17,8)=pcl_data(1,ip)
bufSNW(i17,9)=pcl_data(2,ip)
bufSNW(i17,10)=pcl_data(3,ip)
bufSNW(i17,11)=pcl_data(4,ip)
buf_tagSNW(i17)=tagi(ip)

  IF (bufSNW(i17,1).LT.xmin) THEN
  bufSNW(i17,1)=bufSNW(i17,1)+(xmax-xmin)
  END IF

  IF (bufSNW(i17,2).GT.ymax) THEN
  bufSNW(i17,2)=bufSNW(i17,2)-(ymax-ymin)
  END IF

END IF

! Case 18: the particle escapes Front-North
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.GT.ymaxp.AND.z.GT.zmaxp) THEN
i18=i18+1

bufSFN(i18,1)=pcl(1,ip)
bufSFN(i18,2)=pcl(2,ip)
bufSFN(i18,3)=pcl(3,ip)
bufSFN(i18,4)=pcl(4,ip)
bufSFN(i18,5)=pcl(5,ip)
bufSFN(i18,6)=pcl(6,ip)
bufSFN(i18,7)=pcl(7,ip)
bufSFN(i18,8)=pcl_data(1,ip)
bufSFN(i18,9)=pcl_data(2,ip)
bufSFN(i18,10)=pcl_data(3,ip)
bufSFN(i18,11)=pcl_data(4,ip)
buf_tagSFN(i18)=tagi(ip)

  IF (bufSFN(i18,2).GT.ymax) THEN
  bufSFN(i18,2)=bufSFN(i18,2)-(ymax-ymin)
  END IF

  IF (bufSFN(i18,3).GT.zmax) THEN
  bufSFN(i18,3)=bufSFN(i18,3)-(zmax-zmin)
  END IF

END IF

! Case 19: the particle escapes Front-East
IF (y.LE.ymaxp.AND.y.GE.yminp.AND.x.GT.xmaxp.AND.z.GT.zmaxp) THEN
i19=i19+1

bufSFE(i19,1)=pcl(1,ip)
bufSFE(i19,2)=pcl(2,ip)
bufSFE(i19,3)=pcl(3,ip)
bufSFE(i19,4)=pcl(4,ip)
bufSFE(i19,5)=pcl(5,ip)
bufSFE(i19,6)=pcl(6,ip)
bufSFE(i19,7)=pcl(7,ip)
bufSFE(i19,8)=pcl_data(1,ip)
bufSFE(i19,9)=pcl_data(2,ip)
bufSFE(i19,10)=pcl_data(3,ip)
bufSFE(i19,11)=pcl_data(4,ip)
buf_tagSFE(i19)=tagi(ip)

  IF (bufSFE(i19,1).GT.xmax) THEN
  bufSFE(i19,1)=bufSFE(i19,1)-(xmax-xmin)
  END IF

  IF (bufSFE(i19,3).GT.zmax) THEN
  bufSFE(i19,3)=bufSFE(i19,3)-(zmax-zmin)
  END IF

END IF

! Case 20: the particle escapes Front-South
IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LT.yminp.AND.z.GT.zmaxp) THEN
i20=i20+1

bufSFS(i20,1)=pcl(1,ip)
bufSFS(i20,2)=pcl(2,ip)
bufSFS(i20,3)=pcl(3,ip)
bufSFS(i20,4)=pcl(4,ip)
bufSFS(i20,5)=pcl(5,ip)
bufSFS(i20,6)=pcl(6,ip)
bufSFS(i20,7)=pcl(7,ip)
bufSFS(i20,8)=pcl_data(1,ip)
bufSFS(i20,9)=pcl_data(2,ip)
bufSFS(i20,10)=pcl_data(3,ip)
bufSFS(i20,11)=pcl_data(4,ip)
buf_tagSFS(i20)=tagi(ip)

  IF (bufSFS(i20,2).LT.ymin) THEN
  bufSFS(i20,2)=bufSFS(i20,2)+(ymax-ymin)
  END IF

  IF (bufSFS(i20,3).GT.zmax) THEN
  bufSFS(i20,3)=bufSFS(i20,3)-(zmax-zmin)
  END IF

END IF

! Case 21: the particle escapes Front-West
IF (y.LE.ymaxp.AND.y.GE.yminp.AND.x.LT.xminp.AND.z.GT.zmaxp) THEN
i21=i21+1

bufSFW(i21,1)=pcl(1,ip)
bufSFW(i21,2)=pcl(2,ip)
bufSFW(i21,3)=pcl(3,ip)
bufSFW(i21,4)=pcl(4,ip)
bufSFW(i21,5)=pcl(5,ip)
bufSFW(i21,6)=pcl(6,ip)
bufSFW(i21,7)=pcl(7,ip)
bufSFW(i21,8)=pcl_data(1,ip)
bufSFW(i21,9)=pcl_data(2,ip)
bufSFW(i21,10)=pcl_data(3,ip)
bufSFW(i21,11)=pcl_data(4,ip)
buf_tagSFW(i21)=tagi(ip)

  IF (bufSFW(i21,1).LT.xmin) THEN
  bufSFW(i21,1)=bufSFW(i21,1)+(xmax-xmin)
  END IF

  IF (bufSFW(i21,3).GT.zmax) THEN
  bufSFW(i21,3)=bufSFW(i21,3)-(zmax-zmin)
  END IF

END IF

! Case 22: the particle escapes Front-North-East
IF (x.GT.xmaxp.AND.y.GT.ymaxp.AND.z.GT.zmaxp) THEN
i22=i22+1

bufSFNE(i22,1)=pcl(1,ip)
bufSFNE(i22,2)=pcl(2,ip)
bufSFNE(i22,3)=pcl(3,ip)
bufSFNE(i22,4)=pcl(4,ip)
bufSFNE(i22,5)=pcl(5,ip)
bufSFNE(i22,6)=pcl(6,ip)
bufSFNE(i22,7)=pcl(7,ip)
bufSFNE(i22,8)=pcl_data(1,ip)
bufSFNE(i22,9)=pcl_data(2,ip)
bufSFNE(i22,10)=pcl_data(3,ip)
bufSFNE(i22,11)=pcl_data(4,ip)
buf_tagSFNE(i22)=tagi(ip)

  IF (bufSFNE(i22,1).GT.xmax) THEN
  bufSFNE(i22,1)=bufSFNE(i22,1)-(xmax-xmin)
  END IF

  IF (bufSFNE(i22,2).GT.ymax) THEN
  bufSFNE(i22,2)=bufSFNE(i22,2)-(ymax-ymin)
  END IF

  IF (bufSFNE(i22,3).GT.zmax) THEN
  bufSFNE(i22,3)=bufSFNE(i22,3)-(zmax-zmin)
  END IF

END IF

! Case 23: the particle escapes Front-South-East
IF (x.GT.xmaxp.AND.y.LT.yminp.AND.z.GT.zmaxp) THEN
i23=i23+1

bufSFSE(i23,1)=pcl(1,ip)
bufSFSE(i23,2)=pcl(2,ip)
bufSFSE(i23,3)=pcl(3,ip)
bufSFSE(i23,4)=pcl(4,ip)
bufSFSE(i23,5)=pcl(5,ip)
bufSFSE(i23,6)=pcl(6,ip)
bufSFSE(i23,7)=pcl(7,ip)
bufSFSE(i23,8)=pcl_data(1,ip)
bufSFSE(i23,9)=pcl_data(2,ip)
bufSFSE(i23,10)=pcl_data(3,ip)
bufSFSE(i23,11)=pcl_data(4,ip)
buf_tagSFSE(i23)=tagi(ip)

  IF (bufSFSE(i23,1).GT.xmax) THEN
  bufSFSE(i23,1)=bufSFSE(i23,1)-(xmax-xmin)
  END IF

  IF (bufSFSE(i23,2).LT.ymin) THEN
  bufSFSE(i23,2)=bufSFSE(i23,2)+(ymax-ymin)
  END IF

  IF (bufSFSE(i23,3).GT.zmax) THEN
  bufSFSE(i23,3)=bufSFSE(i23,3)-(zmax-zmin)
  END IF

END IF

! Case 24: the particle escapes Front-South-West
IF (x.LT.xminp.AND.y.LT.yminp.AND.z.GT.zmaxp) THEN
i24=i24+1

bufSFSW(i24,1)=pcl(1,ip)
bufSFSW(i24,2)=pcl(2,ip)
bufSFSW(i24,3)=pcl(3,ip)
bufSFSW(i24,4)=pcl(4,ip)
bufSFSW(i24,5)=pcl(5,ip)
bufSFSW(i24,6)=pcl(6,ip)
bufSFSW(i24,7)=pcl(7,ip)
bufSFSW(i24,8)=pcl_data(1,ip)
bufSFSW(i24,9)=pcl_data(2,ip)
bufSFSW(i24,10)=pcl_data(3,ip)
bufSFSW(i24,11)=pcl_data(4,ip)
buf_tagSFSW(i24)=tagi(ip)

  IF (bufSFSW(i24,1).LT.xmin) THEN
  bufSFSW(i24,1)=bufSFSW(i24,1)+(xmax-xmin)
  END IF

  IF (bufSFSW(i24,2).LT.ymin) THEN
  bufSFSW(i24,2)=bufSFSW(i24,2)+(ymax-ymin)
  END IF

  IF (bufSFSW(i24,3).GT.zmax) THEN
  bufSFSW(i24,3)=bufSFSW(i24,3)-(zmax-zmin)
  END IF

END IF

! Case 25: the particle escapes Front-North-West
IF (x.LT.xminp.AND.y.GT.ymaxp.AND.z.GT.zmaxp) THEN
i25=i25+1

bufSFNW(i25,1)=pcl(1,ip)
bufSFNW(i25,2)=pcl(2,ip)
bufSFNW(i25,3)=pcl(3,ip)
bufSFNW(i25,4)=pcl(4,ip)
bufSFNW(i25,5)=pcl(5,ip)
bufSFNW(i25,6)=pcl(6,ip)
bufSFNW(i25,7)=pcl(7,ip)
bufSFNW(i25,8)=pcl_data(1,ip)
bufSFNW(i25,9)=pcl_data(2,ip)
bufSFNW(i25,10)=pcl_data(3,ip)
bufSFNW(i25,11)=pcl_data(4,ip)
buf_tagSFNW(i25)=tagi(ip)

  IF (bufSFNW(i25,1).LT.xmin) THEN
  bufSFNW(i25,1)=bufSFNW(i25,1)+(xmax-xmin)
  END IF

  IF (bufSFNW(i25,2).GT.ymax) THEN
  bufSFNW(i25,2)=bufSFNW(i25,2)-(ymax-ymin)
  END IF

  IF (bufSFNW(i25,3).GT.zmax) THEN
  bufSFNW(i25,3)=bufSFNW(i25,3)-(zmax-zmin)
  END IF

END IF

! Case 26: the particle escapes Front Center
IF (x.GE.xminp.AND.x.LE.xmaxp.AND.y.GE.yminp.AND.y.LE.ymaxp.AND.z.GT.zmaxp) THEN
i26=i26+1

bufSF(i26,1)=pcl(1,ip)
bufSF(i26,2)=pcl(2,ip)
bufSF(i26,3)=pcl(3,ip)
bufSF(i26,4)=pcl(4,ip)
bufSF(i26,5)=pcl(5,ip)
bufSF(i26,6)=pcl(6,ip)
bufSF(i26,7)=pcl(7,ip)
bufSF(i26,8)=pcl_data(1,ip)
bufSF(i26,9)=pcl_data(2,ip)
bufSF(i26,10)=pcl_data(3,ip)
bufSF(i26,11)=pcl_data(4,ip)
buf_tagSF(i26)=tagi(ip)

  IF (bufSF(i26,3).GT.zmax) THEN
  bufSF(i26,3)=bufSF(i26,3)-(zmax-zmin)
  END IF

END IF

ENDDO

!***********************************************************************
! EXCHANGE OF PARTICLES BETWEEN PROCESSES
!***********************************************************************

! Exchange sizes of the OUTGOING buffered arrays.

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(NESC(1),1,MPI_INTEGER,ngh(1),tag1,&
                  NINC(20),1,MPI_INTEGER,ngh(20),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(2),1,MPI_INTEGER,ngh(2),tag2,&
                  NINC(21),1,MPI_INTEGER,ngh(21),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(3),1,MPI_INTEGER,ngh(3),tag3,&
                  NINC(18),1,MPI_INTEGER,ngh(18),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(4),1,MPI_INTEGER,ngh(4),tag4,&
                  NINC(19),1,MPI_INTEGER,ngh(19),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(5),1,MPI_INTEGER,ngh(5),tag5,&
                  NINC(24),1,MPI_INTEGER,ngh(24),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(6),1,MPI_INTEGER,ngh(6),tag6,&
                  NINC(25),1,MPI_INTEGER,ngh(25),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(7),1,MPI_INTEGER,ngh(7),tag7,&
                  NINC(22),1,MPI_INTEGER,ngh(22),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(8),1,MPI_INTEGER,ngh(8),tag8,&
                  NINC(23),1,MPI_INTEGER,ngh(23),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(9),1,MPI_INTEGER,ngh(9),tag9,&
                  NINC(26),1,MPI_INTEGER,ngh(26),tag9,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(10),1,MPI_INTEGER,ngh(10),tag10,&
                  NINC(12),1,MPI_INTEGER,ngh(12),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(11),1,MPI_INTEGER,ngh(11),tag11,&
                  NINC(13),1,MPI_INTEGER,ngh(13),tag11,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(12),1,MPI_INTEGER,ngh(12),tag12,&
                  NINC(10),1,MPI_INTEGER,ngh(10),tag12,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(13),1,MPI_INTEGER,ngh(13),tag13,&
                  NINC(11),1,MPI_INTEGER,ngh(11),tag13,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(14),1,MPI_INTEGER,ngh(14),tag14,&
                  NINC(16),1,MPI_INTEGER,ngh(16),tag14,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(15),1,MPI_INTEGER,ngh(15),tag15,&
                  NINC(17),1,MPI_INTEGER,ngh(17),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(16),1,MPI_INTEGER,ngh(16),tag16,&
                  NINC(14),1,MPI_INTEGER,ngh(14),tag16,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(17),1,MPI_INTEGER,ngh(17),tag17,&
                  NINC(15),1,MPI_INTEGER,ngh(15),tag17,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(18),1,MPI_INTEGER,ngh(18),tag18,&
                  NINC(3),1,MPI_INTEGER,ngh(3),tag18,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(19),1,MPI_INTEGER,ngh(19),tag19,&
                  NINC(4),1,MPI_INTEGER,ngh(4),tag19,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(20),1,MPI_INTEGER,ngh(20),tag20,&
                  NINC(1),1,MPI_INTEGER,ngh(1),tag20,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(21),1,MPI_INTEGER,ngh(21),tag21,&
                  NINC(2),1,MPI_INTEGER,ngh(2),tag21,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(22),1,MPI_INTEGER,ngh(22),tag22,&
                  NINC(7),1,MPI_INTEGER,ngh(7),tag22,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(23),1,MPI_INTEGER,ngh(23),tag23,&
                  NINC(8),1,MPI_INTEGER,ngh(8),tag23,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(24),1,MPI_INTEGER,ngh(24),tag24,&
                  NINC(5),1,MPI_INTEGER,ngh(5),tag24,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(25),1,MPI_INTEGER,ngh(25),tag25,&
                  NINC(6),1,MPI_INTEGER,ngh(6),tag25,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(26),1,MPI_INTEGER,ngh(26),tag26,&
                  NINC(9),1,MPI_INTEGER,ngh(9),tag26,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(NESC(1),1,MPI_INTEGER,ngh(1),tag1,&
                  NINC(20),1,MPI_INTEGER,ngh(20),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(2),1,MPI_INTEGER,ngh(2),tag2,&
                  NINC(21),1,MPI_INTEGER,ngh(21),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(3),1,MPI_INTEGER,ngh(3),tag3,&
                  NINC(18),1,MPI_INTEGER,ngh(18),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(4),1,MPI_INTEGER,ngh(4),tag4,&
                  NINC(19),1,MPI_INTEGER,ngh(19),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(5),1,MPI_INTEGER,ngh(5),tag5,&
                  NINC(24),1,MPI_INTEGER,ngh(24),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(6),1,MPI_INTEGER,ngh(6),tag6,&
                  NINC(25),1,MPI_INTEGER,ngh(25),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(7),1,MPI_INTEGER,ngh(7),tag7,&
                  NINC(22),1,MPI_INTEGER,ngh(22),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(8),1,MPI_INTEGER,ngh(8),tag8,&
                  NINC(23),1,MPI_INTEGER,ngh(23),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(9),1,MPI_INTEGER,ngh(9),tag9,&
                  NINC(26),1,MPI_INTEGER,ngh(26),tag9,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(10),1,MPI_INTEGER,ngh(10),tag10,&
                  NINC(12),1,MPI_INTEGER,ngh(12),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(11),1,MPI_INTEGER,ngh(11),tag11,&
                  NINC(13),1,MPI_INTEGER,ngh(13),tag11,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(12),1,MPI_INTEGER,ngh(12),tag12,&
                  NINC(10),1,MPI_INTEGER,ngh(10),tag12,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(13),1,MPI_INTEGER,ngh(13),tag13,&
                  NINC(11),1,MPI_INTEGER,ngh(11),tag13,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(14),1,MPI_INTEGER,ngh(14),tag14,&
                  NINC(16),1,MPI_INTEGER,ngh(16),tag14,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(15),1,MPI_INTEGER,ngh(15),tag15,&
                  NINC(17),1,MPI_INTEGER,ngh(17),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(16),1,MPI_INTEGER,ngh(16),tag16,&
                  NINC(14),1,MPI_INTEGER,ngh(14),tag16,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(17),1,MPI_INTEGER,ngh(17),tag17,&
                  NINC(15),1,MPI_INTEGER,ngh(15),tag17,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(18),1,MPI_INTEGER,ngh(18),tag18,&
                  NINC(3),1,MPI_INTEGER,ngh(3),tag18,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(19),1,MPI_INTEGER,ngh(19),tag19,&
                  NINC(4),1,MPI_INTEGER,ngh(4),tag19,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(20),1,MPI_INTEGER,ngh(20),tag20,&
                  NINC(1),1,MPI_INTEGER,ngh(1),tag20,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(21),1,MPI_INTEGER,ngh(21),tag21,&
                  NINC(2),1,MPI_INTEGER,ngh(2),tag21,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(22),1,MPI_INTEGER,ngh(22),tag22,&
                  NINC(7),1,MPI_INTEGER,ngh(7),tag22,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(23),1,MPI_INTEGER,ngh(23),tag23,&
                  NINC(8),1,MPI_INTEGER,ngh(8),tag23,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(24),1,MPI_INTEGER,ngh(24),tag24,&
                  NINC(5),1,MPI_INTEGER,ngh(5),tag24,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(25),1,MPI_INTEGER,ngh(25),tag25,&
                  NINC(6),1,MPI_INTEGER,ngh(6),tag25,COMM,stat,ierr)

CALL MPI_SENDRECV(NESC(26),1,MPI_INTEGER,ngh(26),tag26,&
                  NINC(9),1,MPI_INTEGER,ngh(9),tag26,COMM,stat,ierr)

ENDIF

!***********************************************************************
! Total number of unchanged particles
NPT=NPI-(NESC(1)+NESC(2)+NESC(3)+NESC(4)+NESC(5)+NESC(6)+NESC(7)+NESC(8)+NESC(9)+&
NESC(10)+NESC(11)+NESC(12)+NESC(13)+NESC(14)+NESC(15)+NESC(16)+NESC(17)+&
NESC(18)+NESC(19)+NESC(20)+NESC(21)+NESC(22)+NESC(23)+NESC(24)+NESC(25)+NESC(26))

! Total number of particles (+incoming-outcoming)
NPO=NPT+(NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+&
NINC(18)+NINC(19)+NINC(20)+NINC(21)+NINC(22)+NINC(23)+NINC(24)+NINC(25)+NINC(26))
!***********************************************************************

! Exchange of the OUTCOMING and INCOMING buffered arrays.

ALLOCATE(bufRBN(NINC(1),11),bufRBE(NINC(2),11),bufRBS(NINC(3),11))
ALLOCATE(bufRBW(NINC(4),11))
ALLOCATE(bufRBNE(NINC(5),11),bufRBSE(NINC(6),11),bufRBSW(NINC(7),11))
ALLOCATE(bufRBNW(NINC(8),11),bufRB(NINC(9),11),bufRN(NINC(10),11))
ALLOCATE(bufRE(NINC(11),11))
ALLOCATE(bufRS(NINC(12),11),bufRW(NINC(13),11),bufRNE(NINC(14),11))
ALLOCATE(bufRSE(NINC(15),11))
ALLOCATE(bufRSW(NINC(16),11),bufRNW(NINC(17),11),bufRFN(NINC(18),11))
ALLOCATE(bufRFE(NINC(19),11),bufRFS(NINC(20),11),bufRFW(NINC(21),11))
ALLOCATE(bufRFNE(NINC(22),11),bufRFSE(NINC(23),11),bufRFSW(NINC(24),11))
ALLOCATE(bufRFNW(NINC(25),11),bufRF(NINC(26),11))

ALLOCATE(buf_tagRBN(NINC(1)),buf_tagRBE(NINC(2)),buf_tagRBS(NINC(3)))
ALLOCATE(buf_tagRBW(NINC(4)),buf_tagRBNE(NINC(5)),buf_tagRBSE(NINC(6)))
ALLOCATE(buf_tagRBSW(NINC(7)),buf_tagRBNW(NINC(8)),buf_tagRB(NINC(9)))
ALLOCATE(buf_tagRN(NINC(10)),buf_tagRE(NINC(11)),buf_tagRS(NINC(12)))
ALLOCATE(buf_tagRW(NINC(13)),buf_tagRNE(NINC(14)),buf_tagRSE(NINC(15)))
ALLOCATE(buf_tagRSW(NINC(16)),buf_tagRNW(NINC(17)),buf_tagRFN(NINC(18)))
ALLOCATE(buf_tagRFE(NINC(19)),buf_tagRFS(NINC(20)),buf_tagRFW(NINC(21)))
ALLOCATE(buf_tagRFNE(NINC(22)),buf_tagRFSE(NINC(23)),buf_tagRFSW(NINC(24)))
ALLOCATE(buf_tagRFNW(NINC(25)),buf_tagRF(NINC(26)))

! Exchange particle positions, velocity and weight

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSBN,NESC(1)*11,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
     bufRFS,NINC(20)*11,MPI_DOUBLE_PRECISION,ngh(20),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBE,NESC(2)*11,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
     bufRFW,NINC(21)*11,MPI_DOUBLE_PRECISION,ngh(21),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBS,NESC(3)*11,MPI_DOUBLE_PRECISION,ngh(3),tag3,&
     bufRFN,NINC(18)*11,MPI_DOUBLE_PRECISION,ngh(18),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBW,NESC(4)*11,MPI_DOUBLE_PRECISION,ngh(4),tag4,&
     bufRFE,NINC(19)*11,MPI_DOUBLE_PRECISION,ngh(19),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBNE,NESC(5)*11,MPI_DOUBLE_PRECISION,ngh(5),tag5,&
     bufRFSW,NINC(24)*11,MPI_DOUBLE_PRECISION,ngh(24),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBSE,NESC(6)*11,MPI_DOUBLE_PRECISION,ngh(6),tag6,&
     bufRFNW,NINC(25)*11,MPI_DOUBLE_PRECISION,ngh(25),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBSW,NESC(7)*11,MPI_DOUBLE_PRECISION,ngh(7),tag7,&
     bufRFNE,NINC(22)*11,MPI_DOUBLE_PRECISION,ngh(22),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBNW,NESC(8)*11,MPI_DOUBLE_PRECISION,ngh(8),tag8,&
     bufRFSE,NINC(23)*11,MPI_DOUBLE_PRECISION,ngh(23),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSB,NESC(9)*11,MPI_DOUBLE_PRECISION,ngh(9),tag9,&
     bufRF,NINC(26)*11,MPI_DOUBLE_PRECISION,ngh(26),tag9,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN,NESC(10)*11,MPI_DOUBLE_PRECISION,ngh(10),tag10,&
     bufRS,NINC(12)*11,MPI_DOUBLE_PRECISION,ngh(12),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSE,NESC(11)*11,MPI_DOUBLE_PRECISION,ngh(11),tag11,&
     bufRW,NINC(13)*11,MPI_DOUBLE_PRECISION,ngh(13),tag11,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS,NESC(12)*11,MPI_DOUBLE_PRECISION,ngh(12),tag12,&
     bufRN,NINC(10)*11,MPI_DOUBLE_PRECISION,ngh(10),tag12,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NESC(13)*11,MPI_DOUBLE_PRECISION,ngh(13),tag13,&
     bufRE,NINC(11)*11,MPI_DOUBLE_PRECISION,ngh(11),tag13,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNE,NESC(14)*11,MPI_DOUBLE_PRECISION,ngh(14),tag14,&
     bufRSW,NINC(16)*11,MPI_DOUBLE_PRECISION,ngh(16),tag14,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSE,NESC(15)*11,MPI_DOUBLE_PRECISION,ngh(15),tag15,&
     bufRNW,NINC(17)*11,MPI_DOUBLE_PRECISION,ngh(17),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSW,NESC(16)*11,MPI_DOUBLE_PRECISION,ngh(16),tag16,&
     bufRNE,NINC(14)*11,MPI_DOUBLE_PRECISION,ngh(14),tag16,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNW,NESC(17)*11,MPI_DOUBLE_PRECISION,ngh(17),tag17,&
     bufRSE,NINC(15)*11,MPI_DOUBLE_PRECISION,ngh(15),tag17,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFN,NESC(18)*11,MPI_DOUBLE_PRECISION,ngh(18),tag18,&
     bufRBS,NINC(3)*11,MPI_DOUBLE_PRECISION,ngh(3),tag18,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFE,NESC(19)*11,MPI_DOUBLE_PRECISION,ngh(19),tag19,&
     bufRBW,NINC(4)*11,MPI_DOUBLE_PRECISION,ngh(4),tag19,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFS,NESC(20)*11,MPI_DOUBLE_PRECISION,ngh(20),tag20,&
     bufRBN,NINC(1)*11,MPI_DOUBLE_PRECISION,ngh(1),tag20,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFW,NESC(21)*11,MPI_DOUBLE_PRECISION,ngh(21),tag21,&
     bufRBE,NINC(2)*11,MPI_DOUBLE_PRECISION,ngh(2),tag21,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFNE,NESC(22)*11,MPI_DOUBLE_PRECISION,ngh(22),tag22,&
     bufRBSW,NINC(7)*11,MPI_DOUBLE_PRECISION,ngh(7),tag22,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFSE,NESC(23)*11,MPI_DOUBLE_PRECISION,ngh(23),tag23,&
     bufRBNW,NINC(8)*11,MPI_DOUBLE_PRECISION,ngh(8),tag23,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFSW,NESC(24)*11,MPI_DOUBLE_PRECISION,ngh(24),tag24,&
     bufRBNE,NINC(5)*11,MPI_DOUBLE_PRECISION,ngh(5),tag24,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFNW,NESC(25)*11,MPI_DOUBLE_PRECISION,ngh(25),tag25,&
     bufRBSE,NINC(6)*11,MPI_DOUBLE_PRECISION,ngh(6),tag25,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF,NESC(26)*11,MPI_DOUBLE_PRECISION,ngh(26),tag26,&
     bufRB,NINC(9)*11,MPI_DOUBLE_PRECISION,ngh(9),tag26,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(bufSBN,NESC(1)*11,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
     bufRFS,NINC(20)*11,MPI_DOUBLE_PRECISION,ngh(20),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBE,NESC(2)*11,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
     bufRFW,NINC(21)*11,MPI_DOUBLE_PRECISION,ngh(21),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBS,NESC(3)*11,MPI_DOUBLE_PRECISION,ngh(3),tag3,&
     bufRFN,NINC(18)*11,MPI_DOUBLE_PRECISION,ngh(18),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBW,NESC(4)*11,MPI_DOUBLE_PRECISION,ngh(4),tag4,&
     bufRFE,NINC(19)*11,MPI_DOUBLE_PRECISION,ngh(19),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBNE,NESC(5)*11,MPI_DOUBLE_PRECISION,ngh(5),tag5,&
     bufRFSW,NINC(24)*11,MPI_DOUBLE_PRECISION,ngh(24),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBSE,NESC(6)*11,MPI_DOUBLE_PRECISION,ngh(6),tag6,&
     bufRFNW,NINC(25)*11,MPI_DOUBLE_PRECISION,ngh(25),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBSW,NESC(7)*11,MPI_DOUBLE_PRECISION,ngh(7),tag7,&
     bufRFNE,NINC(22)*11,MPI_DOUBLE_PRECISION,ngh(22),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBNW,NESC(8)*11,MPI_DOUBLE_PRECISION,ngh(8),tag8,&
     bufRFSE,NINC(23)*11,MPI_DOUBLE_PRECISION,ngh(23),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSB,NESC(9)*11,MPI_DOUBLE_PRECISION,ngh(9),tag9,&
     bufRF,NINC(26)*11,MPI_DOUBLE_PRECISION,ngh(26),tag9,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN,NESC(10)*11,MPI_DOUBLE_PRECISION,ngh(10),tag10,&
     bufRS,NINC(12)*11,MPI_DOUBLE_PRECISION,ngh(12),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSE,NESC(11)*11,MPI_DOUBLE_PRECISION,ngh(11),tag11,&
     bufRW,NINC(13)*11,MPI_DOUBLE_PRECISION,ngh(13),tag11,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS,NESC(12)*11,MPI_DOUBLE_PRECISION,ngh(12),tag12,&
     bufRN,NINC(10)*11,MPI_DOUBLE_PRECISION,ngh(10),tag12,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NESC(13)*11,MPI_DOUBLE_PRECISION,ngh(13),tag13,&
     bufRE,NINC(11)*11,MPI_DOUBLE_PRECISION,ngh(11),tag13,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNE,NESC(14)*11,MPI_DOUBLE_PRECISION,ngh(14),tag14,&
     bufRSW,NINC(16)*11,MPI_DOUBLE_PRECISION,ngh(16),tag14,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSE,NESC(15)*11,MPI_DOUBLE_PRECISION,ngh(15),tag15,&
     bufRNW,NINC(17)*11,MPI_DOUBLE_PRECISION,ngh(17),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSW,NESC(16)*11,MPI_DOUBLE_PRECISION,ngh(16),tag16,&
     bufRNE,NINC(14)*11,MPI_DOUBLE_PRECISION,ngh(14),tag16,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNW,NESC(17)*11,MPI_DOUBLE_PRECISION,ngh(17),tag17,&
     bufRSE,NINC(15)*11,MPI_DOUBLE_PRECISION,ngh(15),tag17,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFN,NESC(18)*11,MPI_DOUBLE_PRECISION,ngh(18),tag18,&
     bufRBS,NINC(3)*11,MPI_DOUBLE_PRECISION,ngh(3),tag18,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFE,NESC(19)*11,MPI_DOUBLE_PRECISION,ngh(19),tag19,&
     bufRBW,NINC(4)*11,MPI_DOUBLE_PRECISION,ngh(4),tag19,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFS,NESC(20)*11,MPI_DOUBLE_PRECISION,ngh(20),tag20,&
     bufRBN,NINC(1)*11,MPI_DOUBLE_PRECISION,ngh(1),tag20,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFW,NESC(21)*11,MPI_DOUBLE_PRECISION,ngh(21),tag21,&
     bufRBE,NINC(2)*11,MPI_DOUBLE_PRECISION,ngh(2),tag21,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFNE,NESC(22)*11,MPI_DOUBLE_PRECISION,ngh(22),tag22,&
     bufRBSW,NINC(7)*11,MPI_DOUBLE_PRECISION,ngh(7),tag22,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFSE,NESC(23)*11,MPI_DOUBLE_PRECISION,ngh(23),tag23,&
     bufRBNW,NINC(8)*11,MPI_DOUBLE_PRECISION,ngh(8),tag23,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFSW,NESC(24)*11,MPI_DOUBLE_PRECISION,ngh(24),tag24,&
     bufRBNE,NINC(5)*11,MPI_DOUBLE_PRECISION,ngh(5),tag24,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFNW,NESC(25)*11,MPI_DOUBLE_PRECISION,ngh(25),tag25,&
     bufRBSE,NINC(6)*11,MPI_DOUBLE_PRECISION,ngh(6),tag25,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF,NESC(26)*11,MPI_DOUBLE_PRECISION,ngh(26),tag26,&
     bufRB,NINC(9)*11,MPI_DOUBLE_PRECISION,ngh(9),tag26,COMM,stat,ierr)

ENDIF

! Exchange particle tags

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(buf_tagSBN,NESC(1),MPI_INTEGER8,ngh(1),tag1,&
                  buf_tagRFS,NINC(20),MPI_INTEGER8,ngh(20),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBE,NESC(2),MPI_INTEGER8,ngh(2),tag2,&
                  buf_tagRFW,NINC(21),MPI_INTEGER8,ngh(21),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBS,NESC(3),MPI_INTEGER8,ngh(3),tag3,&
                  buf_tagRFN,NINC(18),MPI_INTEGER8,ngh(18),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBW,NESC(4),MPI_INTEGER8,ngh(4),tag4,&
                  buf_tagRFE,NINC(19),MPI_INTEGER8,ngh(19),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBNE,NESC(5),MPI_INTEGER8,ngh(5),tag5,&
                  buf_tagRFSW,NINC(24),MPI_INTEGER8,ngh(24),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBSE,NESC(6),MPI_INTEGER8,ngh(6),tag6,&
                  buf_tagRFNW,NINC(25),MPI_INTEGER8,ngh(25),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBSW,NESC(7),MPI_INTEGER8,ngh(7),tag7,&
                  buf_tagRFNE,NINC(22),MPI_INTEGER8,ngh(22),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBNW,NESC(8),MPI_INTEGER8,ngh(8),tag8,&
                  buf_tagRFSE,NINC(23),MPI_INTEGER8,ngh(23),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSB,NESC(9),MPI_INTEGER8,ngh(9),tag9,&
                  buf_tagRF,NINC(26),MPI_INTEGER8,ngh(26),tag9,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSN,NESC(10),MPI_INTEGER8,ngh(10),tag10,&
                  buf_tagRS,NINC(12),MPI_INTEGER8,ngh(12),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSE,NESC(11),MPI_INTEGER8,ngh(11),tag11,&
                  buf_tagRW,NINC(13),MPI_INTEGER8,ngh(13),tag11,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSS,NESC(12),MPI_INTEGER8,ngh(12),tag12,&
                  buf_tagRN,NINC(10),MPI_INTEGER8,ngh(10),tag12,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSW,NESC(13),MPI_INTEGER8,ngh(13),tag13,&
                  buf_tagRE,NINC(11),MPI_INTEGER8,ngh(11),tag13,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSNE,NESC(14),MPI_INTEGER8,ngh(14),tag14,&
                  buf_tagRSW,NINC(16),MPI_INTEGER8,ngh(16),tag14,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSSE,NESC(15),MPI_INTEGER8,ngh(15),tag15,&
                  buf_tagRNW,NINC(17),MPI_INTEGER8,ngh(17),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSSW,NESC(16),MPI_INTEGER8,ngh(16),tag16,&
                  buf_tagRNE,NINC(14),MPI_INTEGER8,ngh(14),tag16,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSNW,NESC(17),MPI_INTEGER8,ngh(17),tag17,&
                  buf_tagRSE,NINC(15),MPI_INTEGER8,ngh(15),tag17,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFN,NESC(18),MPI_INTEGER8,ngh(18),tag18,&
                  buf_tagRBS,NINC(3),MPI_INTEGER8,ngh(3),tag18,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFE,NESC(19),MPI_INTEGER8,ngh(19),tag19,&
                  buf_tagRBW,NINC(4),MPI_INTEGER8,ngh(4),tag19,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFS,NESC(20),MPI_INTEGER8,ngh(20),tag20,&
                  buf_tagRBN,NINC(1),MPI_INTEGER8,ngh(1),tag20,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFW,NESC(21),MPI_INTEGER8,ngh(21),tag21,&
                  buf_tagRBE,NINC(2),MPI_INTEGER8,ngh(2),tag21,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFNE,NESC(22),MPI_INTEGER8,ngh(22),tag22,&
                  buf_tagRBSW,NINC(7),MPI_INTEGER8,ngh(7),tag22,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFSE,NESC(23),MPI_INTEGER8,ngh(23),tag23,&
                  buf_tagRBNW,NINC(8),MPI_INTEGER8,ngh(8),tag23,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFSW,NESC(24),MPI_INTEGER8,ngh(24),tag24,&
                  buf_tagRBNE,NINC(5),MPI_INTEGER8,ngh(5),tag24,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFNW,NESC(25),MPI_INTEGER8,ngh(25),tag25,&
                  buf_tagRBSE,NINC(6),MPI_INTEGER8,ngh(6),tag25,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSF,NESC(26),MPI_INTEGER8,ngh(26),tag26,&
                  buf_tagRB,NINC(9),MPI_INTEGER8,ngh(9),tag26,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(buf_tagSBN,NESC(1),MPI_INTEGER8,ngh(1),tag1,&
                  buf_tagRFS,NINC(20),MPI_INTEGER8,ngh(20),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBE,NESC(2),MPI_INTEGER8,ngh(2),tag2,&
                  buf_tagRFW,NINC(21),MPI_INTEGER8,ngh(21),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBS,NESC(3),MPI_INTEGER8,ngh(3),tag3,&
                  buf_tagRFN,NINC(18),MPI_INTEGER8,ngh(18),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBW,NESC(4),MPI_INTEGER8,ngh(4),tag4,&
                  buf_tagRFE,NINC(19),MPI_INTEGER8,ngh(19),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBNE,NESC(5),MPI_INTEGER8,ngh(5),tag5,&
                  buf_tagRFSW,NINC(24),MPI_INTEGER8,ngh(24),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBSE,NESC(6),MPI_INTEGER8,ngh(6),tag6,&
                  buf_tagRFNW,NINC(25),MPI_INTEGER8,ngh(25),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBSW,NESC(7),MPI_INTEGER8,ngh(7),tag7,&
                  buf_tagRFNE,NINC(22),MPI_INTEGER8,ngh(22),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSBNW,NESC(8),MPI_INTEGER8,ngh(8),tag8,&
                  buf_tagRFSE,NINC(23),MPI_INTEGER8,ngh(23),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSB,NESC(9),MPI_INTEGER8,ngh(9),tag9,&
                  buf_tagRF,NINC(26),MPI_INTEGER8,ngh(26),tag9,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSN,NESC(10),MPI_INTEGER8,ngh(10),tag10,&
                  buf_tagRS,NINC(12),MPI_INTEGER8,ngh(12),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSE,NESC(11),MPI_INTEGER8,ngh(11),tag11,&
                  buf_tagRW,NINC(13),MPI_INTEGER8,ngh(13),tag11,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSS,NESC(12),MPI_INTEGER8,ngh(12),tag12,&
                  buf_tagRN,NINC(10),MPI_INTEGER8,ngh(10),tag12,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSW,NESC(13),MPI_INTEGER8,ngh(13),tag13,&
                  buf_tagRE,NINC(11),MPI_INTEGER8,ngh(11),tag13,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSNE,NESC(14),MPI_INTEGER8,ngh(14),tag14,&
                  buf_tagRSW,NINC(16),MPI_INTEGER8,ngh(16),tag14,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSSE,NESC(15),MPI_INTEGER8,ngh(15),tag15,&
                  buf_tagRNW,NINC(17),MPI_INTEGER8,ngh(17),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSSW,NESC(16),MPI_INTEGER8,ngh(16),tag16,&
                  buf_tagRNE,NINC(14),MPI_INTEGER8,ngh(14),tag16,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSNW,NESC(17),MPI_INTEGER8,ngh(17),tag17,&
                  buf_tagRSE,NINC(15),MPI_INTEGER8,ngh(15),tag17,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFN,NESC(18),MPI_INTEGER8,ngh(18),tag18,&
                  buf_tagRBS,NINC(3),MPI_INTEGER8,ngh(3),tag18,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFE,NESC(19),MPI_INTEGER8,ngh(19),tag19,&
                  buf_tagRBW,NINC(4),MPI_INTEGER8,ngh(4),tag19,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFS,NESC(20),MPI_INTEGER8,ngh(20),tag20,&
                  buf_tagRBN,NINC(1),MPI_INTEGER8,ngh(1),tag20,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFW,NESC(21),MPI_INTEGER8,ngh(21),tag21,&
                  buf_tagRBE,NINC(2),MPI_INTEGER8,ngh(2),tag21,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFNE,NESC(22),MPI_INTEGER8,ngh(22),tag22,&
                  buf_tagRBSW,NINC(7),MPI_INTEGER8,ngh(7),tag22,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFSE,NESC(23),MPI_INTEGER8,ngh(23),tag23,&
                  buf_tagRBNW,NINC(8),MPI_INTEGER8,ngh(8),tag23,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFSW,NESC(24),MPI_INTEGER8,ngh(24),tag24,&
                  buf_tagRBNE,NINC(5),MPI_INTEGER8,ngh(5),tag24,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSFNW,NESC(25),MPI_INTEGER8,ngh(25),tag25,&
                  buf_tagRBSE,NINC(6),MPI_INTEGER8,ngh(6),tag25,COMM,stat,ierr)

CALL MPI_SENDRECV(buf_tagSF,NESC(26),MPI_INTEGER8,ngh(26),tag26,&
                  buf_tagRB,NINC(9),MPI_INTEGER8,ngh(9),tag26,COMM,stat,ierr)

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
z=pcl(3,ip)

IF (x.LE.xmaxp.AND.x.GE.xminp.AND.y.LE.ymaxp.AND.y.GE.yminp.AND.&
    z.LE.zmaxp.AND.z.GE.zminp) THEN
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
pcl_f(1,ip+Ntemp)=bufRBN(ip,1)
pcl_f(2,ip+Ntemp)=bufRBN(ip,2)
pcl_f(3,ip+Ntemp)=bufRBN(ip,3)
pcl_f(4,ip+Ntemp)=bufRBN(ip,4)
pcl_f(5,ip+Ntemp)=bufRBN(ip,5)
pcl_f(6,ip+Ntemp)=bufRBN(ip,6)
pcl_f(7,ip+Ntemp)=bufRBN(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRBN(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRBN(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRBN(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRBN(ip,11)
tagf(ip+Ntemp)=buf_tagRBN(ip)
ENDDO

Ntemp=NPT+NINC(1)

DO ip=1,NINC(2)
pcl_f(1,ip+Ntemp)=bufRBE(ip,1)
pcl_f(2,ip+Ntemp)=bufRBE(ip,2)
pcl_f(3,ip+Ntemp)=bufRBE(ip,3)
pcl_f(4,ip+Ntemp)=bufRBE(ip,4)
pcl_f(5,ip+Ntemp)=bufRBE(ip,5)
pcl_f(6,ip+Ntemp)=bufRBE(ip,6)
pcl_f(7,ip+Ntemp)=bufRBE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRBE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRBE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRBE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRBE(ip,11)
tagf(ip+Ntemp)=buf_tagRBE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)

DO ip=1,NINC(3)
pcl_f(1,ip+Ntemp)=bufRBS(ip,1)
pcl_f(2,ip+Ntemp)=bufRBS(ip,2)
pcl_f(3,ip+Ntemp)=bufRBS(ip,3)
pcl_f(4,ip+Ntemp)=bufRBS(ip,4)
pcl_f(5,ip+Ntemp)=bufRBS(ip,5)
pcl_f(6,ip+Ntemp)=bufRBS(ip,6)
pcl_f(7,ip+Ntemp)=bufRBS(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRBS(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRBS(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRBS(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRBS(ip,11)
tagf(ip+Ntemp)=buf_tagRBS(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)

DO ip=1,NINC(4)
pcl_f(1,ip+Ntemp)=bufRBW(ip,1)
pcl_f(2,ip+Ntemp)=bufRBW(ip,2)
pcl_f(3,ip+Ntemp)=bufRBW(ip,3)
pcl_f(4,ip+Ntemp)=bufRBW(ip,4)
pcl_f(5,ip+Ntemp)=bufRBW(ip,5)
pcl_f(6,ip+Ntemp)=bufRBW(ip,6)
pcl_f(7,ip+Ntemp)=bufRBW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRBW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRBW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRBW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRBW(ip,11)
tagf(ip+Ntemp)=buf_tagRBW(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)

DO ip=1,NINC(5)
pcl_f(1,ip+Ntemp)=bufRBNE(ip,1)
pcl_f(2,ip+Ntemp)=bufRBNE(ip,2)
pcl_f(3,ip+Ntemp)=bufRBNE(ip,3)
pcl_f(4,ip+Ntemp)=bufRBNE(ip,4)
pcl_f(5,ip+Ntemp)=bufRBNE(ip,5)
pcl_f(6,ip+Ntemp)=bufRBNE(ip,6)
pcl_f(7,ip+Ntemp)=bufRBNE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRBNE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRBNE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRBNE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRBNE(ip,11)
tagf(ip+Ntemp)=buf_tagRBNE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)

DO ip=1,NINC(6)
pcl_f(1,ip+Ntemp)=bufRBSE(ip,1)
pcl_f(2,ip+Ntemp)=bufRBSE(ip,2)
pcl_f(3,ip+Ntemp)=bufRBSE(ip,3)
pcl_f(4,ip+Ntemp)=bufRBSE(ip,4)
pcl_f(5,ip+Ntemp)=bufRBSE(ip,5)
pcl_f(6,ip+Ntemp)=bufRBSE(ip,6)
pcl_f(7,ip+Ntemp)=bufRBSE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRBSE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRBSE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRBSE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRBSE(ip,11)
tagf(ip+Ntemp)=buf_tagRBSE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)

DO ip=1,NINC(7)
pcl_f(1,ip+Ntemp)=bufRBSW(ip,1)
pcl_f(2,ip+Ntemp)=bufRBSW(ip,2)
pcl_f(3,ip+Ntemp)=bufRBSW(ip,3)
pcl_f(4,ip+Ntemp)=bufRBSW(ip,4)
pcl_f(5,ip+Ntemp)=bufRBSW(ip,5)
pcl_f(6,ip+Ntemp)=bufRBSW(ip,6)
pcl_f(7,ip+Ntemp)=bufRBSW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRBSW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRBSW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRBSW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRBSW(ip,11)
tagf(ip+Ntemp)=buf_tagRBSW(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)

DO ip=1,NINC(8)
pcl_f(1,ip+Ntemp)=bufRBNW(ip,1)
pcl_f(2,ip+Ntemp)=bufRBNW(ip,2)
pcl_f(3,ip+Ntemp)=bufRBNW(ip,3)
pcl_f(4,ip+Ntemp)=bufRBNW(ip,4)
pcl_f(5,ip+Ntemp)=bufRBNW(ip,5)
pcl_f(6,ip+Ntemp)=bufRBNW(ip,6)
pcl_f(7,ip+Ntemp)=bufRBNW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRBNW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRBNW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRBNW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRBNW(ip,11)
tagf(ip+Ntemp)=buf_tagRBNW(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)

DO ip=1,NINC(9)
pcl_f(1,ip+Ntemp)=bufRB(ip,1)
pcl_f(2,ip+Ntemp)=bufRB(ip,2)
pcl_f(3,ip+Ntemp)=bufRB(ip,3)
pcl_f(4,ip+Ntemp)=bufRB(ip,4)
pcl_f(5,ip+Ntemp)=bufRB(ip,5)
pcl_f(6,ip+Ntemp)=bufRB(ip,6)
pcl_f(7,ip+Ntemp)=bufRB(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRB(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRB(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRB(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRB(ip,11)
tagf(ip+Ntemp)=buf_tagRB(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)

DO ip=1,NINC(10)
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

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)

DO ip=1,NINC(11)
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

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)

DO ip=1,NINC(12)
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

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)

DO ip=1,NINC(13)
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

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)

DO ip=1,NINC(14)
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

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)

DO ip=1,NINC(15)
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

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)

DO ip=1,NINC(16)
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

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)

DO ip=1,NINC(17)
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

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)

DO ip=1,NINC(18)
pcl_f(1,ip+Ntemp)=bufRFN(ip,1)
pcl_f(2,ip+Ntemp)=bufRFN(ip,2)
pcl_f(3,ip+Ntemp)=bufRFN(ip,3)
pcl_f(4,ip+Ntemp)=bufRFN(ip,4)
pcl_f(5,ip+Ntemp)=bufRFN(ip,5)
pcl_f(6,ip+Ntemp)=bufRFN(ip,6)
pcl_f(7,ip+Ntemp)=bufRFN(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRFN(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRFN(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRFN(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRFN(ip,11)
tagf(ip+Ntemp)=buf_tagRFN(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+NINC(18)

DO ip=1,NINC(19)
pcl_f(1,ip+Ntemp)=bufRFE(ip,1)
pcl_f(2,ip+Ntemp)=bufRFE(ip,2)
pcl_f(3,ip+Ntemp)=bufRFE(ip,3)
pcl_f(4,ip+Ntemp)=bufRFE(ip,4)
pcl_f(5,ip+Ntemp)=bufRFE(ip,5)
pcl_f(6,ip+Ntemp)=bufRFE(ip,6)
pcl_f(7,ip+Ntemp)=bufRFE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRFE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRFE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRFE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRFE(ip,11)
tagf(ip+Ntemp)=buf_tagRFE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+NINC(18)+&
NINC(19)

DO ip=1,NINC(20)
pcl_f(1,ip+Ntemp)=bufRFS(ip,1)
pcl_f(2,ip+Ntemp)=bufRFS(ip,2)
pcl_f(3,ip+Ntemp)=bufRFS(ip,3)
pcl_f(4,ip+Ntemp)=bufRFS(ip,4)
pcl_f(5,ip+Ntemp)=bufRFS(ip,5)
pcl_f(6,ip+Ntemp)=bufRFS(ip,6)
pcl_f(7,ip+Ntemp)=bufRFS(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRFS(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRFS(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRFS(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRFS(ip,11)
tagf(ip+Ntemp)=buf_tagRFS(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+NINC(18)+&
NINC(19)+NINC(20)

DO ip=1,NINC(21)
pcl_f(1,ip+Ntemp)=bufRFW(ip,1)
pcl_f(2,ip+Ntemp)=bufRFW(ip,2)
pcl_f(3,ip+Ntemp)=bufRFW(ip,3)
pcl_f(4,ip+Ntemp)=bufRFW(ip,4)
pcl_f(5,ip+Ntemp)=bufRFW(ip,5)
pcl_f(6,ip+Ntemp)=bufRFW(ip,6)
pcl_f(7,ip+Ntemp)=bufRFW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRFW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRFW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRFW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRFW(ip,11)
tagf(ip+Ntemp)=buf_tagRFW(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+NINC(18)+&
NINC(19)+NINC(20)+NINC(21)

DO ip=1,NINC(22)
pcl_f(1,ip+Ntemp)=bufRFNE(ip,1)
pcl_f(2,ip+Ntemp)=bufRFNE(ip,2)
pcl_f(3,ip+Ntemp)=bufRFNE(ip,3)
pcl_f(4,ip+Ntemp)=bufRFNE(ip,4)
pcl_f(5,ip+Ntemp)=bufRFNE(ip,5)
pcl_f(6,ip+Ntemp)=bufRFNE(ip,6)
pcl_f(7,ip+Ntemp)=bufRFNE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRFNE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRFNE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRFNE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRFNE(ip,11)
tagf(ip+Ntemp)=buf_tagRFNE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+NINC(18)+&
NINC(19)+NINC(20)+NINC(21)+NINC(22)

DO ip=1,NINC(23)
pcl_f(1,ip+Ntemp)=bufRFSE(ip,1)
pcl_f(2,ip+Ntemp)=bufRFSE(ip,2)
pcl_f(3,ip+Ntemp)=bufRFSE(ip,3)
pcl_f(4,ip+Ntemp)=bufRFSE(ip,4)
pcl_f(5,ip+Ntemp)=bufRFSE(ip,5)
pcl_f(6,ip+Ntemp)=bufRFSE(ip,6)
pcl_f(7,ip+Ntemp)=bufRFSE(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRFSE(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRFSE(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRFSE(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRFSE(ip,11)
tagf(ip+Ntemp)=buf_tagRFSE(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+NINC(18)+&
NINC(19)+NINC(20)+NINC(21)+NINC(22)+NINC(23)

DO ip=1,NINC(24)
pcl_f(1,ip+Ntemp)=bufRFSW(ip,1)
pcl_f(2,ip+Ntemp)=bufRFSW(ip,2)
pcl_f(3,ip+Ntemp)=bufRFSW(ip,3)
pcl_f(4,ip+Ntemp)=bufRFSW(ip,4)
pcl_f(5,ip+Ntemp)=bufRFSW(ip,5)
pcl_f(6,ip+Ntemp)=bufRFSW(ip,6)
pcl_f(7,ip+Ntemp)=bufRFSW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRFSW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRFSW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRFSW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRFSW(ip,11)
tagf(ip+Ntemp)=buf_tagRFSW(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+NINC(18)+&
NINC(19)+NINC(20)+NINC(21)+NINC(22)+NINC(23)+NINC(24)

DO ip=1,NINC(25)
pcl_f(1,ip+Ntemp)=bufRFNW(ip,1)
pcl_f(2,ip+Ntemp)=bufRFNW(ip,2)
pcl_f(3,ip+Ntemp)=bufRFNW(ip,3)
pcl_f(4,ip+Ntemp)=bufRFNW(ip,4)
pcl_f(5,ip+Ntemp)=bufRFNW(ip,5)
pcl_f(6,ip+Ntemp)=bufRFNW(ip,6)
pcl_f(7,ip+Ntemp)=bufRFNW(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRFNW(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRFNW(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRFNW(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRFNW(ip,11)
tagf(ip+Ntemp)=buf_tagRFNW(ip)
ENDDO

Ntemp=NPT+NINC(1)+NINC(2)+NINC(3)+NINC(4)+NINC(5)+NINC(6)+NINC(7)+NINC(8)+NINC(9)+&
NINC(10)+NINC(11)+NINC(12)+NINC(13)+NINC(14)+NINC(15)+NINC(16)+NINC(17)+NINC(18)+&
NINC(19)+NINC(20)+NINC(21)+NINC(22)+NINC(23)+NINC(24)+NINC(25)

DO ip=1,NINC(26)
pcl_f(1,ip+Ntemp)=bufRF(ip,1)
pcl_f(2,ip+Ntemp)=bufRF(ip,2)
pcl_f(3,ip+Ntemp)=bufRF(ip,3)
pcl_f(4,ip+Ntemp)=bufRF(ip,4)
pcl_f(5,ip+Ntemp)=bufRF(ip,5)
pcl_f(6,ip+Ntemp)=bufRF(ip,6)
pcl_f(7,ip+Ntemp)=bufRF(ip,7)
pcl_data_f(1,ip+Ntemp)=bufRF(ip,8)
pcl_data_f(2,ip+Ntemp)=bufRF(ip,9)
pcl_data_f(3,ip+Ntemp)=bufRF(ip,10)
pcl_data_f(4,ip+Ntemp)=bufRF(ip,11)
tagf(ip+Ntemp)=buf_tagRF(ip)
ENDDO

DEALLOCATE(bufSBN,bufSBE,bufSBS,bufSBW)
DEALLOCATE(bufSBNE,bufSBSE,bufSBSW)
DEALLOCATE(bufSBNW,bufSB,bufSN,bufSE)
DEALLOCATE(bufSS,bufSW,bufSNE,bufSSE)
DEALLOCATE(bufSSW,bufSNW,bufSFN)
DEALLOCATE(bufSFE,bufSFS,bufSFW)
DEALLOCATE(bufSFNE,bufSFSE,bufSFSW)
DEALLOCATE(bufSFNW,bufSF)

DEALLOCATE(buf_tagSBN,buf_tagSBE,buf_tagSBS)
DEALLOCATE(buf_tagSBW,buf_tagSBNE,buf_tagSBSE)
DEALLOCATE(buf_tagSBSW,buf_tagSBNW,buf_tagSB)
DEALLOCATE(buf_tagSN,buf_tagSE,buf_tagSS)
DEALLOCATE(buf_tagSW,buf_tagSNE,buf_tagSSE)
DEALLOCATE(buf_tagSSW,buf_tagSNW,buf_tagSFN)
DEALLOCATE(buf_tagSFE,buf_tagSFS,buf_tagSFW)
DEALLOCATE(buf_tagSFNE,buf_tagSFSE,buf_tagSFSW)
DEALLOCATE(buf_tagSFNW,buf_tagSF)

DEALLOCATE(bufRBN,bufRBE,bufRBS,bufRBW)
DEALLOCATE(bufRBNE,bufRBSE,bufRBSW)
DEALLOCATE(bufRBNW,bufRB,bufRN,bufRE)
DEALLOCATE(bufRS,bufRW,bufRNE,bufRSE)
DEALLOCATE(bufRSW,bufRNW,bufRFN)
DEALLOCATE(bufRFE,bufRFS,bufRFW)
DEALLOCATE(bufRFNE,bufRFSE,bufRFSW)
DEALLOCATE(bufRFNW,bufRF)

DEALLOCATE(buf_tagRBN,buf_tagRBE,buf_tagRBS)
DEALLOCATE(buf_tagRBW,buf_tagRBNE,buf_tagRBSE)
DEALLOCATE(buf_tagRBSW,buf_tagRBNW,buf_tagRB)
DEALLOCATE(buf_tagRN,buf_tagRE,buf_tagRS)
DEALLOCATE(buf_tagRW,buf_tagRNE,buf_tagRSE)
DEALLOCATE(buf_tagRSW,buf_tagRNW,buf_tagRFN)
DEALLOCATE(buf_tagRFE,buf_tagRFS,buf_tagRFW)
DEALLOCATE(buf_tagRFNE,buf_tagRFSE,buf_tagRFSW)
DEALLOCATE(buf_tagRFNW,buf_tagRF)

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
