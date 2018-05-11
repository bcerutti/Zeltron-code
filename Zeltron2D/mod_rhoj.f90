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

MODULE MOD_RHOJ

USE MOD_INPUT

IMPLICIT NONE

PRIVATE

PUBLIC :: RHOJ ! computes the charge density and current on grid nodes
PUBLIC :: COM_RHOJ ! computes the charge density and current at the boundaries between domains
PUBLIC :: JYEE ! computes J in the Yee lattice

 CONTAINS

!***********************************************************************
! Subroutine RHOJ
! This subroutine computes the charge density RHO and *half* of the 
! current density J on grid nodes.

! INPUT: 
! - q: sign electric charge
! - pcl: Particle distribution function
!
! OUTPUT: RHO at t=t+dt and Jx, Jy, Jz at t=t+dt/2
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE RHOJ(q,pcl,rho,Jx,Jy,Jz,xgp,ygp,NPP,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER*8                                   :: ip,NPP
DOUBLE PRECISION                            :: q,rhop,gam,vx,vy,vz
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)      :: pcl
DOUBLE PRECISION                            :: x,y,ux,uy,uz,wght
DOUBLE PRECISION, DIMENSION(1:NXP)          :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)          :: ygp
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP)    :: rho,Jx,Jy,Jz
DOUBLE PRECISION                            :: fp,fq
INTEGER                                     :: i,j
INTEGER, DIMENSION(8)                       :: ngh
INTEGER                                     :: id,COMM,ierr

!***********************************************************************

rho=0.0
Jx=0.0
Jy=0.0
Jz=0.0

DO ip=1,NPP

  x=pcl(1,ip)
  y=pcl(2,ip)
  ux=pcl(4,ip)
  uy=pcl(5,ip)
  uz=pcl(6,ip)
  wght=pcl(7,ip)
  
  ! Charge density of 1 super-particle in the cell (i,j)
  rhop=q*e*wght/(dx*dy)
  
  ! Lorentz factor
  gam=sqrt(1.0+ux*ux+uy*uy+uz*uz)

  ! 3-velocity components
  vx=ux*c/gam
  vy=uy*c/gam
  vz=uz*c/gam

  ! Computation of the nearest node index to (xf,yf), for a constant dx,dy
  i=FLOOR((x-xgp(1))/dx)+1
  j=FLOOR((y-ygp(1))/dy)+1
  
  IF (i.EQ.NXP) THEN
  i=i-1
  END IF
  
  IF (j.EQ.NYP) THEN
  j=j-1
  END IF
  
  fp=(x-xgp(i))/dx
  fq=(y-ygp(j))/dy

  !=====================================================================
  ! Charge density

  rho(i,j)=rho(i,j)+rhop*(1.0-fp)*(1.0-fq)
  rho(i+1,j)=rho(i+1,j)+rhop*fp*(1.0-fq)
  rho(i,j+1)=rho(i,j+1)+rhop*(1.0-fp)*fq
  rho(i+1,j+1)=rho(i+1,j+1)+rhop*fp*fq

  !=====================================================================
  ! Current density

  ! x-component:
  Jx(i,j)=Jx(i,j)+0.5*rhop*vx*(1.0-fp)*(1.0-fq)
  Jx(i+1,j)=Jx(i+1,j)+0.5*rhop*vx*fp*(1.0-fq)
  Jx(i,j+1)=Jx(i,j+1)+0.5*rhop*vx*(1.0-fp)*fq
  Jx(i+1,j+1)=Jx(i+1,j+1)+0.5*rhop*vx*fp*fq

  ! y-component:
  Jy(i,j)=Jy(i,j)+0.5*rhop*vy*(1.0-fp)*(1.0-fq)
  Jy(i+1,j)=Jy(i+1,j)+0.5*rhop*vy*fp*(1.0-fq)
  Jy(i,j+1)=Jy(i,j+1)+0.5*rhop*vy*(1.0-fp)*fq
  Jy(i+1,j+1)=Jy(i+1,j+1)+0.5*rhop*vy*fp*fq

  ! z-component:
  Jz(i,j)=Jz(i,j)+0.5*rhop*vz*(1.0-fp)*(1.0-fq)
  Jz(i+1,j)=Jz(i+1,j)+0.5*rhop*vz*fp*(1.0-fq)
  Jz(i,j+1)=Jz(i,j+1)+0.5*rhop*vz*(1.0-fp)*fq
  Jz(i+1,j+1)=Jz(i+1,j+1)+0.5*rhop*vz*fp*fq

ENDDO

CALL COM_RHOJ(rho,xgp,ygp,id,ngh,COMM,ierr)
CALL COM_RHOJ(Jx,xgp,ygp,id,ngh,COMM,ierr)
CALL COM_RHOJ(Jy,xgp,ygp,id,ngh,COMM,ierr)
CALL COM_RHOJ(Jz,xgp,ygp,id,ngh,COMM,ierr)

END SUBROUTINE RHOJ

!***********************************************************************
! Subroutine COM_RHOJ
! This subroutine computes the charge density RHO and current density J
! at the boundaries between sub-domains.

! INPUT: 
! - fun: RHO or Jx,y,z
! - local grid: xgp,ygp
!
! OUTPUT: RHO or Jx,y,z at the boundaries
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE COM_RHOJ(fun,xgp,ygp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: fun
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
INTEGER, DIMENSION(8)                    :: ngh
INTEGER, PARAMETER                       :: tag1=10,tag2=20,tag3=30,tag4=40
INTEGER, PARAMETER                       :: tag5=50,tag6=60,tag7=70,tag8=80
INTEGER                                  :: id,COMM,ierr
DOUBLE PRECISION, DIMENSION(1:NXP)       :: bufSN,bufSS,bufRN,bufRS
DOUBLE PRECISION, DIMENSION(1:NYP)       :: bufSE,bufSW,bufRE,bufRW
DOUBLE PRECISION                         :: bufSNE,bufSSE,bufSSW,bufSNW
DOUBLE PRECISION                         :: bufRNE,bufRSE,bufRSW,bufRNW
!***********************************************************************

bufSS=fun(:,1)
bufSN=fun(:,NYP)

bufSW=fun(1,:)
bufSE=fun(NXP,:)

bufSNE=fun(NXP,NYP)
bufSSE=fun(NXP,1)
bufSSW=fun(1,1)
bufSNW=fun(1,NYP)

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufRS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag2,&
                  bufRN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag3,&
                  bufRW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag4,&
                  bufRE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag5,&
                  bufRSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag6,&
                  bufRNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag7,&
                  bufRNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag8,&
                  bufRSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag8,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(bufSN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufRS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag2,&
                  bufRN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag3,&
                  bufRW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag4,&
                  bufRE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag5,&
                  bufRSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag6,&
                  bufRNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag7,&
                  bufRNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag8,&
                  bufRSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag8,COMM,stat,ierr)

ENDIF

!***********************************************************************
! Check boundary conditions along X
   
IF (xgp(1).EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.NE."PERIODIC") THEN
   bufRW(:)=0.0
   bufRNW=0.0
   bufRSW=0.0
   END IF

END IF

IF (xgp(NXP).EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.NE."PERIODIC") THEN
   bufRE(:)=0.0
   bufRNE=0.0
   bufRSE=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (ygp(1).EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.NE."PERIODIC") THEN
   bufRS(:)=0.0
   bufRSW=0.0
   bufRSE=0.0
   END IF
   
END IF

IF (ygp(NYP).EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.NE."PERIODIC") THEN
   bufRN(:)=0.0
   bufRNW=0.0
   bufRNE=0.0
   END IF

END IF
   
!***********************************************************************

fun(NXP,:)=fun(NXP,:)+bufRE(:)
fun(1,:)=fun(1,:)+bufRW(:)

fun(:,NYP)=fun(:,NYP)+bufRN(:)
fun(:,1)=fun(:,1)+bufRS(:)

fun(1,1)=fun(1,1)+bufRSW
fun(NXP,1)=fun(NXP,1)+bufRSE
fun(1,NYP)=fun(1,NYP)+bufRNW
fun(NXP,NYP)=fun(NXP,NYP)+bufRNE

END SUBROUTINE COM_RHOJ

!***********************************************************************
! Subroutine JYEE
! This subroutine computes the current density in the Yee lattice.

! INPUT: 
! - Jx: x-component of the current density J at the nodes at t=t+dt/2
! - Jy: y-component of the current density J at the nodes at t=t+dt/2
! - Jz: z-component of the current density J at the nodes at t=t+dt/2
!
! OUTPUT: Jx, Jy, Jz at t=t+dt/2 in the Yee lattice
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE JYEE(Jx,Jy,Jz,xgp,ygp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Jx,Jy,Jz
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION                         :: xmaxp,ymaxp
INTEGER                                  :: ix,iy
INTEGER, DIMENSION(8)                    :: ngh
INTEGER, PARAMETER                       :: tag=5
INTEGER                                  :: id,COMM,ierr
DOUBLE PRECISION, ALLOCATABLE            :: bufS(:),bufR(:)

!***********************************************************************

xmaxp=xgp(NXP)
ymaxp=ygp(NYP)

!***********************************************************************
! Jx

ALLOCATE(bufS(1:NYP),bufR(1:NYP))

bufS=Jx(2,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag,&
                  bufR,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag,&
                  bufR,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag,COMM,stat,ierr)
ENDIF

DO ix=1,NXP-1
Jx(ix,:)=(Jx(ix+1,:)+Jx(ix,:))/2.0
ENDDO

!***********************************************************************
! Check boundary conditions along X
   
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   bufR(:)=0.0
   ENDIF

END IF
!***********************************************************************

Jx(NXP,:)=(bufR(:)+Jx(NXP,:))/2.0

DEALLOCATE(bufS,bufR)

!***********************************************************************
! Jy

ALLOCATE(bufS(1:NXP),bufR(1:NXP))

bufS=Jy(:,2)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag,&
                  bufR,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag,&
                  bufR,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag,COMM,stat,ierr)
ENDIF

DO iy=1,NYP-1
Jy(:,iy)=(Jy(:,iy+1)+Jy(:,iy))/2.0
ENDDO

!***********************************************************************
! Check boundary conditions along Y
   
IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   bufR(:)=0.0
   ENDIF

END IF
!***********************************************************************

Jy(:,NYP)=(bufR(:)+Jy(:,NYP))/2.0

DEALLOCATE(bufS,bufR)

!***********************************************************************
Jz=Jz

END SUBROUTINE JYEE

!***********************************************************************

END MODULE MOD_RHOJ
