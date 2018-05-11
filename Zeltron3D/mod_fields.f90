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

MODULE MOD_FIELDS

USE MOD_INPUT

IMPLICIT NONE

PRIVATE

PUBLIC :: PUSH_EFIELD ! Computes E at t+dt
PUBLIC :: PUSH_BHALF ! Computes B from t to t+dt/2
PUBLIC :: FIELDS_NODES ! Computes E and B at nodes
PUBLIC :: CORRECT_EFIELD ! Ensure that div(E)=4*pi*rho
PUBLIC :: FILTER_FIELD ! Isotropic 3D filter using the nearest neighboring cells

 CONTAINS

!***********************************************************************
! Subroutine PUSH_EFIELD
! This subroutine computes the E field vector at time t+dt

! INPUT: 
! - Bx: x-component of B on the Yee grid at time t+dt/2
! - By: y-component of B on the Yee grid at time t+dt/2
! - Bz: z-component of B on the Yee grid at time t+dt/2
! - Ex: x-component of E on the Yee grid at time t
! - Ey: y-component of E on the Yee grid at time t
! - Ez: z-component of E on the Yee grid at time t
! - Jx: x-component of current J on the Yee grid at time t
! - Jy: y-component of current J on the Yee grid at time t
! - Jz: z-component of current J on the Yee grid at time t
! - xgp,ygp,zgp: Local grid
!
! OUTPUT: Ex,Ey,Ez at time t+dt
!
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE PUSH_EFIELD(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,xgp,ygp,zgp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bx,By,Bz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Jx,Jy,Jz
DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp
DOUBLE PRECISION                               :: xminp,xmaxp,yminp,ymaxp,zminp,zmaxp
INTEGER, DIMENSION(26)                         :: ngh
INTEGER, PARAMETER                             :: tag1=1,tag2=2
INTEGER                                        :: id,COMM,ierr
DOUBLE PRECISION, ALLOCATABLE                  :: bufS1(:,:),BufR1(:,:)
DOUBLE PRECISION, ALLOCATABLE                  :: bufS2(:,:),BufR2(:,:)

! Loop indexes
INTEGER :: ix,iy,iz
!***********************************************************************

xminp=xgp(1)
xmaxp=xgp(NXP)
yminp=ygp(1)
ymaxp=ygp(NYP)
zminp=zgp(1)
zmaxp=zgp(NZP)

!***********************************************************************
! Solve Ex at t=t+dt

ALLOCATE(bufS1(1:NXP,1:NZP),bufR1(1:NXP,1:NZP))
ALLOCATE(bufS2(1:NXP,1:NYP),bufR2(1:NXP,1:NYP))

bufS1=Bz(:,NYP-1,:)
bufS2=By(:,:,NZP-1)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag2,&
                  bufR2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag2,&
                  bufR2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag2,COMM,stat,ierr)
ENDIF

DO iy=2,NYP
DO iz=2,NZP
Ex(:,iy,iz)=Ex(:,iy,iz)+(c*dt/dy)*(Bz(:,iy,iz)-Bz(:,iy-1,iz))-&
                        (c*dt/dz)*(By(:,iy,iz)-By(:,iy,iz-1))-4.0*pi*dt*Jx(:,iy,iz)
ENDDO
ENDDO

DO iz=2,NZP
Ex(:,1,iz)=Ex(:,1,iz)+(c*dt/dy)*(Bz(:,1,iz)-bufR1(:,iz))-&
                      (c*dt/dz)*(By(:,1,iz)-By(:,1,iz-1))-4.0*pi*dt*Jx(:,1,iz)
ENDDO

DO iy=2,NYP
Ex(:,iy,1)=Ex(:,iy,1)+(c*dt/dy)*(Bz(:,iy,1)-Bz(:,iy-1,1))-&
                      (c*dt/dz)*(By(:,iy,1)-bufR2(:,iy))-4.0*pi*dt*Jx(:,iy,1)
ENDDO

Ex(:,1,1)=Ex(:,1,1)+(c*dt/dy)*(Bz(:,1,1)-bufR1(:,1))-&
                    (c*dt/dz)*(By(:,1,1)-bufR2(:,1))-4.0*pi*dt*Jx(:,1,1)

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X
  
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Inside the conductor
   Ex(NXP,:,:)=0.0
   END IF
   
END IF

!***********************************************************************
! Check boundary conditions along Y
 
IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ex(:,1,:)=0.0
   END IF
   
END IF
   
IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ex(:,NYP,:)=0.0
   END IF
   
END IF

!***********************************************************************
! Check boundary conditions along Z
 
IF (zminp.EQ.zmin) THEN

   IF (BOUND_FIELD_ZMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ex(:,:,1)=0.0
   END IF
   
END IF
   
IF (zmaxp.EQ.zmax) THEN

   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ex(:,:,NZP)=0.0
   END IF
   
END IF

!***********************************************************************
! Solve Ey at t=t+dt

ALLOCATE(bufS1(1:NXP,1:NYP),bufR1(1:NXP,1:NYP))
ALLOCATE(bufS2(1:NYP,1:NZP),bufR2(1:NYP,1:NZP))

bufS1=Bx(:,:,NZP-1)
bufS2=Bz(NXP-1,:,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag1,&
                  bufR1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag1,&
                  bufR1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
DO iz=2,NZP
Ey(ix,:,iz)=Ey(ix,:,iz)+(c*dt/dz)*(Bx(ix,:,iz)-Bx(ix,:,iz-1))-&
                        (c*dt/dx)*(Bz(ix,:,iz)-Bz(ix-1,:,iz))-4.0*pi*dt*Jy(ix,:,iz)
ENDDO
ENDDO

DO iz=2,NZP
Ey(1,:,iz)=Ey(1,:,iz)+(c*dt/dz)*(Bx(1,:,iz)-Bx(1,:,iz-1))-&
                      (c*dt/dx)*(Bz(1,:,iz)-bufR2(:,iz))-4.0*pi*dt*Jy(1,:,iz)
ENDDO

DO ix=2,NXP
Ey(ix,:,1)=Ey(ix,:,1)+(c*dt/dz)*(Bx(ix,:,1)-bufR1(ix,:))-&
                      (c*dt/dx)*(Bz(ix,:,1)-Bz(ix-1,:,1))-4.0*pi*dt*Jy(ix,:,1)
ENDDO

Ey(1,:,1)=Ey(1,:,1)+(c*dt/dz)*(Bx(1,:,1)-bufR1(1,:))-&
                    (c*dt/dx)*(Bz(1,:,1)-bufR2(:,1))-4.0*pi*dt*Jy(1,:,1)

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X

IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ey(1,:,:)=0.0
   END IF
   
END IF
   
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ey(NXP,:,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Inside the conductor
   Ey(:,NYP,:)=0.0
   END IF
   
END IF

!***********************************************************************
! Check boundary conditions along Z

IF (zminp.EQ.zmin) THEN

   IF (BOUND_FIELD_ZMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ey(:,:,1)=0.0
   END IF
   
END IF
   
IF (zmaxp.EQ.zmax) THEN

   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ey(:,:,NZP)=0.0
   END IF

END IF

!***********************************************************************
! Solve Ez at t=t+dt

ALLOCATE(bufS1(1:NXP,1:NZP),bufR1(1:NXP,1:NZP))
ALLOCATE(bufS2(1:NYP,1:NZP),bufR2(1:NYP,1:NZP))

bufS1=Bx(:,NYP-1,:)
bufS2=By(NXP-1,:,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
DO iy=2,NYP
Ez(ix,iy,:)=Ez(ix,iy,:)+(c*dt/dx)*(By(ix,iy,:)-By(ix-1,iy,:))-&
                        (c*dt/dy)*(Bx(ix,iy,:)-Bx(ix,iy-1,:))-4.0*pi*dt*Jz(ix,iy,:)
ENDDO
ENDDO

DO iy=2,NYP
Ez(1,iy,:)=Ez(1,iy,:)+(c*dt/dx)*(By(1,iy,:)-bufR2(iy,:))-&
                      (c*dt/dy)*(Bx(1,iy,:)-Bx(1,iy-1,:))-4.0*pi*dt*Jz(1,iy,:)
ENDDO

DO ix=2,NXP
Ez(ix,1,:)=Ez(ix,1,:)+(c*dt/dx)*(By(ix,1,:)-By(ix-1,1,:))-&
                      (c*dt/dy)*(Bx(ix,1,:)-bufR1(ix,:))-4.0*pi*dt*Jz(ix,1,:)
ENDDO

Ez(1,1,:)=Ez(1,1,:)+(c*dt/dx)*(By(1,1,:)-bufR2(1,:))-&
                    (c*dt/dy)*(Bx(1,1,:)-bufR1(1,:))-4.0*pi*dt*Jz(1,1,:)

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X

IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ez(1,:,:)=0.0
   END IF
      
END IF

IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ez(NXP,:,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ez(:,1,:)=0.0
   END IF
   
END IF

IF (ymaxp.EQ.ymax) THEN
   
   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ez(:,NYP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Z

IF (zmaxp.EQ.zmax) THEN

   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   ! Inside the conductor
   Ez(:,:,NZP)=0.0
   END IF
   
END IF

!***********************************************************************

END SUBROUTINE PUSH_EFIELD

!***********************************************************************
! Subroutine PUSH_BHALF
! This subroutine computes the B field vector at t=t+dt/2

! INPUT: 
! - Bx: x-component of initial B on the Yee grid at time t
! - By: y-component of initial B on the Yee grid at time t
! - Bz: z-component of initial B on the Yee grid at time t
! - Ex: x-component of initial E on the Yee grid at time t
! - Ey: y-component of initial E on the Yee grid at time t
! - Ez: z-component of initial E on the Yee grid at time t
! - xgp,ygp,zgp: Local grid
!
! OUTPUT: Bx,By,Bz at time t+dt/2
!
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE PUSH_BHALF(Bx,By,Bz,Ex,Ey,Ez,xgp,ygp,zgp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bx,By,Bz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp
DOUBLE PRECISION                               :: xminp,xmaxp,yminp,ymaxp,zminp,zmaxp
INTEGER, DIMENSION(26)                         :: ngh
INTEGER, PARAMETER                             :: tag1=1,tag2=2
INTEGER                                        :: id,COMM,ierr
DOUBLE PRECISION, ALLOCATABLE                  :: bufS1(:,:),BufR1(:,:)
DOUBLE PRECISION, ALLOCATABLE                  :: bufS2(:,:),BufR2(:,:)

! Loop indexes
INTEGER :: ix,iy,iz
!***********************************************************************

xminp=xgp(1)
xmaxp=xgp(NXP)
yminp=ygp(1)
ymaxp=ygp(NYP)
zminp=zgp(1)
zmaxp=zgp(NZP)

!***********************************************************************
! Bx

ALLOCATE(bufS1(1:NXP,1:NZP),bufR1(1:NXP,1:NZP))
ALLOCATE(bufS2(1:NXP,1:NYP),bufR2(1:NXP,1:NYP))

bufS1=Ez(:,2,:)
bufS2=Ey(:,:,2)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag2,&
                  bufR2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag2,&
                  bufR2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag2,COMM,stat,ierr)
ENDIF

DO iy=1,NYP-1
DO iz=1,NZP-1
Bx(:,iy,iz)=Bx(:,iy,iz)-c*dt/(2.0*dy)*(Ez(:,iy+1,iz)-Ez(:,iy,iz))+&
                        c*dt/(2.0*dz)*(Ey(:,iy,iz+1)-Ey(:,iy,iz))
ENDDO
ENDDO

DO iz=1,NZP-1
Bx(:,NYP,iz)=Bx(:,NYP,iz)-c*dt/(2.0*dy)*(bufR1(:,iz)-Ez(:,NYP,iz))+&
                          c*dt/(2.0*dz)*(Ey(:,NYP,iz+1)-Ey(:,NYP,iz))
ENDDO

DO iy=1,NYP-1
Bx(:,iy,NZP)=Bx(:,iy,NZP)-c*dt/(2.0*dy)*(Ez(:,iy+1,NZP)-Ez(:,iy,NZP))+&
                          c*dt/(2.0*dz)*(bufR2(:,iy)-Ey(:,iy,NZP))
ENDDO

Bx(:,NYP,NZP)=Bx(:,NYP,NZP)-c*dt/(2.0*dy)*(bufR1(:,NZP)-Ez(:,NYP,NZP))+&
                            c*dt/(2.0*dz)*(bufR2(:,NYP)-Ey(:,NYP,NZP))

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X
   
IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! Normal to conductor surface
   Bx(1,:,:)=0.0
   END IF
   
END IF
   
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Normal to conductor surface
   Bx(NXP,:,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Inside the conductor
   Bx(:,NYP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Z

IF (zmaxp.EQ.zmax) THEN

   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   ! Inside the conductor
   Bx(:,:,NZP)=0.0
   END IF

END IF

!***********************************************************************
! By

ALLOCATE(bufS1(1:NXP,1:NYP),bufR1(1:NXP,1:NYP))
ALLOCATE(bufS2(1:NYP,1:NZP),bufR2(1:NYP,1:NZP))

bufS1=Ex(:,:,2)
bufS2=Ez(2,:,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag1,&
                  bufR1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag1,&
                  bufR1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,COMM,stat,ierr)
ENDIF

DO ix=1,NXP-1
DO iz=1,NZP-1
By(ix,:,iz)=By(ix,:,iz)-c*dt/(2.0*dz)*(Ex(ix,:,iz+1)-Ex(ix,:,iz))+&
                        c*dt/(2.0*dx)*(Ez(ix+1,:,iz)-Ez(ix,:,iz))
ENDDO
ENDDO

DO iz=1,NZP-1
By(NXP,:,iz)=By(NXP,:,iz)-c*dt/(2.0*dz)*(Ex(NXP,:,iz+1)-Ex(NXP,:,iz))+&
                          c*dt/(2.0*dx)*(bufR2(:,iz)-Ez(NXP,:,iz))
ENDDO

DO ix=1,NXP-1
By(ix,:,NZP)=By(ix,:,NZP)-c*dt/(2.0*dz)*(bufR1(ix,:)-Ex(ix,:,NZP))+&
                          c*dt/(2.0*dx)*(Ez(ix+1,:,NZP)-Ez(ix,:,NZP))
ENDDO

By(NXP,:,NZP)=By(NXP,:,NZP)-c*dt/(2.0*dz)*(bufR1(NXP,:)-Ex(NXP,:,NZP))+&
                            c*dt/(2.0*dx)*(bufR2(:,NZP)-Ez(NXP,:,NZP))

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X

IF (xmaxp.EQ.xmax) THEN
   
   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Inside the conductor
   By(NXP,:,:)=0.0
   END IF

END IF
   
!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! Normal to conductor surface
   By(:,1,:)=0.0
   END IF
   
END IF

IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Normal to conductor surface
   By(:,NYP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Z

IF (zmaxp.EQ.zmax) THEN
   
   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   ! Inside the conductor
   By(:,:,NZP)=0.0
   END IF

END IF

!***********************************************************************
! Bz

ALLOCATE(bufS1(1:NXP,1:NZP),bufR1(1:NXP,1:NZP))
ALLOCATE(bufS2(1:NYP,1:NZP),bufR2(1:NYP,1:NZP))

bufS1=Ex(:,2,:)
bufS2=Ey(2,:,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,COMM,stat,ierr)
ENDIF

DO ix=1,NXP-1
DO iy=1,NYP-1
Bz(ix,iy,:)=Bz(ix,iy,:)-c*dt/(2.0*dx)*(Ey(ix+1,iy,:)-Ey(ix,iy,:))+&
                        c*dt/(2.0*dy)*(Ex(ix,iy+1,:)-Ex(ix,iy,:))
ENDDO
ENDDO

DO iy=1,NYP-1
Bz(NXP,iy,:)=Bz(NXP,iy,:)-c*dt/(2.0*dx)*(bufR2(iy,:)-Ey(NXP,iy,:))+&
                          c*dt/(2.0*dy)*(Ex(NXP,iy+1,:)-Ex(NXP,iy,:))
ENDDO

DO ix=1,NXP-1
Bz(ix,NYP,:)=Bz(ix,NYP,:)-c*dt/(2.0*dx)*(Ey(ix+1,NYP,:)-Ey(ix,NYP,:))+&
                          c*dt/(2.0*dy)*(bufR1(ix,:)-Ex(ix,NYP,:))
ENDDO

Bz(NXP,NYP,:)=Bz(NXP,NYP,:)-c*dt/(2.0*dx)*(bufR2(NYP,:)-Ey(NXP,NYP,:))+&
                            c*dt/(2.0*dy)*(bufR1(NXP,:)-Ex(NXP,NYP,:))

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X
   
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Inside the conductor
   Bz(NXP,:,:)=0.0
   END IF

END IF
   
!***********************************************************************
! Check boundary conditions along Y

IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Inside the conductor
   Bz(:,NYP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Z

IF (zminp.EQ.zmin) THEN

   IF (BOUND_FIELD_ZMIN.EQ."METAL") THEN
   ! Normal to conductor surface
   Bz(:,:,1)=0.0
   END IF
   
END IF

IF (zmaxp.EQ.zmax) THEN

   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   ! Normal to conductor surface
   Bz(:,:,NZP)=0.0
   END IF

END IF

!***********************************************************************

END SUBROUTINE PUSH_BHALF

!***********************************************************************
! Subroutine FIELDS_NODES
! This subroutine computes E and B at the nodes

! INPUT: 
! - Bx: x-component of initial B on the Yee grid at time t
! - By: y-component of initial B on the Yee grid at time t
! - Bz: z-component of initial B on the Yee grid at time t
! - Ex: x-component of initial E on the Yee grid at time t
! - Ey: y-component of initial E on the Yee grid at time t
! - Ez: z-component of initial E on the Yee grid at time t
! - Bxg: x-component of initial B at the nodes at time t
! - Byg: y-component of initial B at the nodes at time t
! - Bzg: z-component of initial B at the nodes at time t
! - Exg: x-component of initial E at the nodes at time t
! - Eyg: y-component of initial E at the nodes at time t
! - Ezg: z-component of initial E at the nodes at time t
! - xgp,ygp,zgp: Local grid
!
! OUTPUT: Ex,Ey,Ez,Bx,By,Bz at time t at the nodes
!
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE FIELDS_NODES(Bx,By,Bz,Ex,Ey,Ez,Bxg,Byg,Bzg,Exg,Eyg,Ezg,&
                        xgp,ygp,zgp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bx,By,Bz,Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bxg,Byg,Bzg,Exg,Eyg,Ezg
DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp
DOUBLE PRECISION                               :: xminp,xmaxp,yminp,ymaxp,zminp,zmaxp
INTEGER, DIMENSION(26)                         :: ngh
INTEGER, PARAMETER                             :: tag1=1,tag2=2,tag3=3
INTEGER                                        :: id,COMM,ierr
DOUBLE PRECISION, ALLOCATABLE                  :: bufS1(:,:),BufR1(:,:)
DOUBLE PRECISION, ALLOCATABLE                  :: bufS2(:,:),BufR2(:,:)
DOUBLE PRECISION, ALLOCATABLE                  :: bufS3(:),BufR3(:)
! Loop indexes
INTEGER :: ix,iy,iz
!***********************************************************************

xminp=xgp(1)
xmaxp=xgp(NXP)
yminp=ygp(1)
ymaxp=ygp(NYP)
zminp=zgp(1)
zmaxp=zgp(NZP)

!***********************************************************************
! Bxg

ALLOCATE(bufS1(1:NXP,1:NZP),bufR1(1:NXP,1:NZP))
ALLOCATE(bufS2(1:NXP,1:NYP),bufR2(1:NXP,1:NYP))
ALLOCATE(bufS3(1:NXP),bufR3(1:NXP))

bufS1=Bx(:,NYP-1,:)
bufS2=Bx(:,:,NZP-1)
bufS3=Bx(:,NYP-1,NZP-1)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag2,&
                  bufR2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag2,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS3,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag3,&
                  bufR3,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag3,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag2,&
                  bufR2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag2,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS3,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag3,&
                  bufR3,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag3,COMM,stat,ierr)
ENDIF

DO iy=2,NYP
DO iz=2,NZP
Bxg(:,iy,iz)=(Bx(:,iy,iz)+Bx(:,iy-1,iz)+Bx(:,iy,iz-1)+Bx(:,iy-1,iz-1))/4.0
ENDDO
ENDDO

DO iy=2,NYP
Bxg(:,iy,1)=(Bx(:,iy,1)+Bx(:,iy-1,1)+bufR2(:,iy)+bufR2(:,iy-1))/4.0
ENDDO

DO iz=2,NZP
Bxg(:,1,iz)=(Bx(:,1,iz)+BufR1(:,iz)+Bx(:,1,iz-1)+bufR1(:,iz-1))/4.0
ENDDO

Bxg(:,1,1)=(Bx(:,1,1)+bufR1(:,1)+bufR2(:,1)+bufR3(:))/4.0
   
!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN

   DO iz=2,NZP
   Bxg(:,1,iz)=(Bx(:,1,iz)+0.0+Bx(:,1,iz-1)+0.0)/4.0
   ENDDO
   
   Bxg(:,1,1)=(Bx(:,1,1)+0.0+BufR2(:,1)+0.0)/4.0
      
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Z
   
IF (zminp.EQ.zmin) THEN

   IF (BOUND_FIELD_ZMIN.EQ."METAL") THEN
   
   DO iy=2,NYP
   Bxg(:,iy,1)=(Bx(:,iy,1)+Bx(:,iy-1,1)+0.0+0.0)/4.0
   ENDDO
   
   Bxg(:,1,1)=(Bx(:,1,1)+bufR1(:,1)+0.0+0.0)/4.0
   
   END IF
   
END IF

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)
DEALLOCATE(bufS3,bufR3)

!***********************************************************************
! Byg

ALLOCATE(bufS1(1:NYP,1:NZP),bufR1(1:NYP,1:NZP))
ALLOCATE(bufS2(1:NXP,1:NYP),bufR2(1:NXP,1:NYP))
ALLOCATE(bufS3(1:NYP),bufR3(1:NYP))

bufS1=By(NXP-1,:,:)
bufS2=By(:,:,NZP-1)
bufS3=By(NXP-1,:,NZP-1)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag1,&
                  bufR1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag2,&
                  bufR2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag2,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS3,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag3,&
                  bufR3,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag3,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag1,&
                  bufR1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag2,&
                  bufR2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag2,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS3,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag3,&
                  bufR3,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag3,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
DO iz=2,NZP
Byg(ix,:,iz)=(By(ix,:,iz)+By(ix-1,:,iz)+By(ix,:,iz-1)+By(ix-1,:,iz-1))/4.0
ENDDO
ENDDO

DO ix=2,NXP
Byg(ix,:,1)=(By(ix,:,1)+By(ix-1,:,1)+bufR2(ix,:)+bufR2(ix-1,:))/4.0
ENDDO

DO iz=2,NZP
Byg(1,:,iz)=(By(1,:,iz)+bufR1(:,iz)+By(1,:,iz-1)+bufR1(:,iz-1))/4.0
ENDDO

Byg(1,:,1)=(By(1,:,1)+bufR1(:,1)+bufR2(1,:)+bufR3(:))/4.0

!***********************************************************************
! Check boundary conditions along X
   
IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   
   DO iz=2,NZP
   Byg(1,:,iz)=(By(1,:,iz)+0.0+By(1,:,iz-1)+0.0)/4.0
   ENDDO
   
   Byg(1,:,1)=(By(1,:,1)+0.0+bufR2(1,:)+0.0)/4.0
   
   END IF
   
END IF
   
!***********************************************************************
! Check boundary conditions along Z

IF (zminp.EQ.zmin) THEN

   IF (BOUND_FIELD_ZMIN.EQ."METAL") THEN

   DO ix=2,NXP
   Byg(ix,:,1)=(By(ix,:,1)+By(ix-1,:,1)+0.0+0.0)/4.0
   ENDDO
   
   Byg(1,:,1)=(By(1,:,1)+BufR1(:,1)+0.0+0.0)/4.0
      
   END IF

END IF

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)
DEALLOCATE(bufS3,bufR3)

!***********************************************************************
! Bzg

ALLOCATE(bufS1(1:NXP,1:NZP),bufR1(1:NXP,1:NZP))
ALLOCATE(bufS2(1:NYP,1:NZP),bufR2(1:NYP,1:NZP))
ALLOCATE(bufS3(1:NZP),bufR3(1:NZP))

bufS1=Bz(:,NYP-1,:)
bufS2=Bz(NXP-1,:,:)
bufS3=Bz(NXP-1,NYP-1,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS3,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag3,&
                  bufR3,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag3,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS3,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag3,&
                  bufR3,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag3,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
DO iy=2,NYP
Bzg(ix,iy,:)=(Bz(ix,iy,:)+Bz(ix-1,iy,:)+Bz(ix,iy-1,:)+Bz(ix-1,iy-1,:))/4.0
ENDDO
ENDDO

DO ix=2,NXP
Bzg(ix,1,:)=(Bz(ix,1,:)+Bz(ix-1,1,:)+bufR1(ix,:)+bufR1(ix-1,:))/4.0
ENDDO

DO iy=2,NYP
Bzg(1,iy,:)=(Bz(1,iy,:)+bufR2(iy,:)+Bz(1,iy-1,:)+bufR2(iy-1,:))/4.0
ENDDO

Bzg(1,1,:)=(Bz(1,1,:)+BufR2(1,:)+bufR1(1,:)+bufR3(:))/4.0

!***********************************************************************
! Check boundary conditions along X
   
IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   
   DO iy=2,NYP
   Bzg(1,iy,:)=(Bz(1,iy,:)+0.0+Bz(1,iy-1,:)+0.0)/4.0
   ENDDO
   
   Bzg(1,1,:)=(Bz(1,1,:)+0.0+bufR1(1,:)+0.0)/4.0
   
   END IF
   
END IF
   
!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN

   DO ix=2,NXP
   Bzg(ix,1,:)=(Bz(ix,1,:)+Bz(ix-1,1,:)+0.0+0.0)/4.0
   ENDDO
   
   Bzg(1,1,:)=(Bz(1,1,:)+BufR2(1,:)+0.0+0.0)/4.0
      
   END IF

END IF

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)
DEALLOCATE(bufS3,bufR3)

!***********************************************************************
! Exg

ALLOCATE(bufS2(1:NYP,1:NZP),bufR2(1:NYP,1:NZP))

bufS2=Ex(NXP-1,:,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,&
                  bufR2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
Exg(ix,:,:)=(Ex(ix,:,:)+Ex(ix-1,:,:))/2.0
ENDDO

Exg(1,:,:)=(Ex(1,:,:)+bufR2(:,:))/2.0

!***********************************************************************
! Check boundary conditions along X

IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! At the conductor surface
   Exg(1,:,:)=(Ex(1,:,:)+0.0)/2.0
   END IF

END IF

DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Eyg

ALLOCATE(bufS1(1:NXP,1:NZP),bufR1(1:NXP,1:NZP))

bufS1=Ey(:,NYP-1,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufR1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)
ENDIF

DO iy=2,NYP
Eyg(:,iy,:)=(Ey(:,iy,:)+Ey(:,iy-1,:))/2.0
ENDDO

Eyg(:,1,:)=(Ey(:,1,:)+bufR1(:,:))/2.0

!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! At the conductor surface
   Eyg(:,1,:)=(Ey(:,1,:)+0.0)/2.0
   END IF

END IF

DEALLOCATE(bufS1,bufR1)

!***********************************************************************
!Ezg

ALLOCATE(bufS1(1:NXP,1:NYP),bufR1(1:NXP,1:NYP))

bufS1=Ez(:,:,NZP-1)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag1,&
                  bufR1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag1,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag1,&
                  bufR1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag1,COMM,stat,ierr)
ENDIF

DO iz=2,NZP
Ezg(:,:,iz)=(Ez(:,:,iz)+Ez(:,:,iz-1))/2.0
ENDDO

Ezg(:,:,1)=(Ez(:,:,1)+bufR1(:,:))/2.0

!***********************************************************************
! Check boundary conditions along Z

IF (zminp.EQ.zmin) THEN

   IF (BOUND_FIELD_ZMIN.EQ."METAL") THEN
   ! At the conductor surface
   Ezg(:,:,1)=(Ez(:,:,1)+0.0)/2.0
   END IF

END IF

DEALLOCATE(bufS1,bufR1)

END SUBROUTINE FIELDS_NODES

!***********************************************************************
! Subroutine CORRECT_EFIELD
! This subroutine correct the electric field to ensure the conservation 
! of charge, or to ensure that div(E)=4*pi*rho. Poission equation is solved
! using an iterative method (Gauss-Seidel method), with 7 points.

! INPUT: 
! - Ex: x-component of E on the Yee grid at time t
! - Ey: y-component of E on the Yee grid at time t
! - Ez: z-component of E on the Yee grid at time t
! - rho: charge density
!
! OUTPUT: corrected Ex,Ey,Ez at time t+dt
!
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE CORRECT_EFIELD(Ex,Ey,Ez,xgp,ygp,zgp,rho,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: stat
INTEGER                                        :: id,COMM,ierr
INTEGER, DIMENSION(26)                         :: ngh
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: rho,phi

DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp

DOUBLE PRECISION, DIMENSION(1:NYP,1:NZP)       :: bufSE1,bufSW1,bufRE1,bufRW1
DOUBLE PRECISION, DIMENSION(1:NYP,1:NZP)       :: bufSE2,bufRW2

DOUBLE PRECISION, DIMENSION(1:NXP,1:NZP)       :: bufSN1,bufSS1,bufRN1,bufRS1
DOUBLE PRECISION, DIMENSION(1:NXP,1:NZP)       :: bufSN2,bufRS2

DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP)       :: bufSF1,bufSB1,bufRF1,bufRB1
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP)       :: bufSF2,bufRB2

DOUBLE PRECISION, ALLOCATABLE                  :: bufS(:,:),bufR(:,:)

INTEGER, PARAMETER                       :: tag1=1,tag2=2,tag3=3,tag4=4
INTEGER, PARAMETER                       :: tag5=5,tag6=6,tag7=7,tag8=8,tag9=9

DOUBLE PRECISION                         :: denom,xminp,yminp,xmaxp,ymaxp,zminp,zmaxp

! Loop indexes
INTEGER :: ix,iy,iz,iit

!***********************************************************************

xminp=xgp(1)
xmaxp=xgp(NXP)
yminp=ygp(1)
ymaxp=ygp(NYP)
zminp=zgp(1)
zmaxp=zgp(NZP)

!***********************************************************************
! denominator
denom=dx*dx*dy*dy+dx*dx*dz*dz+dy*dy*dz*dz

bufSE2=Ex(NXP-1,:,:)
bufSN2=Ey(:,NYP-1,:)
bufSF2=Ez(:,:,NZP-1)

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSE2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag1,&
                  bufRW2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN2,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag2,&
                  bufRS2,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag3,&
                  bufRB2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag3,COMM,stat,ierr)
                  
ELSE

CALL MPI_SENDRECV(bufSE2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag1,&
                  bufRW2,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN2,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag2,&
                  bufRS2,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag3,&
                  bufRB2,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag3,COMM,stat,ierr)
                  
END IF

! Electric potential phi
phi=0.0

!***********************************************************************
! Beginning iteration
!***********************************************************************

DO iit=1,NIT

bufSE1=phi(NXP-1,:,:)
bufSW1=phi(2,:,:)
bufSN1=phi(:,NYP-1,:)
bufSS1=phi(:,2,:)
bufSF1=phi(:,:,NZP-1)
bufSB1=phi(:,:,2)

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSE1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag4,&
                  bufRW1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag5,&
                  bufRE1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag6,&
                  bufRS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag7,&
                  bufRN1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag8,&
                  bufRB1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSB1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag9,&
                  bufRF1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag9,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(bufSE1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag4,&
                  bufRW1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag5,&
                  bufRE1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag6,&
                  bufRS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag6,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag7,&
                  bufRN1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag8,&
                  bufRB1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSB1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag9,&
                  bufRF1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag9,COMM,stat,ierr)

ENDIF

  DO ix=2,NXP-1
    DO iy=2,NYP-1
      DO iz=2,NZP-1
        phi(ix,iy,iz)=0.5/denom*(&
                      (phi(ix+1,iy,iz)+phi(ix-1,iy,iz))*dy*dy*dz*dz+&
                      (phi(ix,iy+1,iz)+phi(ix,iy-1,iz))*dx*dx*dz*dz+&
                      (phi(ix,iy,iz+1)+phi(ix,iy,iz-1))*dx*dx*dy*dy+&
                      (4.0*pi*rho(ix,iy,iz)-&
                      ((Ex(ix,iy,iz)-Ex(ix-1,iy,iz))/dx+&
                       (Ey(ix,iy,iz)-Ey(ix,iy-1,iz))/dy+&
                       (Ez(ix,iy,iz)-Ez(ix,iy,iz-1))/dz))*&
                      dx*dx*dy*dy*dz*dz)
      ENDDO
    ENDDO
  ENDDO

!***********************************************************************
! Faces of the cube

  ! Surface ix=1
  DO iy=2,NYP-1
    DO iz=2,NZP-1
      phi(1,iy,iz)=0.5/denom*(&
                   (phi(2,iy,iz)+bufRW1(iy,iz))*dy*dy*dz*dz+&
                   (phi(1,iy+1,iz)+phi(1,iy-1,iz))*dx*dx*dz*dz+&
                   (phi(1,iy,iz+1)+phi(1,iy,iz-1))*dx*dx*dy*dy+&
                   (4.0*pi*rho(1,iy,iz)-&
                   ((Ex(1,iy,iz)-bufRW2(iy,iz))/dx+&
                    (Ey(1,iy,iz)-Ey(1,iy-1,iz))/dy+&
                    (Ez(1,iy,iz)-Ez(1,iy,iz-1))/dz))*&
                   dx*dx*dy*dy*dz*dz)
    ENDDO
  ENDDO

  ! Surface iy=1
  DO ix=2,NXP-1
    DO iz=2,NZP-1
      phi(ix,1,iz)=0.5/denom*(&
                   (phi(ix+1,1,iz)+phi(ix-1,1,iz))*dy*dy*dz*dz+&
                   (phi(ix,2,iz)+bufRS1(ix,iz))*dx*dx*dz*dz+&
                   (phi(ix,1,iz+1)+phi(ix,1,iz-1))*dx*dx*dy*dy+&
                   (4.0*pi*rho(ix,1,iz)-&
                   ((Ex(ix,1,iz)-Ex(ix-1,1,iz))/dx+&
                    (Ey(ix,1,iz)-bufRS2(ix,iz))/dy+&
                    (Ez(ix,1,iz)-Ez(ix,1,iz-1))/dz))*&
                   dx*dx*dy*dy*dz*dz)
    ENDDO
  ENDDO

  ! Surface iz=1
  DO ix=2,NXP-1
    DO iy=2,NYP-1
      phi(ix,iy,1)=0.5/denom*(&
                   (phi(ix+1,iy,1)+phi(ix-1,iy,1))*dy*dy*dz*dz+&
                   (phi(ix,iy+1,1)+phi(ix,iy-1,1))*dx*dx*dz*dz+&
                   (phi(ix,iy,2)+bufRB1(ix,iy))*dx*dx*dy*dy+&
                   (4.0*pi*rho(ix,iy,1)-&
                   ((Ex(ix,iy,1)-Ex(ix-1,iy,1))/dx+&
                    (Ey(ix,iy,1)-Ey(ix,iy-1,1))/dy+&
                    (Ez(ix,iy,1)-bufRB2(ix,iy))/dz))*&
                   dx*dx*dy*dy*dz*dz)
    ENDDO
  ENDDO

  ! Surface ix=NXP
  DO iy=2,NYP-1
    DO iz=2,NZP-1
      phi(NXP,iy,iz)=0.5/denom*(&
                    (bufRE1(iy,iz)+phi(NXP-1,iy,iz))*dy*dy*dz*dz+&
                    (phi(NXP,iy+1,iz)+phi(NXP,iy-1,iz))*dx*dx*dz*dz+&
                    (phi(NXP,iy,iz+1)+phi(NXP,iy,iz-1))*dx*dx*dy*dy+&
                    (4.0*pi*rho(NXP,iy,iz)-&
                    ((Ex(NXP,iy,iz)-Ex(NXP-1,iy,iz))/dx+&
                     (Ey(NXP,iy,iz)-Ey(NXP,iy-1,iz))/dy+&
                     (Ez(NXP,iy,iz)-Ez(NXP,iy,iz-1))/dz))*&
                    dx*dx*dy*dy*dz*dz)
    ENDDO
  ENDDO

  ! Surface iy=NYP
  DO ix=2,NXP-1
    DO iz=2,NZP-1
      phi(ix,NYP,iz)=0.5/denom*(&
                     (phi(ix+1,NYP,iz)+phi(ix-1,NYP,iz))*dy*dy*dz*dz+&
                     (bufRN1(ix,iz)+phi(ix,NYP-1,iz))*dx*dx*dz*dz+&
                     (phi(ix,NYP,iz+1)+phi(ix,NYP,iz-1))*dx*dx*dy*dy+&
                     (4.0*pi*rho(ix,NYP,iz)-&
                      ((Ex(ix,NYP,iz)-Ex(ix-1,NYP,iz))/dx+&
                       (Ey(ix,NYP,iz)-Ey(ix,NYP-1,iz))/dy+&
                       (Ez(ix,NYP,iz)-Ez(ix,NYP,iz-1))/dz))*&
                     dx*dx*dy*dy*dz*dz)
    ENDDO
  ENDDO

  ! Surface iz=NZP
  DO ix=2,NXP-1
    DO iy=2,NYP-1
      phi(ix,iy,NZP)=0.5/denom*(&
                     (phi(ix+1,iy,NZP)+phi(ix-1,iy,NZP))*dy*dy*dz*dz+&
                     (phi(ix,iy+1,NZP)+phi(ix,iy-1,NZP))*dx*dx*dz*dz+&
                     (bufRF1(ix,iy)+phi(ix,iy,NZP-1))*dx*dx*dy*dy+&
                     (4.0*pi*rho(ix,iy,NZP)-&
                     ((Ex(ix,iy,NZP)-Ex(ix-1,iy,NZP))/dx+&
                      (Ey(ix,iy,NZP)-Ey(ix,iy-1,NZP))/dy+&
                      (Ez(ix,iy,NZP)-Ez(ix,iy,NZP-1))/dz))*&
                     dx*dx*dy*dy*dz*dz)
    ENDDO
  ENDDO

!***********************************************************************
! Aretes of the cube

  DO ix=2,NXP-1
    phi(ix,1,1)=0.5/denom*(&
                (phi(ix+1,1,1)+phi(ix-1,1,1))*dy*dy*dz*dz+&
                (phi(ix,2,1)+bufRS1(ix,1))*dx*dx*dz*dz+&
                (phi(ix,1,2)+bufRB1(ix,1))*dx*dx*dy*dy+&
                (4.0*pi*rho(ix,1,1)-&
                ((Ex(ix,1,1)-Ex(ix-1,1,1))/dx+&
                 (Ey(ix,1,1)-bufRS2(ix,1))/dy+&
                 (Ez(ix,1,1)-bufRB2(ix,1))/dz))*&
                dx*dx*dy*dy*dz*dz)
  ENDDO

  DO ix=2,NXP-1
    phi(ix,NYP,NZP)=0.5/denom*(&
                    (phi(ix+1,NYP,NZP)+phi(ix-1,NYP,NZP))*dy*dy*dz*dz+&
                    (bufRN1(ix,NZP)+phi(ix,NYP-1,NZP))*dx*dx*dz*dz+&
                    (bufRF1(ix,NYP)+phi(ix,NYP,NZP-1))*dx*dx*dy*dy+&
                    (4.0*pi*rho(ix,NYP,NZP)-&
                    ((Ex(ix,NYP,NZP)-Ex(ix-1,NYP,NZP))/dx+&
                     (Ey(ix,NYP,NZP)-Ey(ix,NYP-1,NZP))/dy+&
                     (Ez(ix,NYP,NZP)-Ez(ix,NYP,NZP-1))/dz))*&
                    dx*dx*dy*dy*dz*dz)
  ENDDO

  DO ix=2,NXP-1
    phi(ix,NYP,1)=0.5/denom*(&
                  (phi(ix+1,NYP,1)+phi(ix-1,NYP,1))*dy*dy*dz*dz+&
                  (bufRN1(ix,1)+phi(ix,NYP-1,1))*dx*dx*dz*dz+&
                  (phi(ix,NYP,2)+bufRB1(ix,NYP))*dx*dx*dy*dy+&
                  (4.0*pi*rho(ix,NYP,1)-&
                  ((Ex(ix,NYP,1)-Ex(ix-1,NYP,1))/dx+&
                   (Ey(ix,NYP,1)-Ey(ix,NYP-1,1))/dy+&
                   (Ez(ix,NYP,1)-bufRB2(ix,NYP))/dz))*&
                  dx*dx*dy*dy*dz*dz)
  ENDDO

  DO ix=2,NXP-1
    phi(ix,1,NZP)=0.5/denom*(&
                  (phi(ix+1,1,NZP)+phi(ix-1,1,NZP))*dy*dy*dz*dz+&
                  (phi(ix,2,NZP)+bufRS1(ix,NZP))*dx*dx*dz*dz+&
                  (bufRF1(ix,1)+phi(ix,1,NZP-1))*dx*dx*dy*dy+&
                  (4.0*pi*rho(ix,1,NZP)-&
                  ((Ex(ix,1,NZP)-Ex(ix-1,1,NZP))/dx+&
                   (Ey(ix,1,NZP)-bufRS2(ix,NZP))/dy+&
                   (Ez(ix,1,NZP)-Ez(ix,1,NZP-1))/dz))*&
                  dx*dx*dy*dy*dz*dz)
  ENDDO

  DO iz=2,NZP-1
    phi(1,1,iz)=0.5/denom*(&
                (phi(2,1,iz)+bufRW1(1,iz))*dy*dy*dz*dz+&
                (phi(1,2,iz)+bufRS1(1,iz))*dx*dx*dz*dz+&
                (phi(1,1,iz+1)+phi(1,1,iz-1))*dx*dx*dy*dy+&
                (4.0*pi*rho(1,1,iz)-&
                ((Ex(1,1,iz)-bufRW2(1,iz))/dx+&
                 (Ey(1,1,iz)-bufRS2(1,iz))/dy+&
                 (Ez(1,1,iz)-Ez(1,1,iz-1))/dz))*&
                dx*dx*dy*dy*dz*dz)
  ENDDO

  DO iz=2,NZP-1
    phi(1,NYP,iz)=0.5/denom*(&
                  (phi(2,NYP,iz)+bufRW1(NYP,iz))*dy*dy*dz*dz+&
                  (bufRN1(1,iz)+phi(1,NYP-1,iz))*dx*dx*dz*dz+&
                  (phi(1,NYP,iz+1)+phi(1,NYP,iz-1))*dx*dx*dy*dy+&
                  (4.0*pi*rho(1,NYP,iz)-&
                  ((Ex(1,NYP,iz)-bufRW2(NYP,iz))/dx+&
                   (Ey(1,NYP,iz)-Ey(1,NYP-1,iz))/dy+&
                   (Ez(1,NYP,iz)-Ez(1,NYP,iz-1))/dz))*&
                  dx*dx*dy*dy*dz*dz)
  ENDDO

  DO iz=2,NZP-1
    phi(NXP,1,iz)=0.5/denom*(&
                  (bufRE1(1,iz)+phi(NXP-1,1,iz))*dy*dy*dz*dz+&
                  (phi(NXP,2,iz)+bufRS1(NXP,iz))*dx*dx*dz*dz+&
                  (phi(NXP,1,iz+1)+phi(NXP,1,iz-1))*dx*dx*dy*dy+&
                  (4.0*pi*rho(NXP,1,iz)-&
                  ((Ex(NXP,1,iz)-Ex(NXP-1,1,iz))/dx+&
                   (Ey(NXP,1,iz)-bufRS2(NXP,iz))/dy+&
                   (Ez(NXP,1,iz)-Ez(NXP,1,iz-1))/dz))*&
                  dx*dx*dy*dy*dz*dz)
  ENDDO

  DO iz=2,NZP-1
    phi(NXP,NYP,iz)=0.5/denom*(&
                    (bufRE1(NYP,iz)+phi(NXP-1,NYP,iz))*dy*dy*dz*dz+&
                    (bufRN1(NXP,iz)+phi(NXP,NYP-1,iz))*dx*dx*dz*dz+&
                    (phi(NXP,NYP,iz+1)+phi(NXP,NYP,iz-1))*dx*dx*dy*dy+&
                    (4.0*pi*rho(NXP,NYP,iz)-&
                    ((Ex(NXP,NYP,iz)-Ex(NXP-1,NYP,iz))/dx+&
                     (Ey(NXP,NYP,iz)-Ey(NXP,NYP-1,iz))/dy+&
                     (Ez(NXP,NYP,iz)-Ez(NXP,NYP,iz-1))/dz))*&
                    dx*dx*dy*dy*dz*dz)
  ENDDO

  DO iy=2,NYP-1
    phi(1,iy,1)=0.5/denom*(&
                (phi(2,iy,1)+bufRW1(iy,1))*dy*dy*dz*dz+&
                (phi(1,iy+1,1)+phi(1,iy-1,1))*dx*dx*dz*dz+&
                (phi(1,iy,2)+bufRB1(1,iy))*dx*dx*dy*dy+&
                (4.0*pi*rho(1,iy,1)-&
                ((Ex(2,iy,1)-bufRW2(iy,1))/dx+&
                 (Ey(1,iy+1,1)-Ey(1,iy-1,1))/dy+&
                 (Ez(1,iy,2)-bufRB2(1,iy))/dz))*&
                dx*dx*dy*dy*dz*dz)
  ENDDO

  DO iy=2,NYP-1
    phi(1,iy,NZP)=0.5/denom*(&
                  (phi(2,iy,NZP)+bufRW1(iy,NZP))*dy*dy*dz*dz+&
                  (phi(1,iy+1,NZP)+phi(1,iy-1,NZP))*dx*dx*dz*dz+&
                  (bufRF1(1,iy)+phi(1,iy,NZP-1))*dx*dx*dy*dy+&
                  (4.0*pi*rho(1,iy,NZP)-&
                  ((Ex(1,iy,NZP)-bufRW2(iy,NZP))/dx+&
                   (Ey(1,iy,NZP)-Ey(1,iy-1,NZP))/dy+&
                   (Ez(1,iy,NZP)-Ez(1,iy,NZP-1))/dz))*&
                  dx*dx*dy*dy*dz*dz)
  ENDDO

  DO iy=2,NYP-1
    phi(NXP,iy,1)=0.5/denom*(&
                  (bufRE1(iy,1)+phi(NXP-1,iy,1))*dy*dy*dz*dz+&
                  (phi(NXP,iy+1,1)+phi(NXP,iy-1,1))*dx*dx*dz*dz+&
                  (phi(NXP,iy,2)+bufRB1(NXP,iy))*dx*dx*dy*dy+&
                  (4.0*pi*rho(NXP,iy,1)-&
                  ((Ex(NXP,iy,1)-Ex(NXP-1,iy,1))/dx+&
                   (Ey(NXP,iy,1)-Ey(NXP,iy-1,1))/dy+&
                   (Ez(NXP,iy,1)-bufRB2(NXP,iy))/dz))*&
                  dx*dx*dy*dy*dz*dz)
  ENDDO

  DO iy=2,NYP-1
    phi(NXP,iy,NZP)=0.5/denom*(&
                    (bufRE1(iy,NZP)+phi(NXP-1,iy,NZP))*dy*dy*dz*dz+&
                    (phi(NXP,iy+1,NZP)+phi(NXP,iy-1,NZP))*dx*dx*dz*dz+&
                    (bufRF1(NXP,iy)+phi(NXP,iy,NZP-1))*dx*dx*dy*dy+&
                    (4.0*pi*rho(NXP,iy,NZP)-&
                    ((Ex(NXP,iy,NZP)-Ex(NXP-1,iy,NZP))/dx+&
                     (Ey(NXP,iy,NZP)-Ey(NXP,iy-1,NZP))/dy+&
                     (Ex(NXP,iy,NZP)-Ez(NXP,iy,NZP-1))/dz))*&
                    dx*dx*dy*dy*dz*dz)
  ENDDO
  
!***********************************************************************
! Corners

  phi(1,1,1)=0.5/denom*(&
             (phi(2,1,1)+bufRW1(1,1))*dy*dy*dz*dz+&
             (phi(1,2,1)+bufRS1(1,1))*dx*dx*dz*dz+&
             (phi(1,1,2)+bufRB1(1,1))*dx*dx*dy*dy+&
             (4.0*pi*rho(1,1,1)-&
             ((Ex(1,1,1)-bufRW2(1,1))/dx+&
              (Ey(1,1,1)-bufRS2(1,1))/dy+&
              (Ez(1,1,1)-bufRB2(1,1))/dz))*&
             dx*dx*dy*dy*dz*dz)

  phi(NXP,1,1)=0.5/denom*(&
               (bufRE1(1,1)+phi(NXP-1,1,1))*dy*dy*dz*dz+&
               (phi(NXP,2,1)+bufRS1(NXP,1))*dx*dx*dz*dz+&
               (phi(NXP,1,2)+bufRB1(NXP,1))*dx*dx*dy*dy+&
               (4.0*pi*rho(NXP,1,1)-&
               ((Ex(NXP,1,1)-Ex(NXP-1,1,1))/dx+&
                (Ey(NXP,1,1)-bufRS2(NXP,1))/dy+&
                (Ez(NXP,1,1)-bufRB2(NXP,1))/dz))*&
               dx*dx*dy*dy*dz*dz)

  phi(1,NYP,1)=0.5/denom*(&
               (phi(2,NYP,1)+bufRW1(NYP,1))*dy*dy*dz*dz+&
               (bufRN1(1,1)+phi(1,NYP-1,1))*dx*dx*dz*dz+&
               (phi(1,NYP,2)+bufRB1(1,NYP))*dx*dx*dy*dy+&
               (4.0*pi*rho(1,NYP,1)-&
               ((Ex(1,NYP,1)-bufRW2(NYP,1))/dx+&
                (Ey(1,NYP,1)-Ey(1,NYP-1,1))/dy+&
                (Ez(1,NYP,1)-bufRB2(1,NYP))/dz))*&
               dx*dx*dy*dy*dz*dz)

  phi(1,1,NZP)=0.5/denom*(&
               (phi(2,1,NZP)+bufRW1(1,NZP))*dy*dy*dz*dz+&
               (phi(1,2,NZP)+bufRS1(1,NZP))*dx*dx*dz*dz+&
               (bufRF1(1,1)+phi(1,1,NZP-1))*dx*dx*dy*dy+&
               (4.0*pi*rho(1,1,NZP)-&
               ((Ex(1,1,NZP)-bufRW2(1,NZP))/dx+&
                (Ey(1,1,NZP)-bufRS2(1,NZP))/dy+&
                (Ez(1,1,NZP)-Ez(1,1,NZP-1))/dz))*&
               dx*dx*dy*dy*dz*dz)

  phi(NXP,NYP,1)=0.5/denom*(&
                 (bufRE1(NYP,1)+phi(NXP-1,NYP,1))*dy*dy*dz*dz+&
                 (bufRN1(NXP,1)+phi(NXP,NYP-1,1))*dx*dx*dz*dz+&
                 (phi(NXP,NYP,2)+bufRB1(NXP,NYP))*dx*dx*dy*dy+&
                 (4.0*pi*rho(NXP,NYP,1)-&
                 ((Ex(NXP,NYP,1)-Ex(NXP-1,NYP,1))/dx+&
                  (Ey(NXP,NYP,1)-Ey(NXP,NYP-1,1))/dy+&
                  (Ez(NXP,NYP,1)-bufRB2(NXP,NYP))/dz))*&
                 dx*dx*dy*dy*dz*dz)

  phi(1,NYP,NZP)=0.5/denom*(&
                 (phi(2,NYP,NZP)+bufRW1(NYP,NZP))*dy*dy*dz*dz+&
                 (bufRN1(1,NZP)+phi(1,NYP-1,NZP))*dx*dx*dz*dz+&
                 (bufRF1(1,NYP)+phi(1,NYP,NZP-1))*dx*dx*dy*dy+&
                 (4.0*pi*rho(1,NYP,NZP)-&
                 ((Ex(1,NYP,NZP)-bufRW2(NYP,NZP))/dx+&
                  (Ey(1,NYP,NZP)-Ey(1,NYP-1,NZP))/dy+&
                  (Ez(1,NYP,NZP)-Ez(1,NYP,NZP-1))/dz))*&
                 dx*dx*dy*dy*dz*dz)

  phi(NXP,1,NZP)=0.5/denom*(&
                 (bufRE1(1,NZP)+phi(NXP-1,1,NZP))*dy*dy*dz*dz+&
                 (phi(NXP,2,NZP)+bufRS1(NXP,NZP))*dx*dx*dz*dz+&
                 (bufRF1(NXP,1)+phi(NXP,1,NZP-1))*dx*dx*dy*dy+&
                 (4.0*pi*rho(NXP,1,NZP)-&
                 ((Ex(NXP,1,NZP)-Ex(NXP-1,1,NZP))/dx+&
                  (Ey(NXP,1,NZP)-bufRS2(NXP,NZP))/dy+&
                  (Ez(NXP,1,NZP)-Ez(NXP,1,NZP-1))/dz))*&
                 dx*dx*dy*dy*dz*dz)

  phi(NXP,NYP,NZP)=0.5/denom*(&
                   (bufRE1(NYP,NZP)+phi(NXP-1,NYP,NZP))*dy*dy*dz*dz+&
                   (bufRN1(NXP,NZP)+phi(NXP,NYP-1,NZP))*dx*dx*dz*dz+&
                   (bufRF1(NXP,NYP)+phi(NXP,NYP,NZP-1))*dx*dx*dy*dy+&
                   (4.0*pi*rho(NXP,NYP,NZP)-&
                   ((Ex(NXP,NYP,NZP)-Ex(NXP-1,NYP,NZP))/dx+&
                    (Ey(NXP,NYP,NZP)-Ey(NXP,NYP-1,NZP))/dy+&
                    (Ez(NXP,NYP,NZP)-Ez(NXP,NYP,NZP-1))/dz))*&
                   dx*dx*dy*dy*dz*dz)

!***********************************************************************
! Check boundary conditions along X

IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   phi(1,:,:)=0.0
   END IF

END IF

IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   phi(NXP,:,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   phi(:,1,:)=0.0
   END IF

END IF
   
IF (ymaxp.EQ.ymax) THEN
   
   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   phi(:,NYP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Z

IF (zminp.EQ.zmin) THEN

   IF (BOUND_FIELD_ZMIN.EQ."METAL") THEN
   phi(:,:,1)=0.0
   END IF

END IF
   
IF (zmaxp.EQ.zmax) THEN
   
   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   phi(:,:,NZP)=0.0
   END IF

END IF

!***********************************************************************

ENDDO

!***********************************************************************
! Corrected nodal electric field
!***********************************************************************

bufSW1(:,:)=phi(2,:,:)
bufSS1(:,:)=phi(:,2,:)
bufSB1(:,:)=phi(:,:,2)

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSW1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag1,&
                  bufRE1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag2,&
                  bufRN1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag3,&
                  bufRB1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag3,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(bufSW1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag1,&
                  bufRE1,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag2,&
                  bufRN1,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag3,&
                  bufRB1,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag3,COMM,stat,ierr)

END IF

!***********************************************************************
! Ex

DO ix=1,NXP-1
Ex(ix,:,:)=Ex(ix,:,:)-(phi(ix+1,:,:)-phi(ix,:,:))/dx
ENDDO

!***********************************************************************
! Check boundary conditions along X

IF (xmaxp.EQ.xmax.AND.BOUND_FIELD_XMAX.NE."PERIODIC") THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   Ex(NXP,:,:)=Ex(NXP,:,:)
   END IF

ELSE

! ix=NXP
Ex(NXP,:,:)=Ex(NXP,:,:)-(bufRE1(:,:)-phi(NXP,:,:))/dx

END IF

!***********************************************************************
! Ey

DO iy=1,NYP-1
Ey(:,iy,:)=Ey(:,iy,:)-(phi(:,iy+1,:)-phi(:,iy,:))/dy
ENDDO

!***********************************************************************
! Check boundary conditions along Y

IF (ymaxp.EQ.ymax.AND.BOUND_FIELD_YMAX.NE."PERIODIC") THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   Ey(:,NYP,:)=Ey(:,NYP,:)
   END IF

ELSE
   
! iy=NYP
Ey(:,NYP,:)=Ey(:,NYP,:)-(bufRN1(:,:)-phi(:,NYP,:))/dy

END IF

!***********************************************************************
! Ez

DO iz=1,NZP-1
Ez(:,:,iz)=Ez(:,:,iz)-(phi(:,:,iz+1)-phi(:,:,iz))/dz
ENDDO

!***********************************************************************
! Check boundary conditions along Z

IF (zmaxp.EQ.zmax.AND.BOUND_FIELD_ZMAX.NE."PERIODIC") THEN

   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   Ez(:,:,NZP)=Ez(:,:,NZP)
   END IF

ELSE

! iz=NZP
Ez(:,:,NZP)=Ez(:,:,NZP)-(bufRF1(:,:)-phi(:,:,NZP))/dz

END IF

END SUBROUTINE CORRECT_EFIELD

!***********************************************************************
! Subroutine FILTER_FIELD
! This subroutine operates an isotropic 3D digital filtering of the field, using
! 27 points.
!
! The 27-points Filter matrix expresses as [Birdsall & Langdon, p.441]
!
!      |a    2a     a|      |2a    4a    2a|        |a    2a     a|
!  M = |2a   4a    2a|      |4a  1-56a   4a|        |2a   4a    2a|
!      |a    2a     a|      |2a    4a    2a|        |a    2a     a|
!           BACK                 CENTER                 FRONT
!
! INPUT: 
! - Field: Field component to filter
!
! OUTPUT: Filtered field Field
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE FILTER_FIELD(Field,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: stat
INTEGER                                        :: id,COMM,ierr
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Field,Field_temp
INTEGER, DIMENSION(26)                         :: ngh
INTEGER, PARAMETER :: tag1=1,tag2=2,tag3=3,tag4=4
INTEGER, PARAMETER :: tag5=5,tag6=6,tag7=7,tag8=8
INTEGER, PARAMETER :: tag9=9,tag10=10,tag11=11,tag12=12
INTEGER, PARAMETER :: tag13=13,tag14=14,tag15=15,tag16=16
INTEGER, PARAMETER :: tag17=17,tag18=18,tag19=19,tag20=20
INTEGER, PARAMETER :: tag21=21,tag22=22,tag23=23,tag24=24
INTEGER, PARAMETER :: tag25=25,tag26=26

DOUBLE PRECISION, DIMENSION(1:NYP,1:NZP) :: bufSE,bufSW,bufRE,bufRW
DOUBLE PRECISION, DIMENSION(1:NXP,1:NZP) :: bufSN,bufSS,bufRN,bufRS
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: bufSF,bufSB,bufRF,bufRB

DOUBLE PRECISION, DIMENSION(1:NZP)       :: bufSNE,bufSSE,bufRNE,bufRSE
DOUBLE PRECISION, DIMENSION(1:NZP)       :: bufSSW,bufSNW,bufRSW,bufRNW

DOUBLE PRECISION, DIMENSION(1:NYP)       :: bufSBE,bufSBW,bufRBE,bufRBW
DOUBLE PRECISION, DIMENSION(1:NYP)       :: bufSFE,bufSFW,bufRFE,bufRFW

DOUBLE PRECISION, DIMENSION(1:NXP)       :: bufSBN,bufSBS,bufRBN,bufRBS
DOUBLE PRECISION, DIMENSION(1:NXP)       :: bufSFN,bufSFS,bufRFN,bufRFS

DOUBLE PRECISION                         ::bufSBNE,bufSBSE,bufRBNE,bufRBSE
DOUBLE PRECISION                         ::bufSBSW,bufSBNW,bufRBSW,bufRBNW
DOUBLE PRECISION                         ::bufSFNE,bufSFSE,bufRFNE,bufRFSE
DOUBLE PRECISION                         ::bufSFSW,bufSFNW,bufRFSW,bufRFNW

INTEGER :: ix,iy,iz

!***********************************************************************

! Faces of the cube
bufSE(:,:)=Field(NXP-1,:,:)
bufSW(:,:)=Field(2,:,:)
bufSN(:,:)=Field(:,NYP-1,:)
bufSS(:,:)=Field(:,2,:)
bufSF(:,:)=Field(:,:,NZP-1)
bufSB(:,:)=Field(:,:,2)

! Aretes

bufSNE(:)=Field(NXP-1,NYP-1,:)
bufSSE(:)=Field(NXP-1,2,:)
bufSSW(:)=Field(2,2,:)
bufSNW(:)=Field(2,NYP-1,:)

bufSBE(:)=Field(NXP-1,:,2)
bufSBW(:)=Field(2,:,2)
bufSFE(:)=Field(NXP-1,:,NZP-1)
bufSFW(:)=Field(2,:,NZP-1)

bufSBN(:)=Field(:,NYP-1,2)
bufSBS(:)=Field(:,2,2)
bufSFN(:)=Field(:,NYP-1,NZP-1)
bufSFS(:)=Field(:,2,NZP-1)

! Corners
bufSBNE=Field(NXP-1,NYP-1,2)
bufSBSE=Field(NXP-1,2,2)
bufSBSW=Field(2,2,2)
bufSBNW=Field(2,NYP-1,2)
bufSFNE=Field(NXP-1,NYP-1,NZP-1)
bufSFSE=Field(NXP-1,2,NZP-1)
bufSFSW=Field(2,2,NZP-1)
bufSFNW=Field(2,NYP-1,NZP-1)

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSE,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag1,&
                  bufRW,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,&
                  bufRE,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSN,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag3,&
                  bufRS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag3,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag4,&
                  bufRN,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag4,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSF,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag5,&
                  bufRB,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag5,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSB,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag6,&
                  bufRF,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag6,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSNE,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag7,&
                  bufRSW,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag7,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSSW,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag8,&
                  bufRNE,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag8,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSSE,NZP,MPI_DOUBLE_PRECISION,ngh(15),tag9,&
                  bufRNW,NZP,MPI_DOUBLE_PRECISION,ngh(17),tag9,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSNW,NZP,MPI_DOUBLE_PRECISION,ngh(17),tag10,&
                  bufRSE,NZP,MPI_DOUBLE_PRECISION,ngh(15),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag11,&
                  bufRFW,NYP,MPI_DOUBLE_PRECISION,ngh(21),tag11,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSFW,NYP,MPI_DOUBLE_PRECISION,ngh(21),tag12,&
                  bufRBE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag12,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSBW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag13,&
                  bufRFE,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag13,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSFE,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag14,&
                  bufRBW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag14,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSBN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag15,&
                  bufRFS,NXP,MPI_DOUBLE_PRECISION,ngh(20),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFS,NXP,MPI_DOUBLE_PRECISION,ngh(20),tag16,&
                  bufRBN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag16,COMM,stat,ierr)                  
                  
CALL MPI_SENDRECV(bufSBS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag17,&
                  bufRFN,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag17,COMM,stat,ierr)  
                  
CALL MPI_SENDRECV(bufSFN,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag18,&
                  bufRBS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag18,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSBNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag19,&
                  bufRFSW,1,MPI_DOUBLE_PRECISION,ngh(24),tag19,COMM,stat,ierr)         
                  
CALL MPI_SENDRECV(bufSFSW,1,MPI_DOUBLE_PRECISION,ngh(24),tag20,&
                  bufRBNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag20,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSBSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag21,&
                  bufRFNW,1,MPI_DOUBLE_PRECISION,ngh(25),tag21,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSFNW,1,MPI_DOUBLE_PRECISION,ngh(25),tag22,&
                  bufRBSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag22,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSBNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag23,&
                  bufRFSE,1,MPI_DOUBLE_PRECISION,ngh(23),tag23,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSFSE,1,MPI_DOUBLE_PRECISION,ngh(23),tag24,&
                  bufRBNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag24,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSBSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag25,&
                  bufRFNE,1,MPI_DOUBLE_PRECISION,ngh(22),tag25,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSFNE,1,MPI_DOUBLE_PRECISION,ngh(22),tag26,&
                  bufRBSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag26,COMM,stat,ierr) 
                   
ELSE

CALL MPI_SENDRECV(bufSE,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag1,&
                  bufRW,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag2,&
                  bufRE,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag2,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSN,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag3,&
                  bufRS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag3,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag4,&
                  bufRN,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag4,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSF,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag5,&
                  bufRB,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag5,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSB,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag6,&
                  bufRF,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag6,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSNE,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag7,&
                  bufRSW,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag7,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSSW,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag8,&
                  bufRNE,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag8,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSSE,NZP,MPI_DOUBLE_PRECISION,ngh(15),tag9,&
                  bufRNW,NZP,MPI_DOUBLE_PRECISION,ngh(17),tag9,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSNW,NZP,MPI_DOUBLE_PRECISION,ngh(17),tag10,&
                  bufRSE,NZP,MPI_DOUBLE_PRECISION,ngh(15),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag11,&
                  bufRFW,NYP,MPI_DOUBLE_PRECISION,ngh(21),tag11,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSFW,NYP,MPI_DOUBLE_PRECISION,ngh(21),tag12,&
                  bufRBE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag12,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSBW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag13,&
                  bufRFE,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag13,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSFE,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag14,&
                  bufRBW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag14,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSBN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag15,&
                  bufRFS,NXP,MPI_DOUBLE_PRECISION,ngh(20),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFS,NXP,MPI_DOUBLE_PRECISION,ngh(20),tag16,&
                  bufRBN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag16,COMM,stat,ierr)                  
                  
CALL MPI_SENDRECV(bufSBS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag17,&
                  bufRFN,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag17,COMM,stat,ierr)  
                  
CALL MPI_SENDRECV(bufSFN,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag18,&
                  bufRBS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag18,COMM,stat,ierr)
                  
CALL MPI_SENDRECV(bufSBNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag19,&
                  bufRFSW,1,MPI_DOUBLE_PRECISION,ngh(24),tag19,COMM,stat,ierr)         
                  
CALL MPI_SENDRECV(bufSFSW,1,MPI_DOUBLE_PRECISION,ngh(24),tag20,&
                  bufRBNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag20,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSBSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag21,&
                  bufRFNW,1,MPI_DOUBLE_PRECISION,ngh(25),tag21,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSFNW,1,MPI_DOUBLE_PRECISION,ngh(25),tag22,&
                  bufRBSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag22,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSBNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag23,&
                  bufRFSE,1,MPI_DOUBLE_PRECISION,ngh(23),tag23,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSFSE,1,MPI_DOUBLE_PRECISION,ngh(23),tag24,&
                  bufRBNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag24,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSBSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag25,&
                  bufRFNE,1,MPI_DOUBLE_PRECISION,ngh(22),tag25,COMM,stat,ierr) 
                  
CALL MPI_SENDRECV(bufSFNE,1,MPI_DOUBLE_PRECISION,ngh(22),tag26,&
                  bufRBSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag26,COMM,stat,ierr) 

END IF

! Field smoothing in each domain

DO ix=2,NXP-1
  DO iy=2,NYP-1
    DO iz=2,NZP-1

    Field_temp(ix,iy,iz)=alpha*(&
    1.0*Field(ix-1,iy+1,iz-1)+2.0*Field(ix,iy+1,iz-1)+1.0*Field(ix+1,iy+1,iz-1)+&
    2.0*Field(ix-1,iy,iz-1)+4.0*Field(ix,iy,iz-1)+2.0*Field(ix+1,iy,iz-1)+&
    1.0*Field(ix-1,iy-1,iz-1)+2.0*Field(ix,iy-1,iz-1)+1.0*Field(ix+1,iy-1,iz-1)+&
    2.0*Field(ix-1,iy+1,iz)+4.0*Field(ix,iy+1,iz)+2.0*Field(ix+1,iy+1,iz)+&
    4.0*Field(ix-1,iy,iz)+(1.0/alpha-56.0)*Field(ix,iy,iz)+4.0*Field(ix+1,iy,iz)+&
    2.0*Field(ix-1,iy-1,iz)+4.0*Field(ix,iy-1,iz)+2.0*Field(ix+1,iy-1,iz)+&
    1.0*Field(ix-1,iy+1,iz+1)+2.0*Field(ix,iy+1,iz+1)+1.0*Field(ix+1,iy+1,iz+1)+&
    2.0*Field(ix-1,iy,iz+1)+4.0*Field(ix,iy,iz+1)+2.0*Field(ix+1,iy,iz+1)+&
    1.0*Field(ix-1,iy-1,iz+1)+2.0*Field(ix,iy-1,iz+1)+1.0*Field(ix+1,iy-1,iz+1))
    
    ENDDO
  ENDDO
ENDDO

! Field smoothing at boundaries
! Faces
! ix=NXP

DO iy=2,NYP-1
  DO iz=2,NZP-1

  Field_temp(NXP,iy,iz)=alpha*(&
  1.0*Field(NXP-1,iy+1,iz-1)+2.0*Field(NXP,iy+1,iz-1)+1.0*bufRE(iy+1,iz-1)+&
  2.0*Field(NXP-1,iy,iz-1)+4.0*Field(NXP,iy,iz-1)+2.0*bufRE(iy,iz-1)+&
  1.0*Field(NXP-1,iy-1,iz-1)+2.0*Field(NXP,iy-1,iz-1)+1.0*bufRE(iy-1,iz-1)+&
  2.0*Field(NXP-1,iy+1,iz)+4.0*Field(NXP,iy+1,iz)+2.0*bufRE(iy+1,iz)+&
  4.0*Field(NXP-1,iy,iz)+(1.0/alpha-56.0)*Field(NXP,iy,iz)+4.0*bufRE(iy,iz)+&
  2.0*Field(NXP-1,iy-1,iz)+4.0*Field(NXP,iy-1,iz)+2.0*bufRE(iy-1,iz)+&
  1.0*Field(NXP-1,iy+1,iz+1)+2.0*Field(NXP,iy+1,iz+1)+1.0*bufRE(iy+1,iz+1)+&
  2.0*Field(NXP-1,iy,iz+1)+4.0*Field(NXP,iy,iz+1)+2.0*bufRE(iy,iz+1)+&
  1.0*Field(NXP-1,iy-1,iz+1)+2.0*Field(NXP,iy-1,iz+1)+1.0*bufRE(iy-1,iz+1))
  
  ENDDO
ENDDO

! ix=1

DO iy=2,NYP-1
  DO iz=2,NZP-1

  Field_temp(1,iy,iz)=alpha*(&
  1.0*bufRW(iy+1,iz-1)+2.0*Field(1,iy+1,iz-1)+1.0*Field(2,iy+1,iz-1)+&
  2.0*bufRW(iy,iz-1)+4.0*Field(1,iy,iz-1)+2.0*Field(2,iy,iz-1)+&
  1.0*bufRW(iy-1,iz-1)+2.0*Field(1,iy-1,iz-1)+1.0*Field(2,iy-1,iz-1)+&
  2.0*bufRW(iy+1,iz)+4.0*Field(1,iy+1,iz)+2.0*Field(2,iy+1,iz)+&
  4.0*bufRW(iy,iz)+(1.0/alpha-56.0)*Field(1,iy,iz)+4.0*Field(2,iy,iz)+&
  2.0*bufRW(iy-1,iz)+4.0*Field(1,iy-1,iz)+2.0*Field(2,iy-1,iz)+&
  1.0*bufRW(iy+1,iz+1)+2.0*Field(1,iy+1,iz+1)+1.0*Field(2,iy+1,iz+1)+&
  2.0*bufRW(iy,iz+1)+4.0*Field(1,iy,iz+1)+2.0*Field(2,iy,iz+1)+&
  1.0*bufRW(iy-1,iz+1)+2.0*Field(1,iy-1,iz+1)+1.0*Field(2,iy-1,iz+1))
  
  ENDDO
ENDDO

! iy=1

DO ix=2,NXP-1
  DO iz=2,NZP-1

  Field_temp(ix,1,iz)=alpha*(&
  1.0*Field(ix-1,2,iz-1)+2.0*Field(ix,2,iz-1)+1.0*Field(ix+1,2,iz-1)+&
  2.0*Field(ix-1,1,iz-1)+4.0*Field(ix,1,iz-1)+2.0*Field(ix+1,1,iz-1)+&
  1.0*bufRS(ix-1,iz-1)+2.0*bufRS(ix,iz-1)+1.0*bufRS(ix+1,iz-1)+&
  2.0*Field(ix-1,2,iz)+4.0*Field(ix,2,iz)+2.0*Field(ix+1,2,iz)+&
  4.0*Field(ix-1,1,iz)+(1.0/alpha-56.0)*Field(ix,1,iz)+4.0*Field(ix+1,1,iz)+&
  2.0*bufRS(ix-1,iz)+4.0*bufRS(ix,iz)+2.0*bufRS(ix+1,iz)+&
  1.0*Field(ix-1,2,iz+1)+2.0*Field(ix,2,iz+1)+1.0*Field(ix+1,2,iz+1)+&
  2.0*Field(ix-1,1,iz+1)+4.0*Field(ix,1,iz+1)+2.0*Field(ix+1,1,iz+1)+&
  1.0*bufRS(ix-1,iz+1)+2.0*bufRS(ix,iz+1)+1.0*bufRS(ix+1,iz+1))

  ENDDO
ENDDO

! iy=NYP

DO ix=2,NXP-1
  DO iz=2,NZP-1

  Field_temp(ix,NYP,iz)=alpha*(&
  1.0*bufRN(ix-1,iz-1)+2.0*bufRN(ix,iz-1)+1.0*bufRN(ix+1,iz-1)+&
  2.0*Field(ix-1,NYP,iz-1)+4.0*Field(ix,NYP,iz-1)+2.0*Field(ix+1,NYP,iz-1)+&
  1.0*Field(ix-1,NYP-1,iz-1)+2.0*Field(ix,NYP-1,iz-1)+1.0*Field(ix+1,NYP-1,iz-1)+&
  2.0*bufRN(ix-1,iz)+4.0*bufRN(ix,iz)+2.0*bufRN(ix+1,iz)+&
  4.0*Field(ix-1,NYP,iz)+(1.0/alpha-56.0)*Field(ix,NYP,iz)+4.0*Field(ix+1,NYP,iz)+&
  2.0*Field(ix-1,NYP-1,iz)+4.0*Field(ix,NYP-1,iz)+2.0*Field(ix+1,NYP-1,iz)+&
  1.0*bufRN(ix-1,iz+1)+2.0*bufRN(ix,iz+1)+1.0*bufRN(ix+1,iz+1)+&
  2.0*Field(ix-1,NYP,iz+1)+4.0*Field(ix,NYP,iz+1)+2.0*Field(ix+1,NYP,iz+1)+&
  1.0*Field(ix-1,NYP-1,iz+1)+2.0*Field(ix,NYP-1,iz+1)+1.0*Field(ix+1,NYP-1,iz+1))
    
  ENDDO
ENDDO

! iz=1

DO ix=2,NXP-1
  DO iy=2,NYP-1

  Field_temp(ix,iy,1)=alpha*(&
  1.0*bufRB(ix-1,iy+1)+2.0*bufRB(ix,iy+1)+1.0*bufRB(ix+1,iy+1)+&
  2.0*bufRB(ix-1,iy)+4.0*bufRB(ix,iy)+2.0*bufRB(ix+1,iy)+&
  1.0*bufRB(ix-1,iy-1)+2.0*bufRB(ix,iy-1)+1.0*bufRB(ix+1,iy-1)+&
  2.0*Field(ix-1,iy+1,1)+4.0*Field(ix,iy+1,1)+2.0*Field(ix+1,iy+1,1)+&
  4.0*Field(ix-1,iy,1)+(1.0/alpha-56.0)*Field(ix,iy,1)+4.0*Field(ix+1,iy,1)+&
  2.0*Field(ix-1,iy-1,1)+4.0*Field(ix,iy-1,1)+2.0*Field(ix+1,iy-1,1)+&
  1.0*Field(ix-1,iy+1,2)+2.0*Field(ix,iy+1,2)+1.0*Field(ix+1,iy+1,2)+&
  2.0*Field(ix-1,iy,2)+4.0*Field(ix,iy,2)+2.0*Field(ix+1,iy,2)+&
  1.0*Field(ix-1,iy-1,2)+2.0*Field(ix,iy-1,2)+1.0*Field(ix+1,iy-1,2))

  ENDDO
ENDDO

! iz=NZP

DO ix=2,NXP-1
  DO iy=2,NYP-1

  Field_temp(ix,iy,NZP)=alpha*(&
  1.0*Field(ix-1,iy+1,NZP-1)+2.0*Field(ix,iy+1,NZP-1)+1.0*Field(ix+1,iy+1,NZP-1)+&
  2.0*Field(ix-1,iy,NZP-1)+4.0*Field(ix,iy,NZP-1)+2.0*Field(ix+1,iy,NZP-1)+&
  1.0*Field(ix-1,iy-1,NZP-1)+2.0*Field(ix,iy-1,NZP-1)+1.0*Field(ix+1,iy-1,NZP-1)+&
  2.0*Field(ix-1,iy+1,NZP)+4.0*Field(ix,iy+1,NZP)+2.0*Field(ix+1,iy+1,NZP)+&
  4.0*Field(ix-1,iy,NZP)+(1.0/alpha-56.0)*Field(ix,iy,NZP)+4.0*Field(ix+1,iy,NZP)+&
  2.0*Field(ix-1,iy-1,NZP)+4.0*Field(ix,iy-1,NZP)+2.0*Field(ix+1,iy-1,NZP)+&
  1.0*bufRF(ix-1,iy+1)+2.0*bufRF(ix,iy+1)+1.0*bufRF(ix+1,iy+1)+&
  2.0*bufRF(ix-1,iy)+4.0*bufRF(ix,iy)+2.0*bufRF(ix+1,iy)+&
  1.0*bufRF(ix-1,iy-1)+2.0*bufRF(ix,iy-1)+1.0*bufRF(ix+1,iy-1))
    
  ENDDO
ENDDO

! Aretes

! ix=NXP,iz=1

DO iy=2,NYP-1

  Field_temp(NXP,iy,1)=alpha*(&
  1.0*bufRB(NXP-1,iy+1)+2.0*bufRB(NXP,iy+1)+1.0*bufRBE(iy+1)+&
  2.0*bufRB(NXP-1,iy)+4.0*bufRB(NXP,iy)+2.0*bufRBE(iy)+&
  1.0*bufRB(NXP-1,iy-1)+2.0*bufRB(NXP,iy-1)+1.0*bufRBE(iy-1)+&
  2.0*Field(NXP-1,iy+1,1)+4.0*Field(NXP,iy+1,1)+2.0*bufRE(iy+1,1)+&
  4.0*Field(NXP-1,iy,1)+(1.0/alpha-56.0)*Field(NXP,iy,1)+4.0*bufRE(iy,1)+&
  2.0*Field(NXP-1,iy-1,1)+4.0*Field(NXP,iy-1,1)+2.0*bufRE(iy-1,1)+&
  1.0*Field(NXP-1,iy+1,2)+2.0*Field(NXP,iy+1,2)+1.0*bufRE(iy+1,2)+&
  2.0*Field(NXP-1,iy,2)+4.0*Field(NXP,iy,2)+2.0*bufRE(iy,2)+&
  1.0*Field(NXP-1,iy-1,2)+2.0*Field(NXP,iy-1,2)+1.0*bufRE(iy-1,2))
    
ENDDO

! ix=NXP,iz=NZP

DO iy=2,NYP-1

  Field_temp(NXP,iy,NZP)=alpha*(&
  1.0*Field(NXP-1,iy+1,NZP-1)+2.0*Field(NXP,iy+1,NZP-1)+1.0*bufRE(iy+1,NZP-1)+&
  2.0*Field(NXP-1,iy,NZP-1)+4.0*Field(NXP,iy,NZP-1)+2.0*bufRE(iy,NZP-1)+&
  1.0*Field(NXP-1,iy-1,NZP-1)+2.0*Field(NXP,iy-1,NZP-1)+1.0*bufRE(iy-1,NZP-1)+&
  2.0*Field(NXP-1,iy+1,NZP)+4.0*Field(NXP,iy+1,NZP)+2.0*bufRE(iy+1,NZP)+&
  4.0*Field(NXP-1,iy,NZP)+(1.0/alpha-56.0)*Field(NXP,iy,NZP)+4.0*bufRE(iy,NZP)+&
  2.0*Field(NXP-1,iy-1,NZP)+4.0*Field(NXP,iy-1,NZP)+2.0*bufRE(iy-1,NZP)+&
  1.0*bufRF(NXP-1,iy+1)+2.0*bufRF(NXP,iy+1)+1.0*bufRFE(iy+1)+&
  2.0*bufRF(NXP-1,iy)+4.0*bufRF(NXP,iy)+2.0*bufRFE(iy)+&
  1.0*bufRF(NXP-1,iy-1)+2.0*bufRF(NXP,iy-1)+1.0*bufRFE(iy-1))
    
ENDDO

! ix=1,iz=1

DO iy=2,NYP-1

  Field_temp(1,iy,1)=alpha*(&
  1.0*bufRBW(iy+1)+2.0*bufRB(1,iy+1)+1.0*bufRB(2,iy+1)+&
  2.0*bufRBW(iy)+4.0*bufRB(1,iy)+2.0*bufRB(2,iy)+&
  1.0*bufRBW(iy-1)+2.0*bufRB(1,iy-1)+1.0*bufRB(2,iy-1)+&
  2.0*bufRW(iy+1,1)+4.0*Field(1,iy+1,1)+2.0*Field(2,iy+1,1)+&
  4.0*bufRW(iy,1)+(1.0/alpha-56.0)*Field(1,iy,1)+4.0*Field(2,iy,1)+&
  2.0*bufRW(iy-1,1)+4.0*Field(1,iy-1,1)+2.0*Field(2,iy-1,1)+&
  1.0*bufRW(iy+1,2)+2.0*Field(1,iy+1,2)+1.0*Field(2,iy+1,2)+&
  2.0*bufRW(iy,2)+4.0*Field(1,iy,2)+2.0*Field(2,iy,2)+&
  1.0*bufRW(iy-1,2)+2.0*Field(1,iy-1,2)+1.0*Field(2,iy-1,2))
    
ENDDO

! ix=1,iz=NZP

DO iy=2,NYP-1

  Field_temp(1,iy,NZP)=alpha*(&
  1.0*bufRW(iy+1,NZP-1)+2.0*Field(1,iy+1,NZP-1)+1.0*Field(2,iy+1,NZP-1)+&
  2.0*bufRW(iy,NZP-1)+4.0*Field(1,iy,NZP-1)+2.0*Field(2,iy,NZP-1)+&
  1.0*bufRW(iy-1,NZP-1)+2.0*Field(1,iy-1,NZP-1)+1.0*Field(2,iy-1,NZP-1)+&
  2.0*bufRW(iy+1,NZP)+4.0*Field(1,iy+1,NZP)+2.0*Field(2,iy+1,NZP)+&
  4.0*bufRW(iy,NZP)+(1.0/alpha-56.0)*Field(1,iy,NZP)+4.0*Field(2,iy,NZP)+&
  2.0*bufRW(iy-1,NZP)+4.0*Field(1,iy-1,NZP)+2.0*Field(2,iy-1,NZP)+&
  1.0*bufRFW(iy+1)+2.0*bufRF(1,iy+1)+1.0*bufRF(2,iy+1)+&
  2.0*bufRFW(iy)+4.0*bufRF(1,iy)+2.0*bufRF(2,iy)+&
  1.0*bufRFW(iy-1)+2.0*bufRF(1,iy-1)+1.0*bufRF(2,iy-1))

ENDDO

! iy=1,iz=1

DO ix=2,NXP-1

  Field_temp(ix,1,1)=alpha*(&
  1.0*bufRB(ix-1,2)+2.0*bufRB(ix,2)+1.0*bufRB(ix+1,2)+&
  2.0*bufRB(ix-1,1)+4.0*bufRB(ix,1)+2.0*bufRB(ix+1,1)+&
  1.0*bufRBS(ix-1)+2.0*bufRBS(ix)+1.0*bufRBS(ix+1)+&
  2.0*Field(ix-1,2,1)+4.0*Field(ix,2,1)+2.0*Field(ix+1,2,1)+&
  4.0*Field(ix-1,1,1)+(1.0/alpha-56.0)*Field(ix,1,1)+4.0*Field(ix+1,1,1)+&
  2.0*bufRS(ix-1,1)+4.0*bufRS(ix,1)+2.0*bufRS(ix+1,1)+&
  1.0*Field(ix-1,2,2)+2.0*Field(ix,2,2)+1.0*Field(ix+1,2,2)+&
  2.0*Field(ix-1,1,2)+4.0*Field(ix,1,2)+2.0*Field(ix+1,1,2)+&
  1.0*bufRS(ix-1,2)+2.0*bufRS(ix,2)+1.0*bufRS(ix+1,2))
    
ENDDO

! iy=NYP,iz=1

DO ix=2,NXP-1

  Field_temp(ix,NYP,1)=alpha*(&
  1.0*bufRBN(ix-1)+2.0*bufRBN(ix)+1.0*bufRBN(ix+1)+&
  2.0*bufRB(ix-1,NYP)+4.0*bufRB(ix,NYP)+2.0*bufRB(ix+1,NYP)+&
  1.0*bufRB(ix-1,NYP-1)+2.0*bufRB(ix,NYP-1)+1.0*bufRB(ix+1,NYP-1)+&
  2.0*bufRN(ix-1,1)+4.0*bufRN(ix,1)+2.0*bufRN(ix+1,1)+&
  4.0*Field(ix-1,NYP,1)+(1.0/alpha-56.0)*Field(ix,NYP,1)+4.0*Field(ix+1,NYP,1)+&
  2.0*Field(ix-1,NYP-1,1)+4.0*Field(ix,NYP-1,1)+2.0*Field(ix+1,NYP-1,1)+&
  1.0*bufRN(ix-1,2)+2.0*bufRN(ix,2)+1.0*bufRN(ix+1,2)+&
  2.0*Field(ix-1,NYP,2)+4.0*Field(ix,NYP,2)+2.0*Field(ix+1,NYP,2)+&
  1.0*Field(ix-1,NYP-1,2)+2.0*Field(ix,NYP-1,2)+1.0*Field(ix+1,NYP-1,2))

ENDDO

! iy=1,iz=NZP

DO ix=2,NXP-1

  Field_temp(ix,1,NZP)=alpha*(&
  1.0*Field(ix-1,2,NZP-1)+2.0*Field(ix,2,NZP-1)+1.0*Field(ix+1,2,NZP-1)+&
  2.0*Field(ix-1,1,NZP-1)+4.0*Field(ix,1,NZP-1)+2.0*Field(ix+1,1,NZP-1)+&
  1.0*bufRS(ix-1,NZP-1)+2.0*bufRS(ix,NZP-1)+1.0*bufRS(ix+1,NZP-1)+&
  2.0*Field(ix-1,2,NZP)+4.0*Field(ix,2,NZP)+2.0*Field(ix+1,2,NZP)+&
  4.0*Field(ix-1,1,NZP)+(1.0/alpha-56.0)*Field(ix,1,NZP)+4.0*Field(ix+1,1,NZP)+&
  2.0*bufRS(ix-1,NZP)+4.0*bufRS(ix,NZP)+2.0*bufRS(ix+1,NZP)+&
  1.0*bufRF(ix-1,2)+2.0*bufRF(ix,2)+1.0*bufRF(ix+1,2)+&
  2.0*bufRF(ix-1,1)+4.0*bufRF(ix,1)+2.0*bufRF(ix+1,1)+&
  1.0*bufRFS(ix-1)+2.0*bufRFS(ix)+1.0*bufRFS(ix+1))

ENDDO

! iy=NYP,iz=NZP

DO ix=2,NXP-1

  Field_temp(ix,NYP,NZP)=alpha*(&
  1.0*bufRN(ix-1,NZP-1)+2.0*bufRN(ix,NZP-1)+1.0*bufRN(ix+1,NZP-1)+&
  2.0*Field(ix-1,NYP,NZP-1)+4.0*Field(ix,NYP,NZP-1)+2.0*Field(ix+1,NYP,NZP-1)+&
  1.0*Field(ix-1,NYP-1,NZP-1)+2.0*Field(ix,NYP-1,NZP-1)+1.0*Field(ix+1,NYP-1,NZP-1)+&
  2.0*bufRN(ix-1,NZP)+4.0*bufRN(ix,NZP)+2.0*bufRN(ix+1,NZP)+&
  4.0*Field(ix-1,NYP,NZP)+(1.0/alpha-56.0)*Field(ix,NYP,NZP)+4.0*Field(ix+1,NYP,NZP)+&
  2.0*Field(ix-1,NYP-1,NZP)+4.0*Field(ix,NYP-1,NZP)+2.0*Field(ix+1,NYP-1,NZP)+&
  1.0*bufRFN(ix-1)+2.0*bufRFN(ix)+1.0*bufRFN(ix+1)+&
  2.0*bufRF(ix-1,NYP)+4.0*bufRF(ix,NYP)+2.0*bufRF(ix+1,NYP)+&
  1.0*bufRF(ix-1,NYP-1)+2.0*bufRF(ix,NYP-1)+1.0*bufRF(ix+1,NYP-1))

ENDDO

! ix=1,iy=1

DO iz=2,NZP-1

  Field_temp(1,1,iz)=alpha*(&
  1.0*bufRW(2,iz-1)+2.0*Field(1,2,iz-1)+1.0*Field(2,2,iz-1)+&
  2.0*bufRW(1,iz-1)+4.0*Field(1,1,iz-1)+2.0*Field(2,1,iz-1)+&
  1.0*bufRSW(iz-1)+2.0*bufRS(1,iz-1)+1.0*bufRS(2,iz-1)+&
  2.0*bufRW(2,iz)+4.0*Field(1,2,iz)+2.0*Field(2,2,iz)+&
  4.0*bufRW(1,iz)+(1.0/alpha-56.0)*Field(1,1,iz)+4.0*Field(2,1,iz)+&
  2.0*bufRSW(iz)+4.0*bufRS(1,iz)+2.0*bufRS(2,iz)+&
  1.0*bufRW(2,iz+1)+2.0*Field(1,2,iz+1)+1.0*Field(2,2,iz+1)+&
  2.0*bufRW(1,iz+1)+4.0*Field(1,1,iz+1)+2.0*Field(2,1,iz+1)+&
  1.0*bufRSW(iz+1)+2.0*bufRS(1,iz+1)+1.0*bufRS(2,iz+1))
    
ENDDO

! ix=1,iy=NYP

DO iz=2,NZP-1

  Field_temp(1,NYP,iz)=alpha*(&
  1.0*bufRNW(iz-1)+2.0*bufRN(1,iz-1)+1.0*bufRN(2,iz-1)+&
  2.0*bufRW(NYP,iz-1)+4.0*Field(1,NYP,iz-1)+2.0*Field(2,NYP,iz-1)+&
  1.0*bufRW(NYP-1,iz-1)+2.0*Field(1,NYP-1,iz-1)+1.0*Field(2,NYP-1,iz-1)+&
  2.0*bufRNW(iz)+4.0*bufRN(1,iz)+2.0*bufRN(2,iz)+&
  4.0*bufRW(NYP,iz)+(1.0/alpha-56.0)*Field(1,NYP,iz)+4.0*Field(2,NYP,iz)+&
  2.0*bufRW(NYP-1,iz)+4.0*Field(1,NYP-1,iz)+2.0*Field(2,NYP-1,iz)+&
  1.0*bufRNW(iz+1)+2.0*bufRN(1,iz+1)+1.0*bufRN(2,iz+1)+&
  2.0*bufRW(NYP,iz+1)+4.0*Field(1,NYP,iz+1)+2.0*Field(2,NYP,iz+1)+&
  1.0*bufRW(NYP-1,iz+1)+2.0*Field(1,NYP-1,iz+1)+1.0*Field(2,NYP-1,iz+1))

ENDDO

! ix=NXP,iy=1

DO iz=2,NZP-1

    Field_temp(NXP,1,iz)=alpha*(&
    1.0*Field(NXP-1,2,iz-1)+2.0*Field(NXP,2,iz-1)+1.0*bufRE(2,iz-1)+&
    2.0*Field(NXP-1,1,iz-1)+4.0*Field(NXP,1,iz-1)+2.0*bufRE(1,iz-1)+&
    1.0*bufRS(NXP-1,iz-1)+2.0*bufRS(NXP,iz-1)+1.0*bufRSE(iz-1)+&
    2.0*Field(NXP-1,2,iz)+4.0*Field(NXP,2,iz)+2.0*bufRE(2,iz)+&
    4.0*Field(NXP-1,1,iz)+(1.0/alpha-56.0)*Field(NXP,1,iz)+4.0*bufRE(1,iz)+&
    2.0*bufRS(NXP-1,iz)+4.0*bufRS(NXP,iz)+2.0*bufRSE(iz)+&
    1.0*Field(NXP-1,2,iz+1)+2.0*Field(NXP,2,iz+1)+1.0*bufRE(2,iz+1)+&
    2.0*Field(NXP-1,1,iz+1)+4.0*Field(NXP,1,iz+1)+2.0*bufRE(1,iz+1)+&
    1.0*bufRS(NXP-1,iz+1)+2.0*bufRS(NXP,iz+1)+1.0*bufRSE(iz+1))

ENDDO

! ix=NXP,iy=NYP

DO iz=2,NZP-1

  Field_temp(NXP,NYP,iz)=alpha*(&
  1.0*bufRN(NXP-1,iz-1)+2.0*bufRN(NXP,iz-1)+1.0*bufRNE(iz-1)+&
  2.0*Field(NXP-1,NYP,iz-1)+4.0*Field(NXP,NYP,iz-1)+2.0*bufRE(NYP,iz-1)+&
  1.0*Field(NXP-1,NYP-1,iz-1)+2.0*Field(NXP,NYP-1,iz-1)+1.0*bufRE(NYP-1,iz-1)+&
  2.0*bufRN(NXP-1,iz)+4.0*bufRN(NXP,iz)+2.0*bufRNE(iz)+&
  4.0*Field(NXP-1,NYP,iz)+(1.0/alpha-56.0)*Field(NXP,NYP,iz)+4.0*bufRE(NYP,iz)+&
  2.0*Field(NXP-1,NYP-1,iz)+4.0*Field(NXP,NYP-1,iz)+2.0*bufRE(NYP-1,iz)+&
  1.0*bufRN(NXP-1,iz+1)+2.0*bufRN(NXP,iz+1)+1.0*bufRNE(iz+1)+&
  2.0*Field(NXP-1,NYP,iz+1)+4.0*Field(NXP,NYP,iz+1)+2.0*bufRE(NYP,iz+1)+&
  1.0*Field(NXP-1,NYP-1,iz+1)+2.0*Field(NXP,NYP-1,iz+1)+1.0*bufRE(NYP-1,iz+1))

ENDDO

! Corners

! 1,1,1

Field_temp(1,1,1)=alpha*(&
1.0*bufRBW(2)+2.0*bufRB(1,2)+1.0*bufRB(2,2)+&
2.0*bufRBW(1)+4.0*bufRB(1,1)+2.0*bufRB(2,1)+&
1.0*bufRBSW+2.0*bufRBS(1)+1.0*bufRBS(2)+&
2.0*bufRW(2,1)+4.0*Field(1,2,1)+2.0*Field(2,2,1)+&
4.0*bufRW(1,1)+(1.0/alpha-56.0)*Field(1,1,1)+4.0*Field(2,1,1)+&
2.0*bufRSW(1)+4.0*bufRS(1,1)+2.0*bufRS(2,1)+&
1.0*bufRW(2,2)+2.0*Field(1,2,2)+1.0*Field(2,2,2)+&
2.0*bufRW(1,2)+4.0*Field(1,1,2)+2.0*Field(2,1,2)+&
1.0*bufRSW(2)+2.0*bufRS(1,2)+1.0*bufRS(2,2))

! NXP,1,1

Field_temp(NXP,1,1)=alpha*(&
1.0*bufRB(NXP-1,2)+2.0*bufRB(NXP,2)+1.0*bufRBE(2)+&
2.0*bufRB(NXP-1,1)+4.0*bufRB(NXP,1)+2.0*bufRBE(1)+&
1.0*bufRBS(NXP-1)+2.0*bufRBS(NXP)+1.0*bufRBSE+&
2.0*Field(NXP-1,2,1)+4.0*Field(NXP,2,1)+2.0*bufRE(2,1)+&
4.0*Field(NXP-1,1,1)+(1.0/alpha-56.0)*Field(NXP,1,1)+4.0*bufRE(1,1)+&
2.0*bufRS(NXP-1,1)+4.0*bufRS(NXP,1)+2.0*bufRSE(1)+&
1.0*Field(NXP-1,2,2)+2.0*Field(NXP,2,2)+1.0*bufRE(2,2)+&
2.0*Field(NXP-1,1,2)+4.0*Field(NXP,1,2)+2.0*bufRE(1,2)+&
1.0*bufRS(NXP-1,2)+2.0*bufRS(NXP,2)+1.0*bufRSE(2))

! 1,1,NZP

Field_temp(1,1,NZP)=alpha*(&
1.0*bufRW(2,NZP-1)+2.0*Field(1,2,NZP-1)+1.0*Field(2,2,NZP-1)+&
2.0*bufRW(1,NZP-1)+4.0*Field(1,1,NZP-1)+2.0*Field(2,1,NZP-1)+&
1.0*bufRSW(NZP-1)+2.0*bufRS(1,NZP-1)+1.0*bufRS(2,NZP-1)+&
2.0*bufRW(2,NZP)+4.0*Field(1,2,NZP)+2.0*Field(2,2,NZP)+&
4.0*bufRW(1,NZP)+(1.0/alpha-56.0)*Field(1,1,NZP)+4.0*Field(2,1,NZP)+&
2.0*bufRSW(NZP)+4.0*bufRS(1,NZP)+2.0*bufRS(2,NZP)+&
1.0*bufRFW(2)+2.0*bufRF(1,2)+1.0*bufRF(2,2)+&
2.0*bufRFW(1)+4.0*bufRF(1,1)+2.0*bufRF(2,1)+&
1.0*bufRFSW+2.0*bufRFS(1)+1.0*bufRFS(2))

! 1,NYP,1

Field_temp(1,NYP,1)=alpha*(&
1.0*bufRBNW+2.0*bufRBN(1)+1.0*bufRBN(2)+&
2.0*bufRBW(NYP)+4.0*bufRB(1,NYP)+2.0*bufRB(2,NYP)+&
1.0*bufRBW(NYP-1)+2.0*bufRB(1,NYP-1)+1.0*bufRB(2,NYP-1)+&
2.0*bufRNW(1)+4.0*bufRN(1,1)+2.0*bufRN(2,1)+&
4.0*bufRW(NYP,1)+(1.0/alpha-56.0)*Field(1,NYP,1)+4.0*Field(2,NYP,1)+&
2.0*bufRW(NYP-1,1)+4.0*Field(1,NYP-1,1)+2.0*Field(2,NYP-1,1)+&
1.0*bufRNW(2)+2.0*bufRN(1,2)+1.0*bufRN(2,2)+&
2.0*bufRW(NYP,2)+4.0*Field(1,NYP,2)+2.0*Field(2,NYP,2)+&
1.0*bufRW(NYP-1,2)+2.0*Field(1,NYP-1,2)+1.0*Field(2,NYP-1,2))

! NXP,1,NZP

Field_temp(NXP,1,NZP)=alpha*(&
1.0*Field(NXP-1,2,NZP-1)+2.0*Field(NXP,2,NZP-1)+1.0*bufRE(2,NZP-1)+&
2.0*Field(NXP-1,1,NZP-1)+4.0*Field(NXP,1,NZP-1)+2.0*bufRE(1,NZP-1)+&
1.0*bufRS(NXP-1,NZP-1)+2.0*bufRS(NXP,NZP-1)+1.0*bufRSE(NZP-1)+&
2.0*Field(NXP-1,2,NZP)+4.0*Field(NXP,2,NZP)+2.0*bufRE(2,NZP)+&
4.0*Field(NXP-1,1,NZP)+(1.0/alpha-56.0)*Field(NXP,1,NZP)+4.0*bufRE(1,NZP)+&
2.0*bufRS(NXP-1,NZP)+4.0*bufRS(NXP,NZP)+2.0*bufRSE(NZP)+&
1.0*bufRF(NXP-1,2)+2.0*bufRF(NXP,2)+1.0*bufRFE(2)+&
2.0*bufRF(NXP-1,1)+4.0*bufRF(NXP,1)+2.0*bufRFE(1)+&
1.0*bufRFS(NXP-1)+2.0*bufRFS(NXP)+1.0*bufRFSE)
        
! NXP,NYP,1
    
Field_temp(NXP,NYP,1)=alpha*(&
1.0*bufRBN(NXP-1)+2.0*bufRBN(NXP)+1.0*bufRBNE+&
2.0*bufRB(NXP-1,NYP)+4.0*bufRB(NXP,NYP)+2.0*bufRBE(NYP)+&
1.0*bufRB(NXP-1,NYP-1)+2.0*bufRB(NXP,NYP-1)+1.0*bufRBE(NYP-1)+&
2.0*bufRN(NXP-1,1)+4.0*bufRN(NXP,1)+2.0*bufRNE(1)+&
4.0*Field(NXP-1,NYP,1)+(1.0/alpha-56.0)*Field(NXP,NYP,1)+4.0*bufRE(NYP,1)+&
2.0*Field(NXP-1,NYP-1,1)+4.0*Field(NXP,NYP-1,1)+2.0*bufRE(NYP-1,1)+&
1.0*bufRN(NXP-1,2)+2.0*bufRN(NXP,2)+1.0*bufRNE(2)+&
2.0*Field(NXP-1,NYP,2)+4.0*Field(NXP,NYP,2)+2.0*bufRE(NYP,2)+&
1.0*Field(NXP-1,NYP-1,2)+2.0*Field(NXP,NYP-1,2)+1.0*bufRE(NYP-1,2))

! 1,NYP,NZP
    
Field_temp(1,NYP,NZP)=alpha*(&
1.0*bufRNW(NZP-1)+2.0*bufRN(1,NZP-1)+1.0*bufRN(2,NZP-1)+&
2.0*bufRW(NYP,NZP-1)+4.0*Field(1,NYP,NZP-1)+2.0*Field(2,NYP,NZP-1)+&
1.0*bufRW(NYP-1,NZP-1)+2.0*Field(1,NYP-1,NZP-1)+1.0*Field(2,NYP-1,NZP-1)+&
2.0*bufRNW(NZP)+4.0*bufRN(1,NZP)+2.0*bufRN(2,NZP)+&
4.0*bufRW(NYP,NZP)+(1.0/alpha-56.0)*Field(1,NYP,NZP)+4.0*Field(2,NYP,NZP)+&
2.0*bufRW(NYP-1,NZP)+4.0*Field(1,NYP-1,NZP)+2.0*Field(2,NYP-1,NZP)+&
1.0*bufRFNW+2.0*bufRFN(1)+1.0*bufRFN(2)+&
2.0*bufRFW(NYP)+4.0*bufRF(1,NYP)+2.0*bufRF(2,NYP)+&
1.0*bufRFW(NYP-1)+2.0*bufRF(1,NYP-1)+1.0*bufRF(2,NYP-1))
  
! NXP,NYP,NZP
    
Field_temp(NXP,NYP,NZP)=alpha*(&
1.0*bufRN(NXP-1,NZP-1)+2.0*bufRN(NXP,NZP-1)+1.0*bufRNE(NZP-1)+&
2.0*Field(NXP-1,NYP,NZP-1)+4.0*Field(NXP,NYP,NZP-1)+2.0*bufRE(NYP,NZP-1)+&
1.0*Field(NXP-1,NYP-1,NZP-1)+2.0*Field(NXP,NYP-1,NZP-1)+1.0*bufRE(NYP-1,NZP-1)+&
2.0*bufRN(NXP-1,NZP)+4.0*bufRN(NXP,NZP)+2.0*bufRNE(NZP)+&
4.0*Field(NXP-1,NYP,NZP)+(1.0/alpha-56.0)*Field(NXP,NYP,NZP)+4.0*bufRE(NYP,NZP)+&
2.0*Field(NXP-1,NYP-1,NZP)+4.0*Field(NXP,NYP-1,NZP)+2.0*bufRE(NYP-1,NZP)+&
1.0*bufRFN(NXP-1)+2.0*bufRFN(NXP)+1.0*bufRFNE+&
2.0*bufRF(NXP-1,NYP)+4.0*bufRF(NXP,NYP)+2.0*bufRFE(NYP)+&
1.0*bufRF(NXP-1,NYP-1)+2.0*bufRF(NXP,NYP-1)+1.0*bufRFE(NYP-1))

! Updating Field
Field(:,:,:)=Field_temp(:,:,:)

END SUBROUTINE FILTER_FIELD

!***********************************************************************

END MODULE MOD_FIELDS
