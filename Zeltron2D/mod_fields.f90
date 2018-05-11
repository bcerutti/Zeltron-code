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
PUBLIC :: FILTER_FIELD ! Isotropic 2D filter using the nearest neighboring cells

 CONTAINS

!***********************************************************************
! Subroutine PUSH_EFIELD
! This subroutine computes the E field vector at time t+dt

! INPUT: 
! - Bx: x-component of B on the Yee grid at time t
! - By: y-component of B on the Yee grid at time t
! - Bz: z-component of B on the Yee grid at time t
! - Ex: x-component of E on the Yee grid at time t
! - Ey: y-component of E on the Yee grid at time t
! - Ez: z-component of E on the Yee grid at time t
! - Jx: x-component of current J on the Yee grid at time t
! - Jy: y-component of current J on the Yee grid at time t
! - Jz: z-component of current J on the Yee grid at time t
! - xgp,ygp: Local grid
!
! OUTPUT: Ex,Ey,Ez at time t+dt
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE PUSH_EFIELD(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,xgp,ygp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bx,By,Bz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Jx,Jy,Jz
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION                         :: xminp,xmaxp,yminp,ymaxp
INTEGER, DIMENSION(8)                    :: ngh
INTEGER, PARAMETER                       :: tag1=51,tag2=52
INTEGER                                  :: id,COMM,ierr
DOUBLE PRECISION, ALLOCATABLE            :: bufS1(:),BufR1(:)
DOUBLE PRECISION, ALLOCATABLE            :: bufS2(:),BufR2(:)

! Loop indexes
INTEGER :: ix,iy

!***********************************************************************

xminp=xgp(1)
xmaxp=xgp(NXP)
yminp=ygp(1)
ymaxp=ygp(NYP)

!***********************************************************************
! Solve Ex at t=t+dt

ALLOCATE(bufS1(1:NXP),bufR1(1:NXP))

bufS1=Bz(:,NYP-1)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
ENDIF

DO iy=2,NYP
Ex(:,iy)=Ex(:,iy)+(c*dt/dy)*(Bz(:,iy)-Bz(:,iy-1))-4.0*pi*dt*Jx(:,iy)
ENDDO

Ex(:,1)=Ex(:,1)+(c*dt/dy)*(Bz(:,1)-bufR1(:))-4.0*pi*dt*Jx(:,1)

DEALLOCATE(bufS1,bufR1)

!***********************************************************************
! Check boundary conditions along X
  
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Inside the conductor
   Ex(NXP,:)=0.0
   END IF
   
END IF

!***********************************************************************
! Check boundary conditions along Y
 
IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ex(:,1)=0.0
   END IF
   
END IF
   
IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ex(:,NYP)=0.0
   END IF
   
END IF
   
!***********************************************************************
! Solve Ey at t=t+dt

ALLOCATE(bufS2(1:NYP),bufR2(1:NYP))

bufS2=Bz(NXP-1,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
Ey(ix,:)=Ey(ix,:)-(c*dt/dx)*(Bz(ix,:)-Bz(ix-1,:))-4.0*pi*dt*Jy(ix,:)
ENDDO

Ey(1,:)=Ey(1,:)-(c*dt/dx)*(Bz(1,:)-bufR2(:))-4.0*pi*dt*Jy(1,:)

DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X

IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ey(1,:)=0.0
   END IF
   
END IF
   
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ey(NXP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Inside the conductor
   Ey(:,NYP)=0.0
   END IF
   
END IF

!***********************************************************************
! Solve Ez at t=t+dt

ALLOCATE(bufS1(1:NXP),bufR1(1:NXP))
ALLOCATE(bufS2(1:NYP),bufR2(1:NYP))

bufS1=Bx(:,NYP-1)
bufS2=By(NXP-1,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
DO iy=2,NYP
Ez(ix,iy)=Ez(ix,iy)+(c*dt/dx)*(By(ix,iy)-By(ix-1,iy))-&
          (c*dt/dy)*(Bx(ix,iy)-Bx(ix,iy-1))-4.0*pi*dt*Jz(ix,iy)
ENDDO
ENDDO

DO iy=2,NYP
Ez(1,iy)=Ez(1,iy)+(c*dt/dx)*(By(1,iy)-bufR2(iy))-&
         (c*dt/dy)*(Bx(1,iy)-Bx(1,iy-1))-4.0*pi*dt*Jz(1,iy)
ENDDO

DO ix=2,NXP
Ez(ix,1)=Ez(ix,1)+(c*dt/dx)*(By(ix,1)-By(ix-1,1))-&
          (c*dt/dy)*(Bx(ix,1)-bufR1(ix))-4.0*pi*dt*Jz(ix,1)
ENDDO

Ez(1,1)=Ez(1,1)+(c*dt/dx)*(By(1,1)-bufR2(1))-&
          (c*dt/dy)*(Bx(1,1)-bufR1(1))-4.0*pi*dt*Jz(1,1)
 
DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X

IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ez(1,:)=0.0
   END IF
      
END IF

IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ez(NXP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ez(:,1)=0.0
   END IF
   
END IF

IF (ymaxp.EQ.ymax) THEN
   
   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Tangent to conductor surface
   Ez(:,NYP)=0.0
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
! - xgp,ygp: Local grid
!
! OUTPUT: Bx,By,Bz at time t+dt/2
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE PUSH_BHALF(Bx,By,Bz,Ex,Ey,Ez,xgp,ygp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bx,By,Bz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION                         :: xminp,xmaxp,yminp,ymaxp
INTEGER, DIMENSION(8)                    :: ngh
INTEGER, PARAMETER                       :: tag1=51,tag2=52
INTEGER                                  :: id,COMM,ierr
DOUBLE PRECISION, ALLOCATABLE            :: bufS1(:),BufR1(:)
DOUBLE PRECISION, ALLOCATABLE            :: bufS2(:),BufR2(:)

! Loop indexes
INTEGER :: ix,iy
!***********************************************************************

xminp=xgp(1)
xmaxp=xgp(NXP)
yminp=ygp(1)
ymaxp=ygp(NYP)

!***********************************************************************
! Bx

ALLOCATE(bufS1(1:NXP),bufR1(1:NXP))

bufS1=Ez(:,2)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,COMM,stat,ierr)
ENDIF

DO iy=1,NYP-1
Bx(:,iy)=Bx(:,iy)-c*dt/(2.0*dy)*(Ez(:,iy+1)-Ez(:,iy))
ENDDO

Bx(:,NYP)=Bx(:,NYP)-c*dt/(2.0*dy)*(bufR1(:)-Ez(:,NYP))

DEALLOCATE(bufS1,bufR1)

!***********************************************************************
! Check boundary conditions along X
   
IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! Normal to conductor surface
   Bx(1,:)=0.0
   END IF
   
END IF
   
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Normal to conductor surface
   Bx(NXP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Inside the conductor
   Bx(:,NYP)=0.0
   END IF

END IF
   
!***********************************************************************
! By

ALLOCATE(bufS2(1:NYP),bufR2(1:NYP))

bufS2=Ez(2,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,COMM,stat,ierr)
ENDIF

DO ix=1,NXP-1
By(ix,:)=By(ix,:)+c*dt/(2.0*dx)*(Ez(ix+1,:)-Ez(ix,:))
ENDDO

By(NXP,:)=By(NXP,:)+c*dt/(2.0*dx)*(bufR2(:)-Ez(NXP,:))

DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X

IF (xmaxp.EQ.xmax) THEN
   
   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Inside the conductor
   By(NXP,:)=0.0
   END IF

END IF
   
!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! Normal to conductor surface
   By(:,1)=0.0
   END IF
   
END IF

IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Normal to conductor surface
   By(:,NYP)=0.0
   END IF

END IF
   
!***********************************************************************
! Bz

ALLOCATE(bufS1(1:NXP),bufR1(1:NXP))
ALLOCATE(bufS2(1:NYP),bufR2(1:NYP))

bufS1=Ex(:,2)
bufS2=Ey(2,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,COMM,stat,ierr)
ENDIF

DO ix=1,NXP-1
DO iy=1,NYP-1
Bz(ix,iy)=Bz(ix,iy)-c*dt/(2.0*dx)*(Ey(ix+1,iy)-Ey(ix,iy))+&
          c*dt/(2.0*dy)*(Ex(ix,iy+1)-Ex(ix,iy))
ENDDO
ENDDO

DO iy=1,NYP-1
Bz(NXP,iy)=Bz(NXP,iy)-c*dt/(2.0*dx)*(bufR2(iy)-Ey(NXP,iy))+&
           c*dt/(2.0*dy)*(Ex(NXP,iy+1)-Ex(NXP,iy))
ENDDO

DO ix=1,NXP-1
Bz(ix,NYP)=Bz(ix,NYP)-c*dt/(2.0*dx)*(Ey(ix+1,NYP)-Ey(ix,NYP))+&
           c*dt/(2.0*dy)*(bufR1(ix)-Ex(ix,NYP))
ENDDO

Bz(NXP,NYP)=Bz(NXP,NYP)-c*dt/(2.0*dx)*(bufR2(NYP)-Ey(NXP,NYP))+&
          c*dt/(2.0*dy)*(bufR1(NXP)-Ex(NXP,NYP))

DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Check boundary conditions along X
   
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   ! Inside the conductor
   Bz(NXP,:)=0.0
   END IF

END IF
   
!***********************************************************************
! Check boundary conditions along Y

IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   ! Inside the conductor
   Bz(:,NYP)=0.0
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
! - xgp,ygp: Local grid
!
! OUTPUT: Ex,Ey,Ez,Bx,By,Bz at time t at the nodes
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE FIELDS_NODES(Bx,By,Bz,Ex,Ey,Ez,Bxg,Byg,Bzg,Exg,Eyg,Ezg,xgp,ygp,&
                        id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bx,By,Bz,Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bxg,Byg,Bzg,Exg,Eyg,Ezg
DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp
DOUBLE PRECISION                         :: xminp,xmaxp,yminp,ymaxp
INTEGER, DIMENSION(8)                    :: ngh
INTEGER, PARAMETER                       :: tag1=51,tag2=52,tag3=53
INTEGER                                  :: id,COMM,ierr
DOUBLE PRECISION, ALLOCATABLE            :: bufS1(:),BufR1(:)
DOUBLE PRECISION, ALLOCATABLE            :: bufS2(:),BufR2(:)
DOUBLE PRECISION                         :: bufS3,BufR3
! Loop indexes
INTEGER :: ix,iy
!***********************************************************************

xminp=xgp(1)
xmaxp=xgp(NXP)
yminp=ygp(1)
ymaxp=ygp(NYP)

!***********************************************************************
! Bxg

ALLOCATE(bufS1(1:NXP),bufR1(1:NXP))

bufS1=Bx(:,NYP-1)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
ENDIF

DO iy=2,NYP
Bxg(:,iy)=(Bx(:,iy)+Bx(:,iy-1))/2.0
ENDDO

Bxg(:,1)=(Bx(:,1)+bufR1(:))/2.0

!***********************************************************************
! Check boundary conditions along Y
   
IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! At the conductor surface
   Bxg(:,1)=(Bx(:,1)+0.0)/2.0
   END IF

END IF
   
DEALLOCATE(bufS1,bufR1)

!***********************************************************************
! Byg

ALLOCATE(bufS2(1:NYP),bufR2(1:NYP))

bufS2=By(NXP-1,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
Byg(ix,:)=(By(ix,:)+By(ix-1,:))/2.0
ENDDO

Byg(1,:)=(By(1,:)+bufR2(:))/2.0

!***********************************************************************
! Check boundary conditions along X
   
IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! At the conductor surface
   Byg(1,:)=(By(1,:)+0.0)/2.0
   END IF

END IF
   
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Bzg

ALLOCATE(bufS1(1:NXP),bufR1(1:NXP))
ALLOCATE(bufS2(1:NYP),bufR2(1:NYP))

bufS1=Bz(:,NYP-1)
bufS2=Bz(NXP-1,:)
bufS3=Bz(NXP-1,NYP-1)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS3,1,MPI_DOUBLE_PRECISION,ngh(5),tag3,&
                  bufR3,1,MPI_DOUBLE_PRECISION,ngh(7),tag3,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
CALL MPI_SENDRECV(bufS3,1,MPI_DOUBLE_PRECISION,ngh(5),tag3,&
                  bufR3,1,MPI_DOUBLE_PRECISION,ngh(7),tag3,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
DO iy=2,NYP
Bzg(ix,iy)=(Bz(ix,iy)+Bz(ix-1,iy)+Bz(ix,iy-1)+Bz(ix-1,iy-1))/4.0
ENDDO
ENDDO

DO ix=2,NXP
Bzg(ix,1)=(Bz(ix,1)+Bz(ix-1,1)+bufR1(ix)+bufR1(ix-1))/4.0
ENDDO

DO iy=2,NYP
Bzg(1,iy)=(Bz(1,iy)+bufR2(iy)+Bz(1,iy-1)+bufR2(iy-1))/4.0
ENDDO

Bzg(1,1)=(Bz(1,1)+BufR2(1)+bufR1(1)+bufR3)/4.0

!***********************************************************************
! Check boundary conditions along X
   
IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   
   DO iy=2,NYP
   Bzg(1,iy)=(Bz(1,iy)+0.0+Bz(1,iy-1)+0.0)/4.0
   ENDDO
   
   Bzg(1,1)=(Bz(1,1)+0.0+bufR1(1)+0.0)/4.0
   
   END IF
   
END IF
   
!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN

   DO ix=2,NXP
   Bzg(ix,1)=(Bz(ix,1)+Bz(ix-1,1)+0.0+0.0)/4.0
   ENDDO
   
   Bzg(1,1)=(Bz(1,1)+BufR2(1)+0.0+0.0)/4.0
      
   END IF

END IF
   
DEALLOCATE(bufS1,bufR1)
DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Exg

ALLOCATE(bufS2(1:NYP),bufR2(1:NYP))

bufS2=Ex(NXP-1,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,&
                  bufR2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,COMM,stat,ierr)
ENDIF

DO ix=2,NXP
Exg(ix,:)=(Ex(ix,:)+Ex(ix-1,:))/2.0
ENDDO

Exg(1,:)=(Ex(1,:)+bufR2(:))/2.0

!***********************************************************************
! Check boundary conditions along X

IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   ! At the conductor surface
   Exg(1,:)=(Ex(1,:)+0.0)/2.0
   END IF

END IF

DEALLOCATE(bufS2,bufR2)

!***********************************************************************
! Eyg

ALLOCATE(bufS1(1:NXP),bufR1(1:NXP))

bufS1=Ey(:,NYP-1)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,&
                  bufR1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,COMM,stat,ierr)
ENDIF

DO iy=2,NYP
Eyg(:,iy)=(Ey(:,iy)+Ey(:,iy-1))/2.0
ENDDO

Eyg(:,1)=(Ey(:,1)+bufR1(:))/2.0

!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   ! At the conductor surface
   Eyg(:,1)=(Ey(:,1)+0.0)/2.0
   END IF

END IF

DEALLOCATE(bufS1,bufR1)

!***********************************************************************
!Ezg

Ezg=Ez

END SUBROUTINE FIELDS_NODES

!***********************************************************************
! Subroutine CORRECT_EFIELD
! This subroutine correct the electric field to ensure the conservation 
! of charge, or to ensure that div(E)=4*pi*rho. Poission equation is solved
! using an iterative method (Gauss-Seidel method), with 5 points.

! INPUT: 
! - Ex: x-component of E on the Yee grid at time t
! - Ey: y-component of E on the Yee grid at time t
! - rho: charge density
!
! OUTPUT: corrected Ex,Ey at time t+dt
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE CORRECT_EFIELD(Ex,Ey,xgp,ygp,rho,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: stat
INTEGER                                  :: id,COMM,ierr
INTEGER, DIMENSION(8)                    :: ngh
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Ex,Ey
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: rho,phi

DOUBLE PRECISION, DIMENSION(1:NXP)       :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)       :: ygp

DOUBLE PRECISION, DIMENSION(1:NYP)       :: bufSE1,bufSW1,bufRE1,bufRW1
DOUBLE PRECISION, DIMENSION(1:NYP)       :: bufSE2,bufRW2

DOUBLE PRECISION, DIMENSION(1:NXP)       :: bufSN1,bufSS1,bufRN1,bufRS1
DOUBLE PRECISION, DIMENSION(1:NXP)       :: bufSN2,bufRS2

INTEGER, PARAMETER                       :: tag1=1,tag2=2,tag3=3,tag4=4
INTEGER, PARAMETER                       :: tag5=5,tag6=6

DOUBLE PRECISION                         :: denom,xminp,yminp,xmaxp,ymaxp

! Loop indexes
INTEGER :: ix,iy,iit

!***********************************************************************

xminp=xgp(1)
xmaxp=xgp(NXP)
yminp=ygp(1)
ymaxp=ygp(NYP)

!***********************************************************************
! denominator
denom=dx*dx+dy*dy

bufSE2=Ex(NXP-1,:)
bufSN2=Ey(:,NYP-1)

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSE2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag1,&
                  bufRW2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN2,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag2,&
                  bufRS2,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag2,COMM,stat,ierr)
                  
ELSE

CALL MPI_SENDRECV(bufSE2,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag1,&
                  bufRW2,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN2,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag2,&
                  bufRS2,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag2,COMM,stat,ierr)

END IF

! Electric potential phi
phi=0.0

!***********************************************************************
! Beginning iteration
!***********************************************************************

DO iit=1,NIT

bufSE1=phi(NXP-1,:)
bufSW1=phi(2,:)
bufSN1=phi(:,NYP-1)
bufSS1=phi(:,2)

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSE1,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag3,&
                  bufRW1,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW1,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag4,&
                  bufRE1,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag5,&
                  bufRS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag6,&
                  bufRN1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag6,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(bufSE1,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag3,&
                  bufRW1,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW1,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag4,&
                  bufRE1,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSN1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag5,&
                  bufRS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag6,&
                  bufRN1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag6,COMM,stat,ierr)

ENDIF

  DO ix=2,NXP-1
    DO iy=2,NYP-1
      phi(ix,iy)=0.5/denom*((phi(ix+1,iy)+phi(ix-1,iy))*dy*dy+&
      (phi(ix,iy+1)+phi(ix,iy-1))*dx*dx+&
      (4.0*pi*rho(ix,iy)-((Ex(ix,iy)-Ex(ix-1,iy))/dx+&
      (Ey(ix,iy)-Ey(ix,iy-1))/dy))*dx*dx*dy*dy)
    ENDDO
  ENDDO

!***********************************************************************
! Computation of phi(ix,1) and phi(1,iy)

DO ix=2,NXP-1
  phi(ix,1)=0.5/denom*((phi(ix+1,1)+phi(ix-1,1))*dy*dy+&
  (phi(ix,2)+bufRS1(ix))*dx*dx+&
  (4.0*pi*rho(ix,1)-((Ex(ix,1)-Ex(ix-1,1))/dx+&
  (Ey(ix,1)-bufRS2(ix))/dy))*dx*dx*dy*dy)
ENDDO

DO iy=2,NYP-1
  phi(1,iy)=0.5/denom*((phi(2,iy)+bufRW1(iy))*dy*dy+&
  (phi(1,iy+1)+phi(1,iy-1))*dx*dx+&
  (4.0*pi*rho(1,iy)-((Ex(1,iy)-bufRW2(iy))/dx+&
  (Ey(1,iy)-Ey(1,iy-1))/dy))*dx*dx*dy*dy)
ENDDO

!***********************************************************************
! Computation of phi(ix,NYP) and phi(NXP,iy)

DO ix=2,NXP-1
  phi(ix,NYP)=0.5/denom*((phi(ix+1,NYP)+phi(ix-1,NYP))*dy*dy+&
  (bufRN1(ix)+phi(ix,NYP-1))*dx*dx+&
  (4.0*pi*rho(ix,NYP)-((Ex(ix,NYP)-Ex(ix-1,NYP))/dx+&
  (Ey(ix,NYP)-Ey(ix,NYP-1))/dy))*dx*dx*dy*dy)
ENDDO

DO iy=2,NYP-1
  phi(NXP,iy)=0.5/denom*((bufRE1(iy)+phi(NXP-1,iy))*dy*dy+&
  (phi(NXP,iy+1)+phi(NXP,iy-1))*dx*dx+&
  (4.0*pi*rho(NXP,iy)-((Ex(NXP,iy)-Ex(NXP-1,iy))/dx+&
  (Ey(NXP,iy)-Ey(NXP,iy-1))/dy))*dx*dx*dy*dy)
ENDDO

!***********************************************************************
! Corners

! Computation of phi(1,1)
phi(1,1)=0.5/denom*((phi(2,1)+bufRW1(1))*dy*dy+&
(phi(1,2)+bufRS1(1))*dx*dx+&
(4.0*pi*rho(1,1)-((Ex(1,1)-bufRW2(1))/dx+&
(Ey(1,1)-bufRS2(1))/dy))*dx*dx*dy*dy)

! Computation of phi(1,NYP)
phi(1,NYP)=0.5/denom*((phi(2,NYP)+bufRW1(NYP))*dy*dy+&
(bufRN1(1)+phi(1,NYP-1))*dx*dx+&
(4.0*pi*rho(1,NYP)-((Ex(1,NYP)-bufRW2(NYP))/dx+&
(Ey(1,NYP)-Ey(1,NYP-1))/dy))*dx*dx*dy*dy)

! Computation of phi(NXP,1)
phi(NXP,1)=0.5/denom*((bufRE1(1)+phi(NXP-1,1))*dy*dy+&
(phi(NXP,2)+bufRS1(NXP))*dx*dx+&
(4.0*pi*rho(NXP,1)-((Ex(NXP,1)-Ex(NXP-1,1))/dx+&
(Ey(NXP,1)-bufRS2(NXP))/dy))*dx*dx*dy*dy)

! Computation of phi(NXP,NYP)
phi(NXP,NYP)=0.5/denom*((bufRE1(NYP)+phi(ix-1,NYP))*dy*dy+&
(bufRN1(NXP)+phi(NXP,NYP-1))*dx*dx+&
(4.0*pi*rho(NXP,NYP)-((Ex(NXP,NYP)-Ex(NXP-1,NYP))/dx+&
(Ey(NXP,NYP)-Ey(NXP,NYP-1))/dy))*dx*dx*dy*dy)

!***********************************************************************
! Check boundary conditions along X

IF (xminp.EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.EQ."METAL") THEN
   phi(1,:)=0.0
   END IF

END IF

IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   phi(NXP,:)=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y

IF (yminp.EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.EQ."METAL") THEN
   phi(:,1)=0.0
   END IF

END IF
   
IF (ymaxp.EQ.ymax) THEN
   
   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   phi(:,NYP)=0.0
   END IF

END IF

!***********************************************************************

ENDDO

!***********************************************************************
! Corrected Yee electric field
!***********************************************************************

bufSS1(:)=phi(:,2)
bufSW1(:)=phi(2,:)

IF (MOD(id,2).EQ.0) THEN

CALL MPI_SENDRECV(bufSS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,&
                  bufRN1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW1,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,&
                  bufRE1,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,COMM,stat,ierr)

ELSE

CALL MPI_SENDRECV(bufSS1,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag1,&
                  bufRN1,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW1,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag2,&
                  bufRE1,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag2,COMM,stat,ierr)

END IF

!***********************************************************************
! Ex

DO ix=1,NXP-1
Ex(ix,:)=Ex(ix,:)-(phi(ix+1,:)-phi(ix,:))/dx
ENDDO

!***********************************************************************
! Check boundary conditions along X

IF (xmaxp.EQ.xmax.AND.BOUND_FIELD_XMAX.NE."PERIODIC") THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   Ex(NXP,:)=Ex(NXP,:)
   END IF

ELSE

! ix=NXP
Ex(NXP,:)=Ex(NXP,:)-(bufRE1(:)-phi(NXP,:))/dx

END IF

!***********************************************************************
! Ey

DO iy=1,NYP-1
Ey(:,iy)=Ey(:,iy)-(phi(:,iy+1)-phi(:,iy))/dy
ENDDO

!***********************************************************************
! Check boundary conditions along Y

IF (ymaxp.EQ.ymax.AND.BOUND_FIELD_YMAX.NE."PERIODIC") THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   Ey(:,NYP)=Ey(:,NYP)
   END IF

ELSE

! iy=NYP
Ey(:,NYP)=Ey(:,NYP)-(bufRN1(:)-phi(:,NYP))/dy

END IF

END SUBROUTINE CORRECT_EFIELD

!***********************************************************************
! Subroutine FILTER_FIELD
! This subroutine operates an isotropic 2D digital filtering of the field, using
! 9 points.
!
! The 9-points Filter matrix expresses as [Birdsall & Langdon, p.441]
!
!      |a    2a     a|
!  M = |2a  1-12a  2a|
!      |a    2a     a|
!
! INPUT: 
! - Field: Field component to filter
!
! OUTPUT: Filtered field Field
!
! REMINDER: North=1,East=2,South=3,West=4,NEast=5,SEast=6,SWest=7,NWest=8
!***********************************************************************

SUBROUTINE FILTER_FIELD(Field,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Field,Field_temp
INTEGER, DIMENSION(8)                    :: ngh
INTEGER, PARAMETER                       :: tag1=1,tag2=2,tag3=3,tag4=4
INTEGER, PARAMETER                       :: tag5=5,tag6=6,tag7=7,tag8=8
INTEGER                                  :: id,COMM,ierr
DOUBLE PRECISION, DIMENSION(1:NXP)       :: bufSN,bufSS,bufRN,bufRS
DOUBLE PRECISION, DIMENSION(1:NYP)       :: bufSE,bufSW,bufRE,bufRW
DOUBLE PRECISION                         :: bufSNE,bufSSE,bufSSW,bufSNW
DOUBLE PRECISION                         :: bufRNE,bufRSE,bufRSW,bufRNW
INTEGER :: ix,iy

!***********************************************************************

bufSN(:)=Field(:,NYP-1)
bufSS(:)=Field(:,2)
bufSE(:)=Field(NXP-1,:)
bufSW(:)=Field(2,:)

bufSNE=Field(NXP-1,NYP-1)
bufSSE=Field(NXP-1,2)
bufSSW=Field(2,2)
bufSNW=Field(2,NYP-1)

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

END IF

! Field smoothing in each domain

DO ix=2,NXP-1
  DO iy=2,NYP-1
  Field_temp(ix,iy)=alpha*(&
             1.0*Field(ix-1,iy+1)+2.0*Field(ix,iy+1)+1.0*Field(ix+1,iy+1)+&
             2.0*Field(ix-1,iy)+(1.0/alpha-12.0)*Field(ix,iy)+2.0*Field(ix+1,iy)+&
             1.0*Field(ix-1,iy-1)+2.0*Field(ix,iy-1)+1.0*Field(ix+1,iy-1))
  ENDDO
ENDDO

! Field smoothing at boundaries

! iy=1
DO ix=2,NXP-1
Field_temp(ix,1)=alpha*(&
           1.0*Field(ix-1,2)+2.0*Field(ix,2)+1.0*Field(ix+1,2)+&
           2.0*Field(ix-1,1)+(1.0/alpha-12.0)*Field(ix,1)+2.0*Field(ix+1,1)+&
           1.0*bufRS(ix-1)+2.0*bufRS(ix)+1.0*bufRS(ix+1))
ENDDO

! iy=NYP
DO ix=2,NXP-1
Field_temp(ix,NYP)=alpha*(&
           1.0*bufRN(ix-1)+2.0*bufRN(ix)+1.0*bufRN(ix+1)+&
           2.0*Field(ix-1,NYP)+(1.0/alpha-12.0)*Field(ix,NYP)+2.0*Field(ix+1,NYP)+&
           1.0*Field(ix-1,NYP-1)+2.0*Field(ix,NYP-1)+1.0*Field(ix+1,NYP-1))
ENDDO

! ix=1
DO iy=2,NYP-1
Field_temp(1,iy)=alpha*(&
           1.0*bufRW(iy+1)+2.0*Field(1,iy+1)+1.0*Field(2,iy+1)+&
           2.0*bufRW(iy)+(1.0/alpha-12.0)*Field(1,iy)+2.0*Field(2,iy)+&
           1.0*bufRW(iy-1)+2.0*Field(1,iy-1)+1.0*Field(2,iy-1))
ENDDO

! ix=NXP
DO iy=2,NYP-1
Field_temp(NXP,iy)=alpha*(&
           1.0*Field(NXP-1,iy+1)+2.0*Field(NXP,iy+1)+1.0*bufRE(iy+1)+&
           2.0*Field(NXP-1,iy)+(1.0/alpha-12.0)*Field(NXP,iy)+2.0*bufRE(iy)+&
           1.0*Field(NXP-1,iy-1)+2.0*Field(NXP,iy-1)+1.0*bufRE(iy-1))
ENDDO

! Smoothing at corners

! ix=1,iy=1
Field_temp(1,1)=alpha*(&
           1.0*bufRW(2)+2.0*Field(1,2)+1.0*Field(2,2)+&
           2.0*bufRW(1)+(1.0/alpha-12.0)*Field(1,1)+2.0*Field(2,1)+&
           1.0*bufRSW+2.0*bufRS(1)+1.0*bufRS(2))

! ix=NXP,iy=1
Field_temp(NXP,1)=alpha*(&
           1.0*Field(NXP-1,2)+2.0*Field(NXP,2)+1.0*bufRE(2)+&
           2.0*Field(NXP-1,1)+(1.0/alpha-12.0)*Field(NXP,1)+2.0*bufRE(1)+&
           1.0*bufRS(NXP-1)+2.0*bufRS(NXP)+1.0*bufRSE)

! ix=1,iy=NYP
Field_temp(1,NYP)=alpha*(&
           1.0*bufRNW+2.0*bufRN(1)+1.0*bufRN(2)+&
           2.0*bufRW(NYP)+(1.0/alpha-12.0)*Field(1,NYP)+2.0*Field(2,NYP)+&
           1.0*bufRW(NYP-1)+2.0*Field(1,NYP-1)+1.0*Field(2,NYP-1))

! ix=NXP,iy=NYP
Field_temp(NXP,NYP)=alpha*(&
           1.0*bufRN(NXP-1)+2.0*bufRN(NXP)+1.0*bufRNE+&
           2.0*Field(NXP-1,NYP)+(1.0/alpha-12.0)*Field(NXP,NYP)+2.0*bufRE(NYP)+&
           1.0*Field(NXP-1,NYP-1)+2.0*Field(NXP,NYP-1)+1.0*bufRE(NYP-1))

! Updating Field
Field(:,:)=Field_temp(:,:)

END SUBROUTINE FILTER_FIELD

!***********************************************************************

END MODULE MOD_FIELDS
