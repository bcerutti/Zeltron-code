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
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE RHOJ(q,pcl,rho,Jx,Jy,Jz,xgp,ygp,zgp,NPP,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER*8                                      :: ip,NPP
DOUBLE PRECISION                               :: q,rhop,gam,vx,vy,vz
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)         :: pcl
DOUBLE PRECISION                               :: x,y,z,ux,uy,uz,wght
DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: rho,Jx,Jy,Jz
DOUBLE PRECISION                               :: fp,fq,fr
INTEGER                                        :: i,j,l
INTEGER, DIMENSION(26)                         :: ngh
INTEGER                                        :: id,COMM,ierr

!***********************************************************************

rho=0.0
Jx=0.0
Jy=0.0
Jz=0.0

DO ip=1,NPP

  x=pcl(1,ip)
  y=pcl(2,ip)
  z=pcl(3,ip)
  ux=pcl(4,ip)
  uy=pcl(5,ip)
  uz=pcl(6,ip)
  wght=pcl(7,ip)

  ! Charge density of 1 super-particle in the cell (i,j,l)
  rhop=q*e*wght/(dx*dy*dz)
  
  ! Lorentz factor
  gam=sqrt(1.0+ux*ux+uy*uy+uz*uz)

  ! 3-velocity components
  vx=ux*c/gam
  vy=uy*c/gam
  vz=uz*c/gam

  ! Computation of the nearest node index to (xf,yf,zf), for a constant dx,dy,dz
  i=FLOOR((x-xgp(1))/dx)+1
  j=FLOOR((y-ygp(1))/dy)+1
  l=FLOOR((z-zgp(1))/dz)+1
  
  IF (i.EQ.NXP) THEN
  i=i-1
  END IF
  
  IF (j.EQ.NYP) THEN
  j=j-1
  END IF
  
  IF (l.EQ.NZP) THEN
  l=l-1
  END IF
  
  fp=(x-xgp(i))/dx
  fq=(y-ygp(j))/dy
  fr=(z-zgp(l))/dz

  !=====================================================================
  ! Charge density

  rho(i,j,l)=rho(i,j,l)+rhop*(1.0-fp)*(1.0-fq)*(1.0-fr)
  rho(i+1,j,l)=rho(i+1,j,l)+rhop*fp*(1.0-fq)*(1.0-fr)
  rho(i,j+1,l)=rho(i,j+1,l)+rhop*(1.0-fp)*fq*(1.0-fr)
  rho(i,j,l+1)=rho(i,j,l+1)+rhop*(1.0-fp)*(1.0-fq)*fr
  rho(i+1,j+1,l)=rho(i+1,j+1,l)+rhop*fp*fq*(1.0-fr)
  rho(i+1,j,l+1)=rho(i+1,j,l+1)+rhop*fp*(1.0-fq)*fr
  rho(i,j+1,l+1)=rho(i,j+1,l+1)+rhop*(1.0-fp)*fq*fr
  rho(i+1,j+1,l+1)=rho(i+1,j+1,l+1)+rhop*fp*fq*fr

  !=====================================================================
  ! Current density

  ! x-component:
  Jx(i,j,l)=Jx(i,j,l)+0.5*rhop*vx*(1.0-fp)*(1.0-fq)*(1.0-fr)
  Jx(i+1,j,l)=Jx(i+1,j,l)+0.5*rhop*vx*fp*(1.0-fq)*(1.0-fr)
  Jx(i,j+1,l)=Jx(i,j+1,l)+0.5*rhop*vx*(1.0-fp)*fq*(1.0-fr)
  Jx(i,j,l+1)=Jx(i,j,l+1)+0.5*rhop*vx*(1.0-fp)*(1.0-fq)*fr
  Jx(i+1,j+1,l)=Jx(i+1,j+1,l)+0.5*rhop*vx*fp*fq*(1.0-fr)
  Jx(i+1,j,l+1)=Jx(i+1,j,l+1)+0.5*rhop*vx*fp*(1.0-fq)*fr
  Jx(i,j+1,l+1)=Jx(i,j+1,l+1)+0.5*rhop*vx*(1.0-fp)*fq*fr
  Jx(i+1,j+1,l+1)=Jx(i+1,j+1,l+1)+0.5*rhop*vx*fp*fq*fr

  ! y-component:
  Jy(i,j,l)=Jy(i,j,l)+0.5*rhop*vy*(1.0-fp)*(1.0-fq)*(1.0-fr)
  Jy(i+1,j,l)=Jy(i+1,j,l)+0.5*rhop*vy*fp*(1.0-fq)*(1.0-fr)
  Jy(i,j+1,l)=Jy(i,j+1,l)+0.5*rhop*vy*(1.0-fp)*fq*(1.0-fr)
  Jy(i,j,l+1)=Jy(i,j,l+1)+0.5*rhop*vy*(1.0-fp)*(1.0-fq)*fr
  Jy(i+1,j+1,l)=Jy(i+1,j+1,l)+0.5*rhop*vy*fp*fq*(1.0-fr)
  Jy(i+1,j,l+1)=Jy(i+1,j,l+1)+0.5*rhop*vy*fp*(1.0-fq)*fr
  Jy(i,j+1,l+1)=Jy(i,j+1,l+1)+0.5*rhop*vy*(1.0-fp)*fq*fr
  Jy(i+1,j+1,l+1)=Jy(i+1,j+1,l+1)+0.5*rhop*vy*fp*fq*fr

  ! z-component:
  Jz(i,j,l)=Jz(i,j,l)+0.5*rhop*vz*(1.0-fp)*(1.0-fq)*(1.0-fr)
  Jz(i+1,j,l)=Jz(i+1,j,l)+0.5*rhop*vz*fp*(1.0-fq)*(1.0-fr)
  Jz(i,j+1,l)=Jz(i,j+1,l)+0.5*rhop*vz*(1.0-fp)*fq*(1.0-fr)
  Jz(i,j,l+1)=Jz(i,j,l+1)+0.5*rhop*vz*(1.0-fp)*(1.0-fq)*fr
  Jz(i+1,j+1,l)=Jz(i+1,j+1,l)+0.5*rhop*vz*fp*fq*(1.0-fr)
  Jz(i+1,j,l+1)=Jz(i+1,j,l+1)+0.5*rhop*vz*fp*(1.0-fq)*fr
  Jz(i,j+1,l+1)=Jz(i,j+1,l+1)+0.5*rhop*vz*(1.0-fp)*fq*fr
  Jz(i+1,j+1,l+1)=Jz(i+1,j+1,l+1)+0.5*rhop*vz*fp*fq*fr

ENDDO

CALL COM_RHOJ(rho,xgp,ygp,zgp,id,ngh,COMM,ierr)
CALL COM_RHOJ(Jx,xgp,ygp,zgp,id,ngh,COMM,ierr)
CALL COM_RHOJ(Jy,xgp,ygp,zgp,id,ngh,COMM,ierr)
CALL COM_RHOJ(Jz,xgp,ygp,zgp,id,ngh,COMM,ierr)

END SUBROUTINE RHOJ

!***********************************************************************
! Subroutine COM_RHOJ
! This subroutine computes the charge density RHO and current density J
! on grid nodes.

! INPUT: 
! - fun: RHO or Jx,y,z
! - local grid: xgp,ygp,zgp
!
! OUTPUT: RHO or Jx,y,z at the boundaries
!
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE COM_RHOJ(fun,xgp,ygp,zgp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: fun
DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp
INTEGER, DIMENSION(26)                         :: ngh
INTEGER                                        :: id,COMM,ierr

INTEGER, PARAMETER                             :: tag1=1,tag2=2,tag3=3,tag4=4
INTEGER, PARAMETER                             :: tag5=5,tag6=6,tag7=7,tag8=8
INTEGER, PARAMETER                             :: tag9=9,tag10=10,tag11=11,tag12=12
INTEGER, PARAMETER                             :: tag13=13,tag14=14,tag15=15,tag16=16
INTEGER, PARAMETER                             :: tag17=17,tag18=18,tag19=19,tag20=20
INTEGER, PARAMETER                             :: tag21=21,tag22=22,tag23=23,tag24=24
INTEGER, PARAMETER                             :: tag25=25,tag26=26

DOUBLE PRECISION, DIMENSION(1:NXP,1:NZP)       :: bufSN,bufSS,bufRN,bufRS
DOUBLE PRECISION, DIMENSION(1:NYP,1:NZP)       :: bufSE,bufSW,bufRE,bufRW
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP)       :: bufSB,bufSF,bufRB,bufRF

DOUBLE PRECISION, DIMENSION(1:NZP)             :: bufSNE,bufSSE,bufSSW,bufSNW
DOUBLE PRECISION, DIMENSION(1:NZP)             :: bufRNE,bufRSE,bufRSW,bufRNW

DOUBLE PRECISION, DIMENSION(1:NXP)             :: bufSFS,bufSFN,bufSBN,bufSBS
DOUBLE PRECISION, DIMENSION(1:NXP)             :: bufRFS,bufRFN,bufRBN,bufRBS

DOUBLE PRECISION, DIMENSION(1:NYP)             :: bufSBW,bufSBE,bufSFW,bufSFE
DOUBLE PRECISION, DIMENSION(1:NYP)             :: bufRBW,bufRBE,bufRFW,bufRFE

DOUBLE PRECISION                               :: bufSBSW,bufSBSE,bufSBNW,bufSFSW
DOUBLE PRECISION                               :: bufSBNE,bufSFNW,bufSFSE,bufSFNE
DOUBLE PRECISION                               :: bufRBSW,bufRBSE,bufRBNW,bufRFSW
DOUBLE PRECISION                               :: bufRBNE,bufRFNW,bufRFSE,bufRFNE

!***********************************************************************

! Faces of the cube

bufSS=fun(:,1,:)
bufSN=fun(:,NYP,:)

bufSW=fun(1,:,:)
bufSE=fun(NXP,:,:)

bufSB=fun(:,:,1)
bufSF=fun(:,:,NZP)

! Aretes of the cube

bufSNE=fun(NXP,NYP,:)
bufSSE=fun(NXP,1,:)
bufSSW=fun(1,1,:)
bufSNW=fun(1,NYP,:)

bufSFS=fun(:,1,NZP)
bufSFN=fun(:,NYP,NZP)
bufSBN=fun(:,NYP,1)
bufSBS=fun(:,1,1)

bufSBW=fun(1,:,1)
bufSBE=fun(NXP,:,1)
bufSFW=fun(1,:,NZP)
bufSFE=fun(NXP,:,NZP)

! Corners

bufSBSW=fun(1,1,1)
bufSBSE=fun(NXP,1,1)
bufSBNW=fun(1,NYP,1)
bufSFSW=fun(1,1,NZP)
bufSBNE=fun(NXP,NYP,1)
bufSFNW=fun(1,NYP,NZP)
bufSFSE=fun(NXP,1,NZP)
bufSFNE=fun(NXP,NYP,NZP)

IF (MOD(id,2).EQ.0) THEN

! Exchange at the faces of the cube:

CALL MPI_SENDRECV(bufSN,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufRS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag2,&
                  bufRN,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSE,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag3,&
                  bufRW,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag4,&
                  bufRE,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSB,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag5,&
                  bufRF,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag6,&
                  bufRB,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag6,COMM,stat,ierr)

! Exchange at the aretes of the cube:

CALL MPI_SENDRECV(bufSNE,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag7,&
                  bufRSW,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSE,NZP,MPI_DOUBLE_PRECISION,ngh(15),tag8,&
                  bufRNW,NZP,MPI_DOUBLE_PRECISION,ngh(17),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSW,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag9,&
                  bufRNE,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag9,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNW,NZP,MPI_DOUBLE_PRECISION,ngh(17),tag10,&
                  bufRSE,NZP,MPI_DOUBLE_PRECISION,ngh(15),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFS,NXP,MPI_DOUBLE_PRECISION,ngh(20),tag11,&
                  bufRBN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag11,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFN,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag12,&
                  bufRBS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag12,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag13,&
                  bufRFS,NXP,MPI_DOUBLE_PRECISION,ngh(20),tag13,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag14,&
                  bufRFN,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag14,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag15,&
                  bufRFE,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag16,&
                  bufRFW,NYP,MPI_DOUBLE_PRECISION,ngh(21),tag16,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFW,NYP,MPI_DOUBLE_PRECISION,ngh(21),tag17,&
                  bufRBE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag17,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFE,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag18,&
                  bufRBW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag18,COMM,stat,ierr)

! Exchange at the corners of the cube:

CALL MPI_SENDRECV(bufSBSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag19,&
                  bufRFNE,1,MPI_DOUBLE_PRECISION,ngh(22),tag19,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag20,&
                  bufRFNW,1,MPI_DOUBLE_PRECISION,ngh(25),tag20,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag21,&
                  bufRFSE,1,MPI_DOUBLE_PRECISION,ngh(23),tag21,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFSW,1,MPI_DOUBLE_PRECISION,ngh(24),tag22,&
                  bufRBNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag22,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag23,&
                  bufRFSW,1,MPI_DOUBLE_PRECISION,ngh(24),tag23,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFNW,1,MPI_DOUBLE_PRECISION,ngh(25),tag24,&
                  bufRBSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag24,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFSE,1,MPI_DOUBLE_PRECISION,ngh(23),tag25,&
                  bufRBNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag25,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFNE,1,MPI_DOUBLE_PRECISION,ngh(22),tag26,&
                  bufRBSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag26,COMM,stat,ierr)

ELSE

! Exchange at the faces of the cube:

CALL MPI_SENDRECV(bufSN,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag1,&
                  bufRS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag1,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag2,&
                  bufRN,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag2,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSE,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag3,&
                  bufRW,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag3,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSW,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag4,&
                  bufRE,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag4,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSB,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag5,&
                  bufRF,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag5,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSF,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag6,&
                  bufRB,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag6,COMM,stat,ierr)

! Exchange at the aretes of the cube:

CALL MPI_SENDRECV(bufSNE,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag7,&
                  bufRSW,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag7,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSE,NZP,MPI_DOUBLE_PRECISION,ngh(15),tag8,&
                  bufRNW,NZP,MPI_DOUBLE_PRECISION,ngh(17),tag8,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSSW,NZP,MPI_DOUBLE_PRECISION,ngh(16),tag9,&
                  bufRNE,NZP,MPI_DOUBLE_PRECISION,ngh(14),tag9,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSNW,NZP,MPI_DOUBLE_PRECISION,ngh(17),tag10,&
                  bufRSE,NZP,MPI_DOUBLE_PRECISION,ngh(15),tag10,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFS,NXP,MPI_DOUBLE_PRECISION,ngh(20),tag11,&
                  bufRBN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag11,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFN,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag12,&
                  bufRBS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag12,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBN,NXP,MPI_DOUBLE_PRECISION,ngh(1),tag13,&
                  bufRFS,NXP,MPI_DOUBLE_PRECISION,ngh(20),tag13,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBS,NXP,MPI_DOUBLE_PRECISION,ngh(3),tag14,&
                  bufRFN,NXP,MPI_DOUBLE_PRECISION,ngh(18),tag14,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag15,&
                  bufRFE,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag15,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag16,&
                  bufRFW,NYP,MPI_DOUBLE_PRECISION,ngh(21),tag16,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFW,NYP,MPI_DOUBLE_PRECISION,ngh(21),tag17,&
                  bufRBE,NYP,MPI_DOUBLE_PRECISION,ngh(2),tag17,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFE,NYP,MPI_DOUBLE_PRECISION,ngh(19),tag18,&
                  bufRBW,NYP,MPI_DOUBLE_PRECISION,ngh(4),tag18,COMM,stat,ierr)

! Exchange at the corners of the cube:

CALL MPI_SENDRECV(bufSBSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag19,&
                  bufRFNE,1,MPI_DOUBLE_PRECISION,ngh(22),tag19,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag20,&
                  bufRFNW,1,MPI_DOUBLE_PRECISION,ngh(25),tag20,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag21,&
                  bufRFSE,1,MPI_DOUBLE_PRECISION,ngh(23),tag21,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFSW,1,MPI_DOUBLE_PRECISION,ngh(24),tag22,&
                  bufRBNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag22,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSBNE,1,MPI_DOUBLE_PRECISION,ngh(5),tag23,&
                  bufRFSW,1,MPI_DOUBLE_PRECISION,ngh(24),tag23,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFNW,1,MPI_DOUBLE_PRECISION,ngh(25),tag24,&
                  bufRBSE,1,MPI_DOUBLE_PRECISION,ngh(6),tag24,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFSE,1,MPI_DOUBLE_PRECISION,ngh(23),tag25,&
                  bufRBNW,1,MPI_DOUBLE_PRECISION,ngh(8),tag25,COMM,stat,ierr)

CALL MPI_SENDRECV(bufSFNE,1,MPI_DOUBLE_PRECISION,ngh(22),tag26,&
                  bufRBSW,1,MPI_DOUBLE_PRECISION,ngh(7),tag26,COMM,stat,ierr)

ENDIF

!***********************************************************************
! Check boundary conditions along X
   
IF (xgp(1).EQ.xmin) THEN

   IF (BOUND_FIELD_XMIN.NE."PERIODIC") THEN
   bufRW(:,:)=0.0
   bufRBW(:)=0.0
   bufRFW(:)=0.0
   bufRNW(:)=0.0
   bufRSW(:)=0.0
   bufRBNW=0.0
   bufRBSW=0.0
   bufRFNW=0.0
   bufRFSW=0.0
   END IF

END IF

IF (xgp(NXP).EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.NE."PERIODIC") THEN
   bufRE(:,:)=0.0
   bufRBE(:)=0.0
   bufRFE(:)=0.0
   bufRNE(:)=0.0
   bufRSE(:)=0.0
   bufRBNE=0.0
   bufRBSE=0.0
   bufRFNE=0.0
   bufRFSE=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Y
   
IF (ygp(1).EQ.ymin) THEN

   IF (BOUND_FIELD_YMIN.NE."PERIODIC") THEN
   bufRS(:,:)=0.0
   bufRBS(:)=0.0
   bufRFS(:)=0.0
   bufRSW(:)=0.0
   bufRSE(:)=0.0
   bufRBSW=0.0
   bufRBSE=0.0
   bufRFSW=0.0
   bufRFSE=0.0
   END IF

END IF

IF (ygp(NYP).EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.NE."PERIODIC") THEN
   bufRN(:,:)=0.0
   bufRBN(:)=0.0
   bufRFN(:)=0.0
   bufRNW(:)=0.0
   bufRNE(:)=0.0
   bufRBNW=0.0
   bufRBNE=0.0
   bufRFNW=0.0
   bufRFNE=0.0
   END IF

END IF

!***********************************************************************
! Check boundary conditions along Z
   
IF (zgp(1).EQ.zmin) THEN

   IF (BOUND_FIELD_ZMIN.NE."PERIODIC") THEN
   bufRB(:,:)=0.0
   bufRBN(:)=0.0
   bufRBS(:)=0.0
   bufRBW(:)=0.0
   bufRBE(:)=0.0
   bufRBSW=0.0
   bufRBSE=0.0
   bufRBNW=0.0
   bufRBNE=0.0
   END IF

END IF

IF (zgp(NZP).EQ.zmax) THEN

   IF (BOUND_FIELD_ZMAX.NE."PERIODIC") THEN
   bufRF(:,:)=0.0
   bufRFN(:)=0.0
   bufRFS(:)=0.0
   bufRFW(:)=0.0
   bufRFE(:)=0.0
   bufRFSW=0.0
   bufRFSE=0.0
   bufRFNW=0.0
   bufRFNE=0.0
   END IF

END IF

!***********************************************************************

! Faces of the cube

fun(NXP,:,:)=fun(NXP,:,:)+bufRE(:,:)
fun(1,:,:)=fun(1,:,:)+bufRW(:,:)

fun(:,1,:)=fun(:,1,:)+bufRS(:,:)
fun(:,NYP,:)=fun(:,NYP,:)+bufRN(:,:)

fun(:,:,1)=fun(:,:,1)+bufRB(:,:)
fun(:,:,NZP)=fun(:,:,NZP)+bufRF(:,:)

! Aretes of the cube

fun(1,1,:)=fun(1,1,:)+bufRSW(:)
fun(NXP,1,:)=fun(NXP,1,:)+bufRSE(:)
fun(1,NYP,:)=fun(1,NYP,:)+bufRNW(:)
fun(NXP,NYP,:)=fun(NXP,NYP,:)+bufRNE(:)

fun(:,1,1)=fun(:,1,1)+bufRBS(:)
fun(:,NYP,1)=fun(:,NYP,1)+bufRBN(:)
fun(:,1,NZP)=fun(:,1,NZP)+bufRFS(:)
fun(:,NYP,NZP)=fun(:,NYP,NZP)+bufRFN(:)

fun(1,:,1)=fun(1,:,1)+bufRBW(:)
fun(NXP,:,1)=fun(NXP,:,1)+bufRBE(:)
fun(1,:,NZP)=fun(1,:,NZP)+bufRFW(:)
fun(NXP,:,NZP)=fun(NXP,:,NZP)+bufRFE(:)

! Corners of the cube

fun(1,1,1)=fun(1,1,1)+bufRBSW
fun(NXP,1,1)=fun(NXP,1,1)+bufRBSE
fun(1,NYP,1)=fun(1,NYP,1)+bufRBNW
fun(1,1,NZP)=fun(1,1,NZP)+bufRFSW
fun(NXP,NYP,1)=fun(NXP,NYP,1)+bufRBNE
fun(NXP,1,NZP)=fun(NXP,1,NZP)+bufRFSE
fun(1,NYP,NZP)=fun(1,NYP,NZP)+bufRFNW
fun(NXP,NYP,NZP)=fun(NXP,NYP,NZP)+bufRFNE

END SUBROUTINE COM_RHOJ

!***********************************************************************
! Subroutine JYEE
! This subroutine computes the current density in the Yee lattice.

! INPUT: 
! - Jx: x-component of the current density J at the nodes at t=t+dt/2
! - Jy: y-component of the current density J at the nodes at t=t+dt/2
! - Jz: z-component of the current density J at the nodes at t=t+dt/2
! - grid: xgp,ygp,zgp
!
! OUTPUT: Jx, Jy, Jz at t=t+dt/2 in the Yee lattice
!
! REMINDER:
!
! BNorth=1,BEast=2,BSouth=3,BWest=4,BNEast=5,BSEast=6,BSWest=7,BNWest=8
! BCenter=9,North=10,East=11,South=12,West=13,NEast=14,SEast=15,SWest=16
! NWest=17,FNorth=18,FEast=19,FSouth=20,FWest=21,FNEast=22,FSEast=23,FSWest=24
! FNWest=25,FCenter=26
!
!***********************************************************************

SUBROUTINE JYEE(Jx,Jy,Jz,xgp,ygp,zgp,id,ngh,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: stat
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Jx,Jy,Jz
INTEGER                                        :: ix,iy,iz
DOUBLE PRECISION, DIMENSION(1:NXP)             :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP)             :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP)             :: zgp
DOUBLE PRECISION                               :: xmaxp,ymaxp,zmaxp
INTEGER, DIMENSION(26)                         :: ngh
INTEGER, PARAMETER                             :: tag=1
INTEGER                                        :: id,COMM,ierr
DOUBLE PRECISION, ALLOCATABLE                  :: bufS(:,:),bufR(:,:)

!***********************************************************************

xmaxp=xgp(NXP)
ymaxp=ygp(NYP)
zmaxp=zgp(NZP)

!***********************************************************************
! Jx

ALLOCATE(bufS(1:NYP,1:NZP),bufR(1:NYP,1:NZP))

bufS=Jx(2,:,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag,&
                  bufR,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(13),tag,&
                  bufR,NYP*NZP,MPI_DOUBLE_PRECISION,ngh(11),tag,COMM,stat,ierr)
ENDIF

DO ix=1,NXP-1
Jx(ix,:,:)=(Jx(ix+1,:,:)+Jx(ix,:,:))/2.0
ENDDO

!***********************************************************************
! Check boundary conditions along X
   
IF (xmaxp.EQ.xmax) THEN

   IF (BOUND_FIELD_XMAX.EQ."METAL") THEN
   bufR(:,:)=0.0
   ENDIF

END IF
!***********************************************************************

Jx(NXP,:,:)=(bufR(:,:)+Jx(NXP,:,:))/2.0

DEALLOCATE(bufS,bufR)

!***********************************************************************
! Jy

ALLOCATE(bufS(1:NXP,1:NZP),bufR(1:NXP,1:NZP))

bufS=Jy(:,2,:)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag,&
                  bufR,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(12),tag,&
                  bufR,NXP*NZP,MPI_DOUBLE_PRECISION,ngh(10),tag,COMM,stat,ierr)
ENDIF

DO iy=1,NYP-1
Jy(:,iy,:)=(Jy(:,iy+1,:)+Jy(:,iy,:))/2.0
ENDDO

!***********************************************************************
! Check boundary conditions along Y
   
IF (ymaxp.EQ.ymax) THEN

   IF (BOUND_FIELD_YMAX.EQ."METAL") THEN
   bufR(:,:)=0.0
   ENDIF

END IF
!***********************************************************************

Jy(:,NYP,:)=(bufR(:,:)+Jy(:,NYP,:))/2.0

DEALLOCATE(bufS,bufR)

!***********************************************************************
! Jz

ALLOCATE(bufS(1:NXP,1:NYP),bufR(1:NXP,1:NYP))

bufS=Jz(:,:,2)

IF (MOD(id,2).EQ.0) THEN
CALL MPI_SENDRECV(bufS,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag,&
                  bufR,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag,COMM,stat,ierr)
ELSE
CALL MPI_SENDRECV(bufS,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(9),tag,&
                  bufR,NXP*NYP,MPI_DOUBLE_PRECISION,ngh(26),tag,COMM,stat,ierr)
ENDIF

DO iz=1,NZP-1
Jz(:,:,iz)=(Jz(:,:,iz+1)+Jz(:,:,iz))/2.0
ENDDO

!***********************************************************************
! Check boundary conditions along Z
   
IF (zmaxp.EQ.zmax) THEN

   IF (BOUND_FIELD_ZMAX.EQ."METAL") THEN
   bufR(:,:)=0.0
   ENDIF

END IF
!***********************************************************************

Jz(:,:,NZP)=(bufR(:,:)+Jz(:,:,NZP))/2.0

DEALLOCATE(bufS,bufR)

END SUBROUTINE JYEE

!***********************************************************************

END MODULE MOD_RHOJ
