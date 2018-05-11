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

MODULE MOD_INTERP

USE MOD_INPUT
USE MOD_MPI

IMPLICIT NONE

PRIVATE

PUBLIC :: TRILINEAR_FIELD ! Bilinear interpolation of a field value
PUBLIC :: TRILINEAR_FIELDS ! Bilinear interpolation of fields

 CONTAINS

!***********************************************************************
! Function TRILINEAR_FIELD: bilinear interpolation of the E and B fields
! INPUT: F (2D field to interpolate), 
! (i,j,k) lower corner of cell within which to inteprolate
! xq, yq, zq : weights within the cell for interpolation
!          (xq = 0 means use F at i, xq = 1 means use F at i+1)
!          (yq = 0 means use F at j, yq = 1 means use F at j+1)
!          (zq = 0 means use F at k, zq = 1 means use F at k+1)
! OUTPUT: interpolated value
!***********************************************************************

FUNCTION TRILINEAR_FIELD(F, i, j, k, xq, yq, zq)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: F
INTEGER                                        :: i,j,k
DOUBLE PRECISION                               :: xq,yq,zq
INTEGER                                        :: ip,jp,kp

DOUBLE PRECISION                               :: f000,f001,f010,f011
DOUBLE PRECISION                               :: f100,f101,f110,f111
DOUBLE PRECISION                               :: TRILINEAR_FIELD
!***********************************************************************

! Don't examine i+1 if it has zero weight, in case i+1 is past the array
ip = i
jp = j
kp = k
IF (xq > 0.) ip = ip + 1
IF (yq > 0.) jp = jp + 1
IF (zq > 0.) kp = kp + 1
f000=F(i,j,k)
f100=F(ip,j,k)
f010=F(i,jp,k)
f110=F(ip,jp,k)
f001=F(i,j,kp)
f101=F(ip,j,kp)
f011=F(i,jp,kp)
f111=F(ip,jp,kp)

TRILINEAR_FIELD = (1.-xq) &
 * ((1.-yq)*(1.-zq)*f000 + yq*(1.-zq)*f010 + yq*zq*f011 + (1.-yq)*zq*f001) &
 + xq &
 * ((1.-yq)*(1.-zq)*f100 + yq*(1.-zq)*f110 + yq*zq*f111 + (1.-yq)*zq*f101)

END FUNCTION

!***********************************************************************
! Function TRILINEAR_FIELDS: bilinear interpolation of the E and B fields
! INPUT: Bx,By,Bz,Ex,Ey,Ez (2D functions to interpolate), 
! (xi,yi) (initial grid where fun2d is defined), and (xf,yf) (final grid).
! OUTPUT: bilinear_fields: fun2d calculated in the (xf,yf) grid.
!***********************************************************************

SUBROUTINE TRILINEAR_FIELDS(xi,yi,zi,xf,yf,zf, Finterp, F1,F2,F3,F4,F5,F6)

IMPLICIT NONE

! INPUT
DOUBLE PRECISION, DIMENSION(1:NXP)          :: xi
DOUBLE PRECISION, DIMENSION(1:NYP)          :: yi
DOUBLE PRECISION, DIMENSION(1:NZP)          :: zi
DOUBLE PRECISION                            :: xf,yf,zf
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP)           :: F1
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP), OPTIONAL :: F2,F3,F4,F5,F6
! OUTPUT 
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Finterp

! indices
INTEGER                                     :: i,j,k
! weights
DOUBLE PRECISION                            :: xq,yq,zq
!***********************************************************************

! Computation of the nearest node index to (xf,yf), for a constant dx,dy
i=FLOOR((xf-xi(1))/dx)+1
j=FLOOR((yf-yi(1))/dy)+1
k=FLOOR((zf-zi(1))/dz)+1

IF (i.EQ.NXP) THEN
i=i-1
END IF
  
IF (j.EQ.NYP) THEN
j=j-1
END IF
  
IF (k.EQ.NZP) THEN
k=k-1
END IF

xq=(xf-xi(i))/dx
yq=(yf-yi(j))/dy
zq=(zf-zi(k))/dz

Finterp(1) = TRILINEAR_FIELD(F1, i, j, k, xq, yq, zq)
IF (PRESENT(F2)) Finterp(2) = TRILINEAR_FIELD(F2, i, j, k, xq, yq, zq)
IF (PRESENT(F3)) Finterp(3) = TRILINEAR_FIELD(F3, i, j, k, xq, yq, zq)
IF (PRESENT(F4)) Finterp(4) = TRILINEAR_FIELD(F4, i, j, k, xq, yq, zq)
IF (PRESENT(F5)) Finterp(5) = TRILINEAR_FIELD(F5, i, j, k, xq, yq, zq)
IF (PRESENT(F6)) Finterp(6) = TRILINEAR_FIELD(F6, i, j, k, xq, yq, zq)

END SUBROUTINE TRILINEAR_FIELDS

!***********************************************************************

END MODULE MOD_INTERP
