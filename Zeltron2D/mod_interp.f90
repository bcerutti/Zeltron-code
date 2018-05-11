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

PUBLIC :: BILINEAR_FIELD ! Bilinear interpolation of a field value
PUBLIC :: BILINEAR_FIELDS ! Bilinear interpolation of fields

 CONTAINS

!***********************************************************************
! Function BILINEAR_FIELD: bilinear interpolation of the E and B fields
! INPUT: F (2D field to interpolate), 
! (i,j) lower corner of cell within which to inteprolate
! xp, yq : weights within the cell for interpolation
!          (xp = 0 means use F at i, xp = 1 means use F at i+1)
!          (yq = 0 means use F at j, yq = 1 means use F at j+1)
! OUTPUT: interpolated value
!***********************************************************************

FUNCTION BILINEAR_FIELD(F, i, j, xp, yq)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP)    :: F
INTEGER                                     :: i,j
DOUBLE PRECISION                            :: xp,yq
INTEGER                                     :: ip,jp

DOUBLE PRECISION                            :: f00,f10,f01,f11
DOUBLE PRECISION                            :: BILINEAR_FIELD
!***********************************************************************

! Don't examine i+1 if it has zero weight, in case i+1 is past the array
ip = i
jp = j
IF (xp > 0.) ip = ip + 1
IF (yq > 0.) jp = jp + 1
f00=F(i,j)
f10=F(ip,j)
f01=F(i,jp)
f11=F(ip,jp)

BILINEAR_FIELD = f00*(1.0-xp)*(1.0-yq) + f10*xp*(1.0-yq) &
                +f01*(1.0-xp)*yq       + f11*xp*yq

END FUNCTION

!***********************************************************************
! Function BILINEAR_FIELDS: bilinear interpolation of the E and B fields
! INPUT: Bx,By,Bz,Ex,Ey,Ez (2D functions to interpolate), 
! (xi,yi) (initial grid where fun2d is defined), and (xf,yf) (final grid).
! OUTPUT: bilinear_fields: fun2d calculated in the (xf,yf) grid.
!***********************************************************************

SUBROUTINE BILINEAR_FIELDS(xi,yi,xf,yf, Finterp, F1,F2,F3,F4,F5,F6)

IMPLICIT NONE

! INPUT
DOUBLE PRECISION, DIMENSION(1:NXP)          :: xi
DOUBLE PRECISION, DIMENSION(1:NYP)          :: yi
DOUBLE PRECISION                            :: xf,yf
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP)    :: F1
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP), OPTIONAL :: F2,F3,F4,F5,F6
! OUTPUT 
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Finterp

! indices
INTEGER                                     :: i,j
! weights
DOUBLE PRECISION                            :: xp,yq
!***********************************************************************

! Computation of the nearest node index to (xf,yf), for a constant dx,dy
i=FLOOR((xf-xi(1))/dx)+1
j=FLOOR((yf-yi(1))/dy)+1

IF (i.EQ.NXP) THEN
i=i-1
END IF
  
IF (j.EQ.NYP) THEN
j=j-1
END IF

xp=(xf-xi(i))/dx
yq=(yf-yi(j))/dy

Finterp(1) = BILINEAR_FIELD(F1, i, j, xp, yq)
IF (PRESENT(F2)) Finterp(2) = BILINEAR_FIELD(F2, i, j, xp, yq)
IF (PRESENT(F3)) Finterp(3) = BILINEAR_FIELD(F3, i, j, xp, yq)
IF (PRESENT(F4)) Finterp(4) = BILINEAR_FIELD(F4, i, j, xp, yq)
IF (PRESENT(F5)) Finterp(5) = BILINEAR_FIELD(F5, i, j, xp, yq)
IF (PRESENT(F6)) Finterp(6) = BILINEAR_FIELD(F6, i, j, xp, yq)

END SUBROUTINE BILINEAR_FIELDS

END MODULE MOD_INTERP
