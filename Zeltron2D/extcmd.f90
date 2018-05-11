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

PROGRAM EXTCMD

USE MOD_EXTCMD

IMPLICIT NONE

INTEGER :: cmd, step

PRINT *, "Directions:"
PRINT *, "1. Create a file simCmd.txt (in any directory but data/)"
PRINT *, "2. simCmd.txt should contain one line; at the time of this writing"
PRINT *, "   it contains 4 numbers: timeStep quit dump save"
PRINT *, "   where quit, dump, and save are either 0 (don't do) or 1 (do)."
PRINT *, "3. Run ./testExtcmd to verify that simCmd is what it should be."
PRINT *, "4. Move or copy simCmd.txt to the data/ directory."

PRINT *, " "
PRINT *, "Executing the command to read simCmd.txt:"

CALL READ_COMMAND("simCmd.txt", step, cmd)
CALL PRINT_COMMAND(step, cmd)

END PROGRAM EXTCMD
