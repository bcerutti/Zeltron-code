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

MODULE MOD_EXTCMD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module reads the contents of a (command) file.  The file should
! have a single line (subsequent lines will be ignored) with 4 integers
! separated by a space:
!   timeStep quit dump save
! where timeStep is an integer timeStep at which the command should be
! executed, and quit, dump, and save are either 0 or 1, depending
! on whether those actions should be performed at the time step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE

PRIVATE

INTEGER, PARAMETER, PUBLIC :: EXTCMD_DO_NOTHING = 0
INTEGER, PARAMETER, PUBLIC :: EXTCMD_QUIT = 1
INTEGER, PARAMETER, PUBLIC :: EXTCMD_DUMP = 2
INTEGER, PARAMETER, PUBLIC :: EXTCMD_SAVE = 4

PUBLIC :: READ_COMMAND
PUBLIC :: PRINT_COMMAND

 CONTAINS

FUNCTION BoolToInt(b)

IMPLICIT NONE

LOGICAL :: b
INTEGER :: BoolToInt

IF (b) THEN
  BoolToInt = 1
ELSE
  BoolToInt = 0
ENDIF

END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine READ_COMMAND
!   Reads a command file, and returns a timeStep and a command.
!
! INPUT:
! - fileName the name of the command file
! OUTPUT:
! - timeStep: an integer time step
! - cmd: an integer command, the sum of the EXTCMD_* flags indicating
!        the actions to be performed.
!
! If the command file does not exist, or there are any errors reading
! the command file, this returns timeStep = -1 
! and cmd = 0 = EXTCMD_DO_NOTHING.
!
! Passing a command through a file is very klugey approach; we want
! this subroutine never to fail, never to result in a crash.
! By using the IOSTAT option for file operations, we try to avoid
! failure, and in such cases return innocuous EXTCMD_DO_NOTHING.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE READ_COMMAND(fileName, timeStep, cmd)

IMPLICIT NONE

! INPUT
CHARACTER(len=*)     :: fileName
INTEGER, INTENT(OUT) :: timeStep, cmd

! LOCAL
LOGICAL              :: cmdFileExists, fail
INTEGER              :: step
INTEGER              :: doNothing, doQuit, doDump, doSave
INTEGER              :: ioStatus
!***********************************************************************

! check for command file
fail = .TRUE.
INQUIRE(FILE=fileName, EXIST=cmdFileExists, IOSTAT = ioStatus)
IF (ioStatus == 0) THEN
  IF (cmdFileExists) THEN
    OPEN(10, FILE = fileName, IOSTAT = ioStatus)
    IF (ioStatus == 0) THEN
      READ(10, *, IOSTAT = ioStatus) step, doQuit, doDump, doSave
      IF (ioStatus == 0) THEN
        timeStep = step
        cmd = EXTCMD_QUIT*BoolToInt(doQuit>0) + EXTCMD_DUMP*BoolToInt(doDump > 0) &
          + EXTCMD_SAVE*BoolToInt(doSave > 0)
        fail = .FALSE.
      ELSE
        PRINT *, "Error: IOSTAT=", ioStatus, " while reading command file ", fileName
      ENDIF
      CLOSE(10, IOSTAT = ioStatus)
    ELSE
      PRINT *, "Error: IOSTAT=", ioStatus, " while opening command file ", fileName
    ENDIF
  ENDIF
ELSE
  PRINT *, "Error: IOSTAT=", ioStatus, " while inquiring about existence of file ", fileName
ENDIF

IF (fail) THEN
  timeStep = -1
  cmd = EXTCMD_DO_NOTHING
ENDIF

END SUBROUTINE READ_COMMAND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine PRINT_COMMAND
!   Prints the interpretation of a command from a command file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PRINT_COMMAND(timeStep, cmd)

IMPLICIT NONE

INTEGER              :: timeStep, cmd
!***********************************************************************

PRINT *, "At timestep", timeStep, ":"
IF (cmd == 0) THEN
  PRINT *, "   do nothing"
ELSE 
  IF (IAND(cmd, EXTCMD_QUIT) /= 0) PRINT *, "   quit" 
  IF (IAND(cmd, EXTCMD_DUMP) /= 0) PRINT *, "   dump" 
  IF (IAND(cmd, EXTCMD_SAVE) /= 0) PRINT *, "   save" 
ENDIF

END SUBROUTINE PRINT_COMMAND

!***********************************************************************

END MODULE MOD_EXTCMD
