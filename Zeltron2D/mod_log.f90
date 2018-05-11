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

MODULE MOD_LOG

IMPLICIT NONE

PRIVATE

INTEGER, PARAMETER :: maxFileNameLen = 256

TYPE DataLogObj
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: logEntries
  CHARACTER(len=maxFileNameLen)                 :: fileName
  INTEGER                                       :: numEntries
  DOUBLE PRECISION                              :: waitTime
  DOUBLE PRECISION                              :: lastWriteTime
ENDTYPE DataLogObj

PUBLIC :: DataLogObj

PUBLIC :: START_LOG
PUBLIC :: END_LOG
PUBLIC :: WRITE_LOG_TO_DISK
PUBLIC :: ADD_LOG_DATA

 CONTAINS

!***********************************************************************
! Subroutine START_LOG
! Allocate data log array
!
! IN/OUT:
! - dataLog: the log 
! - numColumns: the number of columns per entry
! - entriesPerDump: the maximum number of entries to store before writing to disk
! - timeBetweenDumps: the desired time between dumps, in seconds
! 
!***********************************************************************

SUBROUTINE START_LOG(dataLog, fileName, numColumns, entriesPerDump, &
  timeBetweenDumps)

IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(DataLogObj), INTENT(INOUT)           :: dataLog
CHARACTER(len=*)                          :: fileName
INTEGER                                   :: numColumns, entriesPerDump
DOUBLE PRECISION                          :: timeBetweenDumps
!***********************************************************************

ALLOCATE(dataLog%logEntries(numColumns, entriesPerDump))
dataLog%numEntries = 0
dataLog%fileName = " "
dataLog%fileName(1:LEN(fileName)) = fileName
dataLog%waitTime = timeBetweenDumps
dataLog%lastWriteTime = MPI_WTIME()

END SUBROUTINE START_LOG

!***********************************************************************
! Subroutine END_LOG
! Finish writing log, never to add any more data
!
! IN/OUT:
! - dataLog: the log 
! 
!***********************************************************************

SUBROUTINE END_LOG(dataLog)

IMPLICIT NONE
TYPE(DataLogObj), INTENT(INOUT)           :: dataLog
!***********************************************************************

CALL WRITE_LOG_TO_DISK(dataLog)
DEALLOCATE(dataLog%logEntries)

END SUBROUTINE END_LOG

!***********************************************************************
! Subroutine WRITE_LOG_TO_DISK
! Write all log data to disk.
!
! IN/OUT:
! - dataLog: the log 
! 
!***********************************************************************

SUBROUTINE WRITE_LOG_TO_DISK(dataLog)

IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(DataLogObj), INTENT(INOUT)           :: dataLog
! LOCAL
CHARACTER(len=280) :: fName
INTEGER :: i,j, numColumns
!***********************************************************************

! Open the file only if there's some data to be written
IF (dataLog%numEntries > 0) THEN
  fname = trim(adjustl(dataLog%fileName))
  OPEN(9,FILE=fname,POSITION='APPEND')
  numColumns = SIZE(dataLog%logEntries, DIM = 1)
  DO j=1,dataLog%numEntries
    WRITE(9,*) (dataLog%logEntries(i,j), i=1, numColumns)
  ENDDO
  CLOSE(9)
  dataLog%numEntries = 0
  dataLog%lastWriteTime = MPI_WTIME()
ENDIF

END SUBROUTINE WRITE_LOG_TO_DISK

!***********************************************************************
! Subroutine ADD_LOG_DATA
! Add a line of data to the log.
!
! INPUT:
! - newLine: a new line of data to be logged
! IN/OUT:
! - dataLog: the to which data should be added
! 
!***********************************************************************

SUBROUTINE ADD_LOG_DATA(dataLog, newLine)

IMPLICIT NONE

INCLUDE 'mpif.h'

TYPE(DataLogObj), INTENT(INOUT)           :: dataLog
DOUBLE PRECISION, DIMENSION(:)            :: newLine
!***********************************************************************

dataLog%numEntries = dataLog%numEntries + 1
dataLog%logEntries(:, dataLog%numEntries) = newLine
IF ((dataLog%numEntries == SIZE(dataLog%logEntries, DIM = 2))  &
    .OR. (MPI_WTIME() - dataLog%lastWriteTime >= dataLog%waitTime)) THEN
  !print *, "Logger", dataLog%fileName, dataLog%numEntries, SIZE(dataLog%logEntries, DIM=2)
  !print *, "      ", MPI_WTIME() - dataLog%lastWriteTime, dataLog%waitTime
  CALL WRITE_LOG_TO_DISK(dataLog)
ENDIF

END SUBROUTINE ADD_LOG_DATA

!***********************************************************************

END MODULE MOD_LOG
