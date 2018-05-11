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

MODULE MOD_IO

USE HDF5
USE MOD_INPUT
USE MOD_MPI
USE MOD_INTERP

IMPLICIT NONE

LOGICAL, PARAMETER, PUBLIC :: USE_HDF5 = .TRUE.

PRIVATE

PUBLIC :: DUMP_FIELD ! Write a field to disk collectively
PUBLIC :: SAVE_FIELDS ! Write fields to disk for later restoration
PUBLIC :: RESTORE_FIELDS ! Read fields from disk collectively
PUBLIC :: DUMP_GLOBAL_FIELD_1D ! Write a field on one rank to disk 
PUBLIC :: DUMP_GLOBAL_FIELD_3D ! Write a field on one rank to disk 
PUBLIC :: DUMP_PARTICLES ! Write particle data to disk
PUBLIC :: SAVE_PARTICLES ! Write particle data to disk for later restoration
PUBLIC :: RESTORE_PARTICLES ! Read particle data from disk
PUBLIC :: DUMP_TRACKED_PARTICLES_SER ! write table of tracked particles to disk
PUBLIC :: GET_TRACK_SEQ_NUM ! get the sequence number for saving tracking data
PUBLIC :: DUMP_SIM_INFO ! write helpful info

!=======================================================================
! Naming:
! READ/WRITE/APPEND access (datasets/attributes) of already-opened files.
! DUMP/SAVE/LOAD create or open files (and then write data).
!
! APPEND writes to extendible datasets.
!
! DATASET indicates that on rank will read/write the entire dataset.
! (In contrast, FIELD and PARTICLES indicate collective data dumping unless
! _SER is present.)  Even though only one rank may ultimately do the
! reading/writing, all ranks need to call the method; this is for cases
! where all ranks opened a file (to perform a collective action), but
! subsequently a single rank wants to write a small dataset.
! 
! _SER means that only one rank should call the function, which means
! that only one rank will open or has opened the relevant file.
!=======================================================================

 CONTAINS

!***********************************************************************
! Subroutine WRITE_STRING_ATTRIB_SER
! This subroutine writes a string attribute to a dataset
! opened by a single rank (this should be called by one rank)
!
! INPUT:
! 
! - dsetId the hdf5 handle to an already-open dataset 
! - attribName the name of the attribute
! - attribStr the string to store in the attribute
! - COMM: mpi comm
!
!***********************************************************************

SUBROUTINE WRITE_STR_ATTRIB_SER(dsetId, attribName, attribStr)

IMPLICIT NONE

! INPUT
INTEGER(HID_T) :: dsetId
CHARACTER(len=*) :: attribName, attribStr
! LOCAL
INTEGER(HSIZE_T) :: strShape(1)
INTEGER(HID_T)   :: attribId, attrSpace
INTEGER          :: h5err
!***********************************************************************

strShape = (/ LEN(attribStr) /)
CALL h5screate_simple_f(1, strShape, attrSpace, h5err)
CALL h5acreate_f(dsetId, attribName, H5T_NATIVE_CHARACTER, attrSpace, &
  attribId, h5err)
CALL h5awrite_f(attribId, H5T_NATIVE_CHARACTER, attribStr, strShape, h5err)
CALL h5sclose_f(attrSpace, h5err)
CALL h5aclose_f(attribId, h5err)

END SUBROUTINE WRITE_STR_ATTRIB_SER

!***********************************************************************
! Subroutine WRITE_DATASET_1DDBL_SER
! This subroutine writes a 1D double-precision dataset array to an 
! already open h5 file (that has been opened by one rank; this function
! is called only by that rank).
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - ra the full 1D array to be written
! - optional description
!
!***********************************************************************

SUBROUTINE WRITE_DATASET_1DDBL_SER(fileId, datasetName, rankToWrite, ra, &
  desc)

IMPLICIT NONE

INTEGER, PARAMETER :: FDIM = 1

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
DOUBLE PRECISION :: ra(:)
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HSIZE_T) :: raShape(FDIM)
INTEGER(HID_T)   :: filespace, dsetId, attribId, attrSpace
INTEGER          :: h5err
!***********************************************************************

raShape = SHAPE(ra)
CALL h5screate_simple_f(FDIM, raShape, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_DOUBLE, filespace, &
  dsetId, h5err)
CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, raShape, h5err)
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB_SER(dsetId, "description", desc)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)

END SUBROUTINE WRITE_DATASET_1DDBL_SER

!***********************************************************************
! Subroutine WRITE_DATASET_1DINT_SER
! This subroutine writes a 1D integer dataset array to an 
! already open h5 file (opened by a single mpi rank; only that rank
! calls this subroutine).
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - ra the full 1D array to be written
! - optional description
!
!***********************************************************************
SUBROUTINE WRITE_DATASET_1DINT_SER(fileId, datasetName, rankToWrite, ra, &
  desc)

IMPLICIT NONE

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
INTEGER          :: ra(:)
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HSIZE_T) :: raLength(1)
INTEGER(HID_T)   :: filespace, dsetId
INTEGER          :: h5err
!***********************************************************************

raLength = SIZE(ra)
CALL h5screate_simple_f(1, raLength, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_INTEGER, filespace, &
  dsetId, h5err)
CALL h5dwrite_f(dsetId, H5T_NATIVE_INTEGER, ra, raLength, h5err)
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB_SER(dsetId, "description", desc)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)

END SUBROUTINE WRITE_DATASET_1DINT_SER

!***********************************************************************
! Subroutine WRITE_STRING_ATTRIB
! This subroutine writes a string attribute to a dataset
! (called by all ranks)
!
! INPUT:
! 
! - dsetId the hdf5 handle to an already-open dataset 
! - attribName the name of the attribute
! - attribStr the string to store in the attribute
! - COMM: mpi comm
!
!***********************************************************************
SUBROUTINE WRITE_STR_ATTRIB(dsetId, attribName, attribStr, COMM)

IMPLICIT NONE

! INPUT
INTEGER(HID_T) :: dsetId
CHARACTER(len=*) :: attribName, attribStr
INTEGER          :: COMM
! LOCAL
INTEGER(HSIZE_T) :: strShape(1)
INTEGER(HID_T)   :: attribId, attrSpace
INTEGER          :: h5err
!***********************************************************************

strShape = (/ LEN(attribStr) /)
CALL h5screate_simple_f(1, strShape, attrSpace, h5err)
CALL h5acreate_f(dsetId, attribName, H5T_NATIVE_CHARACTER, attrSpace, &
  attribId, h5err)
CALL h5awrite_f(attribId, H5T_NATIVE_CHARACTER, attribStr, strShape, h5err)
CALL h5sclose_f(attrSpace, h5err)
CALL h5aclose_f(attribId, h5err)

END SUBROUTINE WRITE_STR_ATTRIB

!***********************************************************************
! Subroutine WRITE_DATASET_STR_LIST
! This subroutine writes a 1D array of same-size strings to an
! already open h5 file.  One mpi rank does the writing, but this
! subroutine must be called by all mpi ranks (until I figure out
! otherwise), for a file that has been opened by all ranks.
!
! While ra will not be accessed except by the rankToWrite,
! raShape needs to be correct on all ranks.
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - raShape = (stringLength, numStrings)
! - ra a 1D character(len=stringLength*numStrings) string
! - COMM: mpi comm
! - optional description
!
!***********************************************************************
SUBROUTINE WRITE_DATASET_STR_LIST(fileId, datasetName, rankToWrite, raShape, &
  ra, COMM, desc)

IMPLICIT NONE

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
INTEGER(HSIZE_T) :: raShape(2)
CHARACTER        :: ra(:)
INTEGER          :: COMM
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HID_T)   :: plistId, filespace, dsetId
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
CALL h5screate_simple_f(2, raShape, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_CHARACTER, filespace, &
  dsetId, h5err)
IF (rank == rankToWrite) THEN
  CALL h5dwrite_f(dsetId, H5T_NATIVE_CHARACTER, ra, raShape, h5err, &
    xfer_prp=plistId)
ENDIF
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB(dsetId, "description", desc, COMM)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)

END SUBROUTINE WRITE_DATASET_STR_LIST

!***********************************************************************
! Subroutine WRITE_DATASET_1DDBL
! This subroutine writes a 1D double-precision dataset array to an 
! already open h5 file.  One mpi rank does the writing, but this
! subroutine must be called by all mpi ranks (until I figure out
! otherwise), for a file that has been opened by all ranks.
!
! While ra will not be accessed except by the rankToWrite,
! raShape needs to be correct for all ranks.
!
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - raShape = SHAPE(ra)
! - ra the full 1D array to be written
! - COMM: mpi comm
! - optional description
!
!***********************************************************************

SUBROUTINE WRITE_DATASET_1DDBL(fileId, datasetName, rankToWrite, raShape,&
  ra, COMM, desc)

IMPLICIT NONE

INTEGER, PARAMETER :: FDIM = 1

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
INTEGER(HSIZE_T) :: raShape(FDIM)
DOUBLE PRECISION :: ra(:)
INTEGER          :: COMM
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HSIZE_T) :: descShape(1)
INTEGER(HID_T)   :: plistId, filespace, dsetId, attribId, attrSpace
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
CALL h5screate_simple_f(FDIM, raShape, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_DOUBLE, filespace, &
  dsetId, h5err)
IF (rank == rankToWrite) THEN
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, raShape, h5err, &
    xfer_prp=plistId)
ENDIF
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB(dsetId, "description", desc, COMM)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)

END SUBROUTINE WRITE_DATASET_1DDBL

!***********************************************************************
! Subroutine WRITE_DATASET_1DINT
! This subroutine writes a 1D integer dataset array to an 
! already open h5 file.  One mpi rank does the writing, but this
! subroutine must be called by all mpi ranks (until I figure out
! otherwise), for a file that has been opened by all ranks.
!
! While ra will not be accessed except by the rankToWrite,
! raShape needs to be correct on all ranks.
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - raShape = SHAPE(ra)
! - ra the full 1D array to be written
! - COMM: mpi comm
! - optional description
!
!***********************************************************************

SUBROUTINE WRITE_DATASET_1DINT(fileId, datasetName, rankToWrite, raShape, &
  ra, COMM, desc)

IMPLICIT NONE

INTEGER, PARAMETER :: FDIM = 1

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
INTEGER(HSIZE_T) :: raShape(FDIM)
INTEGER          :: ra(:)
INTEGER          :: COMM
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HID_T)   :: plistId, filespace, dsetId
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
CALL h5screate_simple_f(FDIM, raShape, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_INTEGER, filespace, &
  dsetId, h5err)
IF (rank == rankToWrite) THEN
  CALL h5dwrite_f(dsetId, H5T_NATIVE_INTEGER, ra, raShape, h5err, &
    xfer_prp=plistId)
ENDIF
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB(dsetId, "description", desc, COMM)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)

END SUBROUTINE WRITE_DATASET_1DINT

!***********************************************************************
! Subroutine WRITE_DATASET_2DDBL
! This subroutine writes a 2D double-precision dataset array to an 
! already open h5 file.  One mpi rank does the writing, but this
! subroutine must be called by all mpi ranks (until I figure out
! otherwise), for a file that has been opened by all ranks.
!
! While ra will not be accessed except by the rankToWrite,
! raShape needs to be correct for all ranks.
!
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - raShape = SHAPE(ra)
! - ra the full 2D array to be written
! - COMM: mpi comm
! - optional description
!
!***********************************************************************

SUBROUTINE WRITE_DATASET_2DDBL(fileId, datasetName, rankToWrite, raShape,&
  ra, COMM, desc)

IMPLICIT NONE

INTEGER, PARAMETER :: FDIM = 2

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
INTEGER(HSIZE_T) :: raShape(FDIM)
DOUBLE PRECISION :: ra(:,:)
INTEGER          :: COMM
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HID_T)   :: plistId, filespace, dsetId, attribId, attrSpace
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
CALL h5screate_simple_f(FDIM, raShape, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_DOUBLE, filespace, &
  dsetId, h5err)
IF (rank == rankToWrite) THEN
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, raShape, h5err, &
    xfer_prp=plistId)
ENDIF
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB(dsetId, "description", desc, COMM)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)

END SUBROUTINE WRITE_DATASET_2DDBL

!***********************************************************************
! Subroutine WRITE_DATASET_2DINT
! This subroutine writes a 2D integer dataset array to an 
! already open h5 file.  One mpi rank does the writing, but this
! subroutine must be called by all mpi ranks (until I figure out
! otherwise), for a file that has been opened by all ranks.
!
! While ra will not be accessed except by the rankToWrite,
! raShape needs to be correct on all ranks.
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - raShape = SHAPE(ra)
! - ra the full 2D array to be written
! - COMM: mpi comm
! - optional description
!
!***********************************************************************

SUBROUTINE WRITE_DATASET_2DINT(fileId, datasetName, rankToWrite, raShape, &
  ra, COMM, desc)

IMPLICIT NONE

INTEGER, PARAMETER :: FDIM = 2

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
INTEGER(HSIZE_T) :: raShape(FDIM)
INTEGER          :: ra(:,:)
INTEGER          :: COMM
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HID_T)   :: plistId, filespace, dsetId
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
CALL h5screate_simple_f(FDIM, raShape, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_INTEGER, filespace, &
  dsetId, h5err)
IF (rank == rankToWrite) THEN
  CALL h5dwrite_f(dsetId, H5T_NATIVE_INTEGER, ra, raShape, h5err, &
    xfer_prp=plistId)
ENDIF
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB(dsetId, "description", desc, COMM)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)

END SUBROUTINE WRITE_DATASET_2DINT

!***********************************************************************
! Subroutine WRITE_DATASET_1DINT8
! This subroutine writes a 1D integer dataset array to an 
! already open h5 file.  One mpi rank does the writing, but this
! subroutine must be called by all mpi ranks (until I figure out
! otherwise), for a file that has been opened by all ranks.
!
! Currently writes an array of doubles.
!
! While ra will not be accessed except by the rankToWrite,
! raShape needs to be correct on all ranks.
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - raShape = SHAPE(ra)
! - ra the full 1D array to be written
! - COMM: mpi comm
! - optional description
!
!***********************************************************************

SUBROUTINE WRITE_DATASET_1DINT8(fileId, datasetName, rankToWrite, raShape,&
  ra, COMM, desc)

IMPLICIT NONE

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
INTEGER(HSIZE_T) :: raShape(1)
INTEGER*8        :: ra(:)
INTEGER          :: COMM
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HSIZE_T) :: raLength(1)
INTEGER(HID_T)   :: plistId, filespace, dsetId
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

raLength = raShape(1)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
CALL h5screate_simple_f(1, raLength, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, &
  H5T_NATIVE_DOUBLE, filespace, &
  dsetId, h5err)
IF (rank == rankToWrite) THEN
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, &
    DBLE(ra), raLength, h5err, xfer_prp=plistId)
ENDIF
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB(dsetId, "description", desc, COMM)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)

END SUBROUTINE WRITE_DATASET_1DINT8

!***********************************************************************
! Subroutine READ_DATASET_1DINT8
! This subroutine reads a 1D integer dataset array to an 
! already open h5 file.  One mpi rank does the reading, and communicates
! the data to all other ranks.
!
! Currently reads an array of doubles.
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToRead the (0-indexed) rank that will actually read the dataset
! - ra the full 1D array to be read (will be filled on every rank)
!   -- will be allocated to correct size
! - COMM: mpi comm
!
!***********************************************************************

SUBROUTINE READ_DATASET_1DINT8(fileId, datasetName, rankToRead, &
  ra, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
INTEGER(HID_T)                      :: fileId ! an already open file
CHARACTER(len=*)                    :: datasetName
INTEGER                             :: rankToRead
INTEGER*8, ALLOCATABLE, INTENT(OUT) :: ra(:)
INTEGER                             :: COMM
! LOCAL
DOUBLE PRECISION, ALLOCATABLE :: dra(:)
INTEGER(HSIZE_T) :: raShape(1), maxShape(1)
INTEGER(HID_T)   :: plistId, filespace, dsetId
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5dopen_f(fileId, datasetName, dsetId, h5err)
IF (rank == rankToRead) THEN
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
  CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
  CALL h5dget_space_f(dsetId, filespace, h5err)
  CALL h5sget_simple_extent_dims_f(filespace, raShape, maxShape, h5err)
ENDIF
CALL MPI_BCAST(raShape, 1, MPI_INTEGER8, rankToRead, COMM, mpiErr)

ALLOCATE(ra(raShape(1)))

IF (rank == rankToRead) THEN
  ALLOCATE(dra(raShape(1)))
  CALL h5dread_f(dsetId, H5T_NATIVE_DOUBLE, &
    dra, raShape, h5err, xfer_prp=plistId)
  CALL h5sclose_f(filespace, h5err)
  CALL h5pclose_f(plistId, h5err)
  ra = dra
  DEALLOCATE(dra)
ENDIF
CALL MPI_BCAST(ra, raShape(1), MPI_INTEGER8, rankToRead, COMM, mpiErr)
CALL h5dclose_f(dsetId, h5err)

END SUBROUTINE READ_DATASET_1DINT8

!***********************************************************************
! Subroutine WRITE_DATASET_3DDBL
! This subroutine writes a 3D double-precision dataset array to an 
! already open h5 file.  One mpi rank does the writing, but this
! subroutine must be called by all mpi ranks (until I figure out
! otherwise) for a file that has been opened by all ranks.
!
! While ra will not be accessed except by the rankToWrite,
! raShape needs to be correct for all ranks.
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToWrite the (0-indexed) rank that will actually write the dataset
! - raShape = SHAPE(ra)
! - ra the full 1D array to be written
! - COMM: mpi comm
! - optional description
!
!***********************************************************************

SUBROUTINE WRITE_DATASET_3DDBL(fileId, datasetName, rankToWrite, raShape,&
  ra, COMM, desc)

IMPLICIT NONE

INTEGER, PARAMETER :: FDIM = 3

! INPUT
INTEGER(HID_T)   :: fileId ! an already open file
CHARACTER(len=*) :: datasetName
INTEGER          :: rankToWrite
INTEGER(HSIZE_T) :: raShape(FDIM)
DOUBLE PRECISION :: ra(:,:,:)
INTEGER          :: COMM
CHARACTER(len=*), OPTIONAL :: desc
! LOCAL
INTEGER(HID_T)   :: plistId, filespace, dsetId
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
CALL h5screate_simple_f(FDIM, raShape, filespace, h5err)
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_DOUBLE, filespace, &
  dsetId, h5err)
IF (rank == rankToWrite) THEN
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, raShape, h5err, &
    xfer_prp=plistId)
ENDIF
IF (PRESENT(desc)) THEN
  CALL WRITE_STR_ATTRIB(dsetId, "description", desc, COMM)
ENDIF
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)

END SUBROUTINE WRITE_DATASET_3DDBL

!***********************************************************************
! Subroutine READ_DATASET_3DDBL
! This subroutine reads a 3D double-precision dataset array from an
! already open h5 file.  One mpi rank does the reading, and 
! communicates the data to all other ranks.
!
! While ra will not be accessed except by the rankToWrite,
! raShape needs to be correct for all ranks.
!
! INPUT:
! 
! - fileId the hdf5 handle to the already-open h5 file
! - datasetName the name of the field dataset
! - rankToRead the (0-indexed) rank that will actually read the dataset
! - ra the full 1D array to be written
!   -- will be allocated to correct size
! - COMM: mpi comm
!
!***********************************************************************

SUBROUTINE READ_DATASET_3DDBL(fileId, datasetName, rankToRead, ra, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER :: FDIM = 3

! INPUT
INTEGER(HID_T)                             :: fileId ! an already open file
CHARACTER(len=*)                           :: datasetName
INTEGER                                    :: rankToRead
DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT) :: ra(:,:,:)
INTEGER                                    :: COMM
! LOCAL
INTEGER(HSIZE_T) :: raShape(FDIM), maxShape(FDIM)
INTEGER(HID_T)   :: plistId, filespace, dsetId
INTEGER          :: rank, h5err, mpiErr
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL h5dopen_f(fileId, datasetName, dsetId, h5err)
IF (rank == rankToRead) THEN
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
  CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_INDEPENDENT_F, h5err)
  CALL h5dget_space_f(dsetId, filespace, h5err)
  CALL h5sget_simple_extent_dims_f(filespace, raShape, maxShape, h5err)
ENDIF
CALL MPI_BCAST(raShape, FDIM, MPI_INTEGER8, rankToRead, COMM, mpiErr)
ALLOCATE(ra(raShape(1), raShape(2), raShape(3)))
IF (rank == rankToRead) THEN
  CALL h5dread_f(dsetId, H5T_NATIVE_DOUBLE, &
    ra, raShape, h5err, xfer_prp=plistId)
  CALL h5sclose_f(filespace, h5err)
  CALL h5pclose_f(plistId, h5err)
ENDIF
CALL MPI_BCAST(ra, raShape, MPI_DOUBLE_PRECISION, rankToRead, COMM, mpiErr)
CALL h5dclose_f(dsetId, h5err)

END SUBROUTINE READ_DATASET_3DDBL

!***********************************************************************
! Subroutine WRITE_TABLE_DBL
! This subroutine writes a table dataset, i.e., a 2D array of shape
! (fixed, arbitrary) to an already-open (in parallel) h5 file, 
! where each processor contributes
! a section of the table (each processer can have a different number
! of entries in the 2nd dimension).
! Each rank must call this.
!
! INPUT:
! 
! - fileId the handle of the hdf5 file (that's already open)
! - datasetName the name of the field dataset
! - table the 2D array (for the local domain)
! - writeRankInfo whether to write info about the number of entries per rank
! - addBoundingBoxes: whether each the bounds of the first NDIM table
!   entries (for each rank) should be added in a separate table
! - COMM: mpi comm
! - it the timestep
! - t the time
!
!***********************************************************************

SUBROUTINE WRITE_TABLE_DBL(fileId, datasetName, table, writeRankInfo, &
  addBoundingBoxes, COMM, it, t)

IMPLICIT NONE

INCLUDE 'mpif.h'

!INPUT
INTEGER(HID_T)            :: fileId
CHARACTER(len=*)          :: datasetName
DOUBLE PRECISION          :: table(:,:)
LOGICAL                   :: writeRankInfo, addBoundingBoxes
INTEGER                   :: COMM
INTEGER, OPTIONAL         :: it
DOUBLE PRECISION, OPTIONAL:: t

! LOCAL
INTEGER                   :: numRanks, rank, mpiErr
! bounds = [[xLower, yLower, zLower], [xUpper, yUpper, zUpper]], 
!   is the bounds of the first NDIM columns of the table
DOUBLE PRECISION          :: localBounds(2,NDIM)
DOUBLE PRECISION, ALLOCATABLE :: allBounds(:,:,:)
INTEGER(HSIZE_T)          :: allBoundsShape(3)
INTEGER*8, ALLOCATABLE    :: entriesPerRank(:)
INTEGER*8, ALLOCATABLE    :: rankOffset(:)
INTEGER                   :: d
INTEGER*8                 :: numLocalEntries

INTEGER(HSIZE_T)          :: localShape(2), globalShape(2), localOffset(2)
INTEGER(HSIZE_T)          :: bLocalShape(3), bGlobalShape(3), bLocalOffset(3)

INTEGER(HID_T) :: dsetId, plistId2
INTEGER(HID_T) :: filespace, subspace, memspace
INTEGER        :: h5err
INTEGER        :: serialWritingRank = 0
!***********************************************************************

localShape = SHAPE(table)
globalShape(1) = SIZE(table, DIM=1)
numLocalEntries = SIZE(table, DIM=2)

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
ALLOCATE(entriesPerRank(numRanks))
ALLOCATE(rankOffset(numRanks+1))

CALL MPI_ALLGATHER(numLocalEntries, 1, MPI_INTEGER8, entriesPerRank, &
  1, MPI_INTEGER8, COMM, mpiErr)

CALL CUMSUM(entriesPerRank, rankOffset(1:numRanks), rankOffset(numRanks+1))

globalShape(2) = rankOffset(numRanks+1)
localOffset(1) = 0
localOffset(2) = rankOffset(rank+1)
  
IF (writePerRankFiles) THEN
  OPEN(9, FILE=trim(perRankFile), POSITION='APPEND')
  WRITE(9, *) "  entriesPerRank =", entriesPerRank
  WRITE(9, *) "  rankOffset =", rankOffset
  WRITE(9, *) "  localShape =", localShape
  WRITE(9, *) "  globalShape =", globalShape
  WRITE(9, *) "  localOffset =", localOffset
  CLOSE(9)
ENDIF

!{
! create data space for dataset
CALL h5screate_simple_f(2, globalShape, filespace, h5err)
CALL h5screate_simple_f(2, localShape, memspace, h5err)
! create dataset with default properties
CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_DOUBLE, filespace, &
  dsetId, h5err)
! select the part that will be written by this domain
CALL h5dget_space_f(dsetId, subspace, h5err)
CALL h5sselect_hyperslab_f(subspace, H5S_SELECT_SET_F, localOffset, &
  localShape, h5err)
! create property list for collective dataset write
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId2, h5err)
CALL h5pset_dxpl_mpio_f(plistId2, H5FD_MPIO_COLLECTIVE_F, h5err)

! write table collectively
CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, table, localShape, h5err, &
  file_space_id = subspace, mem_space_id = memspace, xfer_prp=plistId2)

CALL h5sclose_f(memspace, h5err)
CALL h5sclose_f(subspace, h5err)

CALL h5pclose_f(plistId2, h5err)
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
!}

IF (addBoundingBoxes) THEN !{
  ! figure out the the bounds for all entries on this rank
  ! (this is for restarting with a different number of processors)
  DO d = 1, NDIM
    IF (localShape(2) > 0) THEN
      localBounds(1, d) = MINVAL(table(d,:))
      localBounds(2, d) = MAXVAL(table(d,:))
    ELSE
      localBounds(1, d) = 0.
      localBounds(2, d) = -1.
    ENDIF
  ENDDO

  bLocalShape(1) = 2
  bLocalShape(2) = NDIM
  bLocalShape(3) = 1
  bGlobalShape(1) = 2
  bGlobalShape(2) = NDIM
  bGlobalShape(3) = numRanks
  bLocalOffset(1) = 0
  bLocalOffset(2) = 0
  bLocalOffset(3) = rank
  
! Write collectively
!  ! create data space for dataset
!  CALL h5screate_simple_f(NDIM+1, bGlobalShape, filespace, h5err)
!  CALL h5screate_simple_f(NDIM+1, bLocalShape, memspace, h5err)
!  ! create dataset with default properties
!  CALL h5dcreate_f(fileId, datasetName // "BoundingBoxesPerRank", H5T_NATIVE_DOUBLE, filespace, &
!    dsetId, h5err)
!  ! select the part that will be written by this domain
!  CALL h5dget_space_f(dsetId, subspace, h5err)
!  CALL h5sselect_hyperslab_f(subspace, H5S_SELECT_SET_F, bLocalOffset, &
!    bLocalShape, h5err)
!  ! create property list for collective dataset write
!  CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId2, h5err)
!  CALL h5pset_dxpl_mpio_f(plistId2, H5FD_MPIO_COLLECTIVE_F, h5err)
!
!  ! write table collectively
!  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, localBounds, bLocalShape, &
!    h5err, file_space_id = subspace, mem_space_id = memspace, xfer_prp=plistId2)
!
!  CALL h5pclose_f(plistId2, h5err)
!  CALL h5dclose_f(dsetId, h5err)
!  CALL h5sclose_f(filespace, h5err)
!  CALL h5sclose_f(subspace, h5err)
!  CALL h5sclose_f(memspace, h5err)

! Message and write with one rank
  allBoundsShape(1) = 2
  allBoundsShape(2) = NDIM
  allBoundsShape(3) = numRanks
  IF (rank == serialWritingRank) THEN
    ALLOCATE(allBounds(2,NDIM,numRanks))
  ENDIF
  CALL MPI_GATHER(localBounds, 2*NDIM, MPI_DOUBLE_PRECISION, &
    allBounds, 2*NDIM, MPI_DOUBLE_PRECISION, serialWritingRank, &
    COMM, mpiErr)
  CALL WRITE_DATASET_3DDBL(fileId, datasetName // "BoundingBoxesPerRank", &
    serialWritingRank, allBoundsShape, allBounds, COMM)
  IF (rank == serialWritingRank) THEN
    DEALLOCATE(allBounds)
  ENDIF

ENDIF !}

IF (writeRankInfo) THEN
  ! write the number of entries for each rank
  CALL WRITE_DATASET_1DINT8(fileId, datasetName // 'EntriesPerRank', serialWritingRank, &
    (/ INT8(numRanks) /), entriesPerRank, COMM)
  ! write the offset for each rank
  CALL WRITE_DATASET_1DINT8(fileId, datasetName // 'StartingEntryPerRank', serialWritingRank,&
    (/ INT8(numRanks+1) /), rankOffset, COMM)
ENDIF

! write step and time in serial
IF (PRESENT(it)) THEN
  CALL WRITE_DATASET_1DINT(fileId, 'step', serialWritingRank, (/ INT8(1) /), &
    (/ it /), COMM)
ENDIF
IF (PRESENT(t)) THEN
  CALL WRITE_DATASET_1DDBL(fileId, 'time', serialWritingRank, (/ INT8(1) /), &
    (/ t /), COMM)
ENDIF

DEALLOCATE(rankOffset)
DEALLOCATE(entriesPerRank)

END SUBROUTINE WRITE_TABLE_DBL

!***********************************************************************
! Subroutine READ_TABLE_DBL
! This subroutine reads a table dataset, i.e., a 2D array of shape
! (fixed, arbitrary) from an already-open (in parallel) h5 file, 
! where each processor contributes
! a section of the table (each processer can have a different number
! of entries in the 2nd dimension).
! Each rank must call this.
!
! INPUT:
! 
! - fileId the handle of the hdf5 file (that's already open)
! - datasetName the name of the field dataset
! - table the 2D array (for the local domain)
! - COMM: mpi comm
! - entriesDataset: the name of the dataset specificying how many entries
!     there are from each rank; by default this is
!     datasetName // "EntriesPerRank"
!
!***********************************************************************

SUBROUTINE READ_TABLE_DBL(fileId, datasetName, table, COMM, &
  entriesDataset)

IMPLICIT NONE

INCLUDE 'mpif.h'

!INPUT
INTEGER(HID_T)                              :: fileId
CHARACTER(len=*)                            :: datasetName
DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT)  :: table(:,:)
INTEGER                                     :: COMM
CHARACTER(len=*), OPTIONAL                  :: entriesDataset

! LOCAL
INTEGER                   :: numRanks, rank, mpiErr
INTEGER                   :: numFileRanks
! bounds = [[xLower, yLower, zLower], [xUpper, yUpper, zUpper]], 
!   is the bounds of the first NDIM columns of the table
INTEGER*8, ALLOCATABLE    :: entriesPerRank(:)
INTEGER*8, ALLOCATABLE    :: rankOffset(:)

INTEGER(HSIZE_T)          :: localShape(2), globalShape(2), localOffset(2)
INTEGER(HSIZE_T)          :: maxShape(2)

INTEGER(HID_T) :: dsetId, plistId
INTEGER(HID_T) :: filespace, subspace, memspace
INTEGER        :: h5err
INTEGER        :: serialReadingRank = 0
CHARACTER(len=256) :: edName
!***********************************************************************

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)

IF (PRESENT(entriesDataset)) THEN
  edName = entriesDataset
ELSE
  edName = datasetName // 'EntriesPerRank'
ENDIF

! read the number of entries for each rank
CALL READ_DATASET_1DINT8(fileId, edName, &
  serialReadingRank, entriesPerRank, COMM)

numFileRanks = SIZE(entriesPerRank)
ALLOCATE(rankOffset(numFileRanks+1))
CALL CUMSUM(entriesPerRank, rankOffset(1:numFileRanks), &
  rankOffset(numFileRanks+1))

IF (numFileRanks /= numRanks) THEN
  IF (rank == 0) THEN
    PRINT *, "Error: cannot restore from run on ", numFileRanks, " ranks"
    PRINT *, "   to run on ", numRanks, " ranks"
  ENDIF
  CALL MPI_BARRIER(COMM, mpiErr)
  CALL MPI_ABORT(COMM, 9123, mpiErr)
  STOP 9123
ENDIF

CALL h5dopen_f(fileId, datasetName, dsetId, h5err)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_COLLECTIVE_F, h5err)
CALL h5dget_space_f(dsetId, filespace, h5err)
CALL h5sget_simple_extent_dims_f(filespace, globalShape, maxShape, h5err)

localShape(1) = globalShape(1)
localShape(2) = entriesPerRank(rank+1)
localOffset(1) = 0
localOffset(2) = rankOffset(rank+1)

ALLOCATE(table(localShape(1), localShape(2)))

CALL h5screate_simple_f(2, localShape, memspace, h5err)
! select the part that will be written by this domain
CALL h5dget_space_f(dsetId, subspace, h5err)
CALL h5sselect_hyperslab_f(subspace, H5S_SELECT_SET_F, localOffset, &
  localShape, h5err)

CALL h5dread_f(dsetId, H5T_NATIVE_DOUBLE, &
  table, localShape, h5err, xfer_prp=plistId, &
  mem_space_id = memspace, file_space_id = subspace)
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(memspace, h5err)
CALL h5sclose_f(subspace, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)

DEALLOCATE(entriesPerRank)
DEALLOCATE(rankOffset)

END SUBROUTINE READ_TABLE_DBL

!***********************************************************************
! Subroutine WRITE_LIST_INT8
! This subroutine writes a table dataset, i.e., a 2D array of shape
! (fixed, arbitrary) to an already-open (in parallel) h5 file, 
! where each processor contributes
! a section of the table (each processer can have a different number
! of entries in the 2nd dimension).
! Each rank must call this.
!
! Currently writes an array of doubles.
!
! INPUT:
! 
! - fileId the handle of the hdf5 file (that's already open)
! - datasetName the name of the field dataset
! - table the 2D array (for the local domain)
! - writeRankInfo whether to write info about the number of entries per rank
! - COMM: mpi comm
! - it the timestep
! - t the time
!
!***********************************************************************

SUBROUTINE WRITE_LIST_INT8(fileId, datasetName, table, writeRankInfo, &
  COMM, it, t)

IMPLICIT NONE

INCLUDE 'mpif.h'

!INPUT
INTEGER(HID_T)            :: fileId
CHARACTER(len=*)          :: datasetName
INTEGER*8                 :: table(:)
LOGICAL                   :: writeRankInfo
INTEGER                   :: COMM
INTEGER, OPTIONAL         :: it
DOUBLE PRECISION, OPTIONAL:: t

! LOCAL
INTEGER                   :: numRanks, rank, mpiErr
! bounds = [[xLower, yLower, zLower], [xUpper, yUpper, zUpper]], 
!   is the bounds of the first NDIM columns of the table
INTEGER*8, ALLOCATABLE    :: entriesPerRank(:)
INTEGER*8, ALLOCATABLE    :: rankOffset(:)
INTEGER                   :: d
INTEGER*8                 :: numLocalEntries

INTEGER(HSIZE_T)          :: localShape(1), globalShape(1), localOffset(1)

INTEGER(HID_T) :: dsetId, plistId2
INTEGER(HID_T) :: filespace, subspace, memspace
INTEGER        :: h5err
INTEGER        :: serialWritingRank = 0

INTEGER(HID_T) :: h5int8
!***********************************************************************

! Can't figure out how to write 8-byte integers from fortran
!h5int8 = h5kind_to_type(8, H5_INTEGER_KIND)
h5int8 = H5T_NATIVE_DOUBLE

localShape = SHAPE(table)
globalShape(1) = SIZE(table, DIM=1)
numLocalEntries = SIZE(table, DIM=1)

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
ALLOCATE(entriesPerRank(numRanks))
ALLOCATE(rankOffset(numRanks+1))

CALL MPI_ALLGATHER(numLocalEntries, 1, MPI_INTEGER8, entriesPerRank, &
  1, MPI_INTEGER8, COMM, mpiErr)

CALL CUMSUM(entriesPerRank, rankOffset(1:numRanks), rankOffset(numRanks+1))

globalShape(1) = rankOffset(numRanks+1)
localOffset(1) = rankOffset(rank+1)

!{
! create data space for dataset
CALL h5screate_simple_f(1, globalShape, filespace, h5err)
CALL h5screate_simple_f(1, localShape, memspace, h5err)
! create dataset with default properties
CALL h5dcreate_f(fileId, datasetName, h5int8, &
  filespace, dsetId, h5err)
! select the part that will be written by this domain
CALL h5dget_space_f(dsetId, subspace, h5err)
CALL h5sselect_hyperslab_f(subspace, H5S_SELECT_SET_F, localOffset, &
  localShape, h5err)
! create property list for collective dataset write
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId2, h5err)
CALL h5pset_dxpl_mpio_f(plistId2, H5FD_MPIO_COLLECTIVE_F, h5err)

! write table collectively
CALL h5dwrite_f(dsetId, h5int8, DBLE(table), localShape, h5err, &
  file_space_id = subspace, mem_space_id = memspace, xfer_prp=plistId2)

CALL h5sclose_f(memspace, h5err)
CALL h5sclose_f(subspace, h5err)

CALL h5pclose_f(plistId2, h5err)
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
!}

IF (writeRankInfo) THEN
  ! write the number of entries for each rank
  CALL WRITE_DATASET_1DINT8(fileId, datasetName // 'EntriesPerRank', serialWritingRank, &
    (/ INT8(numRanks) /), entriesPerRank, COMM)
  ! write the offset for each rank
  CALL WRITE_DATASET_1DINT8(fileId, datasetName // 'StartingEntryPerRank', serialWritingRank,&
    (/ INT8(numRanks+1) /), rankOffset, COMM)
ENDIF

! write step and time in serial
IF (PRESENT(it)) THEN
  CALL WRITE_DATASET_1DINT(fileId, 'step', serialWritingRank, (/ INT8(1) /),&
    (/ it /), COMM)
ENDIF
IF (PRESENT(t)) THEN
  CALL WRITE_DATASET_1DDBL(fileId, 'time', serialWritingRank, (/ INT8(1) /),&
    (/ t /), COMM)
ENDIF

DEALLOCATE(rankOffset)
DEALLOCATE(entriesPerRank)

END SUBROUTINE WRITE_LIST_INT8

!***********************************************************************
! Subroutine READ_LIST_INT8
! This subroutine reads a table dataset, i.e., a 1D array of shape
! (fixed, arbitrary) from an already-open (in parallel) h5 file, 
! where each processor contributes
! a section of the table (each processer can have a different number
! of entries in the 2nd dimension).
! Each rank must call this.
!
! INPUT:
! 
! - fileId the handle of the hdf5 file (that's already open)
! - datasetName the name of the field dataset
! - table the 2D array (for the local domain)
! - COMM: mpi comm
! - entriesDataset: the name of the dataset specificying how many entries
!     there are from each rank; by default this is
!     datasetName // "EntriesPerRank"
!
!***********************************************************************

SUBROUTINE READ_LIST_INT8(fileId, datasetName, table, COMM, &
  entriesDataset)

IMPLICIT NONE

INCLUDE 'mpif.h'

!INPUT
INTEGER(HID_T)                              :: fileId
CHARACTER(len=*)                            :: datasetName
INTEGER*8, ALLOCATABLE, INTENT(OUT)         :: table(:)
INTEGER                                     :: COMM
CHARACTER(len=*), OPTIONAL                  :: entriesDataset

! LOCAL
DOUBLE PRECISION, ALLOCATABLE :: dtable(:)
INTEGER                   :: numRanks, rank, mpiErr
INTEGER                   :: numFileRanks
! bounds = [[xLower, yLower, zLower], [xUpper, yUpper, zUpper]], 
!   is the bounds of the first NDIM columns of the table
INTEGER*8, ALLOCATABLE    :: entriesPerRank(:)
INTEGER*8, ALLOCATABLE    :: rankOffset(:)

INTEGER(HSIZE_T)          :: localShape(2), globalShape(2), localOffset(2)
INTEGER(HSIZE_T)          :: maxShape(2)

INTEGER(HID_T) :: dsetId, plistId
INTEGER(HID_T) :: filespace, subspace, memspace
INTEGER        :: h5err
INTEGER        :: serialReadingRank = 0

CHARACTER(len=256) :: edName
!***********************************************************************

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)

IF (PRESENT(entriesDataset)) THEN
  edName = entriesDataset
ELSE
  edName = datasetName // 'EntriesPerRank'
ENDIF

! read the number of entries for each rank
CALL READ_DATASET_1DINT8(fileId, trim(edName), &
  serialReadingRank, entriesPerRank, COMM)
numFileRanks = SIZE(entriesPerRank)
ALLOCATE(rankOffset(numFileRanks+1))
CALL CUMSUM(entriesPerRank, rankOffset(1:numFileRanks), &
  rankOffset(numFileRanks+1))

IF (numFileRanks /= numRanks) THEN
  IF (rank == 0) THEN
    PRINT *, "Error: cannot restore from run on ", numFileRanks, " ranks"
    PRINT *, "   to run on ", numRanks, " ranks"
  ENDIF
  CALL MPI_BARRIER(COMM, mpiErr)
  CALL MPI_ABORT(COMM, 9124, mpiErr)
  STOP 9124
ENDIF

CALL h5dopen_f(fileId, datasetName, dsetId, h5err)
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId, h5err)
CALL h5pset_dxpl_mpio_f(plistId, H5FD_MPIO_COLLECTIVE_F, h5err)
CALL h5dget_space_f(dsetId, filespace, h5err)
CALL h5sget_simple_extent_dims_f(filespace, globalShape, maxShape, h5err)

localShape(1) = entriesPerRank(rank+1)
localOffset(1) = rankOffset(rank+1)

ALLOCATE(table(localShape(1)))
ALLOCATE(dtable(localShape(1)))

CALL h5screate_simple_f(1, localShape, memspace, h5err)
! select the part that will be written by this domain
CALL h5dget_space_f(dsetId, subspace, h5err)
CALL h5sselect_hyperslab_f(subspace, H5S_SELECT_SET_F, localOffset, &
  localShape, h5err)

CALL h5dread_f(dsetId, H5T_NATIVE_DOUBLE, &
  dtable, localShape, h5err, xfer_prp=plistId, &
  mem_space_id = memspace, file_space_id = subspace)
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(memspace, h5err)
CALL h5sclose_f(subspace, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5pclose_f(plistId, h5err)
table = dtable

DEALLOCATE(dtable)
DEALLOCATE(entriesPerRank)
DEALLOCATE(rankOffset)

END SUBROUTINE READ_LIST_INT8

!***********************************************************************
! Subroutine APPEND_DATASET_1DDBL_SER
! This subroutine appends to (or creates) an extendible 1D dataset in 
! an h5 file opened (in serial) by one rank.
! This should be called by a single rank.
!
! INPUT:
! 
! - fileId the handle to the opened file
! - datasetName the name of the field dataset
! - ra the 1D array to append
! - datasetExists: whether the dataset exists, or (FALSE) must be created
!
!***********************************************************************

SUBROUTINE APPEND_DATASET_1DDBL_SER(fileId, datasetName, ra, datasetExists)

IMPLICIT NONE

INTEGER, PARAMETER        :: FDIM = 1

!INPUT
INTEGER(HID_T)            :: fileId
CHARACTER(len=*)          :: datasetName
DOUBLE PRECISION          :: ra(:)
LOGICAL                   :: datasetExists

! LOCAL
INTEGER(HSIZE_T)          :: localShape(FDIM), maxShape(FDIM), chunkShape(FDIM)
INTEGER(HSIZE_T)          :: offset(FDIM), newShape(FDIM)
INTEGER(HSIZE_T)          :: entries

INTEGER(HID_T)            :: dsetId, praId, filespace, memspace
INTEGER(HID_T)            :: stepDsetId, timeDsetId
INTEGER                   :: h5err
!***********************************************************************

localShape = SHAPE(ra)
entries = localShape(1)

maxShape(1) = H5S_UNLIMITED_F
chunkShape(1) = MAX(1, entries)

IF (.NOT.datasetExists) THEN
  ! Create file with default properties
  CALL h5screate_simple_f(FDIM, localShape, filespace, h5err, maxShape)

  ! Enable chunking (needed for extendible datasets)
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, praId, h5err)
  CALL h5pset_chunk_f(praId, FDIM, chunkShape, h5err)

  ! create dataset
  CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_DOUBLE, filespace, &
    dsetId, h5err, praId)

  ! write
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, localShape, h5err)

  CALL h5pclose_f(praId, h5err)
  CALL h5sclose_f(filespace, h5err)
  CALL h5dclose_f(dsetId, h5err)
ELSE ! file exists already
  ! open dataset (assume dataset exists because file exists
  CALL h5dopen_f(fileId, datasetName, dsetId, h5err)
  ! get filespace
  CALL h5dget_space_f(dsetId, filespace, h5err)
  ! Get current dataset dimensions
  CALL h5sget_simple_extent_dims_f(filespace, offset, maxShape, h5err)
  CALL h5sclose_f(filespace, h5err)
  
  newShape = offset
  newShape(1) = newShape(1) + entries
  ! extend dataset
  CALL h5dset_extent_f(dsetId, newShape, h5err)
  CALL h5dget_space_f(dsetId, filespace, h5err)
  ! write to extended part
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
    offset, localShape, h5err)
  CALL h5screate_simple_f(FDIM, localShape, memspace, h5err)

  ! write
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, localShape, h5err, &
    memspace, filespace)

  CALL h5sclose_f(memspace, h5err)
  CALL h5sclose_f(filespace, h5err)
  CALL h5dclose_f(dsetId, h5err)
ENDIF

END SUBROUTINE APPEND_DATASET_1DDBL_SER

!***********************************************************************
! Subroutine APPEND_DATASET_2DDBL_SER
! This subroutine appends to (or creates) an extendible 2D dataset in 
! an h5 file opened (in serial) by one rank; the extendible dimension
! is the second dimension in fortran (which will be transposed to the 
! first dimension in the h5 file)
! This should be called by a single rank.
!
! INPUT:
! 
! - fileId the handle to the opened file
! - datasetName the name of the field dataset
! - ra the 2D array to append
! - datasetExists: whether the dataset exists, or (FALSE) must be created
!
!***********************************************************************

SUBROUTINE APPEND_DATASET_2DDBL_SER(fileId, datasetName, ra, datasetExists)

IMPLICIT NONE

INTEGER, PARAMETER        :: FDIM = 2

!INPUT
INTEGER(HID_T)            :: fileId
CHARACTER(len=*)          :: datasetName
DOUBLE PRECISION          :: ra(:,:)
LOGICAL                   :: datasetExists

! LOCAL
INTEGER(HSIZE_T)          :: localShape(FDIM), maxShape(FDIM), chunkShape(FDIM)
INTEGER(HSIZE_T)          :: offset(FDIM), newShape(FDIM)
INTEGER(HSIZE_T)          :: entries

INTEGER(HID_T)            :: dsetId, plistId, filespace, memspace
INTEGER(HID_T)            :: stepDsetId, timeDsetId
INTEGER                   :: h5err
!***********************************************************************

localShape = SHAPE(ra)
entries = localShape(FDIM)

maxShape = localShape
maxShape(FDIM) = H5S_UNLIMITED_F
  
chunkShape(1) = localShape(1)
chunkShape(FDIM) = MAX(1, entries)

IF (.NOT.datasetExists) THEN
  CALL h5screate_simple_f(FDIM, localShape, filespace, h5err, maxShape)
  ! Enable chunking (needed for extendible datasets)
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plistId, h5err)
  CALL h5pset_chunk_f(plistId, FDIM, chunkShape, h5err)

  ! create dataset
  CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_DOUBLE, filespace, &
    dsetId, h5err, plistId)

  ! write
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, localShape, h5err)

  CALL h5pclose_f(plistId, h5err)
  CALL h5sclose_f(filespace, h5err)
  CALL h5dclose_f(dsetId, h5err)
ELSE ! file exists already
  ! open dataset (assume dataset exists because file exists
  CALL h5dopen_f(fileId, datasetName, dsetId, h5err)
  ! get filespace
  CALL h5dget_space_f(dsetId, filespace, h5err)
  ! CALL h5sget_simple_extent_ndims_f(filespace, datasetRank, h5err)
  ! Get current dataset dimensions
  CALL h5sget_simple_extent_dims_f(filespace, offset, maxShape, h5err)
  CALL h5sclose_f(filespace, h5err)
  
  newShape = offset
  offset(1) = 0
  newShape(FDIM) = newShape(FDIM) + entries
  ! extend dataset
  CALL h5dset_extent_f(dsetId, newShape, h5err)
  CALL h5dget_space_f(dsetId, filespace, h5err)
  ! write to extended part
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
    offset, localShape, h5err)
  CALL h5screate_simple_f(FDIM, localShape, memspace, h5err)

  ! write
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, localShape, h5err, &
    memspace, filespace)

  CALL h5sclose_f(memspace, h5err)
  CALL h5sclose_f(filespace, h5err)
  CALL h5dclose_f(dsetId, h5err)
ENDIF

END SUBROUTINE APPEND_DATASET_2DDBL_SER

!***********************************************************************
! Subroutine APPEND_TABLE_3DDBL_SER
! This subroutine appends to (or creates) an extendible 3D dataset in 
! an h5 file opened (in serial) by one rank, from a 2D table that is
! reshaped to 3D; the extendible dimension
! is the third dimension in fortran (which will be transposed to the 
! first dimension in the h5 file)
! This should be called by a single rank.
!
! A 2D ra of shape (M, N) will be reshaped to a 3D array of shape
! (M, nMidDim, N/nMidDim).
!
! INPUT:
! 
! - fileId the handle to the opened file
! - datasetName the name of the field dataset
! - ra the 2D array to append
! - datasetExists: whether the dataset exists, or (FALSE) must be created
! - nMidDim the middle dimension of the 3D array
!
!***********************************************************************

SUBROUTINE APPEND_TABLE_3DDBL_SER(fileId, datasetName, ra, datasetExists, &
  nMidDim)

IMPLICIT NONE

INTEGER, PARAMETER        :: FDIM = 3

!INPUT
INTEGER(HID_T)            :: fileId
CHARACTER(len=*)          :: datasetName
DOUBLE PRECISION          :: ra(:,:)
INTEGER                   :: nMidDim
LOGICAL                   :: datasetExists

! LOCAL
INTEGER(HSIZE_T)          :: localShape(FDIM), maxShape(FDIM), chunkShape(FDIM)
INTEGER(HSIZE_T)          :: offset(FDIM), newShape(FDIM)

INTEGER(HID_T)            :: dsetId, plistId, filespace, memspace
INTEGER(HID_T)            :: stepDsetId, timeDsetId
INTEGER                   :: h5err
!***********************************************************************

localShape(1) = SIZE(ra,DIM=1)
localShape(2) = nMidDim
localShape(3) = SIZE(ra,DIM=2) / nMidDim

maxShape(1) = localShape(1)
maxShape(2) = nMidDim
maxShape(FDIM) = H5S_UNLIMITED_F
  
chunkShape(1) = localShape(1)
! A few tests on verus suggest chunkShape(2)=localShape(2) could decrease
! file-writing-time by about 25% ; increasing chunkChape(2) to 2 or 4
! decreases file-writing-time a little bit (3% and 8%...I'm not sure
! what is statistically significant).  But 1 is convenient for reading
! the trajectory of a single particle.
chunkShape(2) = 1
chunkShape(FDIM) = localShape(FDIM)

IF (.NOT.datasetExists) THEN
  CALL h5screate_simple_f(FDIM, localShape, filespace, h5err, maxShape)
  ! Enable chunking (needed for extendible datasets)
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plistId, h5err)
  CALL h5pset_chunk_f(plistId, FDIM, chunkShape, h5err)

  ! create dataset
  CALL h5dcreate_f(fileId, datasetName, H5T_NATIVE_DOUBLE, filespace, &
    dsetId, h5err, plistId)

  ! write
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, localShape, h5err)

  CALL h5pclose_f(plistId, h5err)
  CALL h5sclose_f(filespace, h5err)
  CALL h5dclose_f(dsetId, h5err)
ELSE ! file exists already
  ! open dataset (assume dataset exists because file exists
  CALL h5dopen_f(fileId, datasetName, dsetId, h5err)
  ! get filespace
  CALL h5dget_space_f(dsetId, filespace, h5err)
  ! CALL h5sget_simple_extent_ndims_f(filespace, datasetRank, h5err)
  ! Get current dataset dimensions
  CALL h5sget_simple_extent_dims_f(filespace, offset, maxShape, h5err)
  CALL h5sclose_f(filespace, h5err)
  
  newShape = offset
  offset(1) = 0
  offset(2) = 0
  newShape(FDIM) = newShape(FDIM) + localShape(FDIM)
  ! extend dataset
  CALL h5dset_extent_f(dsetId, newShape, h5err)
  CALL h5dget_space_f(dsetId, filespace, h5err)
  ! write to extended part
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
    offset, localShape, h5err)
  CALL h5screate_simple_f(FDIM, localShape, memspace, h5err)

  ! write
  CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, ra, localShape, h5err, &
    memspace, filespace)

  CALL h5sclose_f(memspace, h5err)
  CALL h5sclose_f(filespace, h5err)
  CALL h5dclose_f(dsetId, h5err)
ENDIF

END SUBROUTINE APPEND_TABLE_3DDBL_SER

!***********************************************************************
! Subroutine WRITE_GLOBAL_FIELD_EXTRAS
! writes time-step and, if desired, coordinate axes and descriptions
!
! INPUT:
! 
! - fileId a collectively-opened hdf5 file 
! - rankToWrite: the (0-indexed) rank that knows F and will write it to disk.
! - it: Timestep number
! - COMM: mpi comm
! - xgp,ygp,zpg: coordinates of field values
! - xLabel, yLabel, zLabel: names of coordinates, like "x" or "phi"
!
!***********************************************************************

SUBROUTINE WRITE_GLOBAL_FIELD_EXTRAS(fileId, rankToWrite, it, COMM, &
  xgp, xLabel, ygp, yLabel, zgp, zLabel)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
INTEGER(HID_T) :: fileId
INTEGER                                  :: rankToWrite
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xgp, ygp, zgp
CHARACTER(len=*), OPTIONAL               :: xLabel, yLabel, zLabel
INTEGER                                  :: it
INTEGER                                  :: COMM

! LOCAL
INTEGER                                  :: FDIM
CHARACTER(len=10)                        :: cit
CHARACTER(len=1) :: dirStr
!***********************************************************************

write(cit,'(i10)') it
cit = adjustl(cit)

FDIM = 1
IF (PRESENT(ygp)) FDIM = 2
IF (PRESENT(zgp)) FDIM = 3

CALL WRITE_DATASET_1DINT(fileId, 'step', rankToWrite, (/ INT8(1) /), &
  (/ it /), COMM)
! hdf5 transposes fortran arrays, so reverse axis labels
IF (PRESENT(xgp) .AND. PRESENT(xLabel)) THEN
  write(dirStr,'(i1)') (FDIM-1)
  CALL WRITE_DATASET_1DDBL(fileId, 'axis' // dirStr // 'coords', &
    rankToWrite, INT8(SHAPE(xgp)), xgp, COMM, xLabel)
ENDIF
IF (PRESENT(ygp) .AND. PRESENT(yLabel)) THEN
  write(dirStr,'(i1)') (FDIM-2)
  CALL WRITE_DATASET_1DDBL(fileId, 'axis' // dirStr // 'coords', &
    rankToWrite, INT8(SHAPE(ygp)), ygp, COMM, yLabel)
ENDIF
IF (PRESENT(zgp) .AND. PRESENT(zLabel)) THEN
  write(dirStr,'(i1)') (FDIM-3)
  CALL WRITE_DATASET_1DDBL(fileId, 'axis' // dirStr // 'coords', &
    rankToWrite, INT8(SHAPE(zgp)), zgp, COMM, zLabel)
ENDIF

END SUBROUTINE WRITE_GLOBAL_FIELD_EXTRAS

!***********************************************************************
! Subroutine DUMP_GLOBAL_FIELD_1D
! This subroutine writes a 1D array (known entirely to one rank) to disk
! This should be called by all ranks, even though only one rank might
! know the field and write it.
!
! While F will not be accessed except by the rankToWrite,
! Fshape needs to be correct for all ranks.
!
! INPUT:
! 
! - fileBaseName the base name of the file, 
!     to which "_it.h5" will be appended (see argument "it", below)
! - datasetname the name of the field dataset
! - Fshape = SHAPE(F)
! - F: (entire) array of field values
! - rankToWrite: the (0-indexed) rank that knows F and will write it to disk.
! - it: Timestep number
! - COMM: mpi comm
! - xgp,ygp,zpg: coordinates of field values
! - xLabel, yLabel, zLabel: names of coordinates, like "x" or "phi"
!
!***********************************************************************

SUBROUTINE DUMP_GLOBAL_FIELD_1D(fileBaseName, datasetname, Fshape, F, &
  rankToWrite, it, COMM, xgp, xLabel, ygp, yLabel, zgp, zLabel)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER                       :: FDIM = 1

! INPUT
CHARACTER(len=*)                         :: fileBaseName
CHARACTER(len=*)                         :: datasetname
! following must have rank=FDIM
DOUBLE PRECISION, DIMENSION(:)           :: F
INTEGER(HSIZE_T)                         :: Fshape(FDIM)
INTEGER                                  :: rankToWrite
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xgp, ygp, zgp
CHARACTER(len=*), OPTIONAL               :: xLabel, yLabel, zLabel
INTEGER                                  :: it
INTEGER                                  :: COMM

! LOCAL
CHARACTER(len=10)                        :: cit

! hdf5 
INTEGER(HID_T) :: fileId
INTEGER(HID_T) :: plistId
INTEGER        :: h5err
!***********************************************************************

write(cit,'(i10)') it
cit = adjustl(cit)

! initialize fortran interface
CALL h5open_f(h5err)

! setup file access property list for parallel i/o
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistId, h5err)
CALL h5pset_fapl_mpio_f(plistId, COMM, MPI_INFO_NULL, h5err)

! create file collectively (erase any existing data)
CALL h5fcreate_f(fileBaseName // "_" // trim(cit) // ".h5", &
  H5F_ACC_TRUNC_F, fileId, h5err, access_prp = plistId)
CALL h5pclose_f(plistId, h5err)

CALL WRITE_DATASET_1DDBL(fileId, datasetname, rankToWrite, Fshape,&
  F, COMM)

CALL WRITE_GLOBAL_FIELD_EXTRAS(fileId, rankToWrite, it, COMM, &
  xgp, xLabel, ygp, yLabel, zgp, zLabel)
  
CALL h5fclose_f(fileId, h5err)
CALL h5close_f(h5err)

END SUBROUTINE DUMP_GLOBAL_FIELD_1D

!***********************************************************************
! Subroutine DUMP_GLOBAL_FIELD_3D
! This subroutine writes a 3D array (known entirely to one rank) to disk
! This should be called by all ranks, even though only one rank might
! know the field and write it.
!
! While F will not be accessed except by the rankToWrite,
! Fshape needs to be correct for all ranks.
!
! INPUT:
! 
! - fileBaseName the base name of the file, 
!     to which "_it.h5" will be appended (see argument "it", below)
! - datasetname the name of the field dataset
! - Fshape = SHAPE(F) : must be correct on all ranks
! - F: (entire) array of field values
! - rankToWrite: the (0-indexed) rank that knows F and will write it to disk.
! - it: Timestep number
! - transp: whether to transpose the axes
!   Note that hdf5 is row-major and fortran is row-minor, so
!   when preserving memory layout, (transp=.FALSE.) the hdf5 array looks
!   transposed compared to fortran.
!   E.g., given F(x,y,z) and transp=.FALSE.
!    the hdf5 array will have axis 0 = z, axis 1 = y, axis 2 = x.
!    In both fortran and hdf5, F(x,y,z) and F(x+dx,y,z) will be adjacent
!    in memory.
!   If transp=.TRUE., then
!   in hdf5, axis 0 = x, axis 1 = y, axis 2 = z, and in the hdf5 layout,
!   F(x,y,z) and F(x,y,z+dz) will be adjacent.
! - COMM: mpi comm
! - xgp,ygp,zpg: coordinates of field values
! - xLabel, yLabel, zLabel: names of coordinates, like "x" or "phi"
!
!***********************************************************************

RECURSIVE SUBROUTINE DUMP_GLOBAL_FIELD_3D(fileBaseName, datasetname, Fshape,&
  F, rankToWrite, it, transp, COMM, xgp, xLabel, ygp, yLabel, zgp, zLabel)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER                       :: FDIM = 3

! INPUT
CHARACTER(len=*)                         :: fileBaseName
CHARACTER(len=*)                         :: datasetname
! following must have rank=FDIM
DOUBLE PRECISION, DIMENSION(:,:,:)       :: F
INTEGER(HSIZE_T), DIMENSION(FDIM)        :: Fshape
INTEGER                                  :: rankToWrite
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xgp, ygp, zgp
CHARACTER(len=*), OPTIONAL               :: xLabel, yLabel, zLabel
INTEGER                                  :: it
LOGICAL                                  :: transp
INTEGER                                  :: COMM

! LOCAL
CHARACTER(len=10)                        :: cit
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE       :: Ft
INTEGER(HSIZE_T), DIMENSION(FDIM)        :: FshapeT
INTEGER                                  :: i, j, k

! hdf5 
INTEGER(HID_T) :: fileId
INTEGER(HID_T) :: plistId
INTEGER        :: h5err
!***********************************************************************

IF (transp .EQV. .TRUE.) THEN !{

  DO i = 1, FDIM
    FshapeT(i) = Fshape(FDIM+1-i)
  ENDDO
  ALLOCATE(Ft(FshapeT(1), FshapeT(2), FshapeT(3)))
  DO j = 1, Fshape(2)
    Ft(:,j,:) = TRANSPOSE(F(:,j,:))
  ENDDO
  IF (PRESENT(zLabel)) THEN
    CALL DUMP_GLOBAL_FIELD_3D(fileBaseName, datasetname, FshapeT,&
      Ft, rankToWrite, it, .FALSE., COMM, &
      zgp, zLabel, ygp, yLabel, xgp, xLabel)
  ELSE
    CALL DUMP_GLOBAL_FIELD_3D(fileBaseName, datasetname, FshapeT,&
      Ft, rankToWrite, it, .FALSE., COMM)
  ENDIF
  DEALLOCATE(Ft)

ELSE ! }{
! Do not transpose

write(cit,'(i10)') it
cit = adjustl(cit)

! initialize fortran interface
CALL h5open_f(h5err)

! setup file access property list for parallel i/o
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistId, h5err)
CALL h5pset_fapl_mpio_f(plistId, COMM, MPI_INFO_NULL, h5err)

! create file collectively (erase any existing data)
CALL h5fcreate_f(fileBaseName // "_" // trim(cit) // ".h5", &
  H5F_ACC_TRUNC_F, fileId, h5err, access_prp = plistId)
CALL h5pclose_f(plistId, h5err)

CALL WRITE_DATASET_3DDBL(fileId, datasetname, rankToWrite, Fshape, F, COMM)

CALL WRITE_GLOBAL_FIELD_EXTRAS(fileId, rankToWrite, it, COMM, &
  xgp, xLabel, ygp, yLabel, zgp, zLabel)
  
CALL h5fclose_f(fileId, h5err)
CALL h5close_f(h5err)

ENDIF !}

END SUBROUTINE DUMP_GLOBAL_FIELD_3D

!***********************************************************************
! Subroutine DUMP_FIELD
! This subroutine writes a 3D array, decomposed in 3D across all ranks,
! to disk. 
!
! INPUT:
! 
! - fileBaseName the base name of the file, 
!     to which "_DUMPNUM.h5" will be appended
! - datasetname the name of the field dataset
! - F: array of (local domain) field values
! - xgp,ygp: coordinates of (local domain) field values
! - it: Timestep number
! - COMM: mpi comm
!
! OUTPUT:
! - ierr: error code
!
!***********************************************************************

SUBROUTINE DUMP_FIELD(fileBaseName, datasetname, F,it, COMM, &
  xgp, xLabel, ygp, yLabel, zgp, zLabel)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER                       :: FDIM = NDIM

! INPUT
CHARACTER(len=*)                         :: fileBaseName
CHARACTER(len=*)                         :: datasetname
! following must have rank=FDIM=NDIM
DOUBLE PRECISION, DIMENSION(:,:,:)         :: F
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xgp, ygp, zgp
CHARACTER(len=*), OPTIONAL               :: xLabel, yLabel, zLabel
INTEGER                                  :: it
INTEGER                                  :: COMM

! LOCAL
CHARACTER(len=10)                        :: cit

! shape of local and global field array
INTEGER(HSIZE_T), DIMENSION(FDIM)   :: localShape
INTEGER(HSIZE_T), DIMENSION(FDIM)   :: globalShape 
! offset of local array origin within global array
INTEGER(HSIZE_T), DIMENSION(FDIM)   :: offset

! hdf5 
INTEGER(HID_T) :: fileId
INTEGER(HID_T) :: dsetId
INTEGER(HID_T) :: filespace, subspace, memspace
INTEGER(HID_T) :: plistId, plistId2
INTEGER        :: h5err
INTEGER        :: mpiErr, numRanks, rank

INTEGER        :: d
INTEGER        :: coordWritingRank = 0
INTEGER        :: stepWritingRank = 0
CHARACTER(len=1) :: dirStr
!***********************************************************************

write(cit,'(i10)') it
cit = adjustl(cit)

localShape = SHAPE(F) 

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL FILL_GLOBAL_FIELD_SHAPE(localShape, globalShape, offset, COMM)

! initialize fortran interface
CALL h5open_f(h5err)

! setup file access property list for parallel i/o
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistId, h5err)
CALL h5pset_fapl_mpio_f(plistId, COMM, MPI_INFO_NULL, h5err)

! create file collectively (erase any existing data)
CALL h5fcreate_f(fileBaseName // "_" // trim(cit) // ".h5", &
  H5F_ACC_TRUNC_F, fileId, h5err, access_prp = plistId)
CALL h5pclose_f(plistId, h5err)

! create data space for dataset
CALL h5screate_simple_f(FDIM, globalShape, filespace, h5err)
CALL h5screate_simple_f(FDIM, localShape, memspace, h5err)
! create dataset with default properties
CALL h5dcreate_f(fileId, datasetname, H5T_NATIVE_DOUBLE, filespace, &
  dsetId, h5err)
! select the part that will be written by this domain
CALL h5dget_space_f(dsetId, subspace, h5err)
CALL h5sselect_hyperslab_f(subspace, H5S_SELECT_SET_F, offset, localShape, h5err)
! create property list for collective dataset write
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId2, h5err)
CALL h5pset_dxpl_mpio_f(plistId2, H5FD_MPIO_COLLECTIVE_F, h5err)

! write field collectively
CALL h5dwrite_f(dsetId, H5T_NATIVE_DOUBLE, F, localShape, h5err, &
  file_space_id = subspace, mem_space_id = memspace, xfer_prp=plistId2)

CALL h5pclose_f(plistId2, h5err)
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5sclose_f(memspace, h5err)
CALL h5sclose_f(subspace, h5err)

! write step and grid-point coordinates in serial
CALL WRITE_DATASET_1DINT(fileId, 'step', stepWritingRank, (/ INT8(1) /), &
  (/ it /), COMM)
! hdf5 transposes fortran arrays, so reverse axis labels
IF (PRESENT(xgp) .AND. PRESENT(xLabel) .AND. FDIM >= 1) THEN
  write(dirStr,'(i1)') (FDIM-1)
  CALL WRITE_DATASET_1DDBL(fileId, 'axis' // dirStr // 'coords', &
    coordWritingRank, INT8(SHAPE(xgp)), xgp, COMM, xLabel)
ENDIF
IF (PRESENT(ygp) .AND. PRESENT(yLabel) .AND. FDIM >= 2) THEN
  write(dirStr,'(i1)') (FDIM-2)
  CALL WRITE_DATASET_1DDBL(fileId, 'axis' // dirStr // 'coords', &
    coordWritingRank, INT8(SHAPE(ygp)), ygp, COMM, yLabel)
ENDIF
IF (PRESENT(zgp) .AND. PRESENT(zLabel) .AND. FDIM >= 3) THEN
  write(dirStr,'(i1)') (FDIM-3)
  CALL WRITE_DATASET_1DDBL(fileId, 'axis' // dirStr // 'coords', &
    coordWritingRank, INT8(SHAPE(zgp)), zgp, COMM, zLabel)
ENDIF

CALL h5fclose_f(fileId, h5err)
CALL h5close_f(h5err)

END SUBROUTINE DUMP_FIELD

!***********************************************************************
! Subroutine SAVE_FIELDS
! This subroutine saves in parallel the E,B fields array (at the nodes) 
! to disk for future restoration of the simulation.
!
! INPUT:
! 
! - Eg,Bg: Yee fields including guard cells
! - it: Timestep number
! - coords: coordinates of each process in the cartesian topology
!
!***********************************************************************

SUBROUTINE SAVE_FIELDS(Exg,Eyg,Ezg,Bxg,Byg,Bzg,it,coords,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: id,COMM,ierr
INTEGER, DIMENSION(2)                    :: coords
! Following is NDIM-dependent
DOUBLE PRECISION, DIMENSION(:,:,:)       :: Bxg,Byg,Bzg,Exg,Eyg,Ezg
INTEGER                                  :: it

!***********************************************************************

CALL DUMP_FIELD("./data_restore/Bx", "field", Bxg, it, COMM)
CALL DUMP_FIELD("./data_restore/By", "field", Byg, it, COMM)
CALL DUMP_FIELD("./data_restore/Bz", "field", Bzg, it, COMM)
CALL DUMP_FIELD("./data_restore/Ex", "field", Exg, it, COMM)
CALL DUMP_FIELD("./data_restore/Ey", "field", Eyg, it, COMM)
CALL DUMP_FIELD("./data_restore/Ez", "field", Ezg, it, COMM)

END SUBROUTINE SAVE_FIELDS

!***********************************************************************
! Subroutine LOAD_FIELD
! This subroutine reads a 2D array, decomposed in 2D across all ranks,
! from disk. 
!
! INPUT:
! 
! - fileBaseName the base name of the file, 
!     to which "_DUMPNUM.h5" will be appended
! - datasetname the name of the field dataset
! - F: array of (local domain) field values
!      must have correct shape
! - it: Timestep number
! - COMM: mpi comm
!
! OUTPUT:
! - ierr: error code
!
!***********************************************************************

SUBROUTINE LOAD_FIELD(fileBaseName, datasetname, it, F, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER                             :: FDIM = NDIM

! INPUT
CHARACTER(len=*)                               :: fileBaseName
CHARACTER(len=*)                               :: datasetname
! following must have rank=FDIM=NDIM
DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT)  :: F
INTEGER                                        :: it
INTEGER                                        :: COMM

! LOCAL
CHARACTER(len=10)                        :: cit

! shape of local and global field array
INTEGER(HSIZE_T), DIMENSION(FDIM)   :: localShape
INTEGER(HSIZE_T), DIMENSION(FDIM)   :: globalShape 
! offset of local array origin within global array
INTEGER(HSIZE_T), DIMENSION(FDIM)   :: offset

! hdf5 
INTEGER(HID_T) :: fileId
INTEGER(HID_T) :: dsetId
INTEGER(HID_T) :: filespace, subspace, memspace
INTEGER(HID_T) :: plistId, plistId2
INTEGER        :: h5err
INTEGER        :: mpiErr, numRanks, rank

INTEGER        :: d
INTEGER        :: coordWritingRank = 0
INTEGER        :: stepWritingRank = 0
CHARACTER(len=1) :: dirStr
!***********************************************************************

write(cit,'(i10)') it
cit = adjustl(cit)

localShape = SHAPE(F) 

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL FILL_GLOBAL_FIELD_SHAPE(localShape, globalShape, offset, COMM)

! initialize fortran interface
CALL h5open_f(h5err)

! setup file access property list for parallel i/o
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistId, h5err)
CALL h5pset_fapl_mpio_f(plistId, COMM, MPI_INFO_NULL, h5err)

! create file collectively (erase any existing data)
CALL h5fopen_f(fileBaseName // "_" // trim(cit) // ".h5", &
  H5F_ACC_RDONLY_F, fileId, h5err, access_prp = plistId)
CALL h5pclose_f(plistId, h5err)

! create data space for dataset
CALL h5screate_simple_f(FDIM, globalShape, filespace, h5err)
CALL h5screate_simple_f(FDIM, localShape, memspace, h5err)
! open dataset
CALL h5dopen_f(fileId, datasetname, dsetId, h5err)
! select the part that will be written by this domain
CALL h5dget_space_f(dsetId, subspace, h5err)
CALL h5sselect_hyperslab_f(subspace, H5S_SELECT_SET_F, offset, localShape, h5err)
! create property list for collective dataset read
CALL h5pcreate_f(H5P_DATASET_XFER_F, plistId2, h5err)
CALL h5pset_dxpl_mpio_f(plistId2, H5FD_MPIO_COLLECTIVE_F, h5err)

! read field collectively
CALL h5dread_f(dsetId, H5T_NATIVE_DOUBLE, F, localShape, h5err, &
  file_space_id = subspace, mem_space_id = memspace, xfer_prp=plistId2)

CALL h5pclose_f(plistId2, h5err)
CALL h5dclose_f(dsetId, h5err)
CALL h5sclose_f(filespace, h5err)
CALL h5sclose_f(memspace, h5err)
CALL h5sclose_f(subspace, h5err)

CALL h5fclose_f(fileId, h5err)
CALL h5close_f(h5err)

END SUBROUTINE LOAD_FIELD

!***********************************************************************
! Subroutine RESTORE_FIELDS
! This subroutine restores the E,B fields array (at the nodes) to disk from 
! saved data.
!
! INPUT:
! 
! - Eg,Bg: Yee fields, including guard cells
! - it: Timestep number
! - coords: coordinates of each process in the cartesian topology
!
!***********************************************************************

SUBROUTINE RESTORE_FIELDS(Exg,Eyg,Ezg,Bxg,Byg,Bzg,it,coords,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                     :: id,COMM,ierr
INTEGER, DIMENSION(NDIM)                    :: coords
! Following must have dimension=NDIM
DOUBLE PRECISION, DIMENSION(:,:,:)          :: Bxg,Byg,Bzg,Exg,Eyg,Ezg
INTEGER                                     :: it
!***********************************************************************

CALL LOAD_FIELD("./data_restore/Bx", "field", it, Bxg, COMM)
CALL LOAD_FIELD("./data_restore/By", "field", it, Byg, COMM)
CALL LOAD_FIELD("./data_restore/Bz", "field", it, Bzg, COMM)
CALL LOAD_FIELD("./data_restore/Ex", "field", it, Exg, COMM)
CALL LOAD_FIELD("./data_restore/Ey", "field", it, Eyg, COMM)
CALL LOAD_FIELD("./data_restore/Ez", "field", it, Ezg, COMM)

END SUBROUTINE RESTORE_FIELDS

!***********************************************************************
! Subroutine DUMP_PARTICLES
! This subroutine writes all the particle data to disk.
!
! INPUT:
! 
! - pcl: Particle distribution function
! - tag: Particle's tag
! - NPP: Number of particles per domain
! - spec: Particle' species
! - sym: Type of particles
! - it: Timestep number
! - restore: whether this is a restore-dump (affects directory)
!
!***********************************************************************

SUBROUTINE DUMP_PARTICLES(pcl,tag,NPP,it,spec,sym,id,restore,COMM, pcl_data)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: id,COMM
INTEGER                                  :: it,i
INTEGER*8                                :: NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
INTEGER*8, DIMENSION(1:NPP)              :: tag
CHARACTER(LEN=10)                        :: spec,sym
LOGICAL                                  :: restore
DOUBLE PRECISION, DIMENSION(1:4,1:NPP), OPTIONAL :: pcl_data


INTEGER(HID_T) :: fileId, plistId
INTEGER        :: h5err
CHARACTER(LEN=10)                        :: cit
CHARACTER(len=14)                        :: dir
!***********************************************************************

! Convert the integer it and id into a string cit
WRITE(cit,'(i10)') it

! This step left-justify the string
cit=adjustl(cit)
spec=adjustl(spec)
sym=adjustl(sym)

! initialize fortran interface
CALL h5open_f(h5err)

dir = "              "
IF (restore) THEN
  dir(1:12) = "data_restore"
ELSE
  dir(1:14) = "data/particles"
ENDIF

IF (writePerRankFiles) THEN
  OPEN(9, FILE=trim(perRankFile), POSITION='APPEND')
  WRITE(9, *) "Dumping ", NPP, " " // trim(spec) // "_" // trim(sym)
  CLOSE(9)
ENDIF
! setup file access property list for parallel i/o
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistId, h5err)
CALL h5pset_fapl_mpio_f(plistId, COMM, MPI_INFO_NULL, h5err)

! create file collectively (erase any existing data)
CALL h5fcreate_f( "./" // trim(dir) // "/" // trim(spec) &
  // "_" // trim(sym) // "_" // trim(cit) // ".h5", &
  H5F_ACC_TRUNC_F, fileId, h5err, access_prp = plistId)
CALL h5pclose_f(plistId, h5err)

CALL WRITE_TABLE_DBL(fileId, "particles", pcl, .TRUE., .TRUE., COMM, it, it*dt)
CALL WRITE_LIST_INT8(fileId, "tags", tag, .FALSE., COMM)

IF (PRESENT(pcl_data)) THEN
  CALL WRITE_TABLE_DBL(fileId, "particleData", &
    pcl_data, .FALSE., .FALSE., COMM)
ENDIF

CALL h5fclose_f(fileId, h5err)
CALL h5close_f(h5err)
  
END SUBROUTINE DUMP_PARTICLES

!***********************************************************************
! Subroutine DUMP_PARTICLES
! This subroutine writes all the particle data to disk.
!
! INPUT:
! 
! - pcl: Particle distribution function
! - tag: Particle's tag
! - NPP: Number of particles per domain
! - spec: Particle' species
! - sym: Type of particles
! - it: Timestep number
! - restore: whether this is a restore-dump (affects directory)
!
!***********************************************************************

SUBROUTINE SAVE_PARTICLES(pcl,pcl_data,tag,NPP,it,spec,sym,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: id,COMM,ierr
INTEGER                                  :: it
INTEGER*8                                :: NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
DOUBLE PRECISION, DIMENSION(1:4,1:NPP)   :: pcl_data
INTEGER*8, DIMENSION(1:NPP)              :: tag
CHARACTER(LEN=10)                        :: spec,sym
!***********************************************************************

CALL DUMP_PARTICLES(pcl,tag,NPP,it,spec,sym,id,.TRUE.,COMM, pcl_data)

END SUBROUTINE SAVE_PARTICLES

!***********************************************************************
! Subroutine LOAD_PARTICLES
! This subroutine reads all the particle data from disk.
!
! INPUT:
! 
! - pcl: shape = (7,?) Particle distribution function
! - pcl_data: shape = (4,?) other particle info
! - tag: Particle's tag
! - spec: Particle' species
! - sym: Type of particles
! - it: Timestep number
!
! The arrays pcl, tag (and if present, pcl_data) will be allocated to
! the correct size.
!
!***********************************************************************

SUBROUTINE LOAD_PARTICLES(pcl,tag,it,spec,sym,COMM, pcl_data)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                                 :: COMM
INTEGER                                                 :: it
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE           :: pcl
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, OPTIONAL :: pcl_data
INTEGER*8, DIMENSION(:), ALLOCATABLE                    :: tag
CHARACTER(LEN=10)                        :: spec,sym

CHARACTER(LEN=10)                        :: cit
INTEGER(HID_T) :: fileId, plistId
INTEGER        :: h5err
!***********************************************************************

! Convert the integer it and id into a string cit
WRITE(cit,'(i10)') it

! This step left-justifies the string
cit=adjustl(cit)
spec=adjustl(spec)
sym=adjustl(sym)

! initialize fortran interface
CALL h5open_f(h5err)

! setup file access property list for parallel i/o
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistId, h5err)
CALL h5pset_fapl_mpio_f(plistId, COMM, MPI_INFO_NULL, h5err)

! create file collectively (erase any existing data)
CALL h5fopen_f( "./data_restore/" // trim(spec) &
  // "_" // trim(sym) // "_" // trim(cit) // ".h5", &
  H5F_ACC_RDONLY_F, fileId, h5err, access_prp = plistId)
IF (h5err /= 0) STOP 123
CALL h5pclose_f(plistId, h5err)
IF (h5err /= 0) STOP 124

CALL READ_TABLE_DBL(fileId, "particles", pcl, COMM, "particlesEntriesPerRank")
CALL READ_LIST_INT8(fileId, "tags", tag, COMM, "particlesEntriesPerRank")

IF (PRESENT(pcl_data)) THEN
  CALL READ_TABLE_DBL(fileId, "particleData", &
    pcl_data, COMM, "particlesEntriesPerRank")
ENDIF

CALL h5fclose_f(fileId, h5err)
CALL h5close_f(h5err)
  
END SUBROUTINE LOAD_PARTICLES

!***********************************************************************
! Subroutine RESTORE_PARTICLES
! This subroutine restores the particle distribution function from 
! saved data.
!
! INPUT:
! 
! - pcl: Particle distribution function
! - pcl_data: Particle distribution data function
! - tag: Particle's tag
! - NPP: Number of particles per domain
! - spec: Particle' species
! - sym: Type of particles
! - it: Timestep number
!
!***********************************************************************

SUBROUTINE RESTORE_PARTICLES(pcl,pcl_data,tag,NPP,it,spec,sym,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                   :: id,COMM,ierr
INTEGER                                   :: it
INTEGER*8, INTENT(OUT)                    :: NPP
DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT):: pcl(:,:),pcl_data(:,:)
INTEGER*8, ALLOCATABLE, INTENT(OUT)       :: tag(:)
CHARACTER(LEN=10)                         :: spec,sym

!***********************************************************************

CALL LOAD_PARTICLES(pcl,tag,it,spec,sym,COMM, pcl_data)
NPP = SIZE(pcl, 2)

END SUBROUTINE RESTORE_PARTICLES

!***********************************************************************
! Subroutine GET_TRACK_SEQ_NUM(fileBaseName, seqNum)
!
! INPUT:
! 
! - fileBaseName: filename, without the suffix _collectorNum_seq.h5
! OUTPUT:
! - seqNum: the least unused number
!***********************************************************************

SUBROUTINE GET_TRACK_SEQ_NUM(fileBaseName, seqNum)

IMPLICIT NONE

! INPUT
CHARACTER(len=*)                         :: fileBaseName
! OUTPUT
INTEGER, INTENT(OUT)                     :: seqNum

! LOCAL
CHARACTER(len=10)                        :: seqStr
LOGICAL                                  :: fileExists
INTEGER                                  :: ioStatus
!***********************************************************************

DO seqNum = 0, 10000
  write(seqStr,'(i10)') seqNum
  seqStr = adjustl(seqStr)
! Make sure following filename agrees with that in DUMP_TRACKED_PARTICLES_SER
  INQUIRE(FILE=fileBaseName // "_cid0_" // trim(seqStr) // ".h5", &
    EXIST = fileExists, IOSTAT = ioStatus)
! Slightly dangerous: if INQUIRE fails every time, we loop to 10000
  IF (ioStatus == 0 .AND. .NOT. fileExists) THEN
    EXIT
  ENDIF
ENDDO

END SUBROUTINE GET_TRACK_SEQ_NUM

!***********************************************************************
! Subroutine DUMP_TRACKED_PARTICLES_SER
! This subroutine writes all the particle data to disk (called by 
! individual ranks that dump their own files.)
!
! INPUT:
! 
! - fileBaseName: filename, without the suffix _collectorNum_seq.h5
! - ptclTable: Particle data
! - collectorNum: the collector number (rank in collectorComm)
! - seqNum: a sequence number for the file (in case we want to write
!           to a new file after some time)
! - numPtcls: the number of particles
! - timeSteps: a sequence of timesteps
! - times: a sequence of times
!
!***********************************************************************

SUBROUTINE DUMP_TRACKED_PARTICLES_SER(fileBaseName, ptclTable, numPtcls, &
  collectorNum, seqNum, timeSteps, times)

IMPLICIT NONE

! INPUT
CHARACTER(len=*)                         :: fileBaseName
DOUBLE PRECISION, DIMENSION(:,:)         :: ptclTable
INTEGER                                  :: numPtcls, collectorNum, seqNum
INTEGER, DIMENSION(:), OPTIONAL          :: timeSteps
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: times
! LOCAL
CHARACTER(len=256)                       :: fileName
LOGICAL                                  :: fileExists
CHARACTER(len=10)                        :: collectorStr, seqStr
INTEGER(HID_T)                           :: fileId
INTEGER                                  :: h5err
!***********************************************************************

write(collectorStr,'(i10)') collectorNum
write(seqStr,'(i10)') seqNum
collectorStr = adjustl(collectorStr)
seqStr = adjustl(seqStr)

! Make sure following filename agrees with that in GET_TRACK_SEQ_NUM
fileName = fileBaseName // "_cid" // trim(collectorStr) // "_" // trim(seqStr) &
  // ".h5"

! initialize fortran interface
CALL h5open_f(h5err)

INQUIRE(FILE=fileName, EXIST=fileExists)
IF (.NOT.fileExists) THEN
  ! Create file with default properties
  CALL h5fcreate_f(fileName, H5F_ACC_EXCL_F, fileId, h5err)
ELSE ! file exists already
  ! open file
  CALL h5fopen_f(fileName, H5F_ACC_RDWR_F, fileId, h5err)
ENDIF

CALL APPEND_TABLE_3DDBL_SER(fileId, "particles", ptclTable, fileExists, &
  numPtcls)

IF (PRESENT(timeSteps)) THEN
  CALL APPEND_DATASET_1DDBL_SER(fileId, "timeSteps", DBLE(timeSteps), fileExists)
ENDIF

IF (PRESENT(times)) THEN
  CALL APPEND_DATASET_1DDBL_SER(fileId, "times", times, fileExists)
ENDIF

! close file
CALL h5fclose_f(fileId, h5err)
! close fortran interface
CALL h5close_f(h5err)

END SUBROUTINE DUMP_TRACKED_PARTICLES_SER

!***********************************************************************
! Subroutine DUMP_SIM_INFO
!
! Save useful information about simulation
!***********************************************************************

SUBROUTINE DUMP_SIM_INFO(hostname, rankOrder, ranksPerNode, &
  domainLbs, domainUbs, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
CHARACTER*(MPI_MAX_PROCESSOR_NAME)    :: hostname
INTEGER          :: rankOrder, ranksPerNode, COMM
DOUBLE PRECISION, DIMENSION(NDIM)      :: domainLbs, domainUbs
! LOCAL
INTEGER          :: rank, numRanks, mpiErr
INTEGER, DIMENSION(NDIM) :: domainDecomp, domainIndex
INTEGER, DIMENSION(:,:), ALLOCATABLE  :: domainIndices
INTEGER                  :: hostLen, maxHostLen
CHARACTER, DIMENSION(:), ALLOCATABLE  :: allHostNames 
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE  :: domainLowerBounds
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE  :: domainUpperBounds

INTEGER(HID_T) :: fileId, plistId
INTEGER        :: h5err
!***********************************************************************

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL FILL_DOMAIN_DECOMP(COMM, mpiErr, domainDecomp, domainIndex)


! initialize fortran interface
CALL h5open_f(h5err)

! setup file access property list for parallel i/o
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistId, h5err)
CALL h5pset_fapl_mpio_f(plistId, COMM, MPI_INFO_NULL, h5err)

! create file collectively (erase any existing data)
CALL h5fcreate_f("./data/info.h5", &
  H5F_ACC_TRUNC_F, fileId, h5err, access_prp = plistId)
CALL h5pclose_f(plistId, h5err)

CALL WRITE_DATASET_1DINT(fileId, "numMpiRanks", 0, &
  INT8( (/ 1 /) ), (/ numRanks /), COMM)

CALL WRITE_DATASET_1DINT(fileId, "ranksPerNode", 0, &
  INT8( (/ 1 /) ), (/ ranksPerNode /), COMM)

CALL WRITE_DATASET_1DINT(fileId, "rankOrder", 0, &
  INT8( (/ 1 /) ), (/ rankOrder /), COMM, &
  "0=round robin, 1=adj ranks on same node")

CALL WRITE_DATASET_1DINT(fileId, "domainDecomp", 0, &
  INT8( (/ NDIM /) ), domainDecomp, COMM)

! Get domain indices
IF (rank == 0) THEN
  ALLOCATE(domainIndices(NDIM, numRanks))
ENDIF
CALL MPI_GATHER(domainIndex, NDIM, MPI_INTEGER, domainIndices, NDIM, &
  MPI_INTEGER, 0, COMM, mpiErr)
CALL WRITE_DATASET_2DINT(fileId, "domainCoords", 0, &
  INT8( (/ NDIM, numRanks /) ), domainIndices, COMM)
IF (rank == 0) THEN
  DEALLOCATE(domainIndices)
ENDIF

! Get domain bounds
IF (rank == 0) THEN
  ALLOCATE(domainLowerBounds(NDIM, numRanks))
  ALLOCATE(domainUpperBounds(NDIM, numRanks))
ENDIF
CALL MPI_GATHER(domainLbs, NDIM, MPI_DOUBLE_PRECISION, domainLowerBounds, &
  NDIM, MPI_DOUBLE_PRECISION, 0, COMM, mpiErr)
CALL MPI_GATHER(domainUbs, NDIM, MPI_DOUBLE_PRECISION, domainUpperBounds, &
  NDIM, MPI_DOUBLE_PRECISION, 0, COMM, mpiErr)
CALL WRITE_DATASET_2DDBL(fileId, "domainLowerBounds", 0, &
  INT8( (/ NDIM, numRanks /) ), domainLowerBounds, COMM, &
  "xminp, yminp, zminp for each rank")
CALL WRITE_DATASET_2DDBL(fileId, "domainUpperBounds", 0, &
  INT8( (/ NDIM, numRanks /) ), domainUpperBounds, COMM, &
  "xmaxp, ymaxp, zmaxp for each rank")
IF (rank == 0) THEN
  DEALLOCATE(domainLowerBounds, domainUpperBounds)
ENDIF

! Get rank hostnames
hostLen = LEN(trim(adjustl(hostname)))
CALL MPI_ALLREDUCE(hostLen, maxHostLen, 1, MPI_INTEGER, MPI_MAX, COMM, mpiErr)
IF (rank == 0) THEN
  ALLOCATE(allHostNames(maxHostLen*numRanks))
ENDIF

CALL MPI_GATHER(trim(adjustl(hostname)), maxHostLen, MPI_CHARACTER, allHostNames, &
  maxHostLen, MPI_CHARACTER, 0, COMM, mpiErr)

CALL WRITE_DATASET_STR_LIST(fileId, "rankNodes", 0, &
  INT8( (/ maxHostLen, numRanks /) ), allHostNames, COMM)

IF (rank == 0) THEN
  DEALLOCATE(allHostNames)
ENDIF

CALL h5fclose_f(fileId, h5err)

END SUBROUTINE DUMP_SIM_INFO

!***********************************************************************

END MODULE
