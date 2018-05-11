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

USE MOD_INPUT
USE MOD_MPI
USE MOD_INTERP

IMPLICIT NONE

LOGICAL, PARAMETER, PUBLIC :: USE_HDF5 = .FALSE.

PRIVATE

PUBLIC :: DUMP_GLOBAL_FIELD_1D ! Write a field on one rank to disk 
PUBLIC :: DUMP_GLOBAL_FIELD_3D ! Write a field on one rank to disk 
PUBLIC :: DUMP_FIELD ! write field to disk
PUBLIC :: SAVE_FIELDS ! Write fields to disk for later restoration
PUBLIC :: RESTORE_FIELDS ! Read fields from disk collectively
PUBLIC :: DUMP_PARTICLES ! Write all the particle data to disk
PUBLIC :: SAVE_PARTICLES ! Write particle data to disk for later restoration
PUBLIC :: RESTORE_PARTICLES ! Read particle data from disk
PUBLIC :: DUMP_TRACKED_PARTICLES_SER ! write table of tracked particles to disk
PUBLIC :: GET_TRACK_SEQ_NUM ! get the sequence number for saving tracking data
PUBLIC :: DUMP_SIM_INFO ! write helpful info

 CONTAINS

!***********************************************************************
! Subroutine DUMP_GLOBAL_FIELD_1D
! This subroutine writes a 1D array (known entirely to one rank) to disk
! This should be called by all ranks, even though only one rank might
! know the field and write it.
!
! INPUT:
! 
! - fileBaseName the base name of the file, 
!     to which "_it.h5" will be appended (see argument "it", below)
! - Fshape = SHAPE(F)
!   - is INT8 for compatibility with hdf5
! - F: (entire) array of field values
! - rankToWrite: the (0-indexed) rank that knows F and will write it to disk.
! - it: Timestep number
! - COMM: mpi comm
! Following are ignored for text files:
! - datasetname the name of the field dataset
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
INTEGER*8                                :: Fshape(FDIM)
! following must have rank=FDIM
DOUBLE PRECISION, DIMENSION(:)           :: F
INTEGER                                  :: rankToWrite
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xgp, ygp, zgp
CHARACTER(len=*), OPTIONAL               :: xLabel, yLabel, zLabel
INTEGER                                  :: it
INTEGER                                  :: COMM

! LOCAL
CHARACTER(len=10)                        :: cit
INTEGER*8                                :: i,j,k
INTEGER                                  :: rank, mpiErr
!***********************************************************************

 write(cit,'(i10)') it
 cit = adjustl(cit)
  
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
IF (rankToWrite == rank) THEN
  OPEN(9,FILE=fileBaseName // trim(cit) // ".dat")
  DO i=1,Fshape(1)
    WRITE(9,"(1" // FMT_DBL // ")") (F(i))
  ENDDO
  CLOSE(9)
ENDIF

END SUBROUTINE DUMP_GLOBAL_FIELD_1D

!***********************************************************************
! Subroutine DUMP_GLOBAL_FIELD_3D
! This subroutine writes a 3D array (known entirely to one rank) to disk
! This should be called by all ranks, even though only one rank might
! know the field and write it.
!
! INPUT:
! 
! - fileBaseName the base name of the file, 
!     to which "_it.h5" will be appended (see argument "it", below)
! - Fshape = SHAPE(F)
!   - is INT8 for compatibility with hdf5
! - F: (entire) array of field values
! - rankToWrite: the (0-indexed) rank that knows F and will write it to disk.
! - it: Timestep number
! - transp: whether to transpose the axes
!   N.B.  transp = .FALSE. means preserving the memory layout
!     which may make it appear that the fortran and disk arrays are
!     actually transposed.
!   E.g., given F(x,y,z) and transp=.FALSE.
!    the disk array will place F(:,y,z) on the same row.
!    In both fortran and disk, F(x,y,z) and F(x+dx,y,z) will be adjacent
!    in memory.
!   If transp=.TRUE., then
!    the disk array will place F(x,y,:) on the same row.
!    F(x,y,z) and F(x,y,z+dz) will be adjacent on disk.
! - COMM: mpi comm
! Following are ignored for text files:
! - datasetname the name of the field dataset
! - xgp,ygp,zpg: coordinates of field values 
! - xLabel, yLabel, zLabel: names of coordinates, like "x" or "phi"
!
!***********************************************************************

SUBROUTINE DUMP_GLOBAL_FIELD_3D(fileBaseName, datasetname, Fshape, F, &
  rankToWrite, it, transp, COMM, xgp, xLabel, ygp, yLabel, zgp, zLabel)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER                       :: FDIM = 3

! INPUT
CHARACTER(len=*)                         :: fileBaseName
CHARACTER(len=*)                         :: datasetname
INTEGER*8, DIMENSION(FDIM)               :: Fshape
! following must have rank=FDIM
DOUBLE PRECISION, DIMENSION(:,:,:)       :: F
INTEGER                                  :: rankToWrite
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xgp, ygp, zgp
CHARACTER(len=*), OPTIONAL               :: xLabel, yLabel, zLabel
INTEGER                                  :: it
LOGICAL                                  :: transp
INTEGER                                  :: COMM

! LOCAL
CHARACTER(len=10)                        :: cit,CF
INTEGER*8                                :: i,j,k
INTEGER                                  :: rank, mpiErr
!***********************************************************************

 write(cit,'(i10)') it
 cit = adjustl(cit)
 WRITE(CF,'(i10)') Fshape(3)
 CF=adjustl(CF)
 
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
IF (rankToWrite == rank) THEN
  OPEN(9,FILE=fileBaseName // trim(cit) // ".dat")

  IF (transp .EQV. .TRUE.) THEN

    DO i=1,Fshape(1)
      DO j=1,Fshape(2)
        WRITE(9,"(" // trim(CF) // FMT_DBL // ")") (F(i,j,k),k=1,Fshape(3))
      ENDDO
    ENDDO

  ELSE

    DO k=1,Fshape(3)
      DO j=1,Fshape(2)
        WRITE(9,"(" // trim(CF) // FMT_DBL // ")") (F(i,j,k),i=1,Fshape(1))
      ENDDO
    ENDDO

  ENDIF

  CLOSE(9)
ENDIF

END SUBROUTINE DUMP_GLOBAL_FIELD_3D

!***********************************************************************
! Subroutine DUMP_FIELD
! This subroutine writes an array to disk. 
!
! INPUT:
! 
! - fileBaseName the base name of the file, 
!     to which "_DUMPNUM.dat" will be appended
!    (or _idRANK_DUMPNUM.dat if each rank dumps separately)
! - datasetName the name of the field dataset (ignored)
! - F: array of (local domain) field values
! - xgp,ygp: coordinates of (local domain) field values
! - it: Timestep number
! - domIndex: NDIM-dimensional index of this rank in the cartesian topology
! - COMM: mpi comm
!
! OUTPUT:
! - ierr: error code
!
!***********************************************************************

SUBROUTINE DUMP_FIELD(fileBaseName, datasetName, F,it, COMM, &
  xgp, xLabel, ygp, yLabel, zgp, zLabel)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER                       :: FDIM = NDIM
! INPUT
CHARACTER(len=*)                         :: fileBaseName
CHARACTER(len=*)                         :: datasetName
! following must have rank=FDIM
DOUBLE PRECISION, DIMENSION(:,:)         :: F
DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xgp, ygp, zgp
CHARACTER(len=*), OPTIONAL               :: xLabel, yLabel, zLabel
INTEGER                                  :: it
INTEGER                                  :: COMM

! LOCAL
CHARACTER(len=10)                        :: cit, cid, CF

! global field, in case one rank dumps all
DOUBLE PRECISION, ALLOCATABLE         :: Ft(:,:)

! shape of local and global field array
INTEGER, DIMENSION(FDIM)   :: localShape
INTEGER, DIMENSION(FDIM)   :: globalShape 
! offset of local array origin within global array
INTEGER, DIMENSION(FDIM)   :: offset

INTEGER        :: mpiErr, numRanks, rank

INTEGER        :: d, i, j
!***********************************************************************

localShape = SHAPE(F) 

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)

 write(cit,'(i10)') it
 cit = adjustl(cit)

IF (DUMP_FIELDS_PARALLEL.EQV..TRUE.) THEN

 WRITE(CF,'(i10)') SIZE(F,1)
 CF=adjustl(CF)

  write(cid,'(i10)') rank
  cid = adjustl(cid)
  OPEN(9,FILE=fileBaseName // "_id" // trim(cid) // '_' // trim(cit) // ".dat")
  DO j=1,SIZE(F,2)
  WRITE(9,"(" // trim(CF) // FMT_DBL // ")") (F(i,j),i=1,SIZE(F,1))
  ENDDO
  CLOSE(9)

ELSE
!===================================================================
! DUMPING OF THE DATA BY ONE PROCESSOR
!===================================================================

  CALL FILL_GLOBAL_FIELD_SHAPE(localShape, globalShape, offset, COMM)
  ALLOCATE(Ft(1:globalShape(1),1:globalShape(2)))
  ! Get all field values on rank 0
  CALL GATHER_GLOBAL_FIELD(F, Ft, offset, 0, COMM)

   WRITE(CF,'(i10)') SIZE(Ft,1)
   CF=adjustl(CF)
  
  IF (rank == 0) THEN
    OPEN(9,FILE=fileBaseName // trim(cit) // ".dat")
    DO j=1,SIZE(Ft,2)
    WRITE(9,"(" // trim(CF) // FMT_DBL // ")") (Ft(i,j),i=1,SIZE(Ft,1))
    ENDDO
    CLOSE(9)
  ENDIF
  
  DEALLOCATE(Ft)

ENDIF

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
DOUBLE PRECISION, DIMENSION(:,:)         :: Bxg,Byg,Bzg,Exg,Eyg,Ezg
CHARACTER(LEN=10)                        :: cit,cid,cx
INTEGER                                  :: it
INTEGER                                  :: ixp,iyp
INTEGER                                  :: nx, ny
!***********************************************************************

! Convert the integer it and id into a string cit and cid
 write(cit,'(i10)') it
 write(cid,'(i10)') id
 
! This step left-justify the string
 cit=adjustl(cit)
 cid=adjustl(cid)
 
nx = SIZE(Bxg, DIM=1)
ny = SIZE(Byg, DIM=2)

 WRITE(CX,'(i10)') nx
 CX=adjustl(CX)

OPEN(9,FILE="./data_restore/Bx_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,ny
WRITE(9,"(" // trim(CX) // FMT_DBL // ")") (Bxg(ixp,iyp),ixp=1,nx)
ENDDO
CLOSE(9)

OPEN(9,FILE="./data_restore/By_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,ny
WRITE(9,"(" // trim(CX) // FMT_DBL // ")") (Byg(ixp,iyp),ixp=1,nx)
ENDDO
CLOSE(9)

OPEN(9,FILE="./data_restore/Bz_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,ny
WRITE(9,"(" // trim(CX) // FMT_DBL // ")") (Bzg(ixp,iyp),ixp=1,nx)
ENDDO
CLOSE(9)

OPEN(9,FILE="./data_restore/Ex_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,ny
WRITE(9,"(" // trim(CX) // FMT_DBL // ")") (Exg(ixp,iyp),ixp=1,nx)
ENDDO
CLOSE(9)

OPEN(9,FILE="./data_restore/Ey_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,ny
WRITE(9,"(" // trim(CX) // FMT_DBL // ")") (Eyg(ixp,iyp),ixp=1,nx)
ENDDO
CLOSE(9)

OPEN(9,FILE="./data_restore/Ez_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,ny
WRITE(9,"(" // trim(CX) // FMT_DBL // ")") (Ezg(ixp,iyp),ixp=1,nx)
ENDDO
CLOSE(9)

END SUBROUTINE SAVE_FIELDS

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

INTEGER                                  :: id,COMM,ierr
INTEGER, DIMENSION(2)                    :: coords
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bxg,Byg,Bzg,Exg,Eyg,Ezg
CHARACTER(LEN=10)                        :: cit,cid,CX
INTEGER                                  :: it
INTEGER                                  :: ixp,iyp
!***********************************************************************

! Convert the integer it and id into a string cit and cid
 write(cit,'(i10)') it
 write(cid,'(i10)') id

 ! This step left-justify the string
 cit=adjustl(cit)
 cid=adjustl(cid)

 WRITE(CX,'(i10)') NXP
 CX=adjustl(CX)
 
OPEN(10,FILE="./data_restore/Bx_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,NYP
READ(10,"(" // trim(CX) // FMT_DBL // ")") (Bxg(ixp,iyp),ixp=1,NXP)
ENDDO
CLOSE(10)

OPEN(10,FILE="./data_restore/By_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,NYP
READ(10,"(" // trim(CX) // FMT_DBL // ")") (Byg(ixp,iyp),ixp=1,NXP)
ENDDO
CLOSE(10)

OPEN(10,FILE="./data_restore/Bz_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,NYP
READ(10,"(" // trim(CX) // FMT_DBL // ")") (Bzg(ixp,iyp),ixp=1,NXP)
ENDDO
CLOSE(10)

OPEN(10,FILE="./data_restore/Ex_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,NYP
READ(10,"(" // trim(CX) // FMT_DBL // ")") (Exg(ixp,iyp),ixp=1,NXP)
ENDDO
CLOSE(10)

OPEN(10,FILE="./data_restore/Ey_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,NYP
READ(10,"(" // trim(CX) // FMT_DBL // ")") (Eyg(ixp,iyp),ixp=1,NXP)
ENDDO
CLOSE(10)

OPEN(10,FILE="./data_restore/Ez_id"// trim(cid) // '_' // trim(cit) // ".dat")
DO iyp=1,NYP
READ(10,"(" // trim(CX) // FMT_DBL // ")") (Ezg(ixp,iyp),ixp=1,NXP)
ENDDO
CLOSE(10)

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
!
!***********************************************************************

SUBROUTINE DUMP_PARTICLES(pcl,tag,NPP,it,spec,sym,id,restore,COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: id,COMM
INTEGER                                  :: it
INTEGER*8                                :: NPP
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
INTEGER*8, DIMENSION(1:NPP)              :: tag
LOGICAL                                  :: restore
CHARACTER(LEN=10)                        :: spec,sym

INTEGER                                  :: i
INTEGER*8                                :: Nids,ip
DOUBLE PRECISION, ALLOCATABLE            :: pcl_print(:,:)
INTEGER*8, ALLOCATABLE                   :: tag_print(:)
CHARACTER(LEN=10)                        :: cit,cid
INTEGER, ALLOCATABLE                     :: Nid(:),displ(:)
INTEGER                                  :: ierr
!***********************************************************************

! Convert the integer it and id into a string cit and cid
WRITE(cit,'(i10)') it
WRITE(cid,'(i10)') id

! This step left-justify the string
 cit=adjustl(cit)
 cid=adjustl(cid)
 spec=adjustl(spec)
 sym=adjustl(sym)

!===================================================================
! PARALLEL DUMPING OF THE DATA
!===================================================================

IF (DUMP_PARTICLES_PARALLEL.EQV..TRUE.) THEN

! Particle distribution function
OPEN(9,FILE="./data/particles/" // trim(spec) // "_id"  // trim(cid) // "_"&
            // trim(sym) // trim(cit) // ".dat")
WRITE(9,"(1" // FMT_INT // ")") NPP
DO ip=1,NPP
WRITE(9,"(7" // FMT_DBL //",1" // FMT_INT // ")") pcl(1,ip),pcl(2,ip),pcl(3,ip),&
                                                  pcl(4,ip),pcl(5,ip),pcl(6,ip),&
                                                  pcl(7,ip),tag(ip)
ENDDO
CLOSE(9)

ELSE

!===================================================================
! DUMPING OF THE DATA BY ONE PROCESSOR
!===================================================================

ALLOCATE(Nid(1:NPX*NPY),displ(1:NPX*NPY))

Nid=0
displ=0

Nids=NPP

CALL MPI_ALLGATHER(Nids,1,MPI_INTEGER,Nid,1,MPI_INTEGER,COMM,ierr)

displ(1)=0
DO i=2,NPX*NPY
displ(i)=displ(i-1)+Nid(i-1)
ENDDO

ALLOCATE(pcl_print(1:7,1:NP))
pcl_print=0d0

CALL MPI_GATHERV(pcl(1,:),NPP,MPI_DOUBLE_PRECISION,pcl_print(1,:),Nid,displ,&
                 MPI_DOUBLE_PRECISION,0,COMM,ierr)

CALL MPI_GATHERV(pcl(2,:),NPP,MPI_DOUBLE_PRECISION,pcl_print(2,:),Nid,displ,&
                 MPI_DOUBLE_PRECISION,0,COMM,ierr)

CALL MPI_GATHERV(pcl(3,:),NPP,MPI_DOUBLE_PRECISION,pcl_print(3,:),Nid,displ,&
                 MPI_DOUBLE_PRECISION,0,COMM,ierr)

CALL MPI_GATHERV(pcl(4,:),NPP,MPI_DOUBLE_PRECISION,pcl_print(4,:),Nid,displ,&
                 MPI_DOUBLE_PRECISION,0,COMM,ierr)

CALL MPI_GATHERV(pcl(5,:),NPP,MPI_DOUBLE_PRECISION,pcl_print(5,:),Nid,displ,&
                 MPI_DOUBLE_PRECISION,0,COMM,ierr)

CALL MPI_GATHERV(pcl(6,:),NPP,MPI_DOUBLE_PRECISION,pcl_print(6,:),Nid,displ,&
                 MPI_DOUBLE_PRECISION,0,COMM,ierr)

CALL MPI_GATHERV(pcl(7,:),NPP,MPI_DOUBLE_PRECISION,pcl_print(7,:),Nid,displ,&
                 MPI_DOUBLE_PRECISION,0,COMM,ierr)

ALLOCATE(tag_print(1:NP))
tag_print=0
                 
CALL MPI_GATHERV(tag,NPP,MPI_INTEGER8,tag_print,Nid,displ,MPI_INTEGER8,0,COMM,ierr)
                 
IF (id.EQ.0) THEN
OPEN(9,FILE="./data/particles/" // trim(spec) // "_"  // trim(sym) // trim(cit) //&
            ".dat")
DO ip=1,NP
WRITE(9,"(7" // FMT_DBL //",1" // FMT_INT // ")") pcl_print(1,ip),pcl_print(2,ip),&
                                                  pcl_print(3,ip),pcl_print(4,ip),&
                                                  pcl_print(5,ip),pcl_print(6,ip),&
                                                  pcl_print(7,ip),tag_print(ip)
ENDDO
CLOSE(9)
END IF

DEALLOCATE(Nid,displ)
DEALLOCATE(pcl_print)
DEALLOCATE(tag_print)

END IF

END SUBROUTINE DUMP_PARTICLES

!***********************************************************************
! Subroutine SAVE_PARTICLES
! This subroutine saves in parallel all the particle data to disk for future
! restoration of the simulation.
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

SUBROUTINE SAVE_PARTICLES(pcl,pcl_data,tag,NPP,it,spec,sym,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER                                  :: id,COMM,ierr
INTEGER                                  :: it
INTEGER*8                                :: NPP,ip
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)   :: pcl
DOUBLE PRECISION, DIMENSION(1:4,1:NPP)   :: pcl_data
INTEGER*8, DIMENSION(1:NPP)              :: tag
CHARACTER(LEN=10)                        :: cit,cid,spec,sym
!***********************************************************************

! Convert the integer it and id into a string cit and cid
WRITE(cit,'(i10)') it
WRITE(cid,'(i10)') id

! This step left-justify the string
 cit=adjustl(cit)
 cid=adjustl(cid)
 spec=adjustl(spec)
 sym=adjustl(sym)

! Particle distribution function
OPEN(9,FILE="./data_restore/" // trim(spec) // "_id"  // trim(cid) // "_"&
            // trim(sym) // trim(cit) // ".dat")
WRITE(9,"(1" // FMT_INT // ")") NPP
DO ip=1,NPP
WRITE(9,"(7" // FMT_DBL //",1" // FMT_INT // ")") pcl(1,ip),pcl(2,ip),pcl(3,ip),&
                                                  pcl(4,ip),pcl(5,ip),pcl(6,ip),&
                                                  pcl(7,ip),tag(ip)
ENDDO
CLOSE(9)

! Particle data distribution function
OPEN(9,FILE="./data_restore/" // trim(spec) // "_data_id"  // trim(cid) // "_"&
            // trim(sym) // trim(cit) // ".dat")
DO ip=1,NPP
WRITE(9,"(4" // FMT_DBL // ")") pcl_data(1,ip),pcl_data(2,ip),pcl_data(3,ip),&
                                pcl_data(4,ip)
ENDDO
CLOSE(9)

END SUBROUTINE SAVE_PARTICLES

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
CHARACTER(LEN=10)                         :: cit,cid,spec,sym
INTEGER*8                                 :: ip
!***********************************************************************

! Convert the integer it and id into a string cit and cid
 write(cit,'(i10)') it
 write(cid,'(i10)') id
 
! This step left-justify the string
 cit=adjustl(cit)
 cid=adjustl(cid)

! Loading the number of particles per domain
OPEN(10,FILE="./data_restore/" // trim(spec) // "_id"  // trim(cid) // "_"&
            // trim(sym) // trim(cit) // ".dat")
READ(10,"(1" // FMT_INT // ")") NPP

ALLOCATE(pcl(1:7,1:NPP))
ALLOCATE(tag(1:NPP))

DO ip=1,NPP
READ(10,"(7" // FMT_DBL //",1" // FMT_INT // ")") pcl(1,ip),pcl(2,ip),pcl(3,ip),&
                                                  pcl(4,ip),pcl(5,ip),pcl(6,ip),&
                                                  pcl(7,ip),tag(ip)
ENDDO
CLOSE(10)

ALLOCATE(pcl_data(1:4,1:NPP))

! Loading the particle distribution function
OPEN(10,FILE="./data_restore/" // trim(spec) // "_data_id"  // trim(cid) // "_"&
            // trim(sym) // trim(cit) // ".dat")
DO ip=1,NPP
READ(10,"(4" // FMT_DBL // ")") pcl_data(1,ip),pcl_data(2,ip),pcl_data(3,ip),&
                                pcl_data(4,ip)
ENDDO
CLOSE(10)

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
  INQUIRE(FILE=fileBaseName // "_cid0_" // trim(seqStr) // ".txt", &
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
! - seqNum: a sequence number for the file, in case we want to 
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
INTEGER, DIMENSION(:)                    :: timeSteps
DOUBLE PRECISION, DIMENSION(:)           :: times
INTEGER                                  :: numPtcls, collectorNum, seqNum
DOUBLE PRECISION, DIMENSION(:,:)         :: ptclTable
! LOCAL
CHARACTER(len=256)                       :: fileName
CHARACTER(len=10)                        :: collectorStr, seqStr
INTEGER :: i,j
!***********************************************************************

write(collectorStr,'(i10)') collectorNum
write(seqStr,'(i10)') seqNum
collectorStr = adjustl(collectorStr)
seqStr = adjustl(seqStr)

! Make sure following filename agrees with that in GET_TRACK_SEQ_NUM
fileName = fileBaseName // "_cid" // trim(collectorStr) // "_" // trim(seqStr) &
  // ".txt"

OPEN(9, FILE=fileName, POSITION='APPEND')
DO i = 1, SIZE(ptclTable,DIM=2)
  WRITE(9,*) (ptclTable(j,i),j=1,SIZE(ptclTable,DIM=1))
ENDDO
CLOSE(9)

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
!***********************************************************************
! Not implemented for text dumping

END SUBROUTINE DUMP_SIM_INFO

END MODULE
