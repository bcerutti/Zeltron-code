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

MODULE MOD_MPI

USE MOD_INPUT

IMPLICIT NONE

! create a global variable to store the mpi rank, to be used for
! debugging output when a function doesn't have immediate access to this
! (except for debugging, this shouldn't be used, global variables being
! deprecated)
INTEGER, PUBLIC :: globalMpiRank

PRIVATE
!PUBLIC

PUBLIC :: CUMSUM ! Form the cumulative sum (running integral) of an array
PUBLIC :: FILL_DOMAIN_INDEX ! Get the NDIM-dimensional domain index in cart. topology
PUBLIC :: FILL_DOMAIN_DECOMP ! Get the NDIM-dimensional domain decomp, with index 
PUBLIC :: ELEMENTS_PER_DOMAIN ! Find the number elements per domain in each 1D direction
PUBLIC :: FILL_GLOBAL_FIELD_SHAPE ! Get the global shape of a field from local shapes
PUBLIC :: GATHER_GLOBAL_FIELD ! Gather the global field on one mpi rank
PUBLIC :: FILL_RANK_ORDER ! Get the #ranks/node and the rank-node ordering
PUBLIC :: ASSERT_DOUBLE ! verify (at compile time) that a variable is a double

INTERFACE CUMSUM
  MODULE PROCEDURE CUMSUM_INT, CUMSUM_INT8  
END INTERFACE 

INTERFACE ELEMENTS_PER_DOMAIN
  MODULE PROCEDURE ELEMENTS_PER_DOMAIN_4, ELEMENTS_PER_DOMAIN_8  
END INTERFACE 

INTERFACE FILL_GLOBAL_FIELD_SHAPE
  MODULE PROCEDURE FILL_GLOBAL_FIELD_SHAPE_4, FILL_GLOBAL_FIELD_SHAPE_8  
END INTERFACE 

 CONTAINS

!***********************************************************************
! Subroutine CUMSUM
! Find the cumulative sum (like a running integral) of an array,
! such that the first element is zero,
! and also find the sum of all elements of raIn
!
! INPUT:
! 
! - raIn: a 1D array
! OUTPUT: 
! - raOut: raOut(n) = Sum_{m=1}^{n-1} raIn(m)  (undefined for n > SIZE(raIn))
!   raOut is the same size as raIn (or longer)
! - total: total = SUM(raIn)
!
!***********************************************************************

SUBROUTINE CUMSUM_INT(raIn, raOut, total)

IMPLICIT NONE

! INPUT
INTEGER, DIMENSION(:), INTENT(IN)  :: raIn
INTEGER, DIMENSION(:), INTENT(OUT) :: raOut
INTEGER, INTENT(OUT)               :: total
! LOCAL
INTEGER :: i
!***********************************************************************

total = 0
IF (SIZE(raIn) > 0) THEN
  raOut(1) = 0
  DO i = 2, SIZE(raIn)
    raOut(i) = raOut(i-1) + raIn(i-1)
  ENDDO
  total = raOut(SIZE(raIn)) + raIn(SIZE(raIn))
ENDIF

END SUBROUTINE CUMSUM_INT

!***********************************************************************
! Subroutine CUMSUM
! Find the cumulative sum (like a running integral) of an array,
! such that the first element is zero,
! and also find the sum of all elements of raIn
!
! INPUT:
! 
! - raIn: a 1D array
! OUTPUT: 
! - raOut: raOut(n) = Sum_{m=1}^{n-1} raIn(m)  (undefined for n > SIZE(raIn))
!   raOut is the same size as raIn (or longer)
! - total: total = SUM(raIn)
!
!***********************************************************************

SUBROUTINE CUMSUM_INT8(raIn, raOut, total)

IMPLICIT NONE

! INPUT
INTEGER*8, DIMENSION(:), INTENT(IN)  :: raIn
INTEGER*8, DIMENSION(:), INTENT(OUT) :: raOut
INTEGER*8, INTENT(OUT)               :: total
! LOCAL
INTEGER*8 :: i
!***********************************************************************

total = 0
IF (SIZE(raIn) > 0) THEN
  raOut(1) = 0
  DO i = 2, SIZE(raIn)
    raOut(i) = raOut(i-1) + raIn(i-1)
  ENDDO
  total = raOut(SIZE(raIn)) + raIn(SIZE(raIn))
ENDIF

END SUBROUTINE CUMSUM_INT8

!***********************************************************************
! Subroutine FILL_DOMAIN_DECOMP
!
! INPUT:
! 
! - COMM: mpi comm
!
! OUTPUT: 
! - ierr: error code
! - domainDecomp the cartesian topology; e.g., in 2D, returns (NXP, NYP)
!   if there are NXP domains in the x-direction, and NYP in y.
! - domainIndex the NDIM-dimensional domain index in cartesian topology
!   (0-indexed)
!
!***********************************************************************

SUBROUTINE FILL_DOMAIN_DECOMP(COMM, err, domainDecomp, domainIndex)

IMPLICIT NONE

INTEGER :: COMM
INTEGER, INTENT(OUT):: err
INTEGER, INTENT(OUT) :: domainDecomp(NDIM)
INTEGER, INTENT(OUT) :: domainIndex(NDIM)
LOGICAL :: periods(NDIM)
INTEGER :: rank
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, err)

CALL MPI_CART_GET(COMM, NDIM, domainDecomp, periods, domainIndex, err)

END SUBROUTINE FILL_DOMAIN_DECOMP

!***********************************************************************
! Subroutine FILL_DOMAIN_INDEX
!
! INPUT:
! 
! - COMM: mpi comm
!
! OUTPUT: 
! - ierr: error code
! - domainIndex the NDIM-dimensional domain index in cartesian topology
!   (0-indexed)
!
!***********************************************************************

SUBROUTINE FILL_DOMAIN_INDEX(COMM, err, domainIndex)

IMPLICIT NONE

INTEGER :: COMM
INTEGER, INTENT(OUT):: err
INTEGER, INTENT(OUT):: domainIndex(NDIM)
INTEGER :: dims(NDIM)
!***********************************************************************

CALL FILL_DOMAIN_DECOMP(COMM, err, dims, domainIndex)

END SUBROUTINE FILL_DOMAIN_INDEX

!***********************************************************************
! Subroutine ELEMENTS_PER_DOMAIN
! This subroutine calculates the number of elements (of a grid), in
! each direction, belonging to each MPI domain.  It assumes a rectangular
! decomposition, so that all domains at the same x-location have the
! same number of elements in the x-direction, etc.
!
! INPUT:
! 
! - localRaShape: the shape of the array of values belonging to this domain
! - domIndex: 0-indexed coordinates of this domain, e.g., (/0,1/)
! - COMM: mpi comm
!
! OUTPUT: the number of elements/proc in each direction
! - xElems: an array of size decompShape(1) 
! - yElems: an array of size decompShape(2)
! E.g., xElems(domIndex(1)) = localRaShape(1)
!       yElems(domIndex(2)) = localRaShape(2)
!
!***********************************************************************

SUBROUTINE ELEMENTS_PER_DOMAIN_8(localRaShape, COMM, xElems, yElems, zElems)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
INTEGER(8), DIMENSION(NDIM)    :: localRaShape
INTEGER                        :: COMM
! OUTPUT
INTEGER, INTENT(OUT)           :: xElems(NPX), yElems(NPY)
INTEGER, INTENT(OUT), OPTIONAL :: zElems(NPZ)
! LOCAL
INTEGER, DIMENSION(NDIM)       :: domIndex
INTEGER                        :: mpiErr, rank
INTEGER                        :: d
INTEGER(8), DIMENSION(NDIM)    :: raShape
INTEGER :: xElemsTmp(NPX), yElemsTmp(NPY), zElemsTmp(NPZ)
!***********************************************************************

xElemsTmp = 0
yElemsTmp = 0
zElemsTmp = 0

CALL FILL_DOMAIN_INDEX(COMM, mpiErr, domIndex)
xElemsTmp(domIndex(1)+1) = localRaShape(1)
yElemsTmp(domIndex(2)+1) = localRaShape(2)

CALL MPI_ALLREDUCE(xElemsTmp, xElems, NPX, MPI_INTEGER, MPI_SUM, COMM, mpiErr)
CALL MPI_ALLREDUCE(yElemsTmp, yElems, NPY, MPI_INTEGER, MPI_SUM, COMM, mpiErr)

IF (NDIM == 2) THEN
  xElems = xElems / NPY
  yElems = yElems / NPX
ELSEIF (NDIM == 3) THEN
  xElems = xElems / (NPY*NPZ)
  yElems = yElems / (NPX*NPY)
  d = MIN(NDIM, 3)
  zElemsTmp(domIndex(d)+1) = localRaShape(d)
  CALL MPI_ALLREDUCE(zElemsTmp, zElems, NPY, MPI_INTEGER, MPI_SUM, COMM, mpiErr)
  ZElems = zElems / (NPX*NPY)
  ! set raShape(3), but use d to avoid error when NDIM=2
  raShape(d) = zElems(domIndex(d)+1)
ENDIF

! check validity
raShape(1) = xElems(domIndex(1)+1)
raShape(2) = yElems(domIndex(2)+1)
DO d = 1, NDIM
  IF (raShape(d) /= localRaShape(d)) THEN
    CALL MPI_COMM_RANK(COMM, rank, mpiErr)
    CALL MPI_COMM_RANK(COMM, rank, mpiErr)
    PRINT *, 'Error (in mod_io 8): rank ', rank, ' should have ', raShape(d), &
      ' elements (like all other procs in same line) in dir ', d, ' but has ', localRaShape(d)
  ENDIF
ENDDO

END SUBROUTINE 

!***********************************************************************

SUBROUTINE ELEMENTS_PER_DOMAIN_4(localRaShape, COMM, xElems, yElems, zElems)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
INTEGER(4), DIMENSION(NDIM)    :: localRaShape
INTEGER                        :: COMM
! OUTPUT
INTEGER, INTENT(OUT)           :: xElems(NPX), yElems(NPY)
INTEGER, INTENT(OUT), OPTIONAL :: zElems(NPZ)
! LOCAL
INTEGER, DIMENSION(NDIM)       :: domIndex
INTEGER                        :: mpiErr, rank
INTEGER                        :: d
INTEGER(4), DIMENSION(NDIM)    :: raShape
INTEGER :: xElemsTmp(NPX), yElemsTmp(NPY), zElemsTmp(NPZ)
!***********************************************************************

xElemsTmp = 0
yElemsTmp = 0
zElemsTmp = 0

CALL FILL_DOMAIN_INDEX(COMM, mpiErr, domIndex)
xElemsTmp(domIndex(1)+1) = localRaShape(1)
yElemsTmp(domIndex(2)+1) = localRaShape(2)

CALL MPI_ALLREDUCE(xElemsTmp, xElems, NPX, MPI_INTEGER, MPI_SUM, COMM, mpiErr)
CALL MPI_ALLREDUCE(yElemsTmp, yElems, NPY, MPI_INTEGER, MPI_SUM, COMM, mpiErr)

IF (NDIM == 2) THEN
  xElems = xElems / NPY
  yElems = yElems / NPX
ELSEIF (NDIM == 3) THEN
  xElems = xElems / (NPY*NPZ)
  yElems = yElems / (NPX*NPY)
  d = MIN(NDIM, 3)
  zElemsTmp(domIndex(d)+1) = localRaShape(d)
  CALL MPI_ALLREDUCE(zElemsTmp, zElems, NPY, MPI_INTEGER, MPI_SUM, COMM, mpiErr)
  ZElems = zElems / (NPX*NPY)
  ! set raShape(3), but use d to avoid error when NDIM=2
  raShape(d) = zElems(domIndex(d)+1)
ENDIF

! check validity
raShape(1) = xElems(domIndex(1)+1)
raShape(2) = yElems(domIndex(2)+1)
DO d = 1, NDIM
  IF (raShape(d) /= localRaShape(d)) THEN
    CALL MPI_COMM_RANK(COMM, rank, mpiErr)
    CALL MPI_COMM_RANK(COMM, rank, mpiErr)
    PRINT *, 'Error (in mod_io 4): rank ', rank, ' should have ', raShape(d), &
      ' elements (like all other procs in same line) in dir ', d, ' but has ', localRaShape(d)
  ENDIF
ENDDO

END SUBROUTINE

!***********************************************************************
! Subroutine FILL_GLOBAL_FIELD_SHAPE
!   Given the local field shape on each rank,
!   find the global field shape, and the offset, in global indices,
!   of the field value with local index (0,0)
!
! INPUT:
! - localShape
! - COMM: mpi comm
! OUTPUT:
! - globalShape
! - offset
!***********************************************************************

SUBROUTINE FILL_GLOBAL_FIELD_SHAPE_4(localShape, globalShape, offset, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
! following must have rank=NDIM
INTEGER*4, DIMENSION(NDIM)              :: localShape
INTEGER*4, DIMENSION(NDIM), INTENT(OUT) :: globalShape, offset
INTEGER                                       :: COMM

! LOCAL
INTEGER, DIMENSION(NDIM)                 :: domIndex
INTEGER, DIMENSION(NPX)                  :: xElemPerProc
INTEGER, DIMENSION(NPY)                  :: yElemPerProc
INTEGER, DIMENSION(NPZ)                  :: zElemPerProc
INTEGER        :: mpiErr, numRanks, rank
!***********************************************************************

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL FILL_DOMAIN_INDEX(COMM, mpiErr, domIndex)

! find shape of local array on each domain
IF (NDIM == 2) THEN
  CALL ELEMENTS_PER_DOMAIN(localShape, COMM, xElemPerProc, yElemPerProc)
ELSE
  CALL ELEMENTS_PER_DOMAIN(localShape, COMM, xElemPerProc, yElemPerProc, zElemPerProc)
ENDIF

globalShape(1) = SUM(xElemPerProc)
globalShape(2) = SUM(yElemPerProc)

offset(1) = SUM(xElemPerProc(:domIndex(1)))
offset(2) = SUM(yElemPerProc(:domIndex(2)))

IF (NDIM >= 3) THEN
  globalShape(MIN(NDIM,3)) = SUM(zElemPerProc)
  offset(MIN(NDIM,3)) = SUM(zElemPerProc(:domIndex(MIN(NDIM,3))))
ENDIF

END SUBROUTINE FILL_GLOBAL_FIELD_SHAPE_4

!***********************************************************************
! Subroutine FILL_GLOBAL_FIELD_SHAPE
!   Given the local field shape on each rank,
!   find the global field shape, and the offset, in global indices,
!   of the field value with local index (0,0)
!
! INPUT:
! - localShape
! - COMM: mpi comm
! OUTPUT:
! - globalShape
! - offset
!***********************************************************************

SUBROUTINE FILL_GLOBAL_FIELD_SHAPE_8(localShape, globalShape, offset, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
! following must have rank=NDIM
INTEGER*8, DIMENSION(NDIM)              :: localShape
INTEGER*8, DIMENSION(NDIM), INTENT(OUT) :: globalShape, offset
INTEGER                                       :: COMM

! LOCAL
INTEGER, DIMENSION(NDIM)                 :: domIndex
INTEGER, DIMENSION(NPX)                  :: xElemPerProc
INTEGER, DIMENSION(NPY)                  :: yElemPerProc
INTEGER, DIMENSION(NPZ)                  :: zElemPerProc
INTEGER        :: mpiErr, numRanks, rank
!***********************************************************************

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL FILL_DOMAIN_INDEX(COMM, mpiErr, domIndex)

! find shape of local array on each domain
IF (NDIM == 2) THEN
  CALL ELEMENTS_PER_DOMAIN(localShape, COMM, xElemPerProc, yElemPerProc)
ELSE
  CALL ELEMENTS_PER_DOMAIN(localShape, COMM, xElemPerProc, yElemPerProc, zElemPerProc)
ENDIF

globalShape(1) = SUM(xElemPerProc)
globalShape(2) = SUM(yElemPerProc)

offset(1) = SUM(xElemPerProc(:domIndex(1)))
offset(2) = SUM(yElemPerProc(:domIndex(2)))

IF (NDIM >= 3) THEN
  globalShape(MIN(NDIM,3)) = SUM(zElemPerProc)
  offset(MIN(NDIM,3)) = SUM(zElemPerProc(:domIndex(MIN(NDIM,3))))
ENDIF

END SUBROUTINE FILL_GLOBAL_FIELD_SHAPE_8

!***********************************************************************
! Subroutine GATHER_GLOBAL_FIELD
! Gather field values from all ranks to a global field on one rank.
!
! INPUT:
! 
! - F the local field values for this rank
! - Ft an array (already allocated) with the global field shape
!   - on rank = gatherRank, this will contain the global field
!   - on other ranks, this is undefined
! - offset an array giving the global index of F(0,0[,0])
! - COMM: mpi comm
!
! OUTPUT: the number of elements/proc in each direction
!
!***********************************************************************

SUBROUTINE GATHER_GLOBAL_FIELD(F, Ft, offset, gatherRank, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
! following must have rank=NDIM
DOUBLE PRECISION, DIMENSION(:,:,:)              :: F
DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: Ft
INTEGER, DIMENSION(NDIM)                      :: offset
INTEGER                                       :: gatherRank
INTEGER                                       :: COMM

! LOCAL
DOUBLE PRECISION, ALLOCATABLE            :: Ft2(:,:,:)
INTEGER, DIMENSION(NDIM)                 :: domIndex
INTEGER, DIMENSION(NPX)                  :: xElemPerProc
INTEGER, DIMENSION(NPY)                  :: yElemPerProc
INTEGER, DIMENSION(NPZ)                  :: zElemPerProc
CHARACTER(len=10)                        :: cit

! shape of local and global field array
INTEGER, DIMENSION(NDIM)   :: localShape
INTEGER, DIMENSION(NDIM)   :: globalShape 

INTEGER        :: mpiErr, numRanks, rank

INTEGER        :: d, i, j, k
!***********************************************************************

localShape = SHAPE(F)
globalShape = SHAPE(Ft)

ALLOCATE(Ft2(globalShape(1), globalShape(2), globalShape(3)))

Ft2 = 0.

DO k=1, localShape(3)
  DO j=1, localShape(2)
    DO i=1, localShape(1)
      Ft2(offset(1) + i, offset(2) + j, offset(3) + k) = F(i,j,k)
    ENDDO
  ENDDO
ENDDO

! :TODO: it would be better to use MPI_GATHERV, but it's a lot easier to use MPI_REDUCE
! Get all field values on rank 0
CALL MPI_REDUCE(Ft2,Ft,PRODUCT(globalShape),MPI_DOUBLE_PRECISION,MPI_SUM,&
  gatherRank,COMM,mpiErr)

DEALLOCATE(Ft2)

END SUBROUTINE GATHER_GLOBAL_FIELD

!***********************************************************************
! Subroutine FILL_RANK_ORDER
! Determine (by examing the hostname of the first two ranks) 
! whether the order is round-robin or not.
!
! INPUT:
! -COMM: mpi comm (should be global comm)
! 
! OUTPUT: Assuming M ranks per node, and N nodes
! - rankOrder 
!   =0 if rank r is on node (r%N)
!   =1 if rank r is on node floor(r/M)
!   If there's only one rank, or if ranksPerNode = 1, then rankOrder = 1
! - ranksPerNode: the number of ranksPerNode
!
!***********************************************************************

SUBROUTINE FILL_RANK_ORDER(COMM, rankOrder, ranksPerNode)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
INTEGER              :: COMM
INTEGER, INTENT(OUT) :: rankOrder, ranksPerNode
! LOCAL
CHARACTER*(MPI_MAX_PROCESSOR_NAME) :: hostName, hostNameRank0
INTEGER :: hostnameLen
INTEGER :: sameNameAsRank0
INTEGER :: rank, numRanks, mpiErr
INTEGER :: i, nextOnDifferentNode, nextOnSameNode
INTEGER :: tag
INTEGER, DIMENSION(:), ALLOCATABLE :: onFirstNode
!***********************************************************************

! Determine which other ranks are running on the same node as rank 0.
! From that, infer the number of ranks per node, and whether the ordering
! has adjacent ranks on the same node, or is round robin.

CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
IF (numRanks == 1) THEN
  rankOrder = 1
  ranksPerNode = 1
ELSE
  CALL MPI_COMM_RANK(COMM, rank, mpiErr)
  CALL MPI_GET_PROCESSOR_NAME(hostName, hostnameLen, mpiErr)
  hostName(hostnameLen+1:MPI_MAX_PROCESSOR_NAME) = ' '
  ! Broadcast hostName of rank 0 to everybody
  IF (rank == 0) THEN
    CALL MPI_BCAST(hostName, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, 0, &
      COMM, mpiErr)
    hostNameRank0 = hostName
    ALLOCATE(onFirstNode(numRanks))
  ELSE
    CALL MPI_BCAST(hostNameRank0, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, 0,&
      COMM, mpiErr)
  ENDIF
  sameNameAsRank0 = 0
  IF (hostName == hostNameRank0) sameNameAsRank0 = 1
  CALL MPI_GATHER(sameNameAsRank0, 1, MPI_INTEGER, onFirstNode, 1,&
    MPI_INTEGER, 0, COMM, mpiErr)
  IF (rank == 0) THEN
    ! Find the next rank that's on the same node, and the next that's
    !   on a different node (from rank 0).
    nextOnDifferentNode = numRanks
    DO i = 2, numRanks
      IF (onFirstNode(i) == 0) THEN
        nextOnDifferentNode = i - 1
        EXIT
      ENDIF
    ENDDO
    nextOnSameNode = numRanks
    DO i = 2, numRanks
      IF (onFirstNode(i) == 1) THEN
        nextOnSameNode = i - 1
        EXIT
      ENDIF
    ENDDO
    IF (nextOnSameNode == 1) THEN
      rankOrder = 1
      ranksPerNode = nextOnDifferentNode
    ENDIF
    IF (nextOnDifferentNode == 1) THEN
      rankOrder = 0
      ! Following assumes a constant ranksPerNode
      ranksPerNode = numRanks / nextOnSameNode 
    ENDIF
    IF (ranksPerNode == 1) THEN
      rankOrder = 1
    ENDIF
    DEALLOCATE(onFirstNode)
  ENDIF
  CALL MPI_BCAST(rankOrder, 1, MPI_INTEGER, 0, COMM, mpiErr)
  CALL MPI_BCAST(ranksPerNode, 1, MPI_INTEGER, 0, COMM, mpiErr)
ENDIF

END SUBROUTINE FILL_RANK_ORDER

!***********************************************************************
! Subroutine CHECK_DOUBLE
!
! This subroutine will cause a compilation error if the argument is not
! a double.
!***********************************************************************

SUBROUTINE ASSERT_DOUBLE(dbl)

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: dbl

END SUBROUTINE ASSERT_DOUBLE

!***********************************************************************

END MODULE

