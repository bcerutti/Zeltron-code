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

MODULE MOD_TRACK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module implements data structures and algorithms for tracking
! particles.  Its purpose is to write tracked-particle data to disk
! efficiently, assuming that, at any given time, a tracked particle may
! be on any MPI rank.
!
! For efficiency, it is (probably) important:
! (1) not to write to disk every time step;
! (2) to write to disk in large contiguous chunks;
! (3) (for analysis) to have particle data fairly conveniently ordered---
!     e.g., one shouldn't have to access a lot of random lines
!    (or worse, files) to get the trajectory of one particle;
! (4) not to have any rank writing super huge amounts of data;
! (5) (for runs on more than about 50,000 ranks) not to use
!     alltoallv (with the global communicator): for a toy test,
!       on 50,000 (Hopper) ranks, takes at least 3 seconds (tiny data), 
!       and for more data 6-7 seconds, and occasionally, depending on
!       load and the number of ranks collecting, 15 or so seconds.
!       on 50,000 (Edison) ranks, takes at least 1.5 seconds (tiny data), 
!       and for more data 4-7 seconds, and occasionally, depending on
!       load and the number of ranks collecting, 10 or so seconds.
!       on 100,000 (Hopper) ranks, takes at least 10 seconds (tiny data), 
!       and for more data 30 seconds, perhaps more.
!       on 100,000 (Edison) ranks, takes at least 5 seconds (tiny data), 
!       and for more data, 6-30 seconds, perhaps more.
!
! In addition, we worry about file corruption: we want to minimize the 
!   chances that tracking data from earlier in the simulation should be
!   lost if the simulation should be killed--e.g., possibly in the middle
!   of writing to disk.
! 
! = The approach taken here =
!
! Each rank chooses which particles to track randomly (if fewer particles
! than ranks, the number of ranks with 1 tracked particles is chosen 
! randomly).  These particles are given unique negative tags, which reflect
! which rank collects that particle (and where in the order of collected
! particles it is).  The absolute value of the tag may be the tag of a
! a non-tracked particle. The tracked-particle tags start at -1 and
! decrease by 1.
!
! Each rank locally collects data from any tracked particles on that rank for
! some number of steps (stepsPerTrackedPtclDump).  After that number of
! steps, data for a fixed set of particles is collected on a subset of
! ranks (called collectors)--e.g., a particles may move from rank to
! rank over time, but it's data will also be collected by the same rank.
! Each collector will sort and *append* its data to its own file.
! (Using separate files will hopefully reduce the chances of losing all
! data due to file corruption.)
!
! To avoid a global alltoallv, data will be gathered locally by each
! collector for nearby ranks.  The collectors will then use alltoallv
! to exchange data.  This is much faster on 100,000 ranks: 0.1 to 2 s
! on Hopper, with a few at 7-13 seconds, depending on the number of
! collectors (and maybe communication load); 
! 0.1-0.5 s on Edison, with a few taking 8-18 s.
!
! In case a file is corrupted, particles that begin on the same rank
! will be collected by many different collectors, so each collector has
! a sort of representative sample.
!
! To reduce the chance of corruption, it may be good not to write all
! the files at the same time.  Perhaps even files should be written one
! step later than odd files?
!
! Also, we may choose, in the middle of the simulation, to start a new
! file for each collector, that way avoiding corruption of the entire
! history.
!
! For each particle, the tag is stored in the next-to-last column, and
! the time step when it was collected is stored in the last column.
!
! A particle with tag t < 0 is collected by collector 
!   r_c = (|t|-1) % numCollectors
! Collector r_c collects particles with tags
!   |t| = 1 + r_c + (n-1)*numCollectors
!   for 
!  n = 1, 2, ..., (NSAMPLE/numCollectors + (r_c < nSAMPLE % numCollectors))
!  Thus
!     n = 1 + (|t| - 1 - r_c) / numCollectors
!     N.B. if "/" is integer division, then n = 1 + (|t|-1)/numCollectors. 
!  If calculating r_c and n (from t) proves to be too time-consuming, 
!    it may help to make numCollectors = 2^m.  In that case,
!    division by and modulo 2^m can be performed by simple bitshift
!    and bitwise-AND operations.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE MOD_INPUT
USE MOD_IO
USE MOD_MPI

IMPLICIT NONE

PRIVATE

PUBLIC :: RANDOM_SELECT ! Select random elements from a list
PUBLIC :: ASSIGN_TRACK_TAGS ! re-number tags of tracked particles
PUBLIC :: INIT_TRACKER_COMM ! Initialize communicators for tracking particles
PUBLIC :: END_TRACKER ! Deallocate communicators for tracking particles
PUBLIC :: TRACKER ! record particle data
! PUBLIC :: ADD_TRACKED_PARTICLE ! add a particle's data to the table
! ARRANGE_TABLE ! rearrange rows of a table into a given order
! SORT_AMONG_RANKS ! message particles from all ranks to nearby collector ranks
! GATHER_TABLE_ON_RANK0 ! combine tables from all ranks onto rank 0
! SORT_AMONG_COLLECTORS ! send particle data to the appropriate collector ranks
PUBLIC :: WRITE_TRACK_DATA ! write tracked-particle data to disk
PUBLIC :: INC_TRACK_SEQ_NUM ! call this after each checkpoint
PUBLIC :: TEST_TRACKER ! For debugging and testing

! A type that holds 2 tables (an array with fixed first dimension); 
! one of them is active.  If the active
! array needs to grow, the other is allocated at a larger size, and 
! the data copied over.
TYPE GrowingTable
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ra1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ra2
  INTEGER :: active ! either 1 or 2
  INTEGER :: numFilled ! number of filled rows of table
ENDTYPE GrowingTable
! Related to GrowingTable
!   Number of columns of data (including tag) and the timestep,
!   to store for tracked particles
INTEGER, PARAMETER :: nTrackCols = 13
!   The factor by which the trackPtcls arrays grow 
REAL, PARAMETER :: raGrowthFactor = 1.5

!====================================================================
! Groups/communicators for messaging tracked particle data
!====================================================================
! The list of ranks (in the global communicator) that collect and write
!   tracked-particle dataa
INTEGER, DIMENSION(:), ALLOCATABLE :: NUM_COLLECT_RANKS
! The group of collecting ranks
INTEGER :: collectGroup
! The group including this rank and one collector, of ranks that 
!   send their tracked-particle data to the one collector
INTEGER :: sendGroup
! The communicator among collecting ranks
INTEGER :: collectorComm
! The communicator among ranks that send to one collector
INTEGER :: sendComm
! Number of collector ranks
INTEGER :: numCollectors
! Whether this rank is a collector
LOGICAL :: iAmCollector

!====================================================================
! Tracked particle data, and data to keep track of when tracked particle
! data has been stored and written.
!====================================================================
! arrays to hold table of local tracked particles
! (we need 2 arrays because they have to grow and shrink)
TYPE(GrowingTable) eTrackPtcls
TYPE(GrowingTable) pTrackPtcls
! Number of particles to collect (finally) on this rank per step
INTEGER :: numCollectedPtclsPerStep
! Dump period for tracking data
INTEGER :: stepsPerTrackedPtclDump
! The first step at which tracked particles are recorded since the last dump
INTEGER :: firstTrackedStepSinceDump
! The last step for which tracked particle data was collected
INTEGER :: lastTrackedStep
! The number of steps on which tracked particle data was collected
!   since the last dump
INTEGER :: trackedStepsSinceDump

! The sequence number: every time a checkpoint is reached, this is 
! incremented, and we create a new file for tracked-particle data.
INTEGER :: trackSeqNum

 CONTAINS

!***********************************************************************
! FUNCTION GET_COLLECTOR_RANK_FOR_PTCL(tPtclData)
!
! This is a function that is passed as an argument to sorting functions;
! it returns the collector rank that should ultimately write data for
! a particle with tracking data tPtclData.
!
!***********************************************************************

FUNCTION GET_COLLECTOR_RANK_FOR_PTCL(tPtclData)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tPtclData
INTEGER :: GET_COLLECTOR_RANK_FOR_PTCL
!***********************************************************************

GET_COLLECTOR_RANK_FOR_PTCL = MOD(INT(-tPtclData(nTrackCols-1)) - 1, &
  numCollectors)

END FUNCTION GET_COLLECTOR_RANK_FOR_PTCL

!***********************************************************************
! SUBROUTINE RANDOM_SELECT
! Select elements randomly from a list.
!
! INPUT
! - numAvailable: the total number of elements available to be selected
!   (the number of elements to select is given by the length of the
!   selected array)
! OUTPUT
! - selected: an array of numToSelect integers, each between 1 and numAvailable
!     (inclusive), in strictly increasing order.
!   (The length of this array gives the number of elements to select.)
!
! This method is nice because it's simple and it yields exactly the 
!   number selected, in order (in 1 pass).
!
! N.B. there's a faster way to do this when selecting a small fraction;
! I haven't looked much, but one approach is in J. S. Vitter, 
! "Faster Methods for Random Sampling," Communications of the ACM 27 (7),
! 703--718 (1984).
!***********************************************************************

SUBROUTINE RANDOM_SELECT(numAvailable, selected)

IMPLICIT NONE

!INPUT
INTEGER :: numAvailable
INTEGER, DIMENSION(:), INTENT(INOUT) :: selected
!LOCAL
INTEGER :: numRemaining, numToSelect, n, numSelected
DOUBLE PRECISION :: numR, numS
DOUBLE PRECISION :: r, p
!***********************************************************************

numRemaining = numAvailable
numToSelect = SIZE(selected)

numR = numRemaining
numS = numToSelect

numSelected = 0
DO n = 1, numAvailable
  ! We need to select numToSelect from the numRemaining elements;
  ! the probability of selecting the next is (numToSelect/numRemaining)
  p = numS/numR
  CALL RANDOM_NUMBER(r)
  IF (r <= p) THEN ! select this
    numSelected = numSelected + 1
    selected(numSelected) = n
    numS = numS - 1.
    IF (numS == 0.) EXIT
  ENDIF
  numR = numR - 1.
ENDDO

END SUBROUTINE RANDOM_SELECT

!***********************************************************************
! SUBROUTINE INIT_TRACKER_COMM
! Set up communicators for sending tracked particle data to the rank
! that will dump it.  This needs to know the rank order and ranks per node
! so it doesn't put two collectors on one node.
!
! INPUT
! - COMM: mpi comm
! - rankOrder: 0 if round robin, rank r on node (r mod ranksPerNode)
!              1 if each node contains adjacent ranks
! - ranksPerNode: the number of ranks per node
! 
!***********************************************************************

SUBROUTINE INIT_TRACKER_COMM(COMM, rankOrder, ranksPerNode)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER :: COMM
INTEGER :: rankOrder, ranksPerNode

! Local
INTEGER :: rank, numRanks, mpiErr, numNodes
INTEGER :: tgtStepsPerDump
INTEGER*8 :: tgtPtclsPerCollector
DOUBLE PRECISION :: tgtMemPerCollector
INTEGER, DIMENSION(NDIM) :: domainDecomp, domainIndex
INTEGER :: i, j, jMax, nCols

INTEGER :: myCollector, collectorRank
INTEGER, DIMENSION(:), ALLOCATABLE :: collectors
INTEGER, DIMENSION(:), ALLOCATABLE :: collectorForRank

INTEGER :: worldGroup
INTEGER, DIMENSION(:), ALLOCATABLE :: sendRanks
INTEGER                            :: numSendRanks

LOGICAL :: roundRobin
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)
CALL FILL_DOMAIN_DECOMP(COMM, mpiErr, domainDecomp, domainIndex)

roundRobin = (rankOrder == 0)

numNodes = numRanks / ranksPerNode

IF (rank == 0) THEN
  nCols = nTrackCols
  ! For large jobs, 120 MB per collector seems about right
  ! And let's target 100 steps per dump
  ! (8 bytes per number)
  ! That's about 120e6 / (8*nCols) = 1e6 ptcl-steps per collector
  ! :KLUGE: this is a hasty algorithm; I'm sure it could be improved.
  tgtMemPerCollector = 120e6 ! Bytes
  tgtStepsPerDump = 100
  tgtPtclsPerCollector = FLOOR(tgtMemPerCollector / (nCols*8) / tgtStepsPerDump)
  numCollectors = MAX(MIN(2, numRanks), 1 + (NSAMPLE-1) / tgtPtclsPerCollector)
  numCollectors = MIN(numCollectors, numNodes)
  stepsPerTrackedPtclDump = CEILING(tgtMemPerCollector / &
    (8.*nCols*NSAMPLE/numCollectors))
  stepsPerTrackedPtclDump = MIN(200, stepsPerTrackedPtclDump)
  PRINT *, "Particle tracking:"
  PRINT *, "  #collectors (hence number of files) =", numCollectors
  PRINT *, "  approx. num. particles per collector =", tgtPtclsPerCollector
  PRINT *, "  timesteps per dump of tracked particle data =", stepsPerTrackedPtclDump
ENDIF
CALL MPI_BCAST(stepsPerTrackedPtclDump, 1, MPI_INTEGER, 0, COMM, mpiErr)

! initialize variables
! following should work equally well at t=0, or later after a restore.
lastTrackedStep = -1
trackedStepsSinceDump = 0
firstTrackedStepSinceDump = -1

! Tell everyone how many collectors there will be
CALL MPI_BCAST(numCollectors, 1, MPI_INTEGER, 0, COMM, mpiErr)

ALLOCATE(collectors(numCollectors))
ALLOCATE(collectorForRank(numRanks))
  
! :KLUGE: do a better job of choosing collectors
! but never more than one collector per node
IF (rank == 0) THEN
  ! choose collectors
  collectorForRank = 0
  IF (roundRobin) THEN
    DO i = 1, numCollectors
      collectors(i) = ((i-1)*numNodes)/numCollectors
    ENDDO
  ELSE
    DO i = 1, numCollectors
      collectors(i) = (((i-1)*numNodes)/numCollectors) * ranksPerNode
    ENDDO
  ENDIF
  PRINT *, "  collector ranks:", collectors
  
ENDIF

! Broadcast the collectors to everyone
CALL MPI_BCAST(collectors, numCollectors, MPI_INTEGER, 0, COMM, mpiErr)

! collectors(i) must be increasing with i
! and never more than one collector per node

! Figure out the collector for each rank
IF (rank == 0) THEN
  ! choose ranks to each collector
  ! Always choose collector on same node if there is a collector on that node
  IF (roundRobin) THEN
    DO i = 1, numCollectors
      jMax = numNodes - 1
      IF (i < numCollectors) jMax = collectors(i+1) - 1
      DO j = collectors(i)+1, jMax+1
        collectorForRank(j:numRanks:numNodes) = collectors(i)
      ENDDO
    ENDDO
  ELSE
    DO i = 1, numCollectors
      jMax = numRanks - 1
      IF (i < numCollectors) jMax = collectors(i+1) - 1
        collectorForRank(collectors(i)+1:jMax+1) = collectors(i)
    ENDDO
  ENDIF
ENDIF
! Send out the collectorsForRank
CALL MPI_BCAST(collectorForRank, numRanks, MPI_INTEGER, 0, COMM, mpiErr)

myCollector = collectorForRank(rank + 1)

IF (writePerRankFiles) THEN
  OPEN(9,FILE=trim(perRankFile), POSITION='APPEND')
  WRITE(9,*) "Collector for this rank is rank ", myCollector
  CLOSE(9)
ENDIF

! Create groups to which this rank might belong:
!  (1) The group of all ranks that send to myCollector.
!  (2) (if rank == myCollector) The group of all collectors

! Allocate (possibly) much more than we need.
ALLOCATE(sendRanks(numRanks))

! collector is the first rank among sendRanks
sendRanks(1) = myCollector
numSendRanks = 1 ! start with 1 for the collector
DO i = 1, numRanks
  IF (myCollector == collectorForRank(i)) THEN
    IF (myCollector /= i-1) THEN
      numSendRanks = numSendRanks + 1
      sendRanks(numSendRanks) = i-1
    ELSE
      ! don't count collector twice (we started numSendRanks at 1)
    ENDIF
  ENDIF
ENDDO

iAmCollector = (rank == myCollector)

IF (iAmCollector) THEN
  IF (writePerRankFiles) THEN
    OPEN(9,FILE=trim(perRankFile), POSITION='APPEND')
    WRITE(9,*) "  Collector", rank, " collects from subgroup of ranks ", &
      sendRanks(1:numSendRanks)
    CLOSE(9)
  ENDIF
ENDIF

! Create groups and associated communicators

CALL MPI_COMM_GROUP(COMM, worldGroup, mpiErr)
CALL MPI_GROUP_INCL(worldGroup, numSendRanks, sendRanks(1:numSendRanks), &
  sendGroup, mpiErr)
CALL MPI_COMM_CREATE(COMM, sendGroup, sendComm, mpiErr)

! Even if this rank is not a collector (not in the collectGroup),
! it needs to participate in the calls to set up the communicators.
CALL MPI_GROUP_INCL(worldGroup, numCollectors, collectors, collectGroup,&
  mpiErr)
CALL MPI_COMM_CREATE(COMM, collectGroup, collectorComm, mpiErr)

! Figure out how many particles to collect (finally) per step
numCollectedPtclsPerStep = 0
IF (iAmCollector) THEN
  CALL MPI_COMM_RANK(collectorComm, collectorRank, mpiErr)
  numCollectedPtclsPerStep = NSAMPLE / numCollectors
  IF (collectorRank < MOD(NSAMPLE, numCollectors)) THEN
    numCollectedPtclsPerStep = numCollectedPtclsPerStep + 1
  ENDIF
ENDIF

IF (writePerRankFiles) THEN
  OPEN(9,FILE=trim(perRankFile), POSITION='APPEND')
  WRITE(9,*) "  Collector", rank, " will collect ", numCollectedPtclsPerStep, &
    " particles (per step)"
  CLOSE(9)
ENDIF

! Determine sequence number
IF (RESTORE) THEN
  ! Make sure following filename (base) agrees with that passed to
  ! DUMP_TRACKED_PARTICLES_SER
  CALL GET_TRACK_SEQ_NUM("data/orbits/trackedElectrons", trackSeqNum)
ELSE
  trackSeqNum = 0
ENDIF

! Initialize data structures

eTrackPtcls%active = 1
eTrackPtcls%numFilled = 0
ALLOCATE(eTrackPtcls%ra1(nTrackCols,1))
pTrackPtcls%active = 1
pTrackPtcls%numFilled = 0
ALLOCATE(pTrackPtcls%ra1(nTrackCols,1))

DEALLOCATE(sendRanks)

DEALLOCATE(collectors)
DEALLOCATE(collectorForRank)

END SUBROUTINE INIT_TRACKER_COMM

!***********************************************************************
! SUBROUTINE ASSIGN_TRACK_TAGS
! Choose and set the tags of tracked particles.
! Also initializes data structures for tracking.
! 
! INPUT
! - spec: the particles to be tracked
! - tags: the tags for all available particles
! - COMM: mpi comm
! OUTPUT
! - tags: are modified (only) for particles that will be tracked
!***********************************************************************

SUBROUTINE ASSIGN_TRACK_TAGS(spec, tags, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
CHARACTER(len=*)                       :: spec
INTEGER*8, DIMENSION(:), INTENT(INOUT) :: tags
INTEGER                                :: COMM

! LOCAL

! The total number of simulated particles can exceed a 4-byte integer,
!   but the number of tracked particles shouldn't, nor should the number
!   of particles per rank.  After all, each particle takes up more than
!   80 B, so 2e9 particles (max signed 4-byte int) would take 160 GB RAM.
! Similarly, 2e9 tracked particles would take 160 GB of disk space per
!  step, which is simply too much.
! If someone from the future reads this with contempt and frustration
!   on a laptop computer that has TBs of memory, I apologize and 
!   congratulate the future on its progress, but express dismay that you
!   haven't just made all INTEGERs 8 bytes to avoid this problem.
INTEGER :: rank, numRanks, mpiErr
INTEGER :: nTracked
INTEGER :: tStartIndex ! 1-indexed
INTEGER, DIMENSION(:), ALLOCATABLE :: selected, tStartIndices, nPerRank
INTEGER :: i
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)

! number of tracked particles on this rank (have to add one for some ranks)
nTracked = NSAMPLE / numRanks
IF (nTracked < 1) THEN
  ! randomly choose ranks with 1 tracked particle on each chosen rank
  IF (rank == 0) THEN
    ALLOCATE(selected(NSAMPLE))
    CALL RANDOM_SELECT(numRanks, selected)
    !print *, NSAMPLE, " ranks selected for 1 particle each:", selected
    ALLOCATE(nPerRank(numRanks))
    nPerRank = 0
    DO i = 1, NSAMPLE
      ! store tStartIndex in nPerRank
      nPerRank(selected(i)) = i
    ENDDO
    !print *, "particles per rank:", nPerRank
  ENDIF
  nTracked = 0
  CALL MPI_SCATTER(nPerRank, 1, MPI_INTEGER, tStartIndex, 1, MPI_INTEGER,&
    0, COMM, mpiErr)
  IF (tStartIndex > 0) nTracked = 1
  IF (rank == 0) THEN
    DEALLOCATE(selected, nPerRank)
  ENDIF
ELSE
  ! Each rank has either nTracked or nTracked+1 particles
  !   randomly selected within each rank
  ! (tStartIndex is 1-indexed, rank is 0-indexed)
  tStartIndex = 1 + rank * nTracked + MIN(MOD(NSAMPLE, numRanks), rank)
  IF (rank < MOD(NSAMPLE, numRanks)) THEN
    nTracked = nTracked + 1
  ENDIF
ENDIF

ALLOCATE(selected(nTracked))
! select particles...
CALL RANDOM_SELECT(SIZE(tags), selected)

IF (writePerRankFiles) THEN
  OPEN(9,FILE=trim(perRankFile), POSITION='APPEND')
  WRITE(9,*) "Tracking ", nTracked, " out of ", SIZE(tags), " local particles:", selected
  CLOSE(9)
ENDIF

! ...and mark them as tracked, by giving them negative tags
!  but also give the tag the right value so we know where to store
!  that particle's tracked data.
DO i = 1, nTracked
  tags(selected(i)) = -tStartIndex + 1 - i
ENDDO

DEALLOCATE(selected)

! [ep]TrackPtcls%ra[12] are deallocated in END_TRACKER

END SUBROUTINE ASSIGN_TRACK_TAGS

!***********************************************************************
! SUBROUTINE END_TRACKER
! Deallocate groups/communicators
!***********************************************************************
SUBROUTINE END_TRACKER()

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER :: mpiErr

CALL MPI_COMM_FREE(sendComm, mpiErr)
CALL MPI_GROUP_FREE(sendGroup, mpiErr)

IF (iAmCollector) THEN
  CALL MPI_COMM_FREE(collectorComm, mpiErr)
  CALL MPI_GROUP_FREE(collectGroup, mpiErr)
ENDIF

IF (eTrackPtcls%active == 1) DEALLOCATE(eTrackPtcls%ra1)
IF (eTrackPtcls%active == 2) DEALLOCATE(eTrackPtcls%ra2)
eTrackPtcls%numFilled = 0
eTrackPtcls%active = -1
IF (pTrackPtcls%active == 1) DEALLOCATE(pTrackPtcls%ra1)
IF (pTrackPtcls%active == 2) DEALLOCATE(pTrackPtcls%ra2)
pTrackPtcls%numFilled = 0
pTrackPtcls%active = -1

END SUBROUTINE END_TRACKER

!***********************************************************************
! SUBROUTINE COPY_TABLE_TO_LARGER
!   Allocate array ra2 with more space, and copy ra1 to it;
!   then deallocate ra1.
!  These are tables (2D arrays); minSize is the minimum size for ra2
!   (in dimension 2).  ra2 will be larger in dimension 2 only.
!***********************************************************************

SUBROUTINE COPY_TABLE_TO_LARGER(ra1, ra2, minSizeOverride)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: ra1, ra2
INTEGER, OPTIONAL :: minSizeOverride

INTEGER :: minSize = 1
INTEGER :: i, nCols, oldSize, newSize
!***********************************************************************

IF (PRESENT(minSizeOverride)) minSize = minSizeOverride

nCols = SIZE(ra1, DIM=1)
oldSize = SIZE(ra1, DIM=2)
newSize = MAX(minSize, CEILING(oldSize * raGrowthFactor) + 1)
ALLOCATE(ra2(nCols, newSize))
ra2(1:nCols, 1:oldSize) = ra1
DEALLOCATE(ra1)

END SUBROUTINE COPY_TABLE_TO_LARGER

!***********************************************************************
! Subroutine TRACKER
! This subroutine writes a sample of particle orbits to disk
!
! INPUT:
! 
! - xp,yp,zp,uxp,uyp,uzp,wt: position, velocity and weight of all the particles
! - Ell,Bpp: parallel electric field and perpendicular magnetic field at (xp,yp)
! - tag: tag of the all the particles
! - NPP: Number of particles per process
! - sym: type of particles
! - step: the current time step
!
!***********************************************************************

SUBROUTINE TRACKER(pcl,pcl_data,tag,NPP,spec, step)

IMPLICIT NONE

! GRW: NPP, ip should be INTEGER*4
! GRW: I can't remember why I wrote the above comment.
INTEGER*8                                       :: NPP,ip
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)          :: pcl
DOUBLE PRECISION, DIMENSION(1:4,1:NPP)          :: pcl_data
INTEGER*8, DIMENSION(1:NPP)                     :: tag
CHARACTER(LEN=10)                               :: spec
INTEGER                                         :: step

!***********************************************************************

IF (lastTrackedStep < step) THEN
  IF (firstTrackedStepSinceDump < 0) firstTrackedStepSinceDump = step
  trackedStepsSinceDump = trackedStepsSinceDump + 1
  lastTrackedStep = step
ENDIF

! This step left-justify the string
!spec=adjustl(spec)

IF (NPP > 0) THEN
  IF (trim(adjustl(spec)) == "electrons") THEN
    DO ip=1,NPP
      IF (tag(ip) < 0) THEN
        CALL ADD_TRACKED_PARTICLE(pcl(:,ip), pcl_data(:,ip), tag(ip), step, eTrackPtcls)
      ENDIF
    ENDDO
  ELSEIF (trim(adjustl(spec)) == "ions") THEN
    DO ip=1,NPP
      IF (tag(ip) < 0) THEN
        CALL ADD_TRACKED_PARTICLE(pcl(:,ip), pcl_data(:,ip), tag(ip), step, pTrackPtcls)
      ENDIF
    ENDDO
  ELSE
    print *, "Error: spec is neither 'electrons' nor 'ions':", spec, ":"
  ENDIF
ENDIF

END SUBROUTINE TRACKER

!***********************************************************************
! Subroutine ADD_TRACKED_PARTICLE
! Add particle tracking data to the local table for a particle
!
! INPUT
! - ptclData: the array of x,y,z, ux, uy, uz, weight for the particle
! - xtraData: other data stored for the tracked particle
! - tag: the tag for the tracked particle
! - timestep: the current timestep (needed to sort particles after
!     collecting them on different ranks, before saving to disk)
! - gra: either eTrackPtcls or pTrackPtcls
!***********************************************************************

SUBROUTINE ADD_TRACKED_PARTICLE(ptclData, xtraData, tag, timestep, gra)

IMPLICIT NONE

DOUBLE PRECISION                  :: ptclData(7), xtraData(4)
INTEGER*8                         :: tag
INTEGER                           :: timestep
TYPE(GrowingTable), INTENT(INOUT) :: gra

INTEGER                        :: raSize, nf
!***********************************************************************

! grow array if necessary
IF (gra%active == 1) THEN
  raSize = SIZE(gra%ra1, DIM=2)
  !print *, "ra1: shape = ", SHAPE(gra%ra1), ", numFilled = ", gra%numFilled
  IF (raSize <= gra%numFilled) THEN
    CALL COPY_TABLE_TO_LARGER(gra%ra1, gra%ra2)
    !print *, "Enlarged array: shape(gra-ra2) = ", SHAPE(gra%ra2)
    gra%active = 2
  ENDIF
ELSEIF (gra%active == 2) THEN
  raSize = SIZE(gra%ra2, DIM=2)
  !print *, "ra2: shape = ", SHAPE(gra%ra2), ", numFilled = ", gra%numFilled, ', size = ', raSize, (raSize <= gra%numFilled) 
  IF (raSize <= gra%numFilled) THEN
    CALL COPY_TABLE_TO_LARGER(gra%ra2, gra%ra1)
    !print *, "Enlarged array: shape(gra-ra1) = ", SHAPE(gra%ra1)
    gra%active = 1
  ENDIF
ELSE
  gra%active = 1
  ! We don't have the information we need here to make a good
  ! first allocation.  That information is present at the end of
  ! ASSIGN_TRACK_TAGS.  Should we worry about that?  The cost of 
  ! increasing to the size of the first allocation is probably small.
  ALLOCATE(gra%ra1(nTrackCols, 10))
ENDIF

gra%numFilled = gra%numFilled + 1
nf = gra%numFilled
IF (gra%active == 1) THEN
  gra%ra1(1:7,nf) = ptclData
  gra%ra1(8:11,nf) = xtraData
  gra%ra1(12, nf) = tag
  gra%ra1(13, nf) = timestep
ELSE
  gra%ra2(1:7,nf) = ptclData
  gra%ra2(8:11,nf) = xtraData
  gra%ra2(12, nf) = tag
  gra%ra2(13, nf) = timestep
ENDIF

END SUBROUTINE ADD_TRACKED_PARTICLE

!***********************************************************************
! Subroutine ARRANGE_TABLE
!  Sort a table of entries in place according to a given arrangement,
! so that
! table[after](:,rowMap(i)) = table[before](:,i)
!
! INPUTS
! - table:
! - rowMap: row i [table(:,i)] should be moved 
!           to row rowMap(i) [table(:,rowMap(i))]
!   rowMap will be altered by this subroutine in an undefined way
!
!***********************************************************************

SUBROUTINE ARRANGE_TABLE(table, rowMap)

IMPLICIT NONE

! INPUT
DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: table
INTEGER, DIMENSION(:), INTENT(INOUT)            :: rowMap

! LOCAL
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tmpRow
INTEGER :: nRows, nCols, indx, newIndx, tmpIndx
!***********************************************************************

! Sorting is easy when you know where things go

nRows = SIZE(table, DIM=2)
IF (nRows > 1) THEN
  nCols = SIZE(table, DIM=1)
  ALLOCATE(tmpRow(nCols))

  ! N.B. This swapping algorithm approaches 3 row-copies per row of table,
  !   in the limit of a large number of rows.
  ! There is a slightly more complicated method approaching 2 row-copies per row.
  ! It can be done with (approx) 1 row-copy per row, but at the additional
  ! cost of finding (and storing) the inverse of rowMap...or is there a better way?

  indx = 1
  DO
    newIndx = rowMap(indx)
    !print *, "indx: ", indx, "->", newIndx
    IF (newIndx == indx) THEN
      indx = indx + 1
      IF (indx > nRows) EXIT
    ELSE
  ! swap rows indx and newIndx
      tmpRow = table(:,newIndx)
      table(:,newIndx) = table(:,indx)
      table(:,indx) = tmpRow
      rowMap(indx) = rowMap(newIndx)
      rowMap(newIndx) = newIndx
    ENDIF
    
  ENDDO
ENDIF

END SUBROUTINE ARRANGE_TABLE

!***********************************************************************
! Subroutine SORT_AMONG_RANKS
! This function, called by all ranks in COMM, takes entries in a table
! (an M x N array with N entries, each with M elements) and sends each
! entry to the appropriate processor.
!
! INPUT:
! 
! - tableIn: the local table -- this will get re-ordered
! - targetRankFn: a function that takes one entry (e.g., table(:,i))
!     and returns the rank (in COMM) that should receive that entry;
!     this function must be local (no mpi commands).
!   N.B. this should return a 0-indexed rank
! - COMM: mpi comm
!
! OUTPUT: 
! - tableOut: the table comprising all entries from all tables 
!     (on all ranks) that should be sent to this rank; its size
!     must be known in advance (that could be changed, however...)
!
!***********************************************************************

SUBROUTINE SORT_AMONG_RANKS(tableIn, tableOut, targetRankFn, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: tableIn
DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT)   :: tableOut
INTEGER                                       :: COMM
INTERFACE
  INTEGER FUNCTION targetRankFn(tableEntry)
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tableEntry
  END FUNCTION
END INTERFACE

! LOCAL
INTEGER :: rank, numRanks, mpiErr, globalRank
! entriesToRank(r+1) = #entries in tableIn going to rank r
INTEGER, DIMENSION(:), ALLOCATABLE :: entriesToRank
! tgtRanks(i) = the rank (+1, i.e., 1-indexed) 
!   to which tableIn(:,i) should be sent
! (and later, a map from current row to ordered row)
INTEGER, DIMENSION(:), ALLOCATABLE :: tgtRanks
! a temporary array to figure out which table row goes where
INTEGER, DIMENSION(:), ALLOCATABLE :: numFilledForRank

! tableIn, reordered so that entries going to the same rank are adjacent,
!   in order of the rank to which they'll be sent
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ordTableIn
! ordTableIn(:,offsetForRank(r+1)+1) is the first entry in ordTableIn
!   going to to rank r
INTEGER, DIMENSION(:), ALLOCATABLE :: offsetsForRank
! the number of entries to receive from each rank
INTEGER, DIMENSION(:), ALLOCATABLE :: entriesFromRank
! tableOut(:,recvOffsets(r+1)+1) will be the first entry in tableOut
!   from rank r
INTEGER, DIMENSION(:), ALLOCATABLE :: recvOffsets

INTEGER :: nCols, nLocalEntries, recvTotal
INTEGER :: i, rp1
!********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)

nCols = SIZE(tableIn, DIM=1)
nLocalEntries = SIZE(tableIn, DIM = 2)

ALLOCATE(entriesToRank(numRanks), entriesFromRank(numRanks), &
  offsetsForRank(numRanks), recvOffsets(numRanks), &
  numFilledForRank(numRanks), &
  tgtRanks(nLocalEntries))

entriesToRank = 0.
! Figure out which entry goes to which rank
DO i = 1, nLocalEntries
  tgtRanks(i) = 1 + targetRankFn(tableIn(:,i))
  entriesToRank(tgtRanks(i)) = entriesToRank(tgtRanks(i)) + 1
ENDDO
! and get the offsets describing where data for each rank should start
CALL CUMSUM(entriesToRank, offsetsForRank, i)

! Communicate above info.
CALL MPI_ALLTOALL(entriesToRank, 1, MPI_INTEGER, &
  entriesFromRank, 1, MPI_INTEGER, COMM, mpiErr)

CALL CUMSUM(entriesFromRank, recvOffsets, recvTotal)

! Reorder table
! Change tgtRanks to the map from current row to future row
numFilledForRank = 0
DO i = 1, nLocalEntries
  rp1 = tgtRanks(i)
  ! recvOffsets values are 0-indexed
  tgtRanks(i) = 1 + offsetsForRank(rp1) + numFilledForRank(rp1)
  numFilledForRank(rp1) = numFilledForRank(rp1) + 1
ENDDO
DEALLOCATE(numFilledForRank)
CALL ARRANGE_TABLE(tableIn, tgtRanks)
DEALLOCATE(tgtRanks)

! We've been keeping track of data in terms of rows; 
! convert to doubles by multiplying by nCols.
entriesToRank = entriesToRank * nCols
offsetsForRank = offsetsForRank * nCols
entriesFromRank = entriesFromRank * nCols
recvOffsets = recvOffsets * nCols

! Check for errors; if possible never let MPI be the first to see an error.
IF (SIZE(tableOut, DIM=1) /= nCols &
    .OR. SIZE(tableOut, DIM=2) /= recvTotal) THEN
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, globalRank, mpiErr)
  print *, "Error in MOD_TRACK::SORT_AMONG_RANKS:", &
    "On global rank ", globalRank, ", local (collector) rank ", rank, &
    "tableOut has shape", SHAPE(tableOut), &
    "but should have shape (", nCols, recvTotal, ")"
  CALL MPI_ABORT(MPI_COMM_WORLD, 343, mpiErr)
  stop 343
ENDIF

! Finally, exchange data with other ranks
CALL MPI_ALLTOALLV(tableIn, entriesToRank, offsetsForRank, &
  MPI_DOUBLE_PRECISION, tableOut, entriesFromRank, recvOffsets, &
  MPI_DOUBLE_PRECISION, COMM, mpiErr)

DEALLOCATE(entriesToRank, entriesFromRank)
DEALLOCATE(offsetsForRank, recvOffsets)

END SUBROUTINE SORT_AMONG_RANKS

!***********************************************************************
! SUBROUTINE GATHER_TABLE_ON_RANK0
! Gather tables---arrays of shape (ncol,:)---from all ranks onto rank 0
!
! INPUT
! - localTable: the local table on this rank, of shape (ncol, :)
! - globalTable: on rank0, this will be allocated to the right size,
!     and filled with the values from all the local tables in COMM
!     (it doesn't get allocated or touched at all on ranks > 0)
!     The calling function must deallocate globalTable (if local rank 0)
! - COMM: mpi comm
!***********************************************************************

SUBROUTINE GATHER_TABLE_ON_RANK0(localTable, globalTable, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                  :: localTable
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)    :: globalTable
INTEGER                                                       :: COMM

! LOCAL
INTEGER, DIMENSION(:), ALLOCATABLE :: numFromEachRank, offsetForEachRank
INTEGER :: rank, numRanks, mpiErr, nCols, nLocalEntries, nGlobalEntries
INTEGER :: nGlobalRows, i
!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)

nCols = SIZE(localTable, DIM=1)
nLocalEntries = SIZE(localTable, DIM = 2)

IF (rank == 0) THEN
  ALLOCATE(numFromEachRank(numRanks), offsetForEachRank(numRanks))
ENDIF

CALL MPI_GATHER((/ nLocalEntries /), 1, MPI_INTEGER, numFromEachRank, 1, &
  MPI_INTEGER, 0, COMM, mpiErr)

IF (rank == 0) THEN
  CALL CUMSUM(numFromEachRank, offsetForEachRank, nGlobalRows)
  IF (ALLOCATED(globalTable)) DEALLOCATE(globalTable)
  ALLOCATE(globalTable(nCols, nGlobalRows))
ENDIF

! convert from rows to elements
IF (rank == 0) THEN
  numFromEachRank = numFromEachRank * nCols
  offsetForEachRank = offsetForEachRank * nCols
ENDIF

CALL MPI_GATHERV(localTable, SIZE(localTable), MPI_DOUBLE_PRECISION, globalTable, &
  numFromEachRank, offsetForEachRank, MPI_DOUBLE_PRECISION, 0, COMM, mpiErr)

IF (rank == 0) DEALLOCATE(numFromEachRank, offsetForEachRank)

END SUBROUTINE GATHER_TABLE_ON_RANK0

!***********************************************************************
! Subroutine SORT_AMONG_COLLECTORS
! This function, called by all ranks takes entries in a table
! (an M x N array with N entries, each with M elements) and sends each
! entry to the appropriate collector rank, in two steps.
! First, the subgroup (sendGroup) gathers all data on its collector;
! then all the collectors communicate via alltoallv.
!
! INPUT:
! 
! - tableIn: the local table -- this will get re-ordered
! - targetCollectorFn: a function that takes one entry (e.g., table(:,i))
!     and returns the rank (in COMM) that should receive that entry;
!     this function must be local (no mpi commands).
!   N.B. this should return a 0-indexed rank
! - COMM: global mpi comm
!
! OUTPUT: 
! - tableOut: the table comprising all entries from all tables 
!     (on all ranks) that should be sent to this rank; its size
!     must be known in advance (that could be changed, however...)
!     If this rank is not a collector, tableOut will be unchanged.
!
!***********************************************************************

SUBROUTINE SORT_AMONG_COLLECTORS(tableIn, tableOut, targetCollectorFn, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! INPUT
DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: tableIn
DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT)   :: tableOut
INTEGER                                       :: COMM
INTERFACE
  INTEGER FUNCTION targetCollectorFn(tableEntry)
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tableEntry
  END FUNCTION
END INTERFACE

! LOCAL
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gatheredTable

INTEGER :: rank, numRanks, mpiErr

INTEGER :: nCols, nLocalEntries, recvTotal
INTEGER :: i, rp1

!********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
!CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)

CALL GATHER_TABLE_ON_RANK0(tableIn, gatheredTable, sendComm)

nCols = SIZE(tableIn, DIM=1)

IF (iAmCollector .EQV. .TRUE.) THEN
  CALL SORT_AMONG_RANKS(gatheredTable, tableOut, targetCollectorFn, collectorComm)
  DEALLOCATE(gatheredTable)
ENDIF

END SUBROUTINE SORT_AMONG_COLLECTORS

!***********************************************************************
! SUBROUTINE WRITE_TRACK_DATA
! 
! INPUT
! - forceWrite: if .TRUE., all tracked particle data will be written;
!               if .FALSE., data will be written only if enough
!                  data has been collected to make it worthwhile
! - COMM: mpi comm
! 
!***********************************************************************

SUBROUTINE WRITE_TRACK_DATA(forceWrite, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

LOGICAL          :: forceWrite
INTEGER          :: COMM

! LOCAL
INTEGER                                       :: stepsToDump
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: finalTable
INTEGER, DIMENSION(:), ALLOCATABLE            :: sortMap, timeSteps
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: times
INTEGER                                       :: i, numPtclSteps, collectorRank
INTEGER                                       :: mpiErr
!***********************************************************************

IF ((forceWrite .OR. trackedStepsSinceDump > stepsPerTrackedPtclDump) &
  .AND. (trackedStepsSinceDump > 0)) THEN

  stepsToDump = trackedStepsSinceDump
  IF (.FALSE. .AND. iAmCollector) THEN
    print *, "Writing track data after last tracked step", lastTrackedStep
    print *, "Writing data for", stepsToDump, "steps"
    print *, "numCollectedPtclsPerStep =", numCollectedPtclsPerStep
    print *, 'seqNum = ', trackSeqNum
  ENDIF
  IF (iAmCollector) THEN
    numPtclSteps = stepsToDump*numCollectedPtclsPerStep
    ALLOCATE(finalTable(nTrackCols, numPtclSteps))
    ALLOCATE(sortMap(numPtclSteps))
    CALL MPI_COMM_RANK(collectorComm, collectorRank, mpiErr)
  ENDIF
! electrons
  IF (eTrackPtcls%active == 1) THEN
    CALL SORT_AMONG_COLLECTORS(eTrackPtcls%ra1(:,1:eTrackPtcls%numFilled), &
      finalTable, GET_COLLECTOR_RANK_FOR_PTCL ,COMM)
  ELSE
    CALL SORT_AMONG_COLLECTORS(eTrackPtcls%ra2(:,1:eTrackPtcls%numFilled), &
      finalTable, GET_COLLECTOR_RANK_FOR_PTCL ,COMM)
  ENDIF
  IF (iAmCollector) THEN
  ! sort particles before dumping
  !     n = 1 + (t - 1 - r_c) / numCollectors
  ! The time step is the last column, the tag the penultimate.
    sortMap = 1 + (-finalTable(nTrackCols-1,:) - 1)/numCollectors &
      + (finalTable(nTrackCols,:) - firstTrackedStepSinceDump)/FTRACK &
        * numCollectedPtclsPerStep
    CALL ARRANGE_TABLE(finalTable, sortMap)
    ALLOCATE(timeSteps(stepsToDump))
    ALLOCATE(times(stepsToDump))
    timeSteps = finalTable(nTrackCols,1::numCollectedPtclsPerStep)
    times = timeSteps * dt
    ! Make sure following filename (base) agrees with that passed to
    ! GET_TRACK_SEQ_NUM
    CALL DUMP_TRACKED_PARTICLES_SER("data/orbits/trackedElectrons", finalTable, &
      numCollectedPtclsPerStep, collectorRank, trackSeqNum, timeSteps, times)
  ENDIF
  eTrackPtcls%numFilled = 0

! ions
  IF (pTrackPtcls%active == 1) THEN
    CALL SORT_AMONG_COLLECTORS(pTrackPtcls%ra1(:,1:pTrackPtcls%numFilled), &
      finalTable, GET_COLLECTOR_RANK_FOR_PTCL ,COMM)
  ELSE
    CALL SORT_AMONG_COLLECTORS(pTrackPtcls%ra2(:,1:pTrackPtcls%numFilled), &
      finalTable, GET_COLLECTOR_RANK_FOR_PTCL ,COMM)
  ENDIF
  IF (iAmCollector) THEN
  ! sort particles before dumping
  !     n = 1 + (t - 1 - r_c) / numCollectors
    sortMap = 1 + (-finalTable(nTrackCols-1,:) - 1)/numCollectors &
      + (finalTable(nTrackCols,:) - firstTrackedStepSinceDump)/FTRACK &
        * numCollectedPtclsPerStep
    CALL ARRANGE_TABLE(finalTable, sortMap)
    CALL DUMP_TRACKED_PARTICLES_SER("data/orbits/trackedIons", finalTable, &
      numCollectedPtclsPerStep, collectorRank, trackSeqNum, timeSteps, times)

    DEALLOCATE(finalTable, sortMap, timeSteps, times)
  ENDIF
  pTrackPtcls%numFilled = 0

  firstTrackedStepSinceDump = -1
  trackedStepsSinceDump = 0

ENDIF

END SUBROUTINE WRITE_TRACK_DATA

!***********************************************************************
! Subroutine INC_TRACK_SEQ_NUM
!
! This should be called after each checkpoint.
! 
!***********************************************************************
SUBROUTINE INC_TRACK_SEQ_NUM()

IMPLICIT NONE 

trackSeqNum = trackSeqNum + 1

END SUBROUTINE INC_TRACK_SEQ_NUM

!***********************************************************************
! Subroutine TEST_TRACKER
! Purely for test purposes, for testing this module, and also testing
! the dumping of tracked particles.
! 
!***********************************************************************

SUBROUTINE TEST_TRACKER(COMM, rankOrder, ranksPerNode)

IMPLICIT NONE 

INCLUDE 'mpif.h'

INTEGER :: COMM, rankOrder, ranksPerNode

INTEGER :: it, nSteps, n, nPrevEntries
INTEGER*8 :: nEntries
INTEGER :: rank, numRanks, mpiErr
DOUBLE PRECISION :: t
DOUBLE PRECISION, ALLOCATABLE :: ptclData(:,:), xtraData(:,:)
INTEGER*8, ALLOCATABLE :: tag(:)
DOUBLE PRECISION :: t1, t2
INTEGER, DIMENSION(:), ALLOCATABLE :: selection
!***********************************************************************

t1 = MPI_WTIME()
nSteps = NT

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)

PRINT *, "Rank", rank, " of ", numRanks, " starting"
IF (rank == 0) THEN
  PRINT *, "Starting test of tracker, running for ", nSteps, " steps"
  ! test random selection
  !ALLOCATE(selection(700))
  !RANDOM_SELECT(30000, selection)
  !PRINT *, "Random selection of 700 out of 30000"
  !PRINT *, selection
  !DEALLOCATE(selection)
ENDIF

CALL INIT_TRACKER_COMM(COMM, rankOrder, ranksPerNode)
PRINT *, "Rank", rank, " initialized"

nEntries = NSAMPLE/numRanks
nPrevEntries = nEntries * rank
IF (rank < MOD(NSAMPLE, numRanks)) nEntries = nEntries + 1
nPrevEntries = nPrevEntries + MIN(rank, MOD(NSAMPLE, numRanks))

ALLOCATE(ptclData(7,nEntries), xtraData(4,nEntries))
ALLOCATE(tag(nEntries))

PRINT *, "Rank", rank, " entering loop with nEntries =", nEntries
PRINT *, "Rank", rank, " n = ", nPrevEntries+1, nPrevEntries+nEntries
DO it = 1, nSteps
  t = it * 0.1
  IF (rank == 0) print *, "Step", it
  DO n = nPrevEntries+1, nPrevEntries + nEntries
    tag(n-nPrevEntries) = -n
    ptclData(:,n-nPrevEntries) = (/ 1.*n, 0.1*n, 0.01*n, 0.001*n, 1e-4*n, 1e-5*n, 1e-6*n /)
    xtraData(:,n-nPrevEntries) = (/ -1.*n, -2.*n, -4.*n, -8.*n /)
  ENDDO
  CALL TRACKER(ptclData, xtraData, tag, nEntries, " electrons", it)
  CALL TRACKER(ptclData, xtraData, tag, nEntries, "      ions", it)
  CALL MPI_BARRIER(COMM, mpiErr)
  CALL WRITE_TRACK_DATA(.FALSE., COMM)
ENDDO
CALL WRITE_TRACK_DATA(.TRUE., COMM)

CALL MPI_BARRIER(COMM, mpiErr)
IF (rank == 0) THEN
  t2 = MPI_WTIME()
  t2 = t2 - t1
  PRINT *, "Elapsed time: ", t2, " seconds"
  PRINT *, "Halting test of tracker and crudely halting program"

ENDIF

CALL END_TRACKER()
DEALLOCATE(ptclData, xtraData)
DEALLOCATE(tag)


CALL MPI_ABORT(COMM, 3418, mpiErr)
STOP 3418

END SUBROUTINE TEST_TRACKER

!***********************************************************************

END MODULE
