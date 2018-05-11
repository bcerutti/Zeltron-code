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

PROGRAM MAIN

USE MOD_INPUT
USE MOD_EXTCMD
USE MOD_MPI
USE MOD_LOG
USE MOD_IO
USE MOD_INITIAL
USE MOD_FIELDS
USE MOD_MOTION
USE MOD_RHOJ
USE MOD_TRACK
USE MOD_ANALYSIS

IMPLICIT NONE

INCLUDE 'mpif.h'

!***********************************************************************

! MPI variables
INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: stat
INTEGER                               :: TOPO_COMM,ierr,NPROC,id
INTEGER, DIMENSION(2)                 :: coords
INTEGER, DIMENSION(8)                 :: ngh,NESC
CHARACTER*(MPI_MAX_PROCESSOR_NAME)    :: hostname
INTEGER                               :: hostnameLen
INTEGER                               :: rankOrder
INTEGER                               :: ranksPerNode

! LOCAL VARIABLES (VALID IN THE DOMAIN OF EACH PROCESS)
!***********************************************************************
! Initial number of particles per species per process
INTEGER*8        :: NPP

! Number of particles per species per process
INTEGER*8        :: NPB,NPD,NEB,NED

! Spatial boundaries in the X-direction of each domain
DOUBLE PRECISION :: xminp,xmaxp

! Spatial boundaries in the Y-direction of each domain
DOUBLE PRECISION :: yminp,ymaxp

! Global nodal grid
DOUBLE PRECISION, DIMENSION(1:NX) :: xg
DOUBLE PRECISION, DIMENSION(1:NY) :: yg
! Nodal grid in each domain
DOUBLE PRECISION, DIMENSION(1:NXP) :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP) :: ygp

! Yee grid in each domain
DOUBLE PRECISION, DIMENSION(1:NXP) :: xyeep
DOUBLE PRECISION, DIMENSION(1:NYP) :: yyeep

! Magnetic and Electric fields components Yee lattice
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bx,By,Bz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Ex,Ey,Ez

! Magnetic and Electric fields components at nodes
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Bxg,Byg,Bzg
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Exg,Eyg,Ezg

! Total charge and current densities
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: Jx0,Jy0,Jz0
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: rho,Jx,Jy,Jz

! Charge and current densities DRIFTING ELECTRONS
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: rhoed,Jxed,Jyed,Jzed
! Charge and current densities DRIFTING IONS
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: rhopd,Jxpd,Jypd,Jzpd

! Charge and current densities BACKGROUND ELECTRONS
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: rhoeb,Jxeb,Jyeb,Jzeb
! Charge and current densities BACKGROUND IONS
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP) :: rhopb,Jxpb,Jypb,Jzpb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reduced nodal grid for fields;  fields are saved on the reduced grid to save space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION, ALLOCATABLE :: xfrg(:), yfrg(:)

! Local reduced nodal grid for fields
DOUBLE PRECISION, ALLOCATABLE :: xfrgp(:), yfrgp(:)

! Shape of the reduced grid
INTEGER, DIMENSION(NDIM) :: localFRedGridShape
INTEGER, DIMENSION(NDIM) :: globalFRedGridShape

! offset of the first element of the reduced grid on this domain
INTEGER, DIMENSION(NDIM) :: localFRedGridOffset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reduced nodal grid for densities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION, ALLOCATABLE :: xdrg(:), ydrg(:)

! Local reduced nodal grid for densities
DOUBLE PRECISION, ALLOCATABLE :: xdrgp(:), ydrgp(:)

! Shape of the reduced grid
INTEGER, DIMENSION(NDIM) :: localDRedGridShape
INTEGER, DIMENSION(NDIM) :: globalDRedGridShape

! offset of the first element of the reduced grid on this domain
INTEGER, DIMENSION(NDIM) :: localDRedGridOffset

!***********************************************************************
! Other plasma paramters
DOUBLE PRECISION  :: mi,nd0,n0,delta,de,B0,sigma,grad

! External photon energy density
DOUBLE PRECISION  :: Uph

! Synchrotron energy losses
DOUBLE PRECISION  :: Esyned,Esynpd,Esyneb,Esynpb,Esyn_electrons,Esyn_ions

! Inverse Compton energy losses
DOUBLE PRECISION  :: Eicsed,Eicspd,Eicseb,Eicspb,Eics_electrons,Eics_ions

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++ PARTICLE DISTRIBUTION GENERATOR ++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Grid for the pre-calculated cumulative distribution functions
INTEGER, PARAMETER :: ND=1000

! Particle momentum parallel to the drift velocity for electrons, ions
DOUBLE PRECISION, DIMENSION(1:ND) :: upe
DOUBLE PRECISION, DIMENSION(1:ND) :: upp

! Cumulative distribution function for up for electrons, ions
DOUBLE PRECISION, DIMENSION(1:ND) :: gFpe
DOUBLE PRECISION, DIMENSION(1:ND) :: gFpp

! Particle momentum ps=uperp/sqrt(1+up*up) for electrons, ions
DOUBLE PRECISION, DIMENSION(1:ND) :: pse
DOUBLE PRECISION, DIMENSION(1:ND) :: psp

! Cumulative distribution function for ps for electrons, ions
DOUBLE PRECISION, DIMENSION(1:ND,1:ND) :: gFse
DOUBLE PRECISION, DIMENSION(1:ND,1:ND) :: gFsp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++ per-rank files for debugging    ++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CHARACTER(len=10) :: rankStr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Computing times
DOUBLE PRECISION  :: t0,t1,t2,t3
! temporary timing, like t0-t3
DOUBLE PRECISION  :: tCheckpt, tStep
DOUBLE PRECISION  :: nextCheckptTime
! More timing
DOUBLE PRECISION  :: initTime, elapsedTime
DOUBLE PRECISION  :: dumpTime, checkptTime, avgStepTime, avgStepCoef
! Number of checkpoint saves and dumps
INTEGER :: numCheckptSaves, numDumps, doCheckpt, doDump, doQuit

! For reading externally-given commands
INTEGER :: extCmd, cmdStep, doStuffComm(3)

!***********************************************************************

! Symbol
CHARACTER(LEN=10) :: bg,drift,electrons,ions,syn,ics

! Loop indexes
INTEGER :: ix,iy,it,it0,ig

!***********************************************************************

! log for step timing
TYPE(DataLogObj) :: timeLog

!***********************************************************************

! some very temporary arrays
DOUBLE PRECISION, ALLOCATABLE :: tmpRa1(:), tmpRa2(:)

!***********************************************************************

! Ion mass
mi=mass_ratio*me

bg="bg"
drift="drift"
electrons="electrons"
ions="ions"
syn="syn"
ics="ics"


!=======================================================================
! COMMAND INITIALIZATION
!=======================================================================
  
extCmd = EXTCMD_DO_NOTHING
cmdStep = -1

!=======================================================================
! SANITY CHECKS
!=======================================================================

! In Jan 2015 FSAVE was changed from an int to a double...this stops
! the code from compiling if FSAVE is an int.
CALL ASSERT_DOUBLE(FSAVE) ! FSAVE should be DOUBLE PRECISION

IF (NSAMPLE > PPC*NCX*NCY) THEN
  PRINT*, 'Error: NSAMPLE is greater than the number of particles'
  STOP 98271
ENDIF


!=======================================================================
! INITIALIZE THE MPI ENVIRONMENT
!=======================================================================

CALL MPI_INIT(ierr)

! Test whether MPI is working properly or not
IF (ierr.NE.MPI_SUCCESS) THEN
   PRINT*,'Error starting MPI program. Terminating.'
   CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
END IF

initTime = MPI_WTIME()
t0 = initTime

! To obtain the ID number of each process
CALL MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
globalMpiRank = id

IF (id.EQ.0) THEN
PRINT*,'*******************************************************'
PRINT*,'         Welcome to the Zeltron code project!          '
PRINT*,'*******************************************************'
PRINT*,' Copyright (C) 2012-2015 Benoît Cerutti & Greg Werner  '
PRINT*,'                                                       '
PRINT*,' This program is free software: you can redistribute   '
PRINT*,' it and/or modify it under the terms of the GNU        '
PRINT*,' General Public License as published by the Free       '
PRINT*,' Software Foundation, either version 3 of the          '
PRINT*,' License, or (at your option) any later version.       '
PRINT*,'*******************************************************'
END IF

! To obtain the number of process
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ierr)
IF (NPROC /= NPX*NPY) THEN
  PRINT*, 'Error: the number of ranks does not equal NPX * NPY (in mod_input.f90)'
  PRINT*, '          #ranks = ', NPROC
  PRINT*, '             NPX = ', NPX
  PRINT*, '             NPY = ', NPY
  PRINT*, '       NPX * NPY = ', (NPX*NPY)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_ABORT(MPI_COMM_WORLD, 524, ierr)
  STOP 524
ENDIF

! To obtain the hostname
CALL MPI_GET_PROCESSOR_NAME(hostname, hostnameLen, ierr)
CALL FILL_RANK_ORDER(MPI_COMM_WORLD, rankOrder, ranksPerNode)
IF (id == 0) THEN
  PRINT*, 'Rank 0 is on ' // hostname(:hostnameLen)
  PRINT*, 'Running ', NPROC, ' ranks, with ', ranksPerNode, ' ranks per node'
  IF (rankOrder == 0) THEN
    PRINT*,   '  with rankOrder = ', rankOrder, ' (round robin)'
  ELSE
    PRINT*,   '  with rankOrder = ', rankOrder, ' (adjacent ranks on same node)'
  ENDIF
ENDIF

! Test: number of processes is a multiple of the number of grid cells
IF (id.EQ.0) THEN
  IF (MOD(NCX,NPX).NE.0.OR.MOD(NCY,NPY).NE.0) THEN
  PRINT*,'*******************************************************'
  PRINT*,'CAUTION: The number of cells along every direction' 
  PRINT*,'should be a multiple of the number of processes!'
  PRINT*,'*******************************************************'
  CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
  END IF
END IF

!=======================================================================
! INITIALIZE PER-RANK FILES
!=======================================================================
IF (writePerRankFiles) THEN
  write(rankStr, '(i10)') id
  rankStr = adjustl(rankStr)
  perRankFile = "rankLog_" // trim(rankStr) // ".txt"
  OPEN(9, FILE=trim(perRankFile))
  WRITE(9, *) "Rank " // trim(rankStr) // " is running on " // hostname(:hostnameLen)
  CLOSE(9)
ENDIF

!=======================================================================
! INITIALIZE THE VIRTUAL TOPOLOGY OF THE PROCESSES
!=======================================================================

CALL COM_TOPOLOGY(ngh,coords,NPROC,id,TOPO_COMM,ierr)

IF (writePerRankFiles) THEN
  OPEN(9, FILE=trim(perRankFile), POSITION='APPEND')
  WRITE(9, *) "Domain index:", coords
  CLOSE(9)
ENDIF

IF (id == 0) THEN
  t1 = MPI_WTIME()
  PRINT *, "MPI set-up time: ", (t1-t0), "seconds"
  t0 = t1
ENDIF

!=======================================================================
! Test tracker
! For debugging only -- the program will halt at the end of this call
!=======================================================================
!IF (TRACK_PARTICLES) THEN
!  CALL TEST_TRACKER(TOPO_COMM, rankOrder, ranksPerNode)
!ENDIF

!=======================================================================
! INITIALIZE DATA LOGS
!=======================================================================

CALL START_LOGS(id)
! save data every 100 steps or 30 seconds, whichever comes first
IF (id == 0) CALL START_LOG(timeLog, "./data/timestep.dat", 1, 100, 30.d0)

!=======================================================================
! SPATIAL BOUNDARIES FOR EACH SUB-DOMAIN
!=======================================================================

xminp=xmin+coords(1)*NCXP*dx
! Try to ensure that xmaxp = xminp of next domain to finite precision.
! Moreover, if domain is uppermost, make sure its xmaxp = xmax
xmaxp=xmin+(coords(1)+1)*NCXP*dx
IF (coords(1) == NPX-1) THEN
  xmaxp = xmax
ENDIF

yminp=ymin+coords(2)*NCYP*dy
ymaxp=ymin+(coords(2)+1)*NCYP*dy
IF (coords(2) == NPY-1) THEN
  ymaxp = ymax
ENDIF

!=======================================================================
! Set SPATIAL GRID (NODAL AND YEE LATTICE)
!=======================================================================

! Nodal lattice
DO ix=1,NX
xg(ix)=(ix-1)*1d0/((NX-1)*1d0)*(xmax-xmin)+xmin
ENDDO

DO iy=1,NY
yg(iy)=(iy-1)*1d0/((NY-1)*1d0)*(ymax-ymin)+ymin
ENDDO

! Local nodes
xgp = xg(coords(1)*NCXP+1 : (coords(1)+1)*NCXP+1)
ygp = yg(coords(2)*NCYP+1 : (coords(2)+1)*NCYP+1)

! Yee lattice
xyeep=xgp+dx/2.0
yyeep=ygp+dy/2.0

!=======================================================================
! Calculate Reduced Grid 
!=======================================================================

! For fields
! Add an upper grid point, which for periodic BCs will be identical
! to the lower-most point.
ALLOCATE(tmpRa1(NFIELDX+1), tmpRa2(NFIELDY+1))

CALL REDUCE_GRID(xg, yg, (/ xminp, yminp /), (/ xmaxp, ymaxp/), TOPO_COMM, &
  globalFRedGridShape, localFRedGridShape, localFRedGridOffset, &
  tmpRa1, tmpRa2)
  
ALLOCATE(xfrg(globalFRedGridShape(1)))
ALLOCATE(yfrg(globalFRedGridShape(2)))

xfrg = tmpRa1(:globalFRedGridShape(1))
yfrg = tmpRa2(:globalFRedGridShape(2))

ALLOCATE(xfrgp(localFRedGridShape(1)))
ALLOCATE(yfrgp(localFRedGridShape(2)))

xfrgp = xfrg(localFRedGridOffset(1):localFRedGridOffset(1)+localFRedGridShape(1)-1)
yfrgp = yfrg(localFRedGridOffset(2):localFRedGridOffset(2)+localFRedGridShape(2)-1)

DEALLOCATE(tmpRa1, tmpRa2)

! For densities
! Add an upper grid point, which for periodic BCs will be identical
! to the lower-most point.
ALLOCATE(tmpRa1(NDX+1), tmpRa2(NDY+1))
CALL REDUCE_GRID(xg, yg, (/ xminp, yminp /), (/ xmaxp, ymaxp/), TOPO_COMM, &
  globalDRedGridShape, localDRedGridShape, localDRedGridOffset, &
  tmpRa1, tmpRa2)
ALLOCATE(xdrg(globalDRedGridShape(1)))
ALLOCATE(ydrg(globalDRedGridShape(2)))
xdrg = tmpRa1(:globalDRedGridShape(1))
ydrg = tmpRa2(:globalDRedGridShape(2))
ALLOCATE(xdrgp(localDRedGridShape(1)))
ALLOCATE(ydrgp(localDRedGridShape(2)))
xdrgp = xdrg(localDRedGridOffset(1):localDRedGridOffset(1)+localDRedGridShape(1)-1)
ydrgp = ydrg(localDRedGridOffset(2):localDRedGridOffset(2)+localDRedGridShape(2)-1)
DEALLOCATE(tmpRa1, tmpRa2)

IF (writePerRankFiles) THEN
  OPEN(9, FILE=trim(perRankFile), POSITION='APPEND')
  WRITE(9, *) "SIZE(xg) = ", SIZE(xg)
  WRITE(9, *) "xg: ", xg
  WRITE(9, *) "SIZE(yg) = ", SIZE(yg)
  WRITE(9, *) "yg: ", yg
  WRITE(9, *) "xgp: ", xgp
  WRITE(9, *) "ygp: ", ygp
  WRITE(9, *) "Global reduced grid for fields:"
  WRITE(9, *) "xfrg: ", xfrg
  WRITE(9, *) "yfrg: ", yfrg
  WRITE(9, *) "Local reduced grid for fields:"
  WRITE(9, *) "localOffset ", localFRedGridOffset
  WRITE(9, *) "localShape ", localFRedGridShape
  WRITE(9, *) "globalShape ", globalFRedGridShape
  WRITE(9, *) "xfrgp: ", xfrgp
  WRITE(9, *) "yfrgp: ", yfrgp
  WRITE(9, *) "Global reduced grid for densities:"
  WRITE(9, *) "xdrg: ", xdrg
  WRITE(9, *) "ydrg: ", ydrg
  WRITE(9, *) "Local reduced grid for densities:"
  WRITE(9, *) "localOffset ", localDRedGridOffset
  WRITE(9, *) "localShape ", localDRedGridShape
  WRITE(9, *) "globalShape ", globalDRedGridShape
  WRITE(9, *) "xdrgp: ", xdrgp
  WRITE(9, *) "ydrgp: ", ydrgp
  CLOSE(9)
ENDIF


!=======================================================================
! Estimate timing: this is very rough, obviously
!   These values will be replaced with actual values later
!=======================================================================

avgStepTime = (4d0 * NP / NPROC) * 2e-6

! These estimates can be pretty lousy.
! It will be a good thing if the simulation dumps (thereby getting a 
! better estimate of dumpTime) well before its first checkpoint.
dumpTime = MIN(3600d0, 10 + 1d0*NCX*NCY/1d7 * 100)
checkptTime = MIN(FSAVE/2, MAX(3600d0, 10 + 1d0*NCX*NCY/1d7 * 50))

numDumps = 0 
numCheckptSaves = 0

!=======================================================================
! Initialization of the tracker
!=======================================================================
IF (TRACK_PARTICLES) THEN
  CALL INIT_TRACKER_COMM(TOPO_COMM, rankOrder, ranksPerNode)
ENDIF 

!=======================================================================
!  Write a file with pertinent simulation information
!=======================================================================

CALL DUMP_SIM_INFO(hostname, rankOrder, ranksPerNode, &
  (/ xminp, yminp /), (/ xmaxp, ymaxp /), TOPO_COMM)

IF (id == 0) THEN
  t1 = MPI_WTIME()
  PRINT *, "Grid set-up time: ", (t1-t0), "seconds"
  t0 = t1
ENDIF

!=======================================================================
! Test restore to see if it results in the same dump.
!   (Then halt simulation.)
!=======================================================================

IF (TEST_RESTORE) THEN
  IF (id.EQ.0) THEN
    t0=MPI_WTIME()
  END IF

  IF (id.EQ.0) THEN
  PRINT*,'*******************************************************'
  PRINT*,'1. BEGIN TEST OF RESTORATION '
  END IF

  it0=time_ref
  
  IF (id.EQ.0) PRINT*,'   Restoring fields'
  ! Restore the Yee fields
  CALL RESTORE_FIELDS(Ex,Ey,Ez,Bx,By,Bz,it0,coords,id,TOPO_COMM,ierr)


  IF (id.EQ.0) PRINT*,'   Restoring particles'
  ! Restore the particle distribution functions
  CALL RESTORE_PARTICLES(pcl_ed,pcl_data_ed,taged,NED,it0,electrons,drift,&
                         id,TOPO_COMM,ierr)
  IF (id.EQ.0) PRINT*,'   Rank 0 restored ', NED, ' drift electrons'
  CALL RESTORE_PARTICLES(pcl_eb,pcl_data_eb,tageb,NEB,it0,electrons,bg,&
                         id,TOPO_COMM,ierr)
  IF (id.EQ.0) PRINT*,'   Rank 0 restored ', NEB, ' bg electrons'
  CALL RESTORE_PARTICLES(pcl_pd,pcl_data_pd,tagpd,NPD,it0,ions,drift,&
                         id,TOPO_COMM,ierr)
  IF (id.EQ.0) PRINT*,'   Rank 0 restored ', NPD, ' drift ions'
  CALL RESTORE_PARTICLES(pcl_pb,pcl_data_pb,tagpb,NPB,it0,ions,bg,&
                         id,TOPO_COMM,ierr)
  IF (id.EQ.0) PRINT*,'   Rank 0 restored ', NPB, ' bg ions'

  ! Restore the nodal fields
  !CALL FIELDS_NODES(Bx,By,Bz,Ex,Ey,Ez,Bxg,Byg,Bzg,Exg,Eyg,Ezg,xgp,ygp,&
  !                  id,ngh,TOPO_COMM,ierr)

  IF (id.EQ.0) THEN
  PRINT*,'2. SESSION RESTORED SUCCESSFULLY!'
  END IF

  IF (id.EQ.0) THEN
    t3=MPI_WTIME()
    PRINT*,'Restoration time:',t3-t0,'seconds'
    t0=t3
    PRINT*,'*******************************************************'
    PRINT*,'3. NOW DUMPING SAME DATA AT NEXT TIMESTEP'
  END IF

  ! Advance time, and then dump, so we can compare files
  it = it0 + 1
  IF (id.EQ.0) THEN
    t0=MPI_WTIME()
  END IF
  
  ! Save the Yee fields to disk
  IF (id.EQ.0) PRINT*,'   Dumping Yee fields'
  CALL SAVE_FIELDS(Ex,Ey,Ez,Bx,By,Bz,it,coords,id,TOPO_COMM,ierr)
  
  ! Save the particle data to disk
  IF (id.EQ.0) PRINT*,'   Dumping bg electrons'
  CALL SAVE_PARTICLES(pcl_eb,pcl_data_eb,tageb,NEB,it,electrons,bg,&
                      id,TOPO_COMM,ierr)
  IF (id.EQ.0) PRINT*,'   Dumping drift electrons'
  CALL SAVE_PARTICLES(pcl_ed,pcl_data_ed,taged,NED,it,electrons,drift,&
                      id,TOPO_COMM,ierr)
  IF (id.EQ.0) PRINT*,'   Dumping bg ions'
  CALL SAVE_PARTICLES(pcl_pb,pcl_data_pb,tagpb,NPB,it,ions,bg,&
                      id,TOPO_COMM,ierr)
  IF (id.EQ.0) PRINT*,'   Dumping drift ions'
  CALL SAVE_PARTICLES(pcl_pd,pcl_data_pd,tagpd,NPD,it,ions,drift,&
                      id,TOPO_COMM,ierr)
  
  IF (id.EQ.0) THEN
    ! To measure time in seconds
    t3=MPI_WTIME()
    PRINT*,'4. COMPLETED DUMP SUCCESSFULLY!'
    PRINT*,'Checkpoint dump time:',t3-t0,'seconds'
    PRINT*,'*******************************************************'
    PRINT*,'Halting simulation after concluding restoration test.'
    PRINT*,'(Do not be alarmed at the MPI_ABORT error following this.)'
  END IF

  CALL MPI_BARRIER(TOPO_COMM, ierr)
  CALL MPI_ABORT(TOPO_COMM, 1047, ierr)
  STOP 1047
ENDIF


!=======================================================================
! Computation of the cumulative distribution functions
!=======================================================================

upe=0.0
upp=0.0
pse=0.0
psp=0.0

! Cumulative distribution function
gFpe=0.0
gFse=0.0
gFpp=0.0
gFsp=0.0

DO ig=1,ND
upe(ig)=(ig-1)*1d0/((ND-1)*1d0)*(udemax-udemin)+udemin
upp(ig)=(ig-1)*1d0/((ND-1)*1d0)*(udpmax-udemin)+udpmin
ENDDO

DO ig=1,ND
pse(ig)=(ig-1)*1d0/((ND-1)*1d0)*(udemax-0.0)+0.0
psp(ig)=(ig-1)*1d0/((ND-1)*1d0)*(udpmax-0.0)+0.0
ENDDO

! Numerical calculation of the cumulative distribution functions
CALL INIT_DRIFT_MAXWELLIAN(upe,gFpe,pse,gFse,thde,ND, TOPO_COMM)
IF (mass_ratio == 1 .AND. udemin == udpmin .AND. udemax == udpmax &
  .AND. thde == thdi) THEN
  gFpp = gFpe
  gFsp = gFse
ELSE
  CALL INIT_DRIFT_MAXWELLIAN(upp,gFpp,psp,gFsp,thdi,ND, TOPO_COMM)
ENDIF

IF (id == 0) THEN
  t1 = MPI_WTIME()
  PRINT *, "Rel. Maxwellian distribution set-up time: ", (t1-t0), "seconds"
  t0 = t1
ENDIF

!=======================================================================

IF (RESTORE.EQV..FALSE.) THEN

!=======================================================================
! Set initial FIELDS at t=0
!=======================================================================

CALL SET_FIELDS(nd0,delta,de,B0,grad,sigma,Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,&
                xgp,ygp,xyeep,yyeep)
IF (id == 0) THEN
  t1 = MPI_WTIME()
  PRINT *, "Field initialization time: ", (t1-t0), "seconds"
  t0 = t1
ENDIF

! External radiation field energy density
Uph=udens_ratio*(B0*B0/(8.0*pi))
                  
!=======================================================================
! DISPLAY INPUT PARAMETERS
!=======================================================================

IF (id.EQ.0) THEN
PRINT*,'*******************************************************'
PRINT*,'                    INPUT PARAMETERS                   '
PRINT*,'*******************************************************'
PRINT*,'Number of process          =',NPROC
PRINT*,'# cells x-direction        =',NCX
PRINT*,'# cells y-direction        =',NCY
PRINT*,'# Particles per cells      =',4*PPC
PRINT*,'Total # of Particles       =',4*NP
PRINT*,'# of tracked Particles     =',NSAMPLE
PRINT*,'rhoc [cm]                  =',rhoc
PRINT*,'rhoc/dx                    =',rhoc/((xmax-xmin)/(NX-1))
PRINT*,'rhoc/dy                    =',rhoc/((ymax-ymin)/(NY-1))
PRINT*,'# of time steps            =',NT
PRINT*,'Dump frequency             =',FDUMP
PRINT*,'Checkpoint frequency       =',FSAVE
PRINT*,'Drifting temp. e- [mc2]    =',thde
PRINT*,'Drifting temp. ion [mc2]   =',thdi
PRINT*,'Background temp. e- [mc2]  =',thbe
PRINT*,'Background temp. ion [mc2] =',thbi
PRINT*,'Density ratio              =',density_ratio
PRINT*,'Mass ions/electrons        =',mass_ratio
PRINT*,'delta/rhoc                 =',delta/rhoc
PRINT*,'de/rhoc                    =',de/rhoc
PRINT*,'B0 [Gauss]                 =',B0
PRINT*,'1/omegac [s]               =',rhoc/c
PRINT*,'dt*omegac                  =',dt/(rhoc/c)
PRINT*,'gamma_rad                  =',grad
PRINT*,'udens_ratio                =',udens_ratio
PRINT*,'gamma_drift                =',1.0/sqrt(1.0-betad*betad)
PRINT*,'Drift. dens. [cm-3]        =',nd0
PRINT*,'sigma                      =',sigma
PRINT*,'Amplitude Perturbation     =',perturb_amp
PRINT*,'Guide field strength in B0 =',guide_field
PRINT*,'*******************************************************'

OPEN(9,FILE="./data/input_params.dat")
WRITE(9,'(27A26)') 'NPROC','NPROC along X','NPROC along Y','NX',&
                   'NY','NC','NP','NT','FDUMP','FSAVE','xmin [cm]',&
                   'xmax [cm]','ymin [cm]','ymax [cm]','dx [cm]','dy [cm]',&
                   'dt [s]','umin','umax','betad','thde','thbe',&
                   'thdi','thbi','density ratio','mass ratio','energy density ratio'
WRITE(9,'(9I26,18E26.16E3)') NPROC,NPX,NPY,NX,NY,PPC,NP,NT,FDUMP,FSAVE,&
                             xmin,xmax,ymin,ymax,dx,dy,dt,umin,umax,betad,&
                             thde,thbe,thdi,thbi,density_ratio,mass_ratio,udens_ratio
CLOSE(9)

OPEN(9,FILE="./data/phys_params.dat")
WRITE(9,'(8A26)') 'delta [cm]','de [cm]','B0 [G]','omegac [s^-1]',&
                  'Drift. dens. [cm^-3]','Uph [erg/cm^-3]','sigma','gamma_rad'
WRITE(9,'(8E26.16E3)') delta,de,B0,c/rhoc,nd0,Uph,sigma,grad
CLOSE(9)

END IF

!=======================================================================
! Generate initial DRIFTING PARTICLES distribution function at t=0
!=======================================================================

! Initial number of particles per species per process
NPP=NCXP*NCYP*PPC

NED=NPP
NPD=NPP

ALLOCATE(pcl_ed(1:7,1:NPP))
ALLOCATE(pcl_data_ed(1:4,1:NPP))

pcl_data_ed=0.0

! DRIFTING ELECTRONS
CALL SET_DRIFT_MAXWELLIAN(-1d0,pcl_ed,thde,upe,gFpe,pse,gFse,xminp,yminp,ND,NPP)

ALLOCATE(pcl_pd(1:7,1:NPP))
ALLOCATE(pcl_data_pd(1:4,1:NPP))

pcl_data_pd=0.0

! DRIFTING IONS
CALL SET_DRIFT_MAXWELLIAN(1d0,pcl_pd,thdi,upp,gFpp,psp,gFsp,xminp,yminp,ND,NPP)

! Weight of drifting particles
CALL WEIGHT(pcl_ed,delta,NPP)
pcl_ed(7,:)=pcl_ed(7,:)*nd0*(xmax-xmin)*(ymax-ymin)/NP

CALL WEIGHT(pcl_pd,delta,NPP)
pcl_pd(7,:)=pcl_pd(7,:)*nd0*(xmax-xmin)*(ymax-ymin)/NP

! Tag of drifting particles
ALLOCATE(taged(1:NPP),tagpd(1:NPP))

CALL SET_TAG(taged,id,NPP)
CALL SET_TAG(tagpd,id,NPP)

!=======================================================================
! Generate initial BACKGROUND PARTICLES distribution function at t=0
!=======================================================================

n0=density_ratio*nd0

NEB=NPP
NPB=NPP

ALLOCATE(pcl_eb(1:7,1:NPP))
ALLOCATE(pcl_data_eb(1:4,1:NPP))

pcl_data_eb=0.0

! Set initial BACKGROUND ELECTRONS distribution function
CALL SET_MAXWELLIAN(pcl_eb,thbe,xminp,yminp,NPP)

ALLOCATE(pcl_pb(1:7,1:NPP))
ALLOCATE(pcl_data_pb(1:4,1:NPP))

pcl_data_pb=0.0

! Set initial BACKGROUND IONS distribution function
CALL SET_MAXWELLIAN(pcl_pb,thbi,xminp,yminp,NPP)

! Weight of background particles
pcl_eb(7,:)=density_ratio*nd0*(xmax-xmin)*(ymax-ymin)/NP
pcl_pb(7,:)=density_ratio*nd0*(xmax-xmin)*(ymax-ymin)/NP

! Tag background particles
ALLOCATE(tageb(1:NPP),tagpb(1:NPP))

CALL SET_TAG(tageb,id,NPP)
CALL SET_TAG(tagpb,id,NPP)
  
! Tag particles for tracking
IF (TRACK_PARTICLES) THEN
  CALL ASSIGN_TRACK_TAGS(electrons, tageb, TOPO_COMM)
  CALL ASSIGN_TRACK_TAGS(ions, tagpb, TOPO_COMM)
ENDIF

IF (id == 0) THEN
  t1 = MPI_WTIME()
  PRINT *, "Particle/tracking initialization time: ", (t1-t0), "seconds"
  t0 = MPI_WTIME()
ENDIF

!=======================================================================
! CHARGES AND (HALF-)CURRENT DENSITIES at t=0 at THE NODES
!=======================================================================

! DRIFTING ELECTRONS
CALL RHOJ(-1d0,pcl_ed,rhoed,Jxed,Jyed,Jzed,xgp,ygp,NPP,id,ngh,TOPO_COMM,ierr)

! DRIFTING IONS
CALL RHOJ(1d0,pcl_pd,rhopd,Jxpd,Jypd,Jzpd,xgp,ygp,NPP,id,ngh,TOPO_COMM,ierr)

! BACKGROUND ELECTRONS
CALL RHOJ(-1d0,pcl_eb,rhoeb,Jxeb,Jyeb,Jzeb,xgp,ygp,NPP,id,ngh,TOPO_COMM,ierr)

! BACKGROUND IONS
CALL RHOJ(1d0,pcl_pb,rhopb,Jxpb,Jypb,Jzpb,xgp,ygp,NPP,id,ngh,TOPO_COMM,ierr)

! TOTAL CHARGES
rho=rhoed+rhopd+rhoeb+rhopb

! TOTAL CURRENT DENSITIES AT t=0
Jx0=2.0*(Jxed+Jxeb+Jxpd+Jxpb)
Jy0=2.0*(Jyed+Jyeb+Jypd+Jypb)
Jz0=2.0*(Jzed+Jzeb+Jzpd+Jzpb)

!=======================================================================
! CORRECTION OF THE ELECTRIC FIELD
!=======================================================================

! Solve Poisson equation to ensure that div(E)=4*pi*rho
CALL CORRECT_EFIELD(Ex,Ey,xgp,ygp,rho,id,ngh,TOPO_COMM,ierr)

! From Yee to Nodes
CALL FIELDS_NODES(Bx,By,Bz,Ex,Ey,Ez,Bxg,Byg,Bzg,Exg,Eyg,Ezg,xgp,ygp,&
                  id,ngh,TOPO_COMM,ierr)

!=======================================================================
! ANALYZING THE INITIAL DATA at t=0
!=======================================================================

! Total magnetic & electric energies
CALL EM_ENERGY(Bxg,Byg,Bzg,Exg,Eyg,Ezg,xgp,ygp,id,TOPO_COMM,ierr)

! Particles' kinetic energy, electrons and ions, bg and drift. particles
CALL KIN_ENERGY(me,pcl_eb,NEB,electrons,bg,id,TOPO_COMM,ierr)
CALL KIN_ENERGY(me,pcl_ed,NED,electrons,drift,id,TOPO_COMM,ierr)
CALL KIN_ENERGY(mi,pcl_pb,NPB,ions,bg,id,TOPO_COMM,ierr)
CALL KIN_ENERGY(mi,pcl_pd,NPD,ions,drift,id,TOPO_COMM,ierr)

! Radiative energy losses, background and drifting particles
! Initialization (No losses)

! Synchrotron
Esyneb=0.0
Esynpb=0.0
Esyned=0.0
Esynpd=0.0
Esyn_electrons=0.0
Esyn_ions=0.0

CALL RAD_ENERGY(Esyneb,Esyned,Esyn_electrons,electrons,syn,id,TOPO_COMM,ierr)
CALL RAD_ENERGY(Esynpb,Esynpd,Esyn_ions,ions,syn,id,TOPO_COMM,ierr)

! Inverse Compton
Eicseb=0.0
Eicspb=0.0
Eicsed=0.0
Eicspd=0.0
Eics_electrons=0.0
Eics_ions=0.0

CALL RAD_ENERGY(Eicseb,Eicsed,Eics_electrons,electrons,ics,id,TOPO_COMM,ierr)
CALL RAD_ENERGY(Eicspb,Eicspd,Eics_ions,ions,ics,id,TOPO_COMM,ierr)

CALL DUMP_GRID("field", xfrg, yfrg, TOPO_COMM)

IF (WRITE_FIELDS.EQV..TRUE.) THEN
! Write fields to disk
CALL DUMP_REDUCED_FIELDS(Exg,Eyg,Ezg,Bxg,Byg,Bzg,xgp,ygp,xfrg,yfrg, &
  localFRedGridShape, localFRedGridOffset, 0,TOPO_COMM)
END IF

IF (WRITE_RHOJ.EQV..TRUE.) THEN
! Write charge and current densities to disk
CALL DUMP_REDUCED_RHOJ(rhoed+rhoeb,2.0*Jxed+2.0*Jxeb,&
  2.0*Jyed+2.0*Jyeb,2.0*Jzed+2.0*Jzeb,&
               xgp,ygp,xfrg,yfrg, localFRedGridShape, localFRedGridOffset,&
               electrons,0,TOPO_COMM)
CALL DUMP_REDUCED_RHOJ(rhopd+rhopb,2.0*Jxpd+2.0*Jxpb,&
  2.0*Jypd+2.0*Jypb,2.0*Jzpd+2.0*Jzpb,&
               xgp,ygp,xfrg,yfrg, localFRedGridShape, localFRedGridOffset,&
               ions,0,TOPO_COMM)
END IF

IF (WRITE_PARTICLES.EQV..TRUE.) THEN
! Write raw particle data to disk
CALL DUMP_PARTICLES(pcl_eb,tageb,NPP,0,electrons,bg,id, .FALSE., TOPO_COMM)
CALL DUMP_PARTICLES(pcl_ed,taged,NPP,0,electrons,drift,id, .FALSE.,TOPO_COMM)
CALL DUMP_PARTICLES(pcl_pb,tagpb,NPP,0,ions,bg,id, .FALSE.,TOPO_COMM)
CALL DUMP_PARTICLES(pcl_pd,tagpd,NPP,0,ions,drift,id, .FALSE.,TOPO_COMM)
END IF

IF (ANALYZE_DISTRIBUTIONS.EQV..TRUE.) THEN
! Spectrum and angular distributions (bottom layer only)
CALL SPECTRUM_ANGULAR(pcl_ed,NPP,0,electrons,drift,id,TOPO_COMM,ierr)
CALL SPECTRUM_ANGULAR(pcl_pd,NPP,0,ions,drift,id,TOPO_COMM,ierr)
CALL SPECTRUM_ANGULAR(pcl_eb,NPP,0,electrons,bg,id,TOPO_COMM,ierr)
CALL SPECTRUM_ANGULAR(pcl_pb,NPP,0,ions,bg,id,TOPO_COMM,ierr)
END IF

CALL DUMP_GRID("", xdrg, ydrg, TOPO_COMM)

IF (ANALYZE_DENSITIES.EQV..TRUE.) THEN
! Spatial distribution
CALL MAP_XY(pcl_eb,NPP, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
  0,electrons,bg,TOPO_COMM)
CALL MAP_XY(pcl_ed,NPP, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
  0,electrons,drift,TOPO_COMM)
CALL MAP_XY(pcl_pb,NPP, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
  0,ions,bg,TOPO_COMM)
CALL MAP_XY(pcl_pd,NPP, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
  0,ions,drift,TOPO_COMM)
END IF

IF (ANALYZE_FLUID.EQV..TRUE.) THEN
! Relativistic macroscopic fluid quantities
CALL MAP_FLUID(me,pcl_eb,NPP, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
  0,electrons,bg,TOPO_COMM)
CALL MAP_FLUID(me,pcl_ed,NPP, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
  0,electrons,drift,TOPO_COMM)
CALL MAP_FLUID(mi,pcl_pb,NPP, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
  0,ions,bg,TOPO_COMM)
CALL MAP_FLUID(mi,pcl_pd,NPP, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
  0,ions,drift,TOPO_COMM)
END IF

IF (ANALYZE_RADIATION.EQV..TRUE.) THEN
! ELECTRON synchrotron radiation spectrum and angular distribution
CALL ANALYSIS_SYNC(me,pcl_eb,Bxg,Byg,Bzg,xgp,ygp,0,electrons,bg,NPP,&
                   id,TOPO_COMM,ierr)
END IF
               
IF (id.EQ.0) THEN
  t3=MPI_WTIME()
  PRINT*,'Initial dump time:',t3-t0,'seconds'
  t0=t3
END IF

!=======================================================================
! INITIALISATION OF THE PARTICLE MOTION
!=======================================================================

! DRIFTING ELECTRONS 4-velocity at t=0-dt/2
CALL INITIAL_PUSH(-1d0,me,pcl_ed,Bxg,Byg,Bzg,Exg,Eyg,Ezg,Uph,xgp,ygp,NPP)

! DRIFTING IONS 4-velocity at t=0-dt/2
CALL INITIAL_PUSH(1d0,mi,pcl_pd,Bxg,Byg,Bzg,Exg,Eyg,Ezg,Uph,xgp,ygp,NPP)

! BACKGROUND ELECTRONS 4-velocity at t=0-dt/2
CALL INITIAL_PUSH(-1d0,me,pcl_eb,Bxg,Byg,Bzg,Exg,Eyg,Ezg,Uph,xgp,ygp,NPP)

! BACKGROUND IONS 4-velocity at t=0-dt/2
CALL INITIAL_PUSH(1d0,mi,pcl_pb,Bxg,Byg,Bzg,Exg,Eyg,Ezg,Uph,xgp,ygp,NPP)

!=======================================================================

it0=0

IF (id.EQ.0 .AND. .NOT. RESTORE) THEN
PRINT*,it0,'/',NT
END IF

END IF

!=======================================================================
! RESTORATION OF A SIMULATION
!=======================================================================

IF (RESTORE.EQV..TRUE.) THEN

IF (id.EQ.0) THEN
  t0=MPI_WTIME()
END IF

IF (id.EQ.0) THEN
PRINT*,'*******************************************************'
PRINT*,'1. BEGIN RESTORATION OF THE SIMULATION'
END IF

it0=time_ref

! Restore the particle distribution functions
CALL RESTORE_PARTICLES(pcl_ed,pcl_data_ed,taged,NED,it0,electrons,drift,&
                       id,TOPO_COMM,ierr)
CALL RESTORE_PARTICLES(pcl_eb,pcl_data_eb,tageb,NEB,it0,electrons,bg,&
                       id,TOPO_COMM,ierr)
CALL RESTORE_PARTICLES(pcl_pd,pcl_data_pd,tagpd,NPD,it0,ions,drift,&
                       id,TOPO_COMM,ierr)
CALL RESTORE_PARTICLES(pcl_pb,pcl_data_pb,tagpb,NPB,it0,ions,bg,&
                       id,TOPO_COMM,ierr)

! Restore the Yee fields
CALL RESTORE_FIELDS(Ex,Ey,Ez,Bx,By,Bz,it0,coords,id,TOPO_COMM,ierr)

! Restore the nodal fields
CALL FIELDS_NODES(Bx,By,Bz,Ex,Ey,Ez,Bxg,Byg,Bzg,Exg,Eyg,Ezg,xgp,ygp,&
                  id,ngh,TOPO_COMM,ierr)

IF (id.EQ.0) THEN
PRINT*,'2. SESSION RESTORED SUCCESSFULLY!'
PRINT*,'*******************************************************'
PRINT*,it0,'/',it0+NT
END IF

IF (id.EQ.0) THEN
  t3=MPI_WTIME()
  PRINT*,'Restoration time:',t3-t0,'seconds'
  t0=t3
END IF

END IF

!=======================================================================
! BEGINNING OF THE LOOP ON TIME
!=======================================================================

IF (id.EQ.0) THEN
! To measure time in seconds
t1=MPI_WTIME()
END IF



!=======================================================================
! START MAIN LOOP
!=======================================================================

DO it=it0+1,it0+NT

  !=======================================================================
  ! DRIFTING ELECTRONS
  !=======================================================================
  
  ! Update u from t-dt/2 to t+dt/2, particle positions unchanged
  CALL BORIS_PUSH(-1d0,me,pcl_ed,pcl_data_ed,Bxg,Byg,Bzg,Exg,Eyg,Ezg,&
                  Uph,xgp,ygp,Esyned,Eicsed,NED)
  
  ! Computation of rho at t and half of the current density
  CALL RHOJ(-1d0,pcl_ed,rhoed,Jx0,Jy0,Jz0,xgp,ygp,NED,id,ngh,TOPO_COMM,ierr)
  
  ! Push the particles from t to t+dt
  CALL PUSH_PARTICLES(pcl_ed,NED)
  
  ! Applying boundary conditions to the particles
  CALL BOUNDARIES_PARTICLES(pcl_ed,pcl_data_ed,taged,NED)

  ! Counting the particles leaving each subdomain
  CALL COUNT_ESCAPE(pcl_ed,xminp,xmaxp,yminp,ymaxp,NED,NESC)
                      
  ! Exchange of particles at the boundaries between processes
  CALL COM_PARTICLES(pcl_ed,pcl_data_ed,taged,xminp,xmaxp,yminp,ymaxp,&
                     NED,NESC,id,ngh,TOPO_COMM,ierr)
                     
  ! Computation of rho at t+dt and the other half of the current density
  CALL RHOJ(-1d0,pcl_ed,rhoed,Jxed,Jyed,Jzed,xgp,ygp,NED,id,ngh,TOPO_COMM,ierr)

  ! Total current density at t+dt/2
  Jxed=Jxed+Jx0
  Jyed=Jyed+Jy0
  Jzed=Jzed+Jz0
  
  !=======================================================================
  ! DRIFTING IONS
  !=======================================================================
  
  ! Update u from t-dt/2 to t+dt/2, particle positions unchanged
  CALL BORIS_PUSH(1d0,mi,pcl_pd,pcl_data_pd,Bxg,Byg,Bzg,Exg,Eyg,Ezg,&
                  Uph,xgp,ygp,Esynpd,Eicspd,NPD)
  
  ! Computation of rho at t and half of the current density
  CALL RHOJ(1d0,pcl_pd,rhopd,Jx0,Jy0,Jz0,xgp,ygp,NPD,id,ngh,TOPO_COMM,ierr)
  
  ! Push the particles from t to t+dt
  CALL PUSH_PARTICLES(pcl_pd,NPD)
  
  ! Applying boundary conditions to the particles
  CALL BOUNDARIES_PARTICLES(pcl_pd,pcl_data_pd,tagpd,NPD)
  
  ! Counting the particles leaving each subdomain
  CALL COUNT_ESCAPE(pcl_pd,xminp,xmaxp,yminp,ymaxp,NPD,NESC)
  
  ! Exchange of particles at the boundaries between processes                   
  CALL COM_PARTICLES(pcl_pd,pcl_data_pd,tagpd,xminp,xmaxp,yminp,ymaxp,&
                     NPD,NESC,id,ngh,TOPO_COMM,ierr)
                     
  ! Computation of rho at t+dt and the other half of the current density
  CALL RHOJ(1d0,pcl_pd,rhopd,Jxpd,Jypd,Jzpd,xgp,ygp,NPD,id,ngh,TOPO_COMM,ierr)
  
  ! Total current density at t+dt/2
  Jxpd=Jxpd+Jx0
  Jypd=Jypd+Jy0
  Jzpd=Jzpd+Jz0
  
  !=======================================================================
  ! BACKGROUND ELECTRONS
  !=======================================================================
  
  ! Update u from t-dt/2 to t+dt/2, particle positions unchanged
  CALL BORIS_PUSH(-1d0,me,pcl_eb,pcl_data_eb,Bxg,Byg,Bzg,Exg,Eyg,Ezg,&
                  Uph,xgp,ygp,Esyneb,Eicseb,NEB)

  ! Computation of rho at t and half of the current density
  CALL RHOJ(-1d0,pcl_eb,rhoeb,Jx0,Jy0,Jz0,xgp,ygp,NEB,id,ngh,TOPO_COMM,ierr)
  
  ! Push the particles from t to t+dt
  CALL PUSH_PARTICLES(pcl_eb,NEB)
                      
  ! Applying boundary conditions to the particles
  CALL BOUNDARIES_PARTICLES(pcl_eb,pcl_data_eb,tageb,NEB)

  ! Counting the particles leaving each subdomain
  CALL COUNT_ESCAPE(pcl_eb,xminp,xmaxp,yminp,ymaxp,NEB,NESC)
  
  ! Exchange of particles at the boundaries between processes
  CALL COM_PARTICLES(pcl_eb,pcl_data_eb,tageb,xminp,xmaxp,yminp,ymaxp,&
                     NEB,NESC,id,ngh,TOPO_COMM,ierr)
                     
  ! Computation of rho at t+dt and the other half of the current density
  CALL RHOJ(-1d0,pcl_eb,rhoeb,Jxeb,Jyeb,Jzeb,xgp,ygp,NEB,id,ngh,TOPO_COMM,ierr)
  
  ! Total current density at t+dt/2
  Jxeb=Jxeb+Jx0
  Jyeb=Jyeb+Jy0
  Jzeb=Jzeb+Jz0

  !=======================================================================
  ! BACKGROUND IONS
  !=======================================================================
  
  ! Update u from t-dt/2 to t+dt/2, particle positions unchanged
  CALL BORIS_PUSH(1d0,mi,pcl_pb,pcl_data_pb,Bxg,Byg,Bzg,Exg,Eyg,Ezg,&
                  Uph,xgp,ygp,Esynpb,Eicspb,NPB)
  
  ! Computation of rho at t and half of the current density
  CALL RHOJ(1d0,pcl_pb,rhopb,Jx0,Jy0,Jz0,xgp,ygp,NPB,id,ngh,TOPO_COMM,ierr)
  
  ! Push the particles from t to t+dt
  CALL PUSH_PARTICLES(pcl_pb,NPB)
  
  ! Applying boundary conditions to the particles
  CALL BOUNDARIES_PARTICLES(pcl_pb,pcl_data_pb,tagpb,NPB)
                      
  ! Counting the particles leaving each subdomain
  CALL COUNT_ESCAPE(pcl_pb,xminp,xmaxp,yminp,ymaxp,NPB,NESC)
                      
  ! Exchange of particles at the boundaries between processes
  CALL COM_PARTICLES(pcl_pb,pcl_data_pb,tagpb,xminp,xmaxp,yminp,ymaxp,&
                     NPB,NESC,id,ngh,TOPO_COMM,ierr)

  ! Computation of rho at t+dt and the other half of the current density
  CALL RHOJ(1d0,pcl_pb,rhopb,Jxpb,Jypb,Jzpb,xgp,ygp,NPB,id,ngh,TOPO_COMM,ierr)
  
  ! Total current density at t+dt/2
  Jxpb=Jxpb+Jx0
  Jypb=Jypb+Jy0
  Jzpb=Jzpb+Jz0
                     
  !=======================================================================
  ! TOTAL CHARGES AND CURRENT DENSITIES
  !=======================================================================

  ! TOTAL CHARGES AT t=t+dt
  rho=rhoed+rhopd+rhoeb+rhopb

  ! TOTAL CURRENT DENSITIES AT t=t+dt/2
  Jx=Jxed+Jxeb+Jxpd+Jxpb
  Jy=Jyed+Jyeb+Jypd+Jypb
  Jz=Jzed+Jzeb+Jzpd+Jzpb

  ! Transformation on the Yee grid
  CALL JYEE(Jx,Jy,Jz,xgp,ygp,id,ngh,TOPO_COMM,ierr)

  !=======================================================================
  ! B FIELD at t=t+dt/2
  !=======================================================================

  CALL PUSH_BHALF(Bx,By,Bz,Ex,Ey,Ez,xgp,ygp,id,ngh,TOPO_COMM,ierr)

  !=======================================================================
  ! SOLVE MAXWELL'S EQUATIONS at t=t+dt
  !=======================================================================
  
  ! E FIELD at t=t+dt
  CALL PUSH_EFIELD(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,xgp,ygp,id,ngh,TOPO_COMM,ierr)

  IF (MOD(it,FREQ_POISSON).EQ.0) THEN
  ! Solve Poisson equation to ensure that div(E)=4*pi*rho
  CALL CORRECT_EFIELD(Ex,Ey,xgp,ygp,rho,id,ngh,TOPO_COMM,ierr)
  END IF
  
  ! B FIELD at t=t+dt
  CALL PUSH_BHALF(Bx,By,Bz,Ex,Ey,Ez,xgp,ygp,id,ngh,TOPO_COMM,ierr)
    
  CALL FIELDS_NODES(Bx,By,Bz,Ex,Ey,Ez,Bxg,Byg,Bzg,Exg,Eyg,Ezg,xgp,ygp,&
                   id,ngh,TOPO_COMM,ierr)
  
  !=======================================================================
  ! ANALYZING THE DATA
  !=======================================================================

  ! Total magnetic & electric energies
  CALL EM_ENERGY(Bxg,Byg,Bzg,Exg,Eyg,Ezg,xgp,ygp,id,TOPO_COMM,ierr)

  ! Particles' kinetic energy, electrons and ions, bg and drift. particle
  CALL KIN_ENERGY(me,pcl_eb,NEB,electrons,bg,id,TOPO_COMM,ierr)
  CALL KIN_ENERGY(me,pcl_ed,NED,electrons,drift,id,TOPO_COMM,ierr)
  CALL KIN_ENERGY(mi,pcl_pb,NPB,ions,bg,id,TOPO_COMM,ierr)
  CALL KIN_ENERGY(mi,pcl_pd,NPD,ions,drift,id,TOPO_COMM,ierr)
  
  ! Radiative energy losses, background and drifting particles
  ! Synchrotron
  CALL RAD_ENERGY(Esyneb,Esyned,Esyn_electrons,electrons,syn,id,TOPO_COMM,ierr)
  CALL RAD_ENERGY(Esynpb,Esynpd,Esyn_ions,ions,syn,id,TOPO_COMM,ierr)
  
  ! Inverse Compton
  CALL RAD_ENERGY(Eicseb,Eicsed,Eics_electrons,electrons,ics,id,TOPO_COMM,ierr)
  CALL RAD_ENERGY(Eicspb,Eicspd,Eics_ions,ions,ics,id,TOPO_COMM,ierr)
  
  !=======================================================================
  ! PARTICLE TRACKER
  !=======================================================================

  IF ((TRACK_PARTICLES) .AND. (it >= time_track) &
    .AND. (MOD(it, FTRACK)==0)) THEN
  
  ! Call the particle tracker to write orbits to disk 
  CALL TRACKER(pcl_eb,pcl_data_eb,tageb,NEB,electrons, it)
  CALL TRACKER(pcl_pb,pcl_data_pb,tagpb,NPB,ions, it)
  
  ! Dump tracker data to disk if it's time to do that.
  CALL WRITE_TRACK_DATA(.FALSE., TOPO_COMM)
  
  END IF

  !=======================================================================
  ! Figure out if we need to dump/save/quit
  !=======================================================================
  
  doQuit = 0
  doDump = 0
  doCheckpt = 0

  IF (MOD(it,FDUMP).EQ.0) THEN
    doDump = 1
  ENDIF
  
  ! Always checkpoint after last step, if CHECKPOINT is .TRUE.
  IF (it == it0 + NT .AND. CHECKPOINT) THEN
    doCheckpt = 1
  ENDIF

  ! Because we have to (MPI) communicate to decide whether to checkpoint,
  ! consider doing it only every FGENCOMM steps.  But: also consider
  ! checkpointing at any dump step, because (1) we'd prefer to checkpoint
  ! at the same time as dumping, and (2) that always offers us the
  ! opportunity to checkpoint before dumping.
  
  IF (MOD(it, FGENCOMM) == 0 .OR. (doDump>0)) THEN
    IF (CHECKPOINT .AND. FSAVE > 0) THEN
    
      ! Because different ranks may have slightly different times, 
      ! rank 0 has to decide and communicate to others.
      IF (id.EQ.0) THEN
      
        tCheckpt=MPI_WTIME()
        elapsedTime = tCheckpt - initTime
        
        ! The checkpoint needs to be completed before nextCheckPtTime
        nextCheckPtTime = (numCheckptSaves + 1) * FSAVE 
        !print *, "nextCheckPtTime =", nextCheckPtTime
        
        ! Try to estimate the start time with a 1.5 safety factor
        nextCheckPtTime = nextCheckPtTime - 1.5*checkptTime &
          - (FGENCOMM-1)*avgStepTime
        !print *, "nextCheckPtTime =", nextCheckPtTime
        
        IF (doDump > 0) THEN
        
          ! Consider that after the dump, elapsedTime maybe be greater
          ! than nextCheckPtTime, so we should checkpt now.
          nextCheckPtTime = nextCheckPtTime - dumpTime 
          !print *, "(at dump) nextCheckPtTime =", nextCheckPtTime
          
          ! Furthermore, we'd rather checkpoint and dump at the same
          ! time, if it doesn't alter the checkpoint time by too much,
          ! and if checkpointing at the next dump will be too late.
          nextCheckPtTime = nextCheckPtTime &
            - MIN(FSAVE_SHIFT_FRAC*FSAVE, FDUMP*avgStepTime)  
          !print *, "(at dump) nextCheckPtTime =", nextCheckPtTime
          
        ENDIF
        
        !print *, "elapsedTime =", elapsedTime
        
        IF (elapsedTime > nextCheckPtTime) THEN
          doCheckpt = 1
        ENDIF
        
      END IF
      
    ENDIF
    
    IF (id == 0) THEN
    
      ! To limit communication, combine communication of whether to
      ! checkpt with any user-issued commands 
      CALL READ_COMMAND("data/simCmd.txt", cmdStep, extCmd)
      
      doStuffComm(1) = doCheckpt
      doStuffComm(2) = cmdStep
      doStuffComm(3) = extCmd
      
    ENDIF
    
    CALL MPI_BCAST(doStuffComm, 3, MPI_INTEGER, 0, TOPO_COMM, ierr)
    
    doCheckpt = doStuffComm(1)
    cmdStep = doStuffComm(2)
    extCmd = doStuffComm(3)
    
  ENDIF

  IF (it == cmdStep .AND. extCmd /= EXTCMD_DO_NOTHING) THEN
  
    IF (IAND(extCmd, EXTCMD_QUIT) /= 0) doQuit = 1
    IF (IAND(extCmd, EXTCMD_DUMP) /= 0) doDump = 1
    IF (IAND(extCmd, EXTCMD_SAVE) /= 0) doCheckpt = 1
    
  ENDIF
  
  !=======================================================================
  ! CHECKPOINT: SAVING THE SIMULATION STATE
  !=======================================================================
  
  IF (doCheckpt > 0) THEN
    
    IF (id.EQ.0) THEN
      t0=MPI_WTIME()
      PRINT *, "Checkpoint: saving state at step", it, &
        "after ", elapsedTime, "seconds"
    END IF

    CALL WRITE_LOGS(id)
    IF (id == 0) CALL WRITE_LOG_TO_DISK(timeLog)
  
    ! Force saving of tracked-particle data (so that we can begin new
    ! file after this)
    IF (TRACK_PARTICLES) THEN
    
      CALL WRITE_TRACK_DATA(.TRUE., TOPO_COMM)
      CALL INC_TRACK_SEQ_NUM()
      
    ENDIF
    
    ! Save the Yee fields to disk
    CALL SAVE_FIELDS(Ex,Ey,Ez,Bx,By,Bz,it,coords,id,TOPO_COMM,ierr)
    
    ! Save the particle data to disk
    CALL SAVE_PARTICLES(pcl_eb,pcl_data_eb,tageb,NEB,it,electrons,bg,&
                        id,TOPO_COMM,ierr)
    CALL SAVE_PARTICLES(pcl_ed,pcl_data_ed,taged,NED,it,electrons,drift,&
                        id,TOPO_COMM,ierr)
    CALL SAVE_PARTICLES(pcl_pb,pcl_data_pb,tagpb,NPB,it,ions,bg,&
                        id,TOPO_COMM,ierr)
    CALL SAVE_PARTICLES(pcl_pd,pcl_data_pd,tagpd,NPD,it,ions,drift,&
                        id,TOPO_COMM,ierr)
    
    IF (id.EQ.0) THEN
      ! To measure time in seconds
      t3=MPI_WTIME()
      IF (numCheckptSaves > 0) THEN
        checkptTime = 0.5 * (t3-t0) + 0.5 * checkptTime
      ELSE
        checkptTime = t3-t0
      ENDIF
      PRINT*,'Checkpoint time:', t3-t0, 'seconds'
      !print*,'Est. checkpoint time:', checkptTime, ' seconds'
    END IF

    ! GRW: I think this is too dangerous...sometimes file system glitches
    ! cause extremely long save times; we might not want to halt a 
    ! simulation just because this happens once.
    !
    ! If FSAVE <= checkptTime, this simulation is going to spend all
    ! its time saving stuff, so we might as well quit.
    ! But make sure we've saved everything before quitting
    !IF (id.EQ.0) THEN
    !  IF (checkptTime > 0.9 * FSAVE) THEN
    !    PRINT *, "Quitting early because the time to save the simulation state"
    !    PRINT *, "is greater than 80% of the time (FSAVE) between checkpoints,"
    !    PRINT *, "so the simulation would spend all its time checkpointing."
    !    PRINT *, "time to save the simulation state =", checkptTime, " s"
    !    PRINT *, "FSAVE =", FSAVE, " s"
    !    doQuit = 1
    !  ENDIF
    !ENDIF
    !CALL MPI_BCAST(doQuit, 1, MPI_INTEGER, 0, TOPO_COMM, ierr)

    numCheckptSaves = numCheckptSaves + 1

  END IF

  !=======================================================================
  ! Dumping and analyzing data every FDUMP timesteps
  !=======================================================================
  
  IF (doDump > 0) THEN

  IF (id.EQ.0) THEN
    t0=MPI_WTIME()
  END IF
  
  CALL WRITE_LOGS(id)
  IF (id == 0) CALL WRITE_LOG_TO_DISK(timeLog)

  IF (WRITE_FIELDS.EQV..TRUE.) THEN
  ! Write fields to disk
  CALL DUMP_REDUCED_FIELDS(Exg,Eyg,Ezg,Bxg,Byg,Bzg,xgp,ygp,xfrg,yfrg, &
    localFRedGridShape, localFRedGridOffset, it,TOPO_COMM)
  END IF

  IF (WRITE_RHOJ.EQV..TRUE.) THEN
  ! Write charge and current densities to disk
  CALL DUMP_REDUCED_RHOJ(rhoed+rhoeb,Jxed+Jxeb,Jyed+Jyeb,Jzed+Jzeb,&
               xgp,ygp,xfrg,yfrg, localFRedGridShape, localFRedGridOffset,&
               electrons,it,TOPO_COMM)
  CALL DUMP_REDUCED_RHOJ(rhopd+rhopb,Jxpd+Jxpb,Jypd+Jypb,Jzpd+Jzpb,&
               xgp,ygp,xfrg,yfrg, localFRedGridShape, localFRedGridOffset,&
               ions,it,TOPO_COMM)
  END IF
  
  IF (WRITE_PARTICLES.EQV..TRUE.) THEN
  ! Write raw particle data to disk
  CALL DUMP_PARTICLES(pcl_eb,tageb,NEB,it,electrons,bg,id, .FALSE.,TOPO_COMM)
  CALL DUMP_PARTICLES(pcl_ed,taged,NED,it,electrons,drift,id, .FALSE.,TOPO_COMM)
  CALL DUMP_PARTICLES(pcl_pb,tagpb,NPB,it,ions,bg,id, .FALSE.,TOPO_COMM)
  CALL DUMP_PARTICLES(pcl_pd,tagpd,NPD,it,ions,drift,id, .FALSE.,TOPO_COMM)
  END IF
  
  IF (ANALYZE_DISTRIBUTIONS.EQV..TRUE.) THEN
  ! Spectrum and angular distribution (bottom layer only)
  CALL SPECTRUM_ANGULAR(pcl_ed,NED,it,electrons,drift,id,TOPO_COMM,ierr)
  CALL SPECTRUM_ANGULAR(pcl_pd,NPD,it,ions,drift,id,TOPO_COMM,ierr)
  CALL SPECTRUM_ANGULAR(pcl_eb,NEB,it,electrons,bg,id,TOPO_COMM,ierr)
  CALL SPECTRUM_ANGULAR(pcl_pb,NPB,it,ions,bg,id,TOPO_COMM,ierr)
  END IF
  
  IF (ANALYZE_DENSITIES.EQV..TRUE.) THEN
  ! Spatial distribution
  CALL MAP_XY(pcl_eb,NEB, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
    it,electrons,bg,TOPO_COMM)
  CALL MAP_XY(pcl_ed,NED, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
    it,electrons,drift,TOPO_COMM)
  CALL MAP_XY(pcl_pb,NPB, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
    it,ions,bg,TOPO_COMM)
  CALL MAP_XY(pcl_pd,NPD, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
    it,ions,drift,TOPO_COMM)
  END IF

  IF (ANALYZE_FLUID.EQV..TRUE.) THEN
  ! Relativistic macroscopic fluid quantities
  CALL MAP_FLUID(me,pcl_eb,NEB, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
    it,electrons,bg,TOPO_COMM)
  CALL MAP_FLUID(me,pcl_ed,NED, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
    it,electrons,drift,TOPO_COMM)
  CALL MAP_FLUID(mi,pcl_pb,NPB, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
    it,ions,bg,TOPO_COMM)
  CALL MAP_FLUID(mi,pcl_pd,NPD, xdrg, ydrg, localDRedGridShape, localDRedGridOffset,&
    it,ions,drift,TOPO_COMM)
  END IF
  
  IF (ANALYZE_RADIATION.EQV..TRUE.) THEN
  ! ELECTRON synchrotron radiation spectrum and angular distribution
  CALL ANALYSIS_SYNC(me,pcl_eb,Bxg,Byg,Bzg,xgp,ygp,it,electrons,bg,NEB,&
                     id,TOPO_COMM,ierr)
  END IF
                         
  IF (id.EQ.0) THEN
    ! To measure time in seconds
    t3=MPI_WTIME()
    IF (numDumps > 0) THEN
      dumpTime = 0.5 * (t3-t0) + 0.5 * dumpTime
    ELSE
      dumpTime = t3-t0
    ENDIF
    PRINT*,'Dumping time:', t3-t0, 'seconds'
    !print*,'Est. dump time:', dumpTime, ' seconds'
  END IF
  numDumps = numDumps + 1
  
  END IF

  !=======================================================================
  
  IF (id.EQ.0) THEN

  PRINT*,it,'/',it0+NT

  ! To measure time in seconds
  t2=MPI_WTIME()
  tStep = t2 - t1
  t1 = t2

  PRINT*,'Computing time:',tStep,'seconds'
  avgStepCoef = 1./MIN(it-it0, FDUMP)
  avgStepTime = avgStepCoef*tStep + (1.-avgStepCoef)*avgStepTime
  ! :TODO: remove following printout
  !PRINT*,'Est. time/step:', avgStepTime, ' seconds'

  ! Write computation time to disk
  IF (id == 0) CALL ADD_LOG_DATA(timeLog, (/ tStep /) )

  END IF
  
  !=======================================================================

  IF (doQuit > 0) THEN
    IF (id == 0) THEN
      PRINT *, "Quitting early as commanded."
    ENDIF
    EXIT
  ENDIF
  
ENDDO

!=======================================================================
! Wind up and deallocate
!=======================================================================

! Dump remaining tracker data to disk 
IF (TRACK_PARTICLES) THEN
  CALL WRITE_TRACK_DATA(.TRUE., TOPO_COMM)
  CALL END_TRACKER()
ENDIF

! finish witing logs
CALL END_LOGS(id)
IF (id == 0) CALL END_LOG(timeLog)

! Deallocate

DEALLOCATE(xfrg, yfrg, xfrgp, yfrgp)
DEALLOCATE(xdrg, ydrg, xdrgp, ydrgp)


!=======================================================================
! TERMINATE THE MPI ENVIRONMENT
!=======================================================================

CALL MPI_FINALIZE(ierr)

END PROGRAM MAIN
