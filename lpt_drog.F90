! The routines in this file have been adapted from drog2dsp_hybrid.F90.

      SUBROUTINE LPT_Drog_TimeStep

      USE LPT_Comm_Module
      USE LPT_Data_Module, ONLY: DiameterEffective,                    &
     &                           DT1 => TimeStepSizes,                 &
     &                           JJDR => PartElemJ,                    &
     &                           LLDR => PartElemL,                    &
     &                           NDR => NumParticlesLocal,             &
     &                           NOMP_MAP => MapOpenMP,                &
     &                           NumParticlesGlobal,                   &
     &                           NumTrackingSnaps,                     &
     &                           PartDiameterLocal,                    &
     &                           PartDomainLocal,                      &
     &                           ParticleInputMethod,                  &
     &                           ParticleSourceX,                      &
     &                           ParticleSourceY,                      &
     &                           ParticleSourceZ,                      &
     &                           SimulationLength,                     &
     &                           T1 => Time1,                          &
     &                           T2 => Time2,                          &
     &                           XDR => PartLonLocal,                  &
     &                           YDR => PartLatLocal,                  &
     &                           ZDR => PartDepthLocal
      USE LPT_Read_Module, ONLY: LPT_Read_From_Particle_File,          &
     &                           LPT_Read_Update_Velocities
      USE LPT_Write_Module, ONLY: LPT_Write_Particle_Snap

      IMPLICIT NONE

      CHARACTER(LEN=100) :: JunkC
      CHARACTER(LEN=100) :: JunkC1

      INTEGER            :: I
      INTEGER            :: II
      INTEGER            :: IMOD
      INTEGER            :: IT
      INTEGER            :: JJ
      INTEGER            :: JNEW
      INTEGER            :: LL
      INTEGER            :: LNEW
      INTEGER            :: IMAP

      REAL(8)            :: DDT
      REAL(8)            :: DIAM
      REAL(8)            :: DNEW
      REAL(8)            :: Random
      REAL(8)            :: UNEWS
      REAL(8)            :: VNEWS
      REAL(8)            :: WNEWS
      REAL(8)            :: XO
      REAL(8)            :: YO
      REAL(8)            :: ZO

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the LPT_Drog_TimeStep routine.")
#endif

!pck 10/04/13 - Write Initial Locations
      CALL LPT_Write_Particle_Snap
      T1 = T1 + SimulationLength/DBLE(NumTrackingSnaps)
      T2 = T2 + SimulationLength/DBLE(NumTrackingSnaps)

      DO I=1,NumTrackingSnaps

#if VERBOSE > 0
         WRITE(JunkC,'(I24)') I
         WRITE(JunkC1,'(I24)') NumTrackingSnaps
         CALL LPT_Print(0,"INFO","Beginning time snap "                &
     &            //TRIM(ADJUSTL(JunkC))//" out of "                   &
     &            //TRIM(ADJUSTL(JunkC1))//".")
#endif

         ! Update the velocities to the correct time.
         CALL LPT_Read_Update_Velocities(I)

         ! Add particles from the source term, if necessary.
         IF(INDEX(ParticleInputMethod,"READFROMFILE").GT.0)THEN
            CALL LPT_Read_From_Particle_File(I)
         ELSEIF(INDEX(ParticleInputMethod,"SOURCE").GT.0)THEN
            DO II=1,10
               CALL RANDOM_NUMBER(Random)
               IMOD = MOD(I,NumTrackerCores)
               IF(IMOD.EQ.0) IMOD = NumTrackerCores
               DO IT=1,NumTrackerCores
                  IF(MyRank.EQ.TrackerRanks(IT).AND.IMOD.EQ.IT)THEN
                     ! Select a particle size at random.
                     DNEW = 0.000050D0+Random*(0.000300D0-0.000050D0)
!                    DNEW = DiameterEffective
                     CALL LPT_Oil_Source(ParticleSourceX,              &
     &                        ParticleSourceY,ParticleSourceZ,DNEW,    &
     &                        NumParticlesGlobal+1,0)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         ! Loop over the particles.
         IF(AmTrackerCore)THEN
#if VERBOSE > 0
            IF(MyRank.EQ.0)THEN
               IF(SUM(PartDomainLocal).GE.100)THEN
                  WRITE(*,'(A,$)') "LPT: INFO: Tracking the particle " &
     &                  //"locations: +"
               ELSE
                  CALL LPT_Print(0,"INFO","Tracking the particle "     &
     &                     //"locations.")
               ENDIF
            ENDIF
#endif
            !$ TIME_TC = TIME_TC - omp_get_wtime()
            !$OMP PARALLEL DEFAULT (SHARED) PRIVATE(IMAP,II,XO,YO,JJ,UNEWS,VNEWS,JNEW,DDT)
            !$OMP  DO
            DO IMAP=1,NDR

!             II = NOMP_MAP(IMAP)
              II = IMAP
              XO = XDR(II)
              YO = YDR(II)
              ZO = ZDR(II)
              JJ = JJDR(II)
              LL = LLDR(II)
              DDT = DT1(II)
              IF(ALLOCATED(PartDiameterLocal))THEN
                DIAM = PartDiameterLocal(II)
              ELSE
                DIAM = -99999.D0
              ENDIF
!             IF (XO<-999998.d0.and.YO<-999999.d0) CYCLE ! cjt
              IF(PartDomainLocal(II).EQ.0) CYCLE

              ! Determine the components of flow at (XO,YO).
              CALL LPT_Drog_VELS(JJ,LL,XO,YO,ZO,UNEWS,VNEWS,WNEWS,DIAM)

              !TRACK THIS PARTICLE FROM TIME T1 TO T2
              CALL LPT_Drog_TRACK(JJ,JNEW,LL,LNEW,XO,YO,ZO,            &
     &                 UNEWS,VNEWS,WNEWS,T1,T2,DDT,II,DIAM)

              XDR(II) = XO
              YDR(II) = YO
              ZDR(II) = ZO
              JJDR(II) = JNEW
              LLDR(II) = LNEW
              DT1(II) = DDT

#if VERBOSE > 0
               IF(MyRank.EQ.0)THEN
                  IF(SUM(PartDomainLocal).GE.100)THEN
                     CALL LPT_Progress(IMAP,SUM(PartDomainLocal))
                  ENDIF
               ENDIF
#endif

            ENDDO
            !$OMP  END DO
            !$OMP END PARALLEL
            !$ TIME_TC = TIME_TC + omp_get_wtime()
         ENDIF

         ! Write the particle locations to the output file.
         CALL LPT_Write_Particle_Snap

         T1 = T1 + SimulationLength/DBLE(NumTrackingSnaps)
         T2 = T2 + SimulationLength/DBLE(NumTrackingSnaps)

      ENDDO

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Drog_TimeStep routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      ! The following routine is called immediately after the
      ! unstructured mesh is read.  We read the element-to-vertex
      ! connectivity table from the ADCIRC mesh.  This routine creates
      ! the element-to-element and vertex-to-element connectivity
      ! tables, and it also determines the boundary segment information
      ! that is needed to prevent the particles from leaving the domain.

      SUBROUTINE LPT_Drog_MAK_NEINFO

      USE LPT_Data_Module, ONLY: NMEL => NumElems,                     &
     &                           NMND => NumVerts,                     &
     &                           ELEMS => MeshConnEV,                  &
     &                           ICEE => MeshConnEE,                   &
     &                           ICNE => MeshConnVE,                   &
     &                           NSEG => NumBoundarySegments,          &
     &                           IBSEG => BoundarySegmentVerts,        &
     &                           ISEGF => BoundarySegmentFollowing,    &
     &                           ISEGP => BoundarySegmentPreceding,    &
     &                           IBSEGEL => BoundarySegmentElement
      USE LPT_Comm_Module, ONLY: MyRank

      IMPLICIT NONE

      INTEGER:: I, J, M, N, I1, I2, M1, M2, N1, N2, N3, ISEG, maxne

      INTEGER,ALLOCATABLE :: NPROP(:)

      INTEGER, PARAMETER:: MAP(4) = (/1, 2, 3, 1 /)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the LPT_MAK_NEINFO routine.")
#endif

      ALLOCATE(NPROP(1:NMND))

      ! Vertex-to-element adjacency information.
      NPROP(1:NMND) = 0
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Computing the number of "         &
     &        //"elements connected to each vertex: +"
      ENDIF
#endif      
      DO M = 1, NMEL
        DO J = 1, 3
          NPROP(ELEMS(M,J)) = NPROP(ELEMS(M,J)) + 1
        ENDDO
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(M,NMEL)
        ENDIF
#endif
      ENDDO
      maxne = maxval(NPROP(1:NMND)) + 1
      ALLOCATE( ICEE(1:NMEL,1:3), ICNE(1:NMND,1:maxne) )
      ICNE(1:NMND,1:maxne) = 0
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Computing the vertex-to-element " &
     &        //"connectivity table: +"
      ENDIF
#endif      
      DO M = 1, NMEL
        DO J = 1, 3
          ICNE(ELEMS(M,J),1) = ICNE(ELEMS(M,J),1) + 1
          ICNE(ELEMS(M,J),ICNE(ELEMS(M,J),1)+1) = M
        ENDDO
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(M,NMEL)
        ENDIF
#endif
      ENDDO

      ! Element-to-element adjacency information. 
      !       1
      !       |\
      !    (1)| \(3)
      !       |__\
      !      2 (2)3
      ICEE(1:NMEL,1:3) = 0
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Computing the element-to-"        &
     &        //"element connectivity table: +"
      ENDIF
#endif      
      DO M = 1, NMEL
        NM_LOOP:DO J = 1, 3
          DO I1 = 1, ICNE(ELEMS(M,MAP(J)),1)
            M1 = ICNE(ELEMS(M,MAP(J)),I1+1)
            IF( M1 == M ) CYCLE
            DO I2 = 1, ICNE(ELEMS(M,MAP(J+1)),1)
              M2 = ICNE(ELEMS(M,MAP(J+1)),I2+1)
              IF( M1 == M2 ) THEN
                ICEE(M,J) = M1
                CYCLE NM_LOOP
              ENDIF
            ENDDO
          ENDDO
        ENDDO NM_LOOP
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(M,NMEL)
        ENDIF
#endif
      ENDDO

#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Determining which vertices are "  &
     &        //"on boundaries: +"
      ENDIF
#endif      
      NPROP(1:NMND)=0
      DO I = 1, NMEL
        N1 = ELEMS(I,1)
        N2 = ELEMS(I,2)
        N3 = ELEMS(I,3)
        NPROP(N1) = NPROP(N1) + N2 - N3
        NPROP(N2) = NPROP(N2) + N3 - N1
        NPROP(N3) = NPROP(N3) + N1 - N2
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(I,NMEL)
        ENDIF
#endif
      ENDDO

      ! Count total number of vertices on boundary.
      I = 0
      DO N = 1, NMND
        IF( NPROP(N) == 0 ) CYCLE
        I = I +1
      ENDDO
      ALLOCATE( IBSEG(1:I*2,1:2),IBSEGEL(1:I*2) )

      ! Compute the boundary segments, which occur in an element
      ! that has a zero in the element-to-element adjacency table.
      ! The bouundary segments are matched with adjacent elements.
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Determining the boundary "        &
     &        //"segments: +"
      ENDIF
#endif      
      ISEG = 0
      DO M = 1, NMEL
        DO J = 1, 3
          IF( ICEE(M,J) == 0 ) then
            ISEG = ISEG + 1
            IBSEG(ISEG,1) = ELEMS(M,MAP(J))
            IBSEG(ISEG,2) = ELEMS(M,MAP(J+1))
            IBSEGEL(ISEG) = M
          ENDIF
        ENDDO
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(M,NMEL)
        ENDIF
#endif
      ENDDO
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Connecting the boundary "         &
     &        //"segments: +"
      ENDIF
#endif      
      NSEG = ISEG
      ALLOCATE( ISEGP(1:NSEG),ISEGF(1:NSEG) )
      ISEGP(:) = 0
      ISEGF(:) = 0
      DO I = 1, NSEG
        DO J = I+1, NSEG
          IF( IBSEG(I,1) == IBSEG(J,2) ) then
             ISEGP(I) = J
             ISEGF(J) = I
          ENDIF
          IF( IBSEG(I,2) == IBSEG(J,1) ) then
             ISEGF(I) = J
             ISEGP(J) = I
          ENDIF
        ENDDO
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(I,NSEG)
        ENDIF
#endif
      ENDDO
 
      DEALLOCATE(NPROP)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Drog_MAK_NEINFO "       &
     &         //"routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_MAK_NEINFO



      ! The following routine contains computations that are performed
      ! after the input files have been read, but before the start
      ! of the time-stepping.

      SUBROUTINE LPT_Drog_Initialize

      USE LPT_Comm_Module
      USE LPT_Data_Module, ONLY: A => ShapeFunc_A,                     &
     &                           A0 => ShapeFunc_A0,                   &
     &                           AR => ElementArea_AR,                 &
     &                           B => ShapeFunc_B,                     &
     &                           DT1 => TimeStepSizes,                 &
     &                           ELEMS => MeshConnEV,                  &
     &                           ICURBS => PartBoundary,               &
     &                           JJDR => PartElemJ,                    &
     &                           KODE => PartBoundaryCode,             &
     &                           LLDR => PartElemL,                    &
     &                           MinimumTimeStep,                      &
     &                           NDR => NumParticlesLocal,             &
     &                           NMEL => NumElems,                     &
     &                           NMND => NumVerts,                     &
     &                           NOMP_MAP => MapOpenMP,                &
     &                           NumTrackingSnaps,                     &
     &                           PartDomainLocal,                      &
     &                           SimulationLength,                     &
     &                           StartingTime,                         &
     &                           T => ElementArea_T,                   &
     &                           T1 => Time1,                          &
     &                           T2 => Time2,                          &
     &                           X => MeshLon,                         &
     &                           XDR => PartLonLocal,                  &
     &                           Y => MeshLat,                         &
     &                           YDR => PartLatLocal,                  &
     &                           ZDR => PartDepthLocal
      USE LPT_Data_Lattice_Table
#ifdef KDTREE
      USE LPT_KDTREE_Module, ONLY: LPT_KDTREE_Search
#endif

      IMPLICIT NONE

      CHARACTER(LEN=100) :: JunkC
      CHARACTER(LEN=100) :: JunkC1

      INTEGER            :: I
      INTEGER            :: ICHECK
      INTEGER            :: III
      INTEGER            :: J
      INTEGER            :: N
      INTEGER            :: N1
      INTEGER            :: N2
      INTEGER            :: N3

      REAL(8)            :: XSTART
      REAL(8)            :: YSTART
      REAL(8)            :: ZSTART

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the LPT_Drog_Initialize "      &
     &         //"routine.")
#endif

      ! Convert other times to seconds.
      SimulationLength = SimulationLength * 3600.D0
      StartingTime = StartingTime
      MinimumTimeStep = MinimumTimeStep * 3600.D0

      ! Establish the limits of integration.
      T1 = StartingTime 
      T2 = StartingTime + SimulationLength/FLOAT(NumTrackingSnaps)

#ifdef MPI
      ! The dedicated reader/writer cores should not execute the code
      ! remaining in this routine.
      IF(UseReaderCore.AND.(MyRank.EQ.0))THEN
         RETURN
      ENDIF
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF
#endif

      ! Compute area coordinates for element interpolation functions.
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Preparing the element "           &
     &        //"interpolation functions: +"
      ENDIF
#endif      
      ALLOCATE( A(1:NMEL,1:3), B(1:NMEL,1:3), A0(1:NMEL,1:2) )
      ALLOCATE( AR(1:NMEL), T(1:NMEL,1:3) )
      DO J=1,NMEL
        N1 = ELEMS(J,1)
        N2 = ELEMS(J,2)
        N3 = ELEMS(J,3)
        A(J,1) = X(N3) - X(N2)
        A(J,2) = X(N1) - X(N3)
        A(J,3) = X(N2) - X(N1)
        B(J,1) = Y(N2) - Y(N3)
        B(J,2) = Y(N3) - Y(N1)
        B(J,3) = Y(N1) - Y(N2)
        A0(J,1) = 0.5* (X(N2)*Y(N3)-X(N3)*Y(N2))
        A0(J,2) = 0.5* (X(N3)*Y(N1)-X(N1)*Y(N3))
        AR(J) = 0.5* (A(J,2)*B(J,1)-A(J,1)*B(J,2))
        T(J,1) = A0(J,1)*2.E0
        T(J,2) = A0(J,2)*2.E0
        T(J,3) = (2.E0*AR(J)-T(J,1)-T(J,2))
#ifdef DEBUG
        IF(AR(J).LE.0)THEN
          WRITE(JunkC,'(I24)') J
          CALL LPT_Print(MyRank,"FATAL ERROR","Element number "        &
     &             //TRIM(ADJUSTL(JunkC))//" has a non-positive area.")
        END IF
#endif
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(J,NMEL)
        ENDIF
#endif
      ENDDO
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The element areas were computed.")
#endif

      ! Make the search table.
      CALL LPT_Drog_MAKE_STAB(NMEL,ELEMS,NMND,X,Y)

      ! Only execute the rest of this routine if an initial set
      ! of particles has also been read.  Otherwise, this information
      ! will be developed in the Oil_Source routine.
      IF(NDR.EQ.0)THEN
         RETURN
      ENDIF

      ! Find the starting element for each particle.
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Finding the starting element "    &
     &        //"for each of the particles: +"
      ENDIF
#endif
      IF(.NOT.ALLOCATED(JJDR)) ALLOCATE(JJDR(NDR))
      IF(.NOT.ALLOCATED(LLDR)) ALLOCATE(LLDR(NDR))
      IF(ALLOCATED(PartDomainLocal)) DEALLOCATE(PartDomainLocal)
      ALLOCATE(PartDomainLocal(1:NDR))
      PartDomainLocal(:) = 0
      DO III=1,NDR
        XSTART = XDR(III)
        YSTART = YDR(III)
        ZSTART = ZDR(III)
        CALL LPT_Drog_Find_Element(XSTART,YSTART,ZSTART,               &
     &           JJDR(III),LLDR(III),ICHECK)
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(III,NDR)
        ENDIF
#endif
        PartDomainLocal(III) = ICHECK
        IF(ICHECK.EQ.1) CYCLE
#ifdef DEBUG
        WRITE(JunkC,'(I24)') III
        WRITE(JunkC1,'(I24)') MyRank
        CALL LPT_Print(MyRank,"WARNING","Tracker core "                &
     &           //TRIM(ADJUSTL(JunkC1))//" could not find the "       &
     &           //"starting element for local particle number "       &
     &           //TRIM(ADJUSTL(JunkC))//".")
#endif
      ENDDO

      IF(.NOT.ALLOCATED(KODE)) ALLOCATE(KODE(NDR))
      IF(.NOT.ALLOCATED(ICURBS)) ALLOCATE(ICURBS(NDR))
      KODE(1:NDR)   = 0
      ICURBS(1:NDR) = 0

      ! Establish the default step size.
      IF(.NOT.ALLOCATED(DT1)) ALLOCATE(DT1(NDR))
      DO III=1,NDR
        DT1(III) = DABS(SimulationLength/FLOAT(NumTrackingSnaps)/10.d0)
      END DO 

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Drog_Initialize "       &
     &         //"routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      ! The following routine finds the horizontal and vertical elements
      ! containing a particle, given its coordinates.

      SUBROUTINE LPT_Drog_Find_Element(XSTART,YSTART,ZSTART,           &
                     JJ,LL,ICHECK)

      USE LPT_Data_Lattice_Table
      USE LPT_Data_Module, ONLY: NDIV => LatticeSearchBins

      IMPLICIT NONE

      INTEGER,INTENT(OUT)   :: ICHECK
      INTEGER,INTENT(INOUT) :: JJ
      INTEGER,INTENT(INOUT) :: LL

      REAL(8),INTENT(IN) :: XSTART
      REAL(8),INTENT(IN) :: YSTART
      REAL(8),INTENT(IN) :: ZSTART

      INTEGER :: I
      INTEGER :: J
      INTEGER :: IND
      INTEGER :: IX
      INTEGER :: IY

#ifdef KDTREE
      ICHECK = 0
      CALL LPT_KDTREE_Search(XSTART,YSTART,I,IND)
      IF(IND.EQ.1)THEN
        JJ = I
        ICHECK = 1
        CALL LPT_Drog_Find_Vertical_Element(I,JJ,XSTART,YSTART,ZSTART)
        LL = I
      ENDIF
#else
      IX = INT( (XSTART-XMIN(1))/DX(1) ) + 1
      IY = INT( (YSTART-XMIN(2))/DX(2) ) + 1
      IX = MAX(0,IX)
      IX = MIN(NDIV+1,IX)
      IY = MAX(0,IY)
      IY = MIN(NDIV+1,IY)
      ICHECK = 0
      DO J=1,NE_PIECE(IX,IY)
        I=NE_PIECE_LIST(J+NE_PIECE_INDEX(IX,IY))
        CALL LPT_Drog_BELEL(I,XSTART,YSTART,IND)
        IF(IND.EQ.1)THEN
          JJ = I
          ICHECK = 1
          CALL LPT_Drog_Find_Vertical_Element(I,JJ,XSTART,YSTART,ZSTART)
          LL = I
          EXIT
        END IF
      ENDDO
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      ! The following routine makes the lattice searching table.
      ! It is called by LPT_Drog_Initialize.

      SUBROUTINE LPT_Drog_MAKE_STAB(NMEL,ELEMS,NMND,X,Y)

      USE LPT_Comm_Module, ONLY: MyRank 
      USE LPT_Data_Lattice_Table

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: ELEMS(NMEL,3)
      INTEGER,INTENT(IN) :: NMEL
      INTEGER,INTENT(IN) :: NMND

      REAL(8),INTENT(IN) :: X(NMND)
      REAL(8),INTENT(IN) :: Y(NMND)

      INTEGER            :: I
      INTEGER            :: ISTART
      INTEGER            :: IX(3)
      INTEGER            :: IY(3)
      INTEGER            :: J
      INTEGER            :: M
      INTEGER            :: Multiplier
      INTEGER            :: N

      REAL(8)            :: XMAX(2)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the LPT_Drog_MAKE_STAB "       &
     &         //"routine.")
#endif

      XMAX(1) = MAXVAL(X(1:NMND))
      XMIN(1) = MINVAL(X(1:NMND))
      XMAX(2) = MAXVAL(Y(1:NMND))
      XMIN(2) = MINVAL(Y(1:NMND))
      DO I = 1, 2
        XMAX(I) = XMAX(I) + 1.D0
        DX(I) = ( XMAX(I) - XMIN(I) ) / DBLE(NDIV)
      ENDDO

      Multiplier = 0
 100  CONTINUE
      Multiplier = Multiplier + 2
      IF(ALLOCATED(NE_PIECE_LIST)) DEALLOCATE(NE_PIECE_LIST)
      ALLOCATE(NE_PIECE_LIST(Multiplier*NMEL))

      ! Search piece index.
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Determining the number of "       &
     &        //"mesh elements in each lattice piece: +"
      ENDIF
#endif
      IF(.NOT.ALLOCATED(NE_PIECE))THEN
         ALLOCATE(NE_PIECE(0:NDIV+1,0:NDIV+1))
      ENDIF
      NE_PIECE(0:NDIV+1,0:NDIV+1) = 0
      DO M = 1, NMEL
        DO J = 1, 3
          N = ELEMS(M,J)
          IX(J) = INT( (X(N)-XMIN(1)) / DX(1) ) + 1
          IY(J) = INT( (Y(N)-XMIN(2)) / DX(2) ) + 1
        ENDDO
        DO I = MINVAL(IX(1:3)), MAXVAL(IX(1:3))
          DO J = MINVAL(IY(1:3)), MAXVAL(IY(1:3))
            NE_PIECE(I,J) = NE_PIECE(I,J) + 1
          ENDDO
        ENDDO
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(M,NMEL)
        ENDIF
#endif
      ENDDO
      IF(.NOT.ALLOCATED(NE_PIECE_INDEX))THEN
         ALLOCATE(NE_PIECE_INDEX(1:NDIV,1:NDIV))
      ENDIF
      ISTART = 0
      DO I = 1, NDIV
        DO J = 1, NDIV
          NE_PIECE_INDEX(I,J) = ISTART
          ISTART = ISTART + NE_PIECE(I,J)
        ENDDO
      ENDDO

      ! Make piece table.
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
        WRITE(*,'(A,$)') "LPT: INFO: Placing the mesh elements in "    &
     &        //"each lattice piece: +"
      ENDIF
#endif      
      NE_PIECE(1:NDIV,1:NDIV) = 0
      DO M = 1, NMEL
        DO J = 1, 3
          N = ELEMS(M,J)
          IX(J) = INT( (X(N)-XMIN(1)) / DX(1) ) + 1
          IY(J) = INT( (Y(N)-XMIN(2)) / DX(2) ) + 1
        ENDDO
        DO I = MINVAL(IX(1:3)), MAXVAL(IX(1:3))
          DO J = MINVAL(IY(1:3)), MAXVAL(IY(1:3))
            ! Count total number of elements in NE_PIECE(I,J).
            NE_PIECE(I,J) = NE_PIECE(I,J) + 1
            IF((NE_PIECE(I,J)+NE_PIECE_INDEX(I,J)).GT.                 &
     &            SIZE(NE_PIECE_LIST))THEN
#if VERBOSE > 2
               IF(MyRank.EQ.0)THEN
                  WRITE(*,'(A)') "+"
               ENDIF
#endif
               GOTO 100
            ENDIF
            ! Store element number in NE_PIECE(I,J).
            NE_PIECE_LIST(NE_PIECE(I,J)+NE_PIECE_INDEX(I,J)) = M
          ENDDO
        ENDDO 
#if VERBOSE > 2
        IF(MyRank.EQ.0)THEN
          CALL LPT_Progress(M,NMEL)
        ENDIF
#endif
      ENDDO

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Drog_MAKE_STAB "        &
     &         //"routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_MAKE_STAB



      ! The following routine determines whether the point (XP,YP)
      ! lies within element J or on its boundaries.  If it does,
      ! then NFLAG=1; otherwise, NFLAG=0.

      SUBROUTINE LPT_Drog_BELEL(J,XP,YP,NFLAG)

      USE LPT_Data_Module, ONLY: ELEMS => MeshConnEV,                  &
     &                           X => MeshLon,                         &
     &                           Y => MeshLat

      IMPLICIT NONE

      INTEGER             :: I
      INTEGER,INTENT(IN)  :: J
      INTEGER             :: K
      INTEGER,INTENT(OUT) :: NFLAG

      REAL(8)             :: CROSSPROD
      REAL(8)             :: D
      REAL(8)             :: DELX
      REAL(8)             :: DELY
      REAL(8)             :: THETA
      REAL(8)             :: VX(3)
      REAL(8)             :: VY(3)
      REAL(8)             :: XLOCAL(3)
      REAL(8),INTENT(IN)  :: XP
      REAL(8)             :: YLOCAL(3)
      REAL(8),INTENT(IN)  :: YP

      ! Extract the local vertex coordinates.
      XLOCAL(1)=X(ELEMS(J,1))
      XLOCAL(2)=X(ELEMS(J,2))
      XLOCAL(3)=X(ELEMS(J,3))
      YLOCAL(1)=Y(ELEMS(J,1))
      YLOCAL(2)=Y(ELEMS(J,2))
      YLOCAL(3)=Y(ELEMS(J,3))

      ! Calculate the x- and y-components of vectors pointing from (XP,YP)
      ! to each vertex on the element.
      DO I=1,3
        DELX=XLOCAL(I)-XP
        DELY=YLOCAL(I)-YP
        D=DSQRT(DELX**2.0+DELY**2.0)
        THETA=DATAN2(DELY,DELX)
        VX(I)=D*DCOS(THETA)
        VY(I)=D*DSIN(THETA)
      ENDDO

      ! Determine if the point is on the element by calculating the
      ! cross-products of neighboring vectors in a direction which will
      ! yield all non-negative numbers if (X,Y) is on the element.
      NFLAG=1
      DO I=1,3
        K=I+1
        IF(I.EQ.3)K=1
        CROSSPROD=VX(I)*VY(K)-VY(I)*VX(K)
        IF(CROSSPROD.LT.0.D0)THEN
          NFLAG=0
          EXIT
        ENDIF
      ENDDO

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_BELEL



      ! The following routine determines the vertical element
      ! in which a particle is located.

      SUBROUTINE LPT_Drog_Find_Vertical_Element(Layer,J,X,Y,Z)

      USE LPT_Data_Module, ONLY: MeshDepth,                            &
     &                           NumLayers,                            &
     &                           Sigma

      IMPLICIT NONE

      ! Argument variables.

      INTEGER,INTENT(IN)    :: J
      INTEGER,INTENT(OUT)   :: Layer

      REAL(8),INTENT(IN)    :: X
      REAL(8),INTENT(IN)    :: Y
      REAL(8),INTENT(INOUT) :: Z

      ! Internal variables.

      INTEGER :: IL

      REAL(8) :: Depths(NumLayers)
      REAL(8) :: LPT_Drog_Interpolate

      ! If we are using a two-dimensional velocity field,
      ! then there can only be one layer.
      IF(NumLayers.EQ.1)THEN
         Layer = 1
         RETURN
      ENDIF

      ! Interpolate the depths for the layers
      ! at the horizontal location.
      Depths(1) = LPT_Drog_Interpolate(J,X,Y,MeshDepth)
      DO IL=2,NumLayers
         Depths(IL) = Depths(1) + (0.D0-Depths(1))                     &
     &              * (Sigma(IL)-Sigma(1))/(Sigma(NumLayers)-Sigma(1))
      ENDDO

      ! If the particle has moved below the sea floor,
      ! then bump it up and return.
      IF(Z.LT.Depths(1))THEN
         Z = Depths(1)
         Layer = 1
         RETURN
      ENDIF

      ! If the particle has moved above the sea surface,
      ! then bump it down and return.
      IF(Z.GT.Depths(NumLayers))THEN
         Z = Depths(NumLayers)
         Layer = NumLayers - 1
         RETURN
      ENDIF

      ! Otherwise, find the element within the water column.
      DO IL=1,NumLayers-1
         IF((Depths(IL).LE.Z).AND.(Z.LE.Depths(IL+1)))THEN
            Layer = IL
            EXIT
         ENDIF
      ENDDO

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_Find_Vertical_Element



      ! The following routine calls related functions to determine
      ! the components of velocity at a point (X,Y,Z).

      SUBROUTINE LPT_Drog_VELS(J,L,X,Y,Z,UF,VF,WF,DIAM)

      USE LPT_Data_Module, ONLY: BuoyancyMethod,                       &
     &                           MeshDepth,                            &
     &                           NumLayers,                            &
     &                           Sigma,                                &
     &                           UNEW => VelU,                         &
     &                           VNEW => VelV,                         &
     &                           WNEW => VelW

      IMPLICIT NONE
 
      INTEGER,INTENT(IN)     :: J
      INTEGER,INTENT(IN)     :: L

      REAL(8), INTENT(IN)    :: DIAM
      REAL(8), INTENT(INOUT) :: UF
      REAL(8), INTENT(INOUT) :: VF
      REAL(8), INTENT(INOUT) :: WF
      REAL(8), INTENT(IN)    :: X
      REAL(8), INTENT(IN)    :: Y
      REAL(8), INTENT(IN)    :: Z

      REAL(8) :: LPT_Drog_Interpolate

      INTEGER :: L1
      INTEGER :: L2

      REAL(8) :: Multiplier
      REAL(8) :: UL1
      REAL(8) :: UL2
      REAL(8) :: VelocityTerminal
      REAL(8) :: VL1
      REAL(8) :: VL2
      REAL(8) :: WL1
      REAL(8) :: WL2
      REAL(8) :: ZDP
      REAL(8) :: ZL1
      REAL(8) :: ZL2

      ! If the velocity field is two-dimensional,
      ! then return the depth-averaged velocities.
      IF(NumLayers.EQ.1)THEN

         UF = LPT_Drog_Interpolate(J,X,Y,UNEW(:,1))
         VF = LPT_Drog_Interpolate(J,X,Y,VNEW(:,1))
         WF = 0.D0

      ! Otherwise we need to interpolate a three-dimensional
      ! velocity field to a specific depth in the water column.
      ELSE

         IF(L.NE.NumLayers)THEN
            L1 = L
            L2 = L+1
         ELSE
            L1 = L-1
            L2 = L
         ENDIF

         ! Find the velocities at the layer below the particle.
         UL1 = LPT_Drog_Interpolate(J,X,Y,UNEW(:,L1))
         VL1 = LPT_Drog_Interpolate(J,X,Y,VNEW(:,L1))
         WL1 = LPT_Drog_Interpolate(J,X,Y,WNEW(:,L1))

         ! Find the velocities at the layer above the particle.
         UL2 = LPT_Drog_Interpolate(J,X,Y,UNEW(:,L2))
         VL2 = LPT_Drog_Interpolate(J,X,Y,VNEW(:,L2))
         WL2 = LPT_Drog_Interpolate(J,X,Y,WNEW(:,L2))

         ! Interpolate vertically to the depth of the particle.
         ZDP = LPT_Drog_Interpolate(J,X,Y,MeshDepth(:))
         ZL1 = ZDP + (0.D0 - ZDP) *                                    &
     &           (Sigma(L1)-Sigma(1))/(Sigma(NumLayers)-Sigma(1))
         ZL2 = ZDP + (0.D0 - ZDP) *                                    &
     &           (Sigma(L2)-Sigma(1))/(Sigma(NumLayers)-Sigma(1))
         Multiplier = (Z-ZL1)/(ZL2-ZL1)
         UF = UL1 + (UL2-UL1) * Multiplier
         VF = VL1 + (VL2-VL1) * Multiplier
         WF = WL1 + (WL2-WL1) * Multiplier

         IF(INDEX(BuoyancyMethod,"NULL").LE.0)THEN
            CALL LPT_Oil_Buoyancy(DIAM,VelocityTerminal)
            WF = WF + VelocityTerminal
         ENDIF

      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_VELS



      ! This function interpolates a 2D field (like velocity or depth)
      ! to the point (X,Y) within horizontal element J.

      REAL(8) FUNCTION LPT_Drog_Interpolate(J,X,Y,Field)

      USE LPT_Data_Module, ONLY: A => ShapeFunc_A,                     &
     &                           A0 => ShapeFunc_A0,                   &
     &                           AR => ElementArea_AR,                 &
     &                           B => ShapeFunc_B,                     &
     &                           ELEMS => MeshConnEV,                  &
     &                           NNMD => NumVerts

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: J

      REAL(8),INTENT(IN) :: Field(NNMD)
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(IN) :: Y

      INTEGER :: N1
      INTEGER :: N2
      INTEGER :: N3

      REAL(8) :: A03
      REAL(8) :: ARI
      REAL(8) :: F1
      REAL(8) :: F2
      REAL(8) :: F3
      REAL(8) :: V1
      REAL(8) :: V2
      REAL(8) :: V3

      LPT_Drog_Interpolate = 0.D0

      N1 = ELEMS(J,1)
      N2 = ELEMS(J,2)
      N3 = ELEMS(J,3)
      A03 = AR(J) - A0(J,1) - A0(J,2)
      ARI = 0.5d0/AR(J)
      F1 = Field(N1)
      F2 = Field(N2)
      F3 = Field(N3)
      V1 = ARI* (B(J,1)*F1+B(J,2)*F2+B(J,3)*F3)
      V2 = ARI* (A(J,1)*F1+A(J,2)*F2+A(J,3)*F3)
      V3 = 2*ARI* (A0(J,1)*F1+A0(J,2)*F2+A03*F3)

      LPT_Drog_Interpolate = LPT_Drog_Interpolate + (V1*X+V2*Y+V3)

      ! We're done here.
      RETURN

      END FUNCTION LPT_Drog_Interpolate



      ! This routine performs a 5th-order Runge-Kutta tracking of (X,Y)
      ! over the time interval (T1 to T2) using adaptive/variable
      ! tracking time step size.  For further details on the general
      ! methodology used in this algorithm, refer to Sections 15.1-15.2
      ! of 'Numerical Recipes' by Press et al.
      !
      ! X and Y are integrated from XSTART and YSTART, respectively,
      ! with a user-specified accuracy EPS.  DT1 is the guessed first
      ! step size and DTMIN is the minimum allowed step size (and can
      ! be zero.  Upon completion of the integrations, XSTART and YSTART
      ! are replaced with the final values for X and Y.
      !
      ! Note that T2 > T1 and T2 < T1 are allowed (i.e., forward and
      ! backward tracking).

      SUBROUTINE LPT_Drog_TRACK(JEL,J,LEL,L,XSTART,YSTART,ZSTART,      &
     &                    USTART,VSTART,WSTART,T1,T2,DT1,IPN,DIAM)

      USE LPT_Data_Module, ONLY: Cx,Cy,                                &
     &                           DiffusionMethod,                      &
     &                           DTMIN => MinimumTimeStep,             &
     &                           Evx,Evy,                              &
     &                           ICURBS => PartBoundary,               &
     &                           KODE => PartBoundaryCode,             &
     &                           NumLayers

      IMPLICIT NONE

      ! Argument variables.

      INTEGER           :: IPN
      INTEGER           :: J
      INTEGER           :: JEL
      INTEGER           :: L
      INTEGER           :: LEL

      REAL(8)           :: DIAM
      REAL(8)           :: DT1
      REAL(8)           :: T1
      REAL(8)           :: T2
      REAL(8)           :: USTART
      REAL(8)           :: VSTART
      REAL(8)           :: WSTART
      REAL(8)           :: XSTART
      REAL(8)           :: YSTART
      REAL(8)           :: ZSTART

      ! Internal parameters.

      INTEGER,PARAMETER :: MAXSTP = 1000

      REAL(8),PARAMETER :: TEST = 1.D-03
      REAL(8),PARAMETER :: ZERO = 0.D0

      ! Internal variables.

      INTEGER           :: ICOUNT
      INTEGER           :: ITEMP
      INTEGER           :: JTEMP
      INTEGER           :: NSTP

      REAL(8)           :: DT
      REAL(8)           :: DTNEXT
      REAL(8)           :: GAMMA
      REAL(8)           :: PDECAY
      REAL(8)           :: R
      REAL(8)           :: Rx
      REAL(8)           :: Ry
      REAL(8)           :: T
      REAL(8)           :: U0
      REAL(8)           :: V0
      REAL(8)           :: W0
      REAL(8)           :: XTEMP
      REAL(8)           :: XX
      REAL(8)           :: YTEMP
      REAL(8)           :: YY
      REAL(8)           :: ZZ

      XX = XSTART
      YY = YSTART
      ZZ = ZSTART
      T = T1
      DT = SIGN(DT1,T2-T1)
      J = JEL
      L = LEL
      ICOUNT = 0
      U0 = USTART
      V0 = VSTART
      W0 = WSTART

      DO NSTP=1,MAXSTP

!        IF((XX<-999998.D0).AND.(YY<-999998.D0))THEN
!           XSTART = XX
!           YSTART = YY
!           ZSTART = ZZ
!           DT1 = DTNEXT
!           RETURN
!        ENDIF

         IF (NSTP.NE.1) CALL LPT_Drog_VELS(J,L,XX,YY,ZZ,U0,V0,W0,DIAM)

         ! Check whether the step size will result in over-shooting
         ! the end of the integration interval T2.  If so, then cut down
         ! the step size.
         IF ((T+DT-T2)* (T+DT-T1).GT.ZERO) DT = T2 - T

         ! If the particle is in the interior of the domain,
         ! then perform a 5th-order Runge-Kutta step.
         IF(KODE(IPN).LE.0) THEN
            CALL LPT_Drog_RKQC(J,L,XX,YY,ZZ,U0,V0,W0,T,DT,DTNEXT,      &
     &               KODE(IPN),ICURBS(IPN),DIAM)

         ! If the particle is on the boundary, then track it there.
         ELSE
            CALL LPT_Drog_BOUNTRK(J,L,XX,YY,ZZ,U0,V0,W0,T,T2,          &
     &               KODE(IPN),ICURBS(IPN),DIAM)
!pck 10/04/13 - Gave DTNEXT a default value - NEEDS TO BE REVISITED
           DTNEXT=DTMIN
         ENDIF

         ! Add horizontal diffusion.
         IF(INDEX(DiffusionMethod,"PROCTOR1994").GT.0)THEN
            IF(Cx*Evx*DT.GT.10E-7 .OR. Cy*Evy*DT.GT.10E-7)THEN
               CALL RANDOM_NUMBER(Rx)
               CALL RANDOM_NUMBER(Ry)
               IF(Cx*Evx*DT.GT.10E-7)THEN
                  XTEMP = XX + (2*Rx-1)*(Cx*Evx*DT)**(1.d0/2.d0)
               ENDIF
               IF(Cy*Evy*DT.GT.10E-7)THEN
                  YTEMP = YY + (2*Ry-1)*(Cy*Evy*DT)**(1.d0/2.d0)
               ENDIF
               CALL LPT_Drog_FNDELE(J,JTEMP,XTEMP,YTEMP,ITEMP)
               IF(ITEMP.NE.0)THEN
                  XX = XTEMP
                  YY = YTEMP
               ENDIF
            ENDIF
         ENDIF

         ! If the integration is complete, then save the position
         ! and recent step size as an estimate for the next tracking.
         IF ((T-T2)* (T2-T1).GE.ZERO) THEN
            XSTART = XX
            YSTART = YY
            ZSTART = ZZ
            DT1 = DTNEXT
            RETURN
         END IF

         IF (ABS(DTNEXT).LT.DTMIN) THEN
           DTNEXT = SIGN(DTMIN,T2-T1)
         END IF

         ! Set the next step size.
         DT = DTNEXT

      ENDDO

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_TRACK



      ! This routine performs a 5th-order Runge-Kutta step with
      ! monitoring of local 4th-order truncation error to ensure
      ! accuracy and adjust step size.
      !
      ! DTTRY is the step size to be attempted, and EPS is the
      ! required accuracy.  Upon finishing the integration step,
      ! XX and YY contain the final tracked location.  DTNEXT is
      ! the estimated next step size.

      SUBROUTINE LPT_Drog_RKQC(JINOUT,LINOUT,XX,YY,ZZ,U,V,W,           &
     &                    T,DTTRY,DTNEXT,KODE,ICURBS,DIAM)

      USE LPT_Comm_Module, ONLY: MyRank
      USE LPT_Data_Module, ONLY: DTMIN => MinimumTimeStep,             &
     &                           EPS => RequiredAccuracy,              &
     &                           IBSEG => BoundarySegmentVerts,        &
     &                           IBSEGEL => BoundarySegmentElement,    &
     &                           NSEG => NumBoundarySegments,          &
     &                           NumLayers,                            &
     &                           X => MeshLon,                         &
     &                           Y => MeshLat

      IMPLICIT NONE

      ! Argument variables.

      INTEGER           :: ICURBS
      INTEGER           :: JINOUT
      INTEGER           :: KODE
      INTEGER           :: LINOUT

      REAL(8)           :: DIAM
      REAL(8)           :: DTNEXT
      REAL(8)           :: DTTRY
      REAL(8)           :: T
      REAL(8)           :: U
      REAL(8)           :: V
      REAL(8)           :: W
      REAL(8)           :: XX
      REAL(8)           :: YY
      REAL(8)           :: ZZ

      ! Internal parameters.

      REAL(8),PARAMETER :: ERRCON = 6.D-04
      REAL(8),PARAMETER :: FCOR   = 1.D0/15.D0
      REAL(8),PARAMETER :: ONE    = 1.D0
      REAL(8),PARAMETER :: PGROW  = -0.2D0
      REAL(8),PARAMETER :: PSHRNK = -0.25D0
      REAL(8),PARAMETER :: SAFETY = 0.9D0

      ! Internal variables.

      INTEGER           :: I
      INTEGER           :: IBNDSLOPEF
      INTEGER           :: ICHECK
      INTEGER           :: IELMIN
      INTEGER           :: IMIN
      INTEGER           :: IND
      INTEGER           :: ITRKSLOPEF
      INTEGER           :: J1HS
      INTEGER           :: J2HS
      INTEGER           :: JBS
      INTEGER           :: JSAV
      INTEGER           :: KODESAV
      INTEGER           :: L1HS
      INTEGER           :: L2HS
      INTEGER           :: LBS
      INTEGER           :: LSAV

      REAL(8)           :: BNDINT
      REAL(8)           :: BNDSLOPE
      REAL(8)           :: C1
      REAL(8)           :: C2
      REAL(8)           :: DIST
      REAL(8)           :: DISTMIN
      REAL(8)           :: DT
      REAL(8)           :: DTDID
      REAL(8)           :: DTH
      REAL(8)           :: ERROR
      REAL(8)           :: TH
      REAL(8)           :: TRKINT
      REAL(8)           :: TRKSLOPE
      REAL(8)           :: TSAV
      REAL(8)           :: USAV
      REAL(8)           :: V_X_VSEG
      REAL(8)           :: V_X_VTRK
      REAL(8)           :: VSAV
      REAL(8)           :: VTRK_X_VSEG
      REAL(8)           :: VX
      REAL(8)           :: VXSEG
      REAL(8)           :: VXTRK
      REAL(8)           :: VY
      REAL(8)           :: VYSEG
      REAL(8)           :: VYTRK
      REAL(8)           :: WSAV
      REAL(8)           :: XB
      REAL(8)           :: XBSMAX
      REAL(8)           :: XBSMIN
      REAL(8)           :: XE
      REAL(8)           :: XI
      REAL(8)           :: XIMIN
      REAL(8)           :: XSAV
      REAL(8)           :: XTEMP
      REAL(8)           :: XTRKMAX
      REAL(8)           :: XTRKMIN
      REAL(8)           :: XXHS
      REAL(8)           :: YB
      REAL(8)           :: YE
      REAL(8)           :: YI
      REAL(8)           :: YIMIN
      REAL(8)           :: YSAV
      REAL(8)           :: YTEMP
      REAL(8)           :: YTRKMAX
      REAL(8)           :: YTRKMIN
      REAL(8)           :: YYHS
      REAL(8)           :: ZSAV
      REAL(8)           :: ZTEMP
      REAL(8)           :: ZZHS

      ! Save the initial values.
      JSAV = JINOUT
      LSAV = LINOUT
      TSAV = T
      USAV = U
      VSAV = V
      WSAV = W
      XSAV = XX
      YSAV = YY
      ZSAV = ZZ
      KODESAV = KODE

      ! Set the step size to the initial value.
      DT = DTTRY

      INFINI_LOOP: DO

         ! Take the first half step.
         DTH = DT/2.d0
         CALL LPT_Drog_RK4(JSAV,J1HS,LSAV,L1HS,XSAV,YSAV,ZSAV,         &
     &            USAV,VSAV,WSAV,TSAV,DTH,XTEMP,YTEMP,ZTEMP,KODE,DIAM)

         ! If KODE.NE.0, then tracking has left the mesh during
         ! the first half of DT.
         !
         ! If KODE.EQ.1, then tracking has left the mesh on the first step.
         ! In this case, depending on KODESAV, either move to the boundary
         ! or just track along the boundary.
         !
         ! If KODESAV.EQ.-1, then the particle is already on the boundary,
         ! so just track it along the boundary.
         !
         ! If KODESAV.EQ.0, then the particle is coming from the interior
         ! of the domain, so move it to the boundary before tracking it there.
         !
         ! The following steps are used to move a particle to the boundary.
         !    1. Check to see which boundary segment(s) was (were) crossed.
         !    2. Determine intersection point(s).
         !    3. If more than one boundary segment was crossed, determine the closest.
         !    4. Put particle at intersection point with closest boundary.
         !    5. Figure out how long it took to get to the boundary and therefore 
         !       how much of DT is left.  This is done using USAV,VSAV.

         ! Extend the time step.
         IF((KODE.NE.0).AND.(KODESAV.EQ.-1))THEN
            RETURN
         ENDIF

         IMIN = 0
         ICHECK = 0

         IF((KODE.EQ.1).AND.(KODESAV.NE.-1))THEN

            DISTMIN=999999999.D0
            XTRKMIN=MIN(XSAV,XTEMP)
            XTRKMAX=MAX(XSAV,XTEMP)
            YTRKMIN=MIN(YSAV,YTEMP)
            YTRKMAX=MAX(YSAV,YTEMP)

            IF((XTEMP-XSAV).NE.0.D0)THEN
               TRKSLOPE=(YTEMP-YSAV)/(XTEMP-XSAV)
               TRKINT=YSAV-TRKSLOPE*XSAV
               ITRKSLOPEF=0
            ELSE
               ITRKSLOPEF=1
!pck 10/04/13 - Gave TRKSLOPE a value - NEEDS TO BE REVISITED
               TRKSLOPE=999999
            ENDIF

            DO I=1,NSEG

               XB=X(IBSEG(I,1))
               YB=Y(IBSEG(I,1))
               XE=X(IBSEG(I,2))
               YE=Y(IBSEG(I,2))

               ! Check if particle is near the boundary segment.
               IF((XB.LT.XTRKMIN).AND.(XE.LT.XTRKMIN)) CYCLE
               IF((XB.GT.XTRKMAX).AND.(XE.GT.XTRKMAX)) CYCLE
               IF((YB.LT.YTRKMIN).AND.(YE.LT.YTRKMIN)) CYCLE
               IF((YB.GT.YTRKMAX).AND.(YE.GT.YTRKMAX)) CYCLE

               ! Determine the slope and intercept of the boundary segment.
               IF((XE-XB).NE.0.D0)THEN
                  BNDSLOPE=(YE-YB)/(XE-XB)
                  BNDINT=YB-BNDSLOPE*XB
                  IBNDSLOPEF=0
               ELSE
                  IBNDSLOPEF=1
!pck 10/04/13 - Gave BNDSLOPE a value - NEEDS TO BE REVISITED
                  BNDSLOPE=999999
               ENDIF

               ! If the track and boundary segment are parallel.
               IF(TRKSLOPE.EQ.BNDSLOPE) CYCLE

               ! Compute the intersection point of track and boundary,
               ! and check if it falls on both the track and boundary.
               ! Keep track of minimum along track distance to the boundary.
               IF(TRKSLOPE.NE.BNDSLOPE)THEN
                  IF((ITRKSLOPEF.NE.1).AND.(IBNDSLOPEF.NE.1))THEN
                     XI=(BNDINT-TRKINT)/(TRKSLOPE-BNDSLOPE)
                     YI=BNDSLOPE*XI+BNDINT
                  ENDIF
                  IF(ITRKSLOPEF.EQ.1)THEN
                     XI=XSAV
                     YI=BNDSLOPE*XI+BNDINT
                  ENDIF
                  IF(IBNDSLOPEF.EQ.1)THEN
                     XI=XB
                     YI=TRKSLOPE*XI+TRKINT
                  ENDIF
                  XBSMIN=MIN(XB,XE)
                  XBSMAX=MAX(XB,XE)
                  IF((XI.GE.XBSMIN).AND.(XI.LE.XBSMAX).AND.            &
     &               (XI.GE.XTRKMIN).AND.(XI.LE.XTRKMAX))THEN
                     DIST=SQRT((XI-XSAV)*(XI-XSAV)+(YI-YSAV)*(YI-YSAV))
                     ICHECK = ICHECK+1
                     IF(DIST.LT.DISTMIN)THEN
                        XIMIN=XI
                        YIMIN=YI
                        DISTMIN=DIST
                        IELMIN=IBSEGEL(I)
                        IMIN=I
                     ENDIF
                  ENDIF
               ENDIF

            ENDDO

            ! The following is a fail-safe calculation in the event that
            ! the boundary segment is not found.  This procedure is the
            ! same as previous, except for the margin of path length.
            IF(IMIN.EQ.0)THEN
               VXTRK=XTEMP-XSAV
               VYTRK=YTEMP-YSAV
               DO I=1,NSEG
                  XB=X(IBSEG(I,1))
                  YB=Y(IBSEG(I,1))
                  XE=X(IBSEG(I,2))
                  YE=Y(IBSEG(I,2))
                  IF((XB.LT.XTRKMIN).AND.(XE.LT.XTRKMIN)) CYCLE
                  IF((XB.GT.XTRKMAX).AND.(XE.GT.XTRKMAX)) CYCLE
                  IF((YB.LT.YTRKMIN).AND.(YE.LT.YTRKMIN)) CYCLE
                  IF((YB.GT.YTRKMAX).AND.(YE.GT.YTRKMAX)) CYCLE
                  VX=XB-XSAV
                  VY=YB-YSAV
                  VXSEG=XE-XB
                  VYSEG=YE-YB
                  VTRK_X_VSEG=VXTRK*VYSEG-VXSEG*VYTRK
                  V_X_VTRK=VX*VYTRK-VXTRK*VY
                  V_X_VSEG=VX*VYSEG-VXSEG*VY
                  IF(DABS(VTRK_X_VSEG).LE.1.D-16) CYCLE
                  C1=V_X_VSEG/VTRK_X_VSEG
                  C2=V_X_VTRK/VTRK_X_VSEG
                  IF((C1.LE.1.1D0).AND.(C2.LE.1.D0))THEN
                     IF((C1.GE.-0.1D0).AND.(C2.GE.0.D0))THEN
                        XIMIN=XB+C2*VXSEG
                        YIMIN=YB+C2*VYSEG
                        DISTMIN=C1*DSQRT(VXTRK*VXTRK+VYTRK*VYTRK)
                        IELMIN=IBSEGEL(I)
                        IMIN=I
                        EXIT
                     ENDIF
                  ENDIF
               ENDDO
               IF(IMIN.EQ.0)THEN
                  IELMIN=0
                  CALL LPT_Print(MyRank,"INFO","Intersection point "   &
     &                     //"was NOT FOUND.")
                  !STOP
               ENDIF
            ENDIF

            ! Move the particle to the boundary.
            XX=XIMIN
            YY=YIMIN
            JINOUT=IELMIN
            ICURBS=IMIN

            DIST=SQRT((XTEMP-XSAV)**2.D0+(YTEMP-YSAV)**2.D0)
            ZZ=ZSAV+(DISTMIN/DIST)*(ZTEMP-ZSAV)
            CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)

            ! Compute the time required for the particle to reach
            ! the boundary, assuming velocity is stationary in time
            ! and space over this interval.
            DTDID=DISTMIN/SQRT(USAV*USAV+VSAV*VSAV+WSAV*WSAV)
            DTNEXT=DT
            T=TSAV+DTDID

            RETURN

         ENDIF

         ! If KODE.EQ.2, then tracking has left the mesh on the 
         ! second step.  In this case, cut down the tracking
         ! time step size.
         IF(KODE.EQ.2)THEN
            DT=DT/4.D0
            CYCLE INFINI_LOOP
         ENDIF

         ! If KODE.EQ.3 or KODE.EQ.4, then tracking has left the mesh
         ! on the third or fourth step.  Cut down the tracking
         ! time step size.
         IF((KODE.EQ.3).OR.(KODE.EQ.4))THEN
            DT=DT/2.D0
            CYCLE INFINI_LOOP
         ENDIF

         ! If the tracking remained in the mesh over the first half step,
         ! then proceed with the second half step.
         TH=TSAV+DTH
         CALL LPT_Drog_VELS(J1HS,L1HS,XTEMP,YTEMP,ZTEMP,U,V,W,DIAM)
         CALL LPT_Drog_RK4(J1HS,J2HS,L1HS,L2HS,XTEMP,YTEMP,ZTEMP,U,V,W,&
     &            TH,DTH,XXHS,YYHS,ZZHS,KODE,DIAM)

         ! If KODE.NE.0, then tracking has left the mesh during
         ! the second half step.  Cut down the tracking step size.
         IF(KODE.NE.0)THEN
            DT=DT/2.D0
            CYCLE INFINI_LOOP
         ENDIF

         ! If the tracking remained in the mesh over both half steps,
         ! then proceed with the large step.
         CALL LPT_Drog_RK4(JSAV,JBS,LSAV,LBS,XSAV,YSAV,ZSAV,           &
     &            USAV,VSAV,WSAV,TSAV,DT,XTEMP,YTEMP,ZTEMP,KODE,DIAM)

         ! If KODE.NE.0, then tracking has left the mesh during
         ! the large step.  Cut down the tracking size.
         IF(KODE.NE.0)THEN
            DT=DT/2.D0
            CYCLE INFINI_LOOP
         ENDIF

         ERROR=DSQRT((ZZHS-ZTEMP)**2+(YYHS-YTEMP)**2+(XXHS-XTEMP)**2)

         ! Scale the error relative to the required accuracy.
         ERROR=ERROR/EPS

         IF((ERROR.GT.ONE).AND.(ABS(DT).GE.DTMIN))THEN

            ! If the truncation error is too large, then reduce
            ! the step size and try again.
            DT=SAFETY*DT*(ERROR**PSHRNK)
            CYCLE INFINI_LOOP

         ELSE

            ! Otherwise, the step succeeded, so estimate the size
            ! of the next step.
            T=TSAV+DT

            IF(ERROR.GT.ERRCON)THEN
               DTNEXT=SAFETY*DT*(ERROR**PGROW)
            ELSE
               DTNEXT=4.D0*DT
            ENDIF
            IF(DT.LT.DTMIN) DT=DTMIN
            IF(DTNEXT.LT.DTMIN) DTNEXT=DTMIN

         ENDIF

         EXIT

      ENDDO INFINI_LOOP

      XX=XXHS
      YY=YYHS
      ZZ=ZZHS

      ! Add the truncation error to increase accuracy from 4th-order
      ! to 5th-order.
      XX=XX+(XX-XTEMP)*FCOR
      YY=YY+(YY-YTEMP)*FCOR
      ZZ=ZZ+(ZZ-ZTEMP)*FCOR

      ! Make sure that the point (XX,YY) is still in the element J2HS.
      ! If not, then don't increase to 5th-order.
      CALL LPT_Drog_BELEL(J2HS,XX,YY,IND)
      IF(IND.EQ.0)THEN
        XX=XXHS
        YY=YYHS
        ZZ=ZZHS
        CALL LPT_Drog_BELEL(J2HS,XX,YY,IND)
      ENDIF
      CALL LPT_Drog_Find_Vertical_Element(L2HS,J2HS,XX,YY,ZZ)

      JINOUT=J2HS
      LINOUT=L2HS

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_RKQC



      ! The following routine uses a 4th-order Runge-Kutta scheme
      ! to advance the solution over a time interval DT and returns
      ! the resulting point (XOUT,YOUT) as the coordinate of the
      ! endpoint of the integration step.

      SUBROUTINE LPT_Drog_RK4(J,JOUT,L,LOUT,XX,YY,ZZ,U,V,W,            &
     &                        T,DT,XOUT,YOUT,ZOUT,KODE,DIAM)

      IMPLICIT NONE

      ! Argument variables.

      INTEGER :: J
      INTEGER :: JOUT
      INTEGER :: KODE
      INTEGER :: L
      INTEGER :: LOUT

      REAL(8) :: DIAM
      REAL(8) :: DT
      REAL(8) :: T
      REAL(8) :: U
      REAL(8) :: V
      REAL(8) :: W
      REAL(8) :: XOUT
      REAL(8) :: XX
      REAL(8) :: YOUT
      REAL(8) :: YY
      REAL(8) :: ZOUT
      REAL(8) :: ZZ

      ! Internal variables.

      INTEGER :: IND
      INTEGER :: INDD
      INTEGER :: JJ
      INTEGER :: JT
      INTEGER :: LT

      REAL(8) :: DT6
      REAL(8) :: DTH
      REAL(8) :: TH
      REAL(8) :: U4
      REAL(8) :: UM
      REAL(8) :: UT
      REAL(8) :: V4
      REAL(8) :: VM
      REAL(8) :: VT
      REAL(8) :: W4
      REAL(8) :: WM
      REAL(8) :: WT
      REAL(8) :: XT
      REAL(8) :: YT
      REAL(8) :: ZT

      KODE=0
      DTH=DT/2.D0
      DT6=DT/6.D0
      TH=T+DTH
      JT=J
      LT=L

      ! First step.
      XT=XX+DTH*U
      YT=YY+DTH*V
      ZT=ZZ+DTH*W
      ! Check whether the point still lies in element J.
      ! If not, then find a new element.
      CALL LPT_Drog_BELEL(JT,XT,YT,IND)
      IF(IND.EQ.0)THEN
         CALL LPT_Drog_FNDELE(JT,JJ,XT,YT,INDD)
         CALL LPT_Drog_FNDELEls(INDD,JJ,XT,YT)
         IF(INDD.EQ.0)THEN
            KODE=1
            XOUT=XT
            YOUT=YT
            ZOUT=ZT
            RETURN
         ENDIF
         JT= J
      ENDIF
      CALL LPT_Drog_Find_Vertical_Element(LT,JT,XT,YT,ZT)

      ! Second step.
      CALL LPT_Drog_VELS(JT,LT,XT,YT,ZT,UT,VT,WT,DIAM)
      XT=XX+DTH*UT
      YT=YY+DTH*VT
      ZT=ZZ+DTH*WT
      CALL LPT_Drog_BELEL(JT,XT,YT,IND)
      IF(IND.EQ.0)THEN
         CALL LPT_Drog_FNDELE(JT,JJ,XT,YT,INDD)
         CALL LPT_Drog_FNDELEls(INDD,JJ,XT,YT)
         IF(INDD.EQ.0)THEN
            KODE=2
            RETURN
         ENDIF
         JT=JJ
      ENDIF
      CALL LPT_Drog_Find_Vertical_Element(LT,JT,XT,YT,ZT)

      ! Third step.
      CALL LPT_Drog_VELS(JT,LT,XT,YT,ZT,UM,VM,WM,DIAM)
      XT=XX+DT*UM
      YT=YY+DT*VM
      ZT=ZZ+DT*WM
      CALL LPT_Drog_BELEL(JT,XT,YT,IND)
      IF(IND.EQ.0)THEN
         CALL LPT_Drog_FNDELE(JT,JJ,XT,YT,INDD)
         CALL LPT_Drog_FNDELEls(INDD,JJ,XT,YT)
         IF(INDD.EQ.0)THEN
            KODE=3
            RETURN
         ENDIF
         JT=JJ
      ENDIF
      CALL LPT_Drog_Find_Vertical_Element(LT,JT,XT,YT,ZT)

      ! Fourth step.
      CALL LPT_Drog_VELS(JT,LT,XT,YT,ZT,U4,V4,W4,DIAM)
      ! Accumulate increments with the proper weights.
      XOUT=XX+(U+2*(UT+UM)+U4)*DT6
      YOUT=YY+(V+2*(VT+VM)+V4)*DT6
      ZOUT=ZZ+(W+2*(WT+WM)+W4)*DT6
      CALL LPT_Drog_BELEL(JT,XOUT,YOUT,IND)
      IF(IND.EQ.0)THEN
         CALL LPT_Drog_FNDELE(JT,JJ,XOUT,YOUT,INDD)
         CALL LPT_Drog_FNDELEls(INDD,JJ,XOUT,YOUT)
         IF(INDD.EQ.0)THEN
            KODE=4
            RETURN
         ENDIF
         JT=JJ
      ENDIF
      CALL LPT_Drog_Find_Vertical_Element(LT,JT,XOUT,YOUT,ZOUT)

      JOUT=JT
      LOUT=LT

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_RK4



      ! The following routine finds new element JJ which contains
      ! the point (XX,YY).  The old element is J.  If the new element
      ! is found, then INDD=1; otherwise, INDD=0.

      SUBROUTINE LPT_Drog_FNDELE(J,JJ,XX,YY,INDD)

      USE LPT_Data_Module, ONLY: ELEMS => MeshConnEV,                  &
     &                           ICEE => MeshConnEE,                   &
     &                           ICNE => MeshConnVE,                   &
     &                           NMEL => NumElems

      IMPLICIT NONE

      ! Argument variables.

      INTEGER,INTENT(OUT)   :: INDD
      INTEGER,INTENT(IN)    :: J
      INTEGER,INTENT(INOUT) :: JJ

      REAL(8),INTENT(IN)    :: XX
      REAL(8),INTENT(IN)    :: YY

      ! Internal variables.

      INTEGER :: I
      INTEGER :: ICOUNT
      INTEGER :: IDONE
      INTEGER :: IND
      INTEGER :: ISUM
      INTEGER :: K
      INTEGER :: KK
      INTEGER :: L1
      INTEGER :: L2
      INTEGER :: N
      INTEGER :: NDONE(1000)
      INTEGER :: NEC(1000)
      INTEGER :: NECN(1000)
      INTEGER :: NL

      INDD = 0

      IDONE = 1
      NDONE(IDONE) = J

      ! Check neighboring elements.
      DO I=1,3
         N=ICEE(J,I)
         IF(N.EQ.0) CYCLE
         CALL LPT_Drog_BELEL(N,XX,YY,IND)
         IF(IND.EQ.1)THEN
            JJ=N
            INDD=1
            RETURN
         ENDIF
         IDONE=IDONE+1
         NDONE(IDONE)=N
      ENDDO
      ISUM=0
      DO L1=1,3
         L2=ELEMS(J,L1)
         KK_LOOP: DO KK=1,ICNE(L2,1)
            N=ICNE(L2,KK+1)
            DO I=1,IDONE
               IF(N.EQ.NDONE(I)) CYCLE KK_LOOP
            ENDDO
            ISUM=ISUM+1
            IDONE=IDONE+1
            NDONE(IDONE)=N
            NEC(ISUM)=N
         ENDDO KK_LOOP
      ENDDO

      ! Progressively search over a larger area for the element
      ! containing (XX,YY).
      ICOUNT=ISUM
      ISUM=0
      IF(ICOUNT.EQ.0)THEN
         DO I=2,IDONE
            ICOUNT=ICOUNT+1
            NEC(ICOUNT)=NDONE(I)
         ENDDO
      ENDIF

      INFINI_LOOP: DO
         DO K=1,ICOUNT
            N=NEC(K)
            CALL LPT_Drog_BELEL(N,XX,YY,IND)
            IF(IND.EQ.1)THEN
               JJ=N
               INDD=1
               RETURN
            ENDIF
         ENDDO
         IF(IDONE.GT.70)THEN
            INDD=0
            RETURN
         END IF
         ! Accumulate list of new elements to be checked on next loop.
         DO K=1,ICOUNT
            N=NEC(K)
            DO L1=1,3
               L2=ELEMS(N,L1)
               KK_LOOP2: DO KK=1,ICNE(L2,1)
                  NL=ICNE(L2,KK+1)
                  DO I=1,IDONE
                     IF(NL.EQ.NDONE(I)) CYCLE KK_LOOP2
                  ENDDO
                  ISUM=ISUM+1
                  IDONE=IDONE+1
                  NDONE(IDONE)=NL
                  NECN(ISUM)=NL
               ENDDO KK_LOOP2
            ENDDO
         ENDDO

         ICOUNT=ISUM
         IF(ICOUNT.EQ.0)THEN
            INDD=0
            RETURN
         ENDIF

         NEC(1:ICOUNT)=NECN(1:ICOUNT)

         ISUM=0

      ENDDO INFINI_LOOP

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_FNDELE



      ! This routine finds the new element JJ that contains
      ! the point (XX,YY), by using a lattice search.

      SUBROUTINE LPT_Drog_FNDELEls(INDD,JJ,XT,YT)

      USE LPT_Data_LATTICE_TABLE

      IMPLICIT NONE

      ! Argument variables.

      INTEGER,INTENT(INOUT) :: INDD
      INTEGER,INTENT(INOUT) :: JJ

      REAL(8),INTENT(IN)    :: XT
      REAL(8),INTENT(IN)    :: YT

      ! Internal variables.

      INTEGER :: I
      INTEGER :: IND
      INTEGER :: IX
      INTEGER :: IY
      INTEGER :: J

      ! If the element was found in LPT_Drog_FNDELE (which is
      ! always called directly before this routine), then
      ! make a fast return.
      IF(INDD.EQ.1) RETURN

      IX=INT((XT-XMIN(1))/DX(1))+1
      IY=INT((YT-XMIN(2))/DX(2))+1
      IX=MAX(0,IX)
      IX=MIN(NDIV+1,IX)
      IY=MAX(0,IY)
      IY=MIN(NDIV+1,IY)
      DO J=1,NE_PIECE(IX,IY)
         I=NE_PIECE_LIST(J+NE_PIECE_INDEX(IX,IY))
         CALL LPT_Drog_BELEL(I,XT,YT,IND)
         IF(IND.EQ.1)THEN
            JJ=I
            INDD=1
            EXIT
         ENDIF
      ENDDO

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_FNDELEls



      ! The following routine tracks a particle along a boundary segment.
      ! The time step size to be attempted is T2-T.  Upon finishing
      ! the integration step, XX and YY contain the final tracked location.

      ! This tracking is done by:
      !  1. Computing the tangential velocity at the particle position.
      !  2. Determining which boundary node is "downstream."
      !  3. Using average of velocities computed in 1 and 2, determine whether 
      !     particle would make it to the node.
      !  4. If the answer in 3 is no, move particle along boundary until DT
      !     is used up.  Remember that particle is on boundary and return.
      !  5. If the answer in 3 is yes, put particle at downstream boundary 
      !     node and then figure how much of DT is left.  Then return and 
      !     see if particle will leave boundary or continue to move along the
      !     next boundary segment.

      SUBROUTINE LPT_Drog_BOUNTRK(JINOUT,LINOUT,XX,YY,ZZ,U,V,W,        &
     &                            T,T2,KODE,ICURBS,DIAM)

      USE LPT_Data_Module, ONLY: IBSEG => BoundarySegmentVerts,        &
     &                           IBSEGEL => BoundarySegmentElement,    &
     &                           ISEGF => BoundarySegmentFollowing,    &
     &                           ISEGP => BoundarySegmentPreceding,    &
     &                           X => MeshLon,                         &
     &                           Y => MeshLat

      IMPLICIT NONE

      ! Argument variables.

      INTEGER :: ICURBS
      INTEGER :: JINOUT
      INTEGER :: KODE
      INTEGER :: LINOUT

      REAL(8) :: DIAM
      REAL(8) :: T
      REAL(8) :: T2
      REAL(8) :: U
      REAL(8) :: V
      REAL(8) :: W
      REAL(8) :: XX
      REAL(8) :: YY
      REAL(8) :: ZZ

      ! Internal variables.

      INTEGER :: JSAV
      INTEGER :: LSAV

      REAL(8) :: COSTHETA
      REAL(8) :: DELTRK
      REAL(8) :: DELXB
      REAL(8) :: DELXE
      REAL(8) :: DELYB
      REAL(8) :: DELYE
      REAL(8) :: DIST
      REAL(8) :: DISTSEG
      REAL(8) :: DTDID
      REAL(8) :: DTTRY
      REAL(8) :: SINTHETA
      REAL(8) :: SPDAVG
      REAL(8) :: SPDTANXM
      REAL(8) :: SPDTANXX
      REAL(8) :: TSAV
      REAL(8) :: UM
      REAL(8) :: USAV
      REAL(8) :: UTANAVG
      REAL(8) :: UTANXX
      REAL(8) :: UTANXM
      REAL(8) :: VM
      REAL(8) :: VSAV
      REAL(8) :: VTANAVG
      REAL(8) :: VTANYY
      REAL(8) :: VTANYM
      REAL(8) :: WM
      REAL(8) :: WSAV
      REAL(8) :: XB
      REAL(8) :: XE
      REAL(8) :: XM
      REAL(8) :: XSAV
      REAL(8) :: YB
      REAL(8) :: YE
      REAL(8) :: YM
      REAL(8) :: YSAV
      REAL(8) :: ZM
      REAL(8) :: ZSAV

      ! Save the initial values.
      JSAV=JINOUT
      LSAV=LINOUT
      TSAV=T
      USAV=U
      VSAV=V
      WSAV=W
      XSAV=XX
      YSAV=YY
      ZSAV=ZZ
      DTTRY=T2-T

      ! Set up the boundary segment.
      XB=X(IBSEG(ICURBS,1))
      YB=Y(IBSEG(ICURBS,1))
      XE=X(IBSEG(ICURBS,2))
      YE=Y(IBSEG(ICURBS,2))
      DELXE=XE-XX
      DELYE=YE-YY
      DELXB=XB-XX
      DELYB=YB-YY

      ! Compute the velocity tangential to boundary at XX,YY and time T.
      DISTSEG=SQRT((XE-XB)*(XE-XB)+(YE-YB)*(YE-YB))
      COSTHETA=(XE-XB)/DISTSEG
      SINTHETA=(YE-YB)/DISTSEG
      UTANXX=USAV*COSTHETA*COSTHETA+VSAV*SINTHETA*COSTHETA
      VTANYY=USAV*COSTHETA*SINTHETA+VSAV*SINTHETA*SINTHETA
      SPDTANXX=SQRT(UTANXX*UTANXX+VTANYY*VTANYY)

      ! If this velocity is zero, then the particle does not move.
      IF((UTANXX.EQ.0.D0).AND.(VTANYY.EQ.0.D0))THEN
         T=T2
         ZZ=DTTRY*WSAV+ZSAV
         CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)
         RETURN
      ENDIF

      ! Determine which end of the boundary that the particle is moving toward. 
      XM=0.D0
      YM=0.D0

      ! If at one end of boundary segment (XE,YE), first check to see if 
      ! particle is moving toward the other end (XB,YB).  IF not, check to 
      ! see if particle is moving to other end of adjacent boundary segment.
      ! If it is moving in neither direction, it is temporarily stuck in 
      ! present position.
      IF((DELXE.EQ.0.D0).AND.(DELYE.EQ.0.D0))THEN
         IF((DELXB*UTANXX.GE.0.D0).AND.(DELYB*VTANYY.GE.0.D0))THEN
            XM=XB
            YM=YB
            GOTO 10
         ENDIF
         XE=X(IBSEG(ISEGF(ICURBS),2))
         YE=Y(IBSEG(ISEGF(ICURBS),2))
         DELXE=XE-XX
         DELYE=YE-YY
         DISTSEG=SQRT(DELXE*DELXE+DELYE*DELYE)
         COSTHETA=DELXE/DISTSEG
         SINTHETA=DELYE/DISTSEG
         UTANXX=USAV*COSTHETA*COSTHETA+VSAV*SINTHETA*COSTHETA
         VTANYY=USAV*COSTHETA*SINTHETA+VSAV*SINTHETA*SINTHETA
         SPDTANXX=SQRT(UTANXX*UTANXX+VTANYY*VTANYY)
         IF((UTANXX.EQ.0.D0).AND.(VTANYY.EQ.0.D0))THEN
            T=T2
            ZZ=DTTRY*WSAV+ZSAV
            CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)
            RETURN
         ENDIF
         IF((DELXE*UTANXX.GE.0.D0).AND.(DELYE*VTANYY.GE.0.D0))THEN
            XM=XE
            YM=YE
            ICURBS=ISEGF(ICURBS)
            JINOUT=IBSEGEL(ICURBS)
            GOTO 10
         ENDIF
         T=T2
         ZZ=DTTRY*WSAV+ZSAV
         CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)
         RETURN
      ENDIF

      ! If at one end of boundary segment (XB,YB), first check to see if 
      ! particle is moving toward the other end (XE,YE).  IF not, check to 
      ! see if particle is moving to other end of adjacent boundary segment.
      ! If it is moving in neither direction, it is temporarily stuck in 
      ! present position.
      IF((DELXB.EQ.0.D0).AND.(DELYB.EQ.0.D0))THEN
         IF((DELXE*UTANXX.GE.0.D0).AND.(DELYE*VTANYY.GE.0.D0))THEN
            XM=XE
            YM=YE
            GOTO 10
         ENDIF
         XB=X(IBSEG(ISEGP(ICURBS),1))
         YB=Y(IBSEG(ISEGP(ICURBS),1))
         DELXB=XB-XX
         DELYB=YB-YY
         DISTSEG=SQRT(DELXB*DELXB+DELYB*DELYB)
         COSTHETA=DELXB/DISTSEG
         SINTHETA=DELYB/DISTSEG
         UTANXX=USAV*COSTHETA*COSTHETA+VSAV*SINTHETA*COSTHETA
         VTANYY=USAV*COSTHETA*SINTHETA+VSAV*SINTHETA*SINTHETA
         SPDTANXX=SQRT(UTANXX*UTANXX+VTANYY*VTANYY)
         IF((UTANXX.EQ.0.D0).AND.(VTANYY.EQ.0.D0))THEN
            T=T2
            ZZ=DTTRY*WSAV+ZSAV
            CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)
            RETURN
         ENDIF
         IF((DELXB*UTANXX.GE.0.D0).AND.(DELYB*VTANYY.GE.0.D0))THEN
            XM=XB
            YM=YB
            ICURBS=ISEGP(ICURBS)
            JINOUT=IBSEGEL(ICURBS)
            GOTO 10
         ENDIF
         T=T2
         ZZ=DTTRY*WSAV+ZSAV
         CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)
         RETURN
      ENDIF

      ! If particle is at neither end of the boundary segment, check to
      ! see if it is moving toward XE,YE.
      IF((DELXE*UTANXX.GE.0.D0).AND.(DELYE*VTANYY.GE.0.D0))THEN
         XM=XE
         YM=YE
         GOTO 10
      ENDIF

      ! If particle is at neither end of the boundary segment, check to
      ! see if it is moving toward XB,YB.
      IF((DELXB*UTANXX.GE.0.D0).AND.(DELYB*VTANYY.GE.0.D0))THEN
         XM=XB
         YM=YB
         GOTO 10
      ENDIF

  10  CONTINUE

      ! If XM,YM have not been set, could not figure which direction particle
      ! is moving.
      IF((XM.EQ.0.D0).AND.(YM.EQ.0.D0))THEN
         T=T2
         ZZ=DTTRY*WSAV+ZSAV
         CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)
         RETURN
      ENDIF

      ! Determine the velocity at the end of the boundary segment where
      ! the particle is moving toward.
      CALL LPT_Drog_VELS(JINOUT,LINOUT,XM,YM,ZM,UM,VM,WM,DIAM)
      UTANXM=UM*COSTHETA*COSTHETA+VM*SINTHETA*COSTHETA
      VTANYM=UM*COSTHETA*SINTHETA+VM*SINTHETA*SINTHETA
      SPDTANXM=SQRT(UTANXM*UTANXM+VTANYM*VTANYM)

      ! Assuming the velocity is stationary in time but varying in space,
      ! compute the track length.
      DIST=SQRT((XM-XX)*(XM-XX)+(YM-YY)*(YM-YY))
      SPDAVG=0.5D0*(SPDTANXX+SPDTANXM)

      ! If the particle overshoots the end of the boundary segment,
      ! then put it at the end.
      IF(DIST.LE.SPDAVG*DTTRY)THEN
         XX=XM
         YY=YM
         DTDID=DIST/SPDAVG
         T=TSAV+DTDID
         KODE=-1
         ZZ=DTDID*WSAV+ZSAV
         CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)

      ! If the particle didn't overshoot the end of the boundary segment,
      ! then move it the appropriate distance along the segment.
      ELSE
         DELTRK=2.D0*SPDTANXX*DTTRY*DIST/(2.D0*DIST-DTTRY*             &
     &              (SPDTANXM-SPDTANXX))
         UTANAVG=(DELTRK*(UTANXM-UTANXX)/DIST+2.D0*UTANXX)/2.D0
         VTANAVG=(DELTRK*(VTANYM-VTANYY)/DIST+2.D0*VTANYY)/2.D0
         XX=DTTRY*UTANAVG+XSAV
         YY=DTTRY*VTANAVG+YSAV
         T=TSAV+DTTRY
         KODE=-1
         ZZ=DTTRY*WSAV   +ZSAV
         CALL LPT_Drog_Find_Vertical_Element(LINOUT,JINOUT,XX,YY,ZZ)
      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE LPT_Drog_BOUNTRK

