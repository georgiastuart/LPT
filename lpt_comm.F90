      MODULE LPT_Comm_Module

      IMPLICIT NONE

#ifdef MPI
      INCLUDE 'mpif.h'
#endif

      INTEGER             :: COMM_LPT_ALL
      INTEGER             :: MyRank
      INTEGER             :: NumRanks
      INTEGER             :: NumTrackerCores
      INTEGER,ALLOCATABLE :: TrackerRanks(:)

      LOGICAL             :: AmTrackerCore
      LOGICAL             :: UseReaderCore
      LOGICAL             :: UseWriterCore

#ifdef MPI
      CHARACTER(LEN=1)     :: DummyC = "0"
      INTEGER,DIMENSION(1) :: DummyI = (/0/)
      REAL(8),DIMENSION(1) :: DummyR = (/0.D0/)
      CHARACTER(LEN=100)   :: MessageC
      INTEGER,ALLOCATABLE  :: MessageI(:)
      REAL(8),ALLOCATABLE  :: MessageR(:)
#endif

      CONTAINS


      SUBROUTINE LPT_Comm_Init

      IMPLICIT NONE

      CHARACTER(LEN=100)  :: JunkC

      INTEGER             :: GROUP_LPT_ALL
      INTEGER             :: IA
      INTEGER             :: IARGC
      INTEGER             :: IERR
      INTEGER             :: IT

      NumRanks = 1
      MyRank = 0

#ifdef MPI
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_DUP(MPI_COMM_WORLD,COMM_LPT_ALL,IERR)
      CALL MPI_COMM_SIZE(COMM_LPT_ALL,NumRanks,IERR)
      CALL MPI_COMM_RANK(COMM_LPT_ALL,MyRank,IERR)
#endif
#if VERBOSE > 0
      WRITE(JunkC,'(I24)') NumRanks
      CALL LPT_Print(0,"INFO","There are a total of "                  &
     &               //TRIM(ADJUSTL(JunkC))//" core(s) available for " &
     &               //"reading, tracking and writing.")
#endif

      UseReaderCore = .FALSE.
      UseWriterCore = .FALSE.
      IF(IARGC().GT.0)THEN
         IA = 0
         DO WHILE (IA.LT.IARGC())
            IA = IA + 1
            CALL GETARG(IA,JunkC)
            IF(JunkC(1:2).EQ."-R")THEN
               UseReaderCore = .TRUE.
            ELSEIF(JunkC(1:2).EQ."-W")THEN
               UseWriterCore = .TRUE.
            ENDIF
         ENDDO
      ENDIF

      IF((NumRanks.EQ.1).AND.UseReaderCore)THEN
#if VERBOSE > 0
         CALL LPT_Print(MyRank,"WARNING","The dedicated reader "       &
     &            //"core cannot be used in a serial run.")
#endif
         UseReaderCore = .FALSE.
      ENDIF

      IF((NumRanks.EQ.1).AND.UseWriterCore)THEN
#if VERBOSE > 0
         CALL LPT_Print(MyRank,"WARNING","The dedicated writer "       &
     &            //"core cannot be used in a serial run.")
#endif
         UseWriterCore = .FALSE.
      ENDIF

      NumTrackerCores = NumRanks
      IF(UseReaderCore)THEN
         NumTrackerCores = NumTrackerCores - 1
      ENDIF
      IF(UseWriterCore)THEN
         NumTrackerCores = NumTrackerCores - 1
      ENDIF
#ifdef DEBUG
      IF(NumTrackerCores.LE.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","There must be at least one "  &
     &            //"core available for tracking.")
      ENDIF
#endif

#if VERBOSE > 0
      IF(UseReaderCore)THEN
         CALL LPT_Print(0,"INFO","The first core will read the input " &
     &            //"files and broadcast to the tracking cores.")
      ENDIF
      IF(UseWriterCore)THEN
         CALL LPT_Print(0,"INFO","The last core will be reserved for " &
     &            //"file output.")
      ENDIF
      WRITE(JunkC,'(I24)') NumTrackerCores
      CALL LPT_Print(0,"INFO","The particles will be tracked by "      &
     &         //TRIM(ADJUSTL(JunkC))//" core(s).")
#endif

#ifdef MPI
      IF(.NOT.ALLOCATED(TrackerRanks))THEN
         ALLOCATE(TrackerRanks(NumTrackerCores))
      ENDIF
      IF(.NOT.UseReaderCore)THEN
         TrackerRanks(1) = 0
      ELSE
         TrackerRanks(1) = 1
      ENDIF
      DO IT=2,NumTrackerCores
         TrackerRanks(IT) = TrackerRanks(IT-1) + 1
      ENDDO
      AmTrackerCore = .FALSE.
      DO IT=1,NumTrackerCores
         IF(MyRank.EQ.TrackerRanks(IT))THEN
            AmTrackerCore = .TRUE.
         ENDIF
      ENDDO
#else
      IF(.NOT.ALLOCATED(TrackerRanks)) ALLOCATE(TrackerRanks(1))
      TrackerRanks(1) = 0
      AmTrackerCore = .TRUE.
#endif

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Comm_Init routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Comm_Final

      IMPLICIT NONE

#ifdef MPI
      INTEGER :: IERR
#endif

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the LPT_Comm_Final routine.")
#endif

#ifdef MPI
      CALL MPI_BARRIER(COMM_LPT_ALL,IERR)
      CALL MPI_FINALIZE(IERR)
#if VERBOSE > 0
      CALL LPT_Print(0,"INFO","The parallel environment was "          &
     &               //"finalized successfully.")
#endif
#endif

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Comm_Final routine.")
#endif

      ! We're done everywhere.
      STOP

      END SUBROUTINE



#ifdef MPI
      SUBROUTINE LPT_Comm_To_Writer(VariableChar,VariableInteger,      &
     &               VariableReal,VariableType,ErrorMessage)

      ! This routine passes a variable from the rank-zero core
      ! to the writer core.  The variable can be of any type,
      ! as controlled by the 'VariableType' input string.
      ! The variable can also be of any length, which is detected
      ! automatically within this routine.

      IMPLICIT NONE

      CHARACTER(*),INTENT(IN)    :: ErrorMessage
      CHARACTER(LEN=100)         :: JunkC
      CHARACTER(LEN=100)         :: Routine
      CHARACTER(*),INTENT(INOUT) :: VariableChar
      CHARACTER(*),INTENT(IN)    :: VariableType

      INTEGER                    :: IC
      INTEGER                    :: IERR
      INTEGER                    :: IREQ
      INTEGER                    :: ISTAT(MPI_STATUS_SIZE)
      INTEGER                    :: ITAG = 0
      INTEGER,INTENT(INOUT)      :: VariableInteger(:)

      REAL(8),INTENT(INOUT)      :: VariableReal(:)

      ! The regular tracker cores should not continue within this routine.
      IF((MyRank.NE.0).AND.AmTrackerCore)THEN
         RETURN
      ENDIF

#ifdef DEBUG_MPI
      WRITE(Routine,'(A)') "LPT_Comm_To_Writer"
#endif

      IF(MyRank.EQ.0)THEN

         IF(INDEX(VariableType,"C").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1a"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_ISEND(VariableChar,LEN(VariableChar),             &
     &               MPI_CHARACTER,NumRanks-1,ITAG,                    &
     &               COMM_LPT_ALL,IREQ,IERR)
         ELSEIF(INDEX(VariableType,"I").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1b"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_ISEND(VariableInteger,SIZE(VariableInteger),      &
     &               MPI_INTEGER,NumRanks-1,ITAG,                      &
     &               COMM_LPT_ALL,IREQ,IERR)
         ELSEIF(INDEX(VariableType,"R").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1c"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_ISEND(VariableReal,SIZE(VariableReal),            &
     &               MPI_REAL8,NumRanks-1,ITAG,                        &
     &               COMM_LPT_ALL,IREQ,IERR)
         ENDIF
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif
#ifdef DEBUG_MPI
         WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 2"
         CALL FLUSH(10000+MyRank)
#endif
         CALL MPI_WAIT(IREQ,ISTAT,IERR)
#ifdef DEBUG
         IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif

      ELSEIF(MyRank+1.EQ.NumRanks)THEN

         IF(INDEX(VariableType,"C").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1a"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableChar,LEN(VariableChar),             &
     &               MPI_CHARACTER,0,MPI_ANY_TAG,                      &
     &               COMM_LPT_ALL,IREQ,IERR)
         ELSEIF(INDEX(VariableType,"I").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1b"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableInteger,SIZE(VariableInteger),      &
     &               MPI_INTEGER,0,MPI_ANY_TAG,                        &
     &               COMM_LPT_ALL,IREQ,IERR)
         ELSEIF(INDEX(VariableType,"R").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1c"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableReal,SIZE(VariableReal),            &
     &               MPI_REAL8,0,MPI_ANY_TAG,                          &
     &               COMM_LPT_ALL,IREQ,IERR)
         ENDIF
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif
#ifdef DEBUG_MPI
         WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 2"
         CALL FLUSH(10000+MyRank)
#endif
         CALL MPI_WAIT(IREQ,ISTAT,IERR)
#ifdef DEBUG
         IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif

      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE
#endif



#ifdef MPI
      SUBROUTINE LPT_Comm_From_Writer(VariableCharacter,               &
     &               VariableInteger,VariableReal,VariableType,        &
     &               ErrorMessage)

      IMPLICIT NONE

      CHARACTER(*),INTENT(IN)    :: ErrorMessage
      CHARACTER(LEN=100)         :: Routine
      CHARACTER(*),INTENT(INOUT) :: VariableCharacter
      CHARACTER(*),INTENT(IN)    :: VariableType

      INTEGER                    :: IC
      INTEGER                    :: IERR
      INTEGER                    :: IREQ(NumRanks)
      INTEGER                    :: ISTAT(MPI_STATUS_SIZE)
      INTEGER                    :: ITAG = 0
      INTEGER,INTENT(INOUT)      :: VariableInteger(:)

      REAL(8),INTENT(INOUT)      :: VariableReal(:)

      ! The dedicated reader core should not continue within this routine.
      IF((MyRank.EQ.0).AND.UseReaderCore)THEN
         RETURN
      ENDIF

#ifdef DEBUG_MPI
      WRITE(Routine,'(A)') "LPT_Comm_From_Writer"
#endif

      IF(MyRank+1.EQ.NumRanks)THEN

         DO IC=1,NumTrackerCores
            IF(TrackerRanks(IC).EQ.MyRank) CYCLE
            IF(INDEX(VariableType,"C").GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1a"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_ISEND(VariableCharacter,LEN(VariableCharacter),&
     &                  MPI_CHARACTER,TrackerRanks(IC),ITAG,           &
     &                  COMM_LPT_ALL,IREQ(TrackerRanks(IC)+1),IERR)
            ELSEIF(INDEX(VariableType,"I").GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1b"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_ISEND(VariableInteger,SIZE(VariableInteger),   &
     &                  MPI_INTEGER,TrackerRanks(IC),ITAG,             &
     &                  COMM_LPT_ALL,IREQ(TrackerRanks(IC)+1),IERR)
            ELSEIF(INDEX(VariableType,"R").GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1c"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_ISEND(VariableReal,SIZE(VariableReal),         &
     &                  MPI_REAL8,TrackerRanks(IC),ITAG,               &
     &                  COMM_LPT_ALL,IREQ(TrackerRanks(IC)+1),IERR)
            ENDIF
#ifdef DEBUG
            IF(IERR.NE.0)THEN
               CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
            ENDIF
#endif
         ENDDO

         DO IC=1,NumTrackerCores
            IF(TrackerRanks(IC).EQ.MyRank) CYCLE
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 2"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_WAIT(IREQ(TrackerRanks(IC)+1),ISTAT,IERR)
#ifdef DEBUG
            IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
               CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
            ENDIF
#endif
         ENDDO

      ELSEIF(AmTrackerCore)THEN

         IF(INDEX(VariableType,"C").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1a"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableCharacter,LEN(VariableCharacter),   &
     &               MPI_CHARACTER,NumRanks-1,MPI_ANY_TAG,COMM_LPT_ALL,&
     &               IREQ(MyRank+1),IERR)
         ELSEIF(INDEX(VariableType,"I").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1b"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableInteger,SIZE(VariableInteger),      &
     &               MPI_INTEGER,NumRanks-1,MPI_ANY_TAG,COMM_LPT_ALL,  &
     &               IREQ(MyRank+1),IERR)
         ELSEIF(INDEX(VariableType,"R").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1c"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableReal,SIZE(VariableReal),            &
     &               MPI_REAL8,NumRanks-1,MPI_ANY_TAG,COMM_LPT_ALL,    &
     &               IREQ(MyRank+1),IERR)
         ENDIF
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif
#ifdef DEBUG_MPI
         WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 2"
         CALL FLUSH(10000+MyRank)
#endif
         CALL MPI_WAIT(IREQ(MyRank+1),ISTAT,IERR)
#ifdef DEBUG
         IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif

      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE
#endif


#ifdef MPI
      SUBROUTINE LPT_Comm_Distribute(VariableChar,VariableInteger,     &
     &               VariableReal,VariableType,ErrorMessage)

      ! This routine distributes a variable from the rank-zero core
      ! to the tracker cores.  The variable can be of any type,
      ! as controlled by the 'VariableType' input string.
      ! The variable can also be of any length, which is detected
      ! automatically within this routine.

      IMPLICIT NONE

      CHARACTER(*),INTENT(IN)    :: ErrorMessage
      CHARACTER(LEN=100)         :: Routine
      CHARACTER(*),INTENT(INOUT) :: VariableChar
      CHARACTER(*),INTENT(IN)    :: VariableType

      INTEGER,SAVE               :: Counter = 0
      INTEGER                    :: IC
      INTEGER                    :: IERR
      INTEGER                    :: IREQ(NumRanks)
      INTEGER                    :: ISTAT(MPI_STATUS_SIZE)
      INTEGER                    :: ITAG = 0
      INTEGER,INTENT(INOUT)      :: VariableInteger(:)

      REAL(8),INTENT(INOUT)      :: VariableReal(:)

      ! The dedicated writer core should not proceed within this routine.
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF

#ifdef DEBUG_MPI
      WRITE(Routine,'(A)') "LPT_Comm_Distribute"
#endif

      IF(MyRank.EQ.0)THEN

         DO IC=1,NumTrackerCores
            IF(TrackerRanks(IC).EQ.MyRank) CYCLE
            IF(INDEX(VariableType,"C").GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1a"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_ISEND(VariableChar,LEN(VariableChar),          &
     &                  MPI_CHARACTER,TrackerRanks(IC),ITAG,           &
     &                  COMM_LPT_ALL,IREQ(TrackerRanks(IC)+1),IERR)
            ELSEIF(INDEX(VariableType,"I").GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1b"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_ISEND(VariableInteger,SIZE(VariableInteger),   &
     &                  MPI_INTEGER,TrackerRanks(IC),ITAG,             &
     &                  COMM_LPT_ALL,IREQ(TrackerRanks(IC)+1),IERR)
            ELSEIF(INDEX(VariableType,"R").GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1c"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_ISEND(VariableReal,SIZE(VariableReal),         &
     &                  MPI_REAL8,TrackerRanks(IC),ITAG,               &
     &                  COMM_LPT_ALL,IREQ(TrackerRanks(IC)+1),IERR)
            ENDIF
#ifdef DEBUG
            IF(IERR.NE.0)THEN
               CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
            ENDIF
#endif
         ENDDO

         DO IC=1,NumTrackerCores
            IF(TrackerRanks(IC).EQ.MyRank) CYCLE
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 2"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_WAIT(IREQ(TrackerRanks(IC)+1),ISTAT,IERR)
#ifdef DEBUG
            IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
               CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
            ENDIF
#endif
         ENDDO

      ELSEIF(AmTrackerCore)THEN

         IF(INDEX(VariableType,"C").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1a"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableChar,LEN(VariableChar),             &
     &               MPI_CHARACTER,0,MPI_ANY_TAG,COMM_LPT_ALL,         &
     &               IREQ(MyRank+1),IERR)
         ELSEIF(INDEX(VariableType,"I").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1b"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableInteger,SIZE(VariableInteger),      &
     &               MPI_INTEGER,0,MPI_ANY_TAG,COMM_LPT_ALL,           &
     &               IREQ(MyRank+1),IERR)
         ELSEIF(INDEX(VariableType,"R").GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1c"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(VariableReal,SIZE(VariableReal),            &
     &               MPI_REAL8,0,MPI_ANY_TAG,COMM_LPT_ALL,             &
     &               IREQ(MyRank+1),IERR)
         ENDIF
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif
#ifdef DEBUG_MPI
         WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 2"
         CALL FLUSH(10000+MyRank)
#endif
         CALL MPI_WAIT(IREQ(MyRank+1),ISTAT,IERR)
#ifdef DEBUG
         IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif

      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE
#endif



#ifdef MPI
      SUBROUTINE LPT_Comm_Collect(VariableChar,VariableInteger,        &
     &               VariableReal,VariableType,ArraySize,ErrorMessage)

      IMPLICIT NONE

      CHARACTER(*),INTENT(IN)    :: ErrorMessage
      CHARACTER(LEN=100)         :: Routine
      CHARACTER(*),INTENT(INOUT) :: VariableChar
      CHARACTER(*),INTENT(IN)    :: VariableType

      INTEGER                    :: ArrayEnd
      INTEGER,INTENT(IN)         :: ArraySize
      INTEGER                    :: ArrayStart
      INTEGER                    :: IC
      INTEGER                    :: IERR
      INTEGER                    :: IREQ(NumRanks)
      INTEGER                    :: ISTAT(MPI_STATUS_SIZE)
      INTEGER                    :: ITAG = 0
      INTEGER                    :: ParcelSize(NumTrackerCores)
      INTEGER,INTENT(INOUT)      :: VariableInteger(:)

      REAL(8),INTENT(INOUT)      :: VariableReal(:)

#ifdef DEBUG_MPI
      WRITE(Routine,'(A)') "LPT_Comm_Collect"
#endif

      ! Unlike in the distribute routine, the sizes of the arrays
      ! in this routine are not the same on each core.
      ! Rather, the tracker cores send their local information,
      ! but the writer core must build a global array.
      ! So the arrays on the writer core are always larger.

      IF(AmTrackerCore.AND.(MyRank+1.NE.NumRanks))THEN

         ! First send the size of the parcel.
#ifdef DEBUG_MPI
         WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1"
         CALL FLUSH(10000+MyRank)
#endif
         CALL MPI_ISEND(ArraySize,1,MPI_INTEGER,NumRanks-1,ITAG,       &
     &            COMM_LPT_ALL,IREQ(MyRank+1),IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif
#ifdef DEBUG_MPI
         WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 2"
         CALL FLUSH(10000+MyRank)
#endif
         CALL MPI_WAIT(IREQ(MyRank+1),ISTAT,IERR)
#ifdef DEBUG
         IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif

         ! Next send the parcel.
         IF(INDEX(VariableType,"I").GT.0)THEN
            IF(ArraySize.GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 3a"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_ISEND(VariableInteger,ArraySize,MPI_INTEGER,   &
     &                  NumRanks-1,ITAG,COMM_LPT_ALL,IREQ(MyRank+1),   &
     &                  IERR)
            ENDIF
         ELSEIF(INDEX(VariableType,"R").GT.0)THEN
            IF(ArraySize.GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 3b"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_ISEND(VariableReal,ArraySize,MPI_REAL8,        &
     &                  NumRanks-1,ITAG,COMM_LPT_ALL,IREQ(MyRank+1),   &
     &                  IERR)
            ENDIF
         ENDIF
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif
         IF(ArraySize.GT.0)THEN
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 4"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_WAIT(IREQ(MyRank+1),ISTAT,IERR)
         ENDIF
#ifdef DEBUG
         IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif

      ELSEIF(MyRank+1.EQ.NumRanks)THEN

         ! First receive the sizes of the parcels.
         ParcelSize(:) = 0
         DO IC=1,NumTrackerCores
            IF(TrackerRanks(IC).EQ.MyRank) CYCLE
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 1"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_IRECV(ParcelSize(IC),1,                           &
     &               MPI_INTEGER,TrackerRanks(IC),                     &
     &               MPI_ANY_TAG,COMM_LPT_ALL,                         &
     &               IREQ(TrackerRanks(IC)+1),IERR)
#ifdef DEBUG
            IF(IERR.NE.0)THEN
               CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
            ENDIF
#endif
         ENDDO
         DO IC=1,NumTrackerCores
            IF(TrackerRanks(IC).EQ.MyRank) CYCLE
#ifdef DEBUG_MPI
            WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 2"
            CALL FLUSH(10000+MyRank)
#endif
            CALL MPI_WAIT(IREQ(TrackerRanks(IC)+1),ISTAT,IERR)
#ifdef DEBUG
            IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
               CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
            ENDIF
#endif
         ENDDO

         ! Next receive the parcels.
         DO IC=1,NumTrackerCores
            IF(TrackerRanks(IC).EQ.MyRank) CYCLE
            IF(IC-1.EQ.0)THEN
               ArrayStart = 1
            ELSE
               ArrayStart = SUM(ParcelSize(1:IC-1))+1
            ENDIF
            ArrayEnd = SUM(ParcelSize(1:IC))
            IF(INDEX(VariableType,"I").GT.0)THEN
               IF(ParcelSize(IC).GT.0)THEN
#ifdef DEBUG_MPI
                  WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))     &
     &                     //" 3a"
                  CALL FLUSH(10000+MyRank)
#endif
                  CALL MPI_IRECV(VariableInteger(ArrayStart:ArrayEnd), &
     &                     ParcelSize(IC),MPI_INTEGER,                 &
     &                     TrackerRanks(IC),MPI_ANY_TAG,               &
     &                     COMM_LPT_ALL,IREQ(TrackerRanks(IC)+1),IERR)
               ENDIF
            ELSEIF(INDEX(VariableType,"R").GT.0)THEN
               IF(ParcelSize(IC).GT.0)THEN
#ifdef DEBUG_MPI
                  WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))     &
     &                     //" 3b"
                  CALL FLUSH(10000+MyRank)
#endif
                  CALL MPI_IRECV(VariableReal(ArrayStart:ArrayEnd),    &
     &                     ParcelSize(IC),MPI_REAL8,                   &
     &                     TrackerRanks(IC),MPI_ANY_TAG,               &
     &                     COMM_LPT_ALL,IREQ(TrackerRanks(IC)+1),IERR)
               ENDIF
            ENDIF
#ifdef DEBUG
            IF(IERR.NE.0)THEN
               CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
            ENDIF
#endif
         ENDDO
         DO IC=1,NumTrackerCores
            IF(TrackerRanks(IC).EQ.MyRank) CYCLE
            IF(ParcelSize(IC).GT.0)THEN
#ifdef DEBUG_MPI
               WRITE(10000+MyRank,'(A)') TRIM(ADJUSTL(Routine))//" 4"
               CALL FLUSH(10000+MyRank)
#endif
               CALL MPI_WAIT(IREQ(TrackerRanks(IC)+1),ISTAT,IERR)
            ENDIF
#ifdef DEBUG
            IF(ISTAT(MPI_ERROR).NE.MPI_SUCCESS.OR.IERR.NE.0)THEN
               CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
            ENDIF
#endif
         ENDDO

      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE
#endif



      END MODULE
