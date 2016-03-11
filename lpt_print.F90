      SUBROUTINE LPT_Print(MessageCore,MessageType,MessageBody)

      USE LPT_Comm_Module, ONLY: MyRank

      IMPLICIT NONE

      CHARACTER(*),INTENT(IN) :: MessageBody
      CHARACTER(*),INTENT(IN) :: MessageType

      INTEGER,INTENT(IN) :: MessageCore

      INTEGER :: IERR

      IF((MessageCore.EQ.-1).OR.(MessageCore.EQ.MyRank))THEN
         WRITE(*,'(A)') "LPT: "//TRIM(ADJUSTL(MessageType))//": "      &
     &         //TRIM(ADJUSTL(MessageBody))
         CALL FLUSH(6)
      ENDIF

      IF(INDEX(MessageType,"FATAL").GT.0)THEN
#ifdef MPI
         CALL MPI_FINALIZE(IERR)
#endif
         STOP
      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Progress(Now,Total)

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: Now
      INTEGER,INTENT(IN) :: Total

      INTEGER            :: C
      INTEGER            :: N

      outer: DO N=1,20
         C = MIN(CEILING(N*0.05*REAL(Total)),Total)
         IF(Now.EQ.C)THEN
            IF(MOD(N,20).EQ.0)THEN
               WRITE(*,'(A)')   "+"
            ELSEIF(MOD(N,5).EQ.0)THEN
               WRITE(*,'(A,$)') "+"
            ELSE
               WRITE(*,'(A,$)') "-"
            ENDIF
            EXIT outer
         ENDIF
      ENDDO outer

      END SUBROUTINE



