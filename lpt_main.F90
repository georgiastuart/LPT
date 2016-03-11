      PROGRAM LPT_Main

      USE LPT_Comm_Module
      USE LPT_Data_Module, ONLY: Deg2Rad,PI
      USE LPT_Read_Module
      USE LPT_Write_Module

      IMPLICIT NONE

      Pi = 2.D0 * DASIN(1.D0)
      Deg2Rad = Pi / 180.D0

      ! Initialize the parallel environment.
      CALL LPT_Comm_Init

      ! Read the input files.
      CALL LPT_Read

      ! Do some initial computations.
      CALL LPT_Drog_Initialize

      ! More like time-trippin', amirite?
      CALL LPT_Drog_TimeStep

      ! Finalize the parallel environment
      CALL LPT_Comm_Final

      ! The program ends in LPT_Comm_Final.
      ! No other instructions should be added here.

      END PROGRAM

