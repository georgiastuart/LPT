#ifdef NETCDF
      MODULE LPT_NetCDF_Module

      USE netcdf

      IMPLICIT NONE

#if NETCDF == 4
      INTEGER :: FillValueI = 0
#endif
      INTEGER :: NC_Counter
      INTEGER :: NC_DIM_Part
      INTEGER :: NC_DIM_Time
      INTEGER :: NC_ID
      INTEGER :: NC_ID_Curr
      INTEGER :: NC_ID_Output
      INTEGER :: NC_ID_Particle
      INTEGER :: NC_ID_Wind
      INTEGER :: NC_Snap
      INTEGER :: NC_SnapCurr
      INTEGER :: NC_SnapWind
      INTEGER :: NC_Status
      INTEGER :: NC_VAR
#if NETCDF == 4
      INTEGER :: NC_VAR_NumParticles
#endif
      INTEGER :: NC_VAR_Time
      INTEGER :: NC_VAR_D
      INTEGER :: NC_VAR_X
      INTEGER :: NC_VAR_Y
      INTEGER :: NC_VAR_Z

      REAL(8)             :: FillValueR = -99999.D0
      REAL(8),ALLOCATABLE :: NC_Temp(:)



      CONTAINS



      SUBROUTINE LPT_NetCDF_Check(NC_Status)

      USE LPT_Comm_Module, ONLY: MyRank
      USE netcdf

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: NC_Status

      IF(NC_Status.NE.NF90_NOERR)THEN
         CALL LPT_Print(MyRank,"FATAL ERROR",                          &
     &            TRIM(NF90_STRERROR(NC_Status)))
      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE



      END MODULE
#endif
