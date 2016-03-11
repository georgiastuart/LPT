      MODULE LPT_Write_Module

      IMPLICIT NONE

      CONTAINS



      SUBROUTINE LPT_Write_Particle_Snap

      USE LPT_Comm_Module
      USE LPT_Data_Module
#ifdef NETCDF
      USE LPT_NetCDF_Module
#endif

      IMPLICIT NONE

      CHARACTER(LEN=100)  :: JunkC
      CHARACTER(LEN=100)  :: JunkC1
      CHARACTER(LEN=100)  :: JunkC2
      CHARACTER(LEN=100)  :: JunkC3
      CHARACTER(LEN=100)  :: JunkC4
      CHARACTER(LEN=100)  :: JunkC5

      INTEGER             :: ArrayEnd
      INTEGER             :: ArrayStart
      INTEGER             :: IA
      INTEGER             :: IP
      INTEGER,ALLOCATABLE :: PartDomainGlobal(:)
      INTEGER,ALLOCATABLE :: PartNumberGlobal(:)
      INTEGER             :: ParticlesOnLocalCore(NumTrackerCores)

      LOGICAL,SAVE        :: FirstCall = .TRUE.

      REAL(8)             :: Conversion
      REAL(8),ALLOCATABLE :: PartDepthGlobal(:)
      REAL(8),ALLOCATABLE :: PartDiameterGlobal(:)
      REAL(8),ALLOCATABLE :: PartLatGlobal(:)
      REAL(8),ALLOCATABLE :: PartLonGlobal(:)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the LPT_Write_Particle_Snap "  &
     &         //"routine.")
#endif

#ifdef MPI
      ! The dedicated reader core should not execute the code
      ! remaining in this routine.
      IF(UseReaderCore.AND.(MyRank.EQ.0))THEN
         RETURN
      ENDIF
#endif

      ! Determine the global number of particles.
#ifdef MPI
      IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
      ! The tracker cores send the local number of particles.
      IF(MyRank+1.NE.NumRanks)THEN
         IF(.NOT.ALLOCATED(MessageI))THEN
            ALLOCATE(MessageI(1))
            MessageI(1) = NumParticlesLocal
         ENDIF
      ! The writer core must be ready to receive the numbers
      ! of particles from all of the tracker cores.
      ELSEIF(MyRank+1.EQ.NumRanks)THEN
         IF(.NOT.ALLOCATED(MessageI))THEN
            ALLOCATE(MessageI(NumTrackerCores))
            MessageI(:) = 0
            IF(.NOT.UseWriterCore)THEN
               MessageI(NumTrackerCores) = NumParticlesLocal
            ENDIF
         ENDIF
      ENDIF
      CALL LPT_Comm_Collect(DummyC,MessageI,DummyR,"I",1,              &
     &         "The local numbers of particles were not received by "  &
     &         //"the writer core.")
#endif
      IF(MyRank+1.EQ.NumRanks)THEN
#ifdef MPI
         ParticlesOnLocalCore(:) = MessageI(:)
         NumParticlesGlobal = SUM(ParticlesOnLocalCore)
#else
         NumParticlesGlobal = NumParticlesLocal
#endif
#if VERBOSE > 1
         WRITE(JunkC,'(I24)') NumParticlesGlobal
         CALL LPT_Print(0,"INFO","There will be "                      &
     &            //TRIM(ADJUSTL(JunkC))//" particles in the next "    &
     &            //"snap of output.")
#endif
      ENDIF

      ! Send the global number of particles back to the tracker cores.
      ! This is necessary for the case when the particles are added
      ! from a source, and we want to keep track of the global numbering
      ! as particles are added on a local core.
#ifdef MPI
      IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
      ALLOCATE(MessageI(1:1))
      MessageI(1) = NumParticlesGlobal
      CALL LPT_Comm_From_Writer(DummyC,MessageI,DummyR,"I",            &
     &         "The global number of particles were not sent to the "  &
     &         //"tracker cores.")
      NumParticlesGlobal = MessageI(1)
#endif

      ! Develop the global arrays of particle locations.
      DO IA=1,2
#ifdef MPI
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         IF(MyRank+1.NE.NumRanks)THEN
            IF(NumParticlesLocal.GT.0)THEN
               ALLOCATE(MessageI(1:NumParticlesLocal))
               IF(IA.EQ.1)THEN
                  MessageI(:) = PartNumberLocal(1:NumParticlesLocal)
               ELSEIF(IA.EQ.2)THEN
                  MessageI(:) = PartDomainLocal(1:NumParticlesLocal)
               ENDIF
            ENDIF
         ELSEIF(MyRank+1.EQ.NumRanks)THEN
            ALLOCATE(MessageI(1:NumParticlesGlobal))
            MessageI(:) = 0
            IF(.NOT.UseWriterCore.AND.NumParticlesLocal.GT.0)THEN
               ArrayStart = SUM(ParticlesOnLocalCore(                  &
     &                          1:NumTrackerCores-1)) + 1
               ArrayEnd = NumParticlesGlobal
               IF(IA.EQ.1)THEN
                  MessageI(ArrayStart:ArrayEnd)                        &
     &                     = PartNumberLocal(1:NumParticlesLocal)
               ELSEIF(IA.EQ.2)THEN
                  MessageI(ArrayStart:ArrayEnd)                        &
     &                     = PartDomainLocal(1:NumParticlesLocal)
               ENDIF
            ENDIF
         ENDIF
         CALL LPT_Comm_Collect(DummyC,MessageI,DummyR,"I",             &
     &            NumParticlesLocal,                                   &
     &            "The local locations of the particles were not "     &
     &         //"received by the writer core.")
#endif
         IF(MyRank+1.EQ.NumRanks)THEN
            IF(IA.EQ.1)THEN
               IF(ALLOCATED(PartNumberGlobal))THEN
                  DEALLOCATE(PartNumberGlobal)
               ENDIF
               IF(.NOT.ALLOCATED(PartNumberGlobal))THEN
                  ALLOCATE(PartNumberGlobal(1:NumParticlesGlobal))
               ENDIF
#ifdef MPI
               PartNumberGlobal(:) = MessageI(:)
#else
               PartNumberGlobal(:)                                     &
     &                   = PartNumberLocal(1:NumParticlesGlobal)
#endif
#if VERBOSE > 1
               CALL LPT_Print(0,"INFO","The global set of particle "   &
     &                  //"numbers is ready for the next snap of "     &
     &                  //"output.")
#endif
            ELSEIF(IA.EQ.2)THEN
               IF(ALLOCATED(PartDomainGlobal))THEN
                  DEALLOCATE(PartDomainGlobal)
               ENDIF
               IF(.NOT.ALLOCATED(PartDomainGlobal))THEN
                  ALLOCATE(PartDomainGlobal(1:NumParticlesGlobal))
               ENDIF
#ifdef MPI
               PartDomainGlobal(:) = MessageI(:)
#else
               PartDomainGlobal(:)                                     &
     &                   = PartDomainLocal(1:NumParticlesGlobal)
#endif
#if VERBOSE > 1
               CALL LPT_Print(0,"INFO","The global set of particle "   &
     &                  //"domains is ready for the next snap of "     &
     &                  //"output.")
#endif
            ENDIF
         ENDIF
      ENDDO
      DO IA=1,4
         IF(IA.EQ.4.AND.INDEX(BuoyancyMethod,"NULL").GT.0) CYCLE
#ifdef MPI
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         IF(MyRank+1.NE.NumRanks)THEN
            IF(NumParticlesLocal.GT.0)THEN
               ALLOCATE(MessageR(NumParticlesLocal))
               IF(IA.EQ.1)THEN
                  MessageR(:) = PartLonLocal(1:NumParticlesLocal)
               ELSEIF(IA.EQ.2)THEN
                  MessageR(:) = PartLatLocal(1:NumParticlesLocal)
               ELSEIF(IA.EQ.3)THEN
                  MessageR(:) = PartDepthLocal(1:NumParticlesLocal)
               ELSEIF(IA.EQ.4)THEN
                  MessageR(:) = PartDiameterLocal(1:NumParticlesLocal)
               ENDIF
            ENDIF
         ELSEIF(MyRank+1.EQ.NumRanks)THEN
            ALLOCATE(MessageR(NumParticlesGlobal))
            MessageR(:) = 0.D0
            IF(.NOT.UseWriterCore.AND.NumParticlesLocal.GT.0)THEN
               ArrayStart = SUM(ParticlesOnLocalCore(                  &
     &                          1:NumTrackerCores-1)) + 1
               ArrayEnd = NumParticlesGlobal
               IF(IA.EQ.1)THEN
                  MessageR(ArrayStart:ArrayEnd)                        &
     &                     = PartLonLocal(1:NumParticlesLocal)
               ELSEIF(IA.EQ.2)THEN
                  MessageR(ArrayStart:ArrayEnd)                        &
     &                     = PartLatLocal(1:NumParticlesLocal)
               ELSEIF(IA.EQ.3)THEN
                  MessageR(ArrayStart:ArrayEnd)                        &
     &                     = PartDepthLocal(1:NumParticlesLocal)
               ELSEIF(IA.EQ.4)THEN
                  MessageR(ArrayStart:ArrayEnd)                        &
     &                     = PartDiameterLocal(1:NumParticlesLocal)
               ENDIF
            ENDIF
         ENDIF
         CALL LPT_Comm_Collect(DummyC,DummyI,MessageR,"R",             &
     &            NumParticlesLocal,                                   &
     &            "The local locations of the particles were not "     &
     &            //"received by the writer core.")
#endif
         IF(MyRank+1.EQ.NumRanks)THEN
            IF(IA.EQ.1)THEN
               IF(ALLOCATED(PartLonGlobal))THEN
                  DEALLOCATE(PartLonGlobal)
               ENDIF
               IF(.NOT.ALLOCATED(PartLonGlobal))THEN
                  ALLOCATE(PartLonGlobal(NumParticlesGlobal))
               ENDIF
#ifdef MPI
               DO IP=1,NumParticlesGlobal
                  PartLonGlobal(PartNumberGlobal(IP))=MessageR(IP)
               ENDDO
#else
               PartLonGlobal(:) = PartLonLocal(1:NumParticlesGlobal)
#endif
#if VERBOSE > 1
               CALL LPT_Print(0,"INFO","The global set of particle "   &
     &                  //"longitudes is ready for the next snap of "  &
     &                  //"output.")
#endif
            ELSEIF(IA.EQ.2)THEN
               IF(ALLOCATED(PartLatGlobal))THEN
                  DEALLOCATE(PartLatGlobal)
               ENDIF
               IF(.NOT.ALLOCATED(PartLatGlobal))THEN
                  ALLOCATE(PartLatGlobal(NumParticlesGlobal))
               ENDIF
#ifdef MPI
               DO IP=1,NumParticlesGlobal
                  PartLatGlobal(PartNumberGlobal(IP))=MessageR(IP)
               ENDDO
#else
               PartLatGlobal(:) = PartLatLocal(1:NumParticlesGlobal)
#endif
#if VERBOSE > 1
               CALL LPT_Print(0,"INFO","The global set of particle "   &
     &                  //"latitudes is ready for the next snap of "   &
     &                  //"output.")
#endif
            ELSEIF(IA.EQ.3)THEN
               IF(ALLOCATED(PartDepthGlobal))THEN
                  DEALLOCATE(PartDepthGlobal)
               ENDIF
               IF(.NOT.ALLOCATED(PartDepthGlobal))THEN
                  ALLOCATE(PartDepthGlobal(NumParticlesGlobal))
               ENDIF
#ifdef MPI
               DO IP=1,NumParticlesGlobal
                  PartDepthGlobal(PartNumberGlobal(IP))=MessageR(IP)
               ENDDO
#else
               PartDepthGlobal(:) = PartDepthLocal(1:NumParticlesGlobal)
#endif
#if VERBOSE > 1
               CALL LPT_Print(0,"INFO","The global set of particle "   &
     &                  //"depths is ready for the next snap of "      &
     &                  //"output.")
#endif
            ELSEIF(IA.EQ.4)THEN
               IF(ALLOCATED(PartDiameterGlobal))THEN
                  DEALLOCATE(PartDiameterGlobal)
               ENDIF
               IF(.NOT.ALLOCATED(PartDiameterGlobal))THEN
                  ALLOCATE(PartDiameterGlobal(NumParticlesGlobal))
               ENDIF
#ifdef MPI
               DO IP=1,NumParticlesGlobal
                  PartDiameterGlobal(PartNumberGlobal(IP))=MessageR(IP)
               ENDDO
#else
               PartDiameterGlobal(:)                                   &
     &                        = PartDiameterLocal(1:NumParticlesGlobal)
#endif
#if VERBOSE > 1
               CALL LPT_Print(0,"INFO","The global set of particle "   &
     &                  //"diameters is ready for the next snap of "   &
     &                  //"output.")
#endif
           ENDIF
         ENDIF
      ENDDO

      ! The tracker cores can go back to, you know, tracking.
      IF(MyRank+1.NE.NumRanks)THEN
         RETURN
      ENDIF

#if VERBOSE > 0
      CALL LPT_Print(MyRank,"INFO","The particle locations are "       &
     &         //"being written.")
#endif

      Conversion = 1.D0
      IF(INDEX(ParticleInputCoordinates,"POLAR").GT.0)THEN
         Conversion = Conversion / EarthRadius
         Conversion = Conversion / Deg2Rad
      ENDIF

      IF(INDEX(OutputFileFormat,"ASCII").GT.0)THEN

         IF(FirstCall)THEN

            OPEN(UNIT=12,FILE=TRIM(ADJUSTL(OutputFileName)),           &
     &           ACTION='WRITE')

            WRITE(12,'(A)') TRIM(ADJUSTL(HeaderParticles))
            WRITE(JunkC1,9001) NumTrackingSnaps+1
            WRITE(JunkC2,9001) NumParticlesGlobal
            WRITE(JunkC3,9002) (SimulationLength/NumTrackingSnaps)*3600
            WRITE(12,'(A)') TRIM(ADJUSTL(JunkC1))//" "                 &
     &            //TRIM(ADJUSTL(JunkC2))//" "//TRIM(ADJUSTL(JunkC3))
#if VERBOSE > 1
            CALL LPT_Print(0,"INFO","The header was written "          &
     &               //"successfully to the output file.")
#endif

         ENDIF

         WRITE(JunkC1,9002) Time1
         WRITE(JunkC2,9001) NumParticlesGlobal
         WRITE(12,'(A)') TRIM(ADJUSTL(JunkC1))//" "                    &
     &         //TRIM(ADJUSTL(JunkC2))
#if VERBOSE > 2
         IF(MyRank.EQ.0)THEN
            IF(NumParticlesGlobal.GE.100)THEN
               WRITE(*,'(A,$)') "LPT: INFO: Writing the particle "     &
     &               //"locations: +"
            ELSE
               CALL LPT_Print(0,"INFO","Writing the particle "         &
     &                  //"locations.")
            ENDIF
         ENDIF
#endif
         DO IP=1,NumParticlesGlobal
            WRITE(JunkC1,9001) IP
            WRITE(JunkC2,9002) PartLonGlobal(IP)*Conversion
            WRITE(JunkC3,9002) PartLatGlobal(IP)*Conversion
            WRITE(JunkC4,9002) PartDepthGlobal(IP)
            IF(INDEX(BuoyancyMethod,"NULL").GT.0)THEN
               WRITE(JunkC5,'(A)') " "
            ELSE
               WRITE(JunkC5,9002) PartDiameterGlobal(IP)
            ENDIF
            WRITE(12,'(A)') TRIM(ADJUSTL(JunkC1))//" "                 &
     &            //TRIM(ADJUSTL(JunkC2))//" "//TRIM(ADJUSTL(JunkC3))  &
     &            //" "//TRIM(ADJUSTL(JunkC4))//" "                    &
     &            //TRIM(ADJUSTL(JunkC5))
#if VERBOSE > 2
            IF(MyRank.EQ.0)THEN
               IF(NumParticlesGlobal.GE.100)THEN
                  CALL LPT_Progress(IP,NumParticlesGlobal)
               ENDIF
            ENDIF
#endif
         ENDDO
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The next snap of particle "          &
     &            //"locations was written successfully to the "       &
     &            //"output file.")
#endif

#ifdef NETCDF
      ELSEIF(INDEX(OutputFileFormat,"NETCDF").GT.0)THEN

         NC_ID = NC_ID_Output

         IF(FirstCall)THEN

            ! Create the output file.
#if NETCDF < 4
            CALL LPT_NetCDF_Check( NF90_CREATE(                        &
     &               TRIM(ADJUSTL(OutputFileName)),                    &
     &               NF90_64BIT_OFFSET,NC_ID) )
#elif NETCDF == 4
!           CALL LPT_NetCDF_Check( NF90_CREATE(                        &
!    &               TRIM(ADJUSTL(OutputFileName)),                    &
!    &               4096,NC_ID) )
!    &               NF90_NETCDF4,NC_ID) )
            CALL LPT_NetCDF_Check( NF90_CREATE(                        &
      &              TRIM(ADJUSTL(OutputFileName)),                    &
      &              IOR(NF90_CLASSIC_MODEL,NF90_HDF5),NC_ID))
#endif

            ! Define the dimensions.
            CALL LPT_NetCDF_Check( NF90_DEF_DIM(NC_ID,"ntimesnap",     &
     &               NF90_UNLIMITED,NC_DIM_Time) )
            ! If the particles are being added from a source,
            ! then the total number of particles will change
            ! in time.  NetCDF version 4 allows for a second
            ! unlimited dimension.
            IF(INDEX(ParticleInputMethod,"READFROMFILE").GT.0.OR.      &
     &         INDEX(ParticleInputMethod,"SOURCE").GT.0)THEN
#if NETCDF < 4
               CALL LPT_NetCDF_Check( NF90_DEF_DIM(NC_ID,"nparticle",  &
     &                  NumTrackingSnaps+1,NC_DIM_Part) )
#elif NETCDF == 4
               CALL LPT_NetCDF_Check( NF90_DEF_DIM(NC_ID,"nparticle",  &
     &                  NF90_UNLIMITED,NC_DIM_Part) )
#endif
            ! Otherwise, the total number of particles will be
            ! constant in time.
            ELSE
               CALL LPT_NetCDF_Check( NF90_DEF_DIM(NC_ID,"nparticle",  &
     &                  NumParticlesGlobal,NC_DIM_Part) )
            ENDIF

            ! Define the variables.
            CALL LPT_NetCDF_Check( NF90_DEF_VAR(NC_ID,"Time_index",    &
     &               NF90_DOUBLE,NC_DIM_Time,NC_VAR_Time) )
            CALL LPT_NetCDF_Check( NF90_DEF_VAR(NC_ID,"X_particle",    &
     &               NF90_DOUBLE,(/NC_DIM_Part,NC_DIM_Time/),NC_VAR_X) )
            CALL LPT_NetCDF_Check( NF90_DEF_VAR(NC_ID,"Y_particle",    &
     &               NF90_DOUBLE,(/NC_DIM_Part,NC_DIM_Time/),NC_VAR_Y) )
            CALL LPT_NetCDF_Check( NF90_DEF_VAR(NC_ID,"Z_particle",    &
     &               NF90_DOUBLE,(/NC_DIM_Part,NC_DIM_Time/),NC_VAR_Z) )
            IF(INDEX(BuoyancyMethod,"NULL").LE.0)THEN
               CALL LPT_NetCDF_Check( NF90_DEF_VAR(NC_ID,              &
     &                  "Particle_diam",NF90_DOUBLE,(/NC_DIM_Part,     &
     &                  NC_DIM_Time/),NC_VAR_D) )
            ENDIF
#if NETCDF == 4
            CALL LPT_NetCDF_Check( NF90_DEF_VAR(NC_ID,                 &
     &               "Particles_per_time_snap",NF90_INT,               &
     &               NC_DIM_Time,NC_Var_NumParticles) )
!...........Cobell - Turn on NetCDF4 compression
            CALL LPT_NetCDF_Check(NF90_DEF_VAR_DEFLATE(NC_ID,          &
     &               NC_VAR_Time,1,1,2))
            CALL LPT_NetCDF_Check( NF90_DEF_VAR_DEFLATE(NC_ID,NC_VAR_X,&
     &               1,1,2))
            CALL LPT_NetCDF_Check( NF90_DEF_VAR_DEFLATE(NC_ID,NC_VAR_Y,&
     &               1,1,2))
            CALL LPT_NetCDF_Check( NF90_DEF_VAR_DEFLATE(NC_ID,NC_VAR_Z,&
     &               1,1,2))
            CALL LPT_NetCDF_Check( NF90_DEF_VAR_DEFLATE(NC_ID,         &
     &               NC_Var_NumParticles,1,1,2))
#endif

            CALL LPT_NetCDF_Check( NF90_PUT_ATT(NC_ID,NC_Var_X,        &
     &               '_FillValue',FillValueR) )
            CALL LPT_NetCDF_Check( NF90_PUT_ATT(NC_ID,NC_Var_Y,        &
     &               '_FillValue',FillValueR) )
            CALL LPT_NetCDF_Check( NF90_PUT_ATT(NC_ID,NC_Var_Z,        &
     &               '_FillValue',FillValueR) )
            CALL LPT_NetCDF_Check( NF90_PUT_ATT(NC_ID,NC_Var_D,        &
     &               '_FillValue',FillValueR) )
#if NETCDF == 4
            CALL LPT_NetCDF_Check( NF90_PUT_ATT(NC_ID,                 &
     &               NC_Var_NumParticles,'_FillValue',FillValueI) )
#endif

            ! End the definition mode.
            CALL LPT_NetCDF_Check( NF90_ENDDEF(NC_ID) )

#if VERBOSE > 1
            CALL LPT_Print(0,"INFO","The header was written "          &
     &               //"successfully to the output file.")
#endif

         ENDIF

         NC_Counter = NC_Counter + 1
         IF(ALLOCATED(NC_Temp)) DEALLOCATE(NC_Temp)
         IF(.NOT.ALLOCATED(NC_Temp))THEN
            ALLOCATE(NC_Temp(NumParticlesGlobal))
         ENDIF

         ! Write the current time.
         CALL LPT_NetCDF_Check( NF90_PUT_VAR(NC_ID,NC_VAR_Time,        &
     &            (/Time1/),START=(/NC_Counter/),COUNT=(/1/)) )

#if NETCDF == 4
         ! Write the number of particles in this snap.
         CALL LPT_NetCDF_Check( NF90_PUT_VAR(NC_ID,                    &
     &            NC_Var_NumParticles,(/NumParticlesGlobal/),          &
     &            START=(/NC_Counter/),COUNT=(/1/)) )
#endif

         ! Write the particle longitudes.
         DO IP=1,NumParticlesGlobal
            NC_Temp(IP) = PartLonGlobal(IP) * Conversion
         ENDDO
         CALL LPT_NetCDF_Check( NF90_PUT_VAR(NC_ID,NC_VAR_X,           &
     &            NC_Temp(:),START=(/1,NC_Counter/),                   &
     &            COUNT=(/NumParticlesGlobal,1/)) )

         ! Write the particle latitudes.
         DO IP=1,NumParticlesGlobal
            NC_Temp(IP) = PartLatGlobal(IP) * Conversion
         ENDDO
         CALL LPT_NetCDF_Check( NF90_PUT_VAR(NC_ID,NC_VAR_Y,           &
     &            NC_Temp(:),START=(/1,NC_Counter/),                   &
     &            COUNT=(/NumParticlesGlobal,1/)) )

         ! Write the particle classes/depths.
         DO IP=1,NumParticlesGlobal
            NC_Temp(IP) = PartDepthGlobal(IP)
         ENDDO
         CALL LPT_NetCDF_Check( NF90_PUT_VAR(NC_ID,NC_VAR_Z,           &
     &            NC_Temp(:),START=(/1,NC_Counter/),                   &
     &            COUNT=(/NumParticlesGlobal,1/)) )

         ! Write the particle diameters.
         IF(INDEX(BuoyancyMethod,"NULL").LE.0)THEN
            DO IP=1,NumParticlesGlobal
               NC_Temp(IP) = PartDiameterGlobal(IP)
            ENDDO
            CALL LPT_NetCDF_Check( NF90_PUT_VAR(NC_ID,NC_VAR_D,        &
     &               NC_Temp(:),START=(/1,NC_Counter/),                &
     &               COUNT=(/NumParticlesGlobal,1/)) )
         ENDIF

#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The next snap of particle "          &
     &            //"locations was written successfully to the "       &
     &            //"output file.")
#endif

         IF(NC_Counter.EQ.NumTrackingSnaps+2)THEN
            CALL LPT_NetCDF_Check( NF90_CLOSE(NC_ID) )
#if VERBOSE > 1
            CALL LPT_Print(0,"INFO","The output file was closed.")
#endif
         ENDIF

         NC_ID_Output = NC_ID
#endif

      ENDIF

      IF(FirstCall)THEN
         FirstCall = .FALSE.
      ENDIF

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Write_Particle_Snap "   &
     &         //"routine.")
#endif

      ! We're done here.
      RETURN

 9001 FORMAT(I24)
 9002 FORMAT(F30.10)

      END SUBROUTINE LPT_Write_Particle_Snap

      END MODULE

