      MODULE LPT_Read_Module

      IMPLICIT NONE

      CHARACTER(LEN=500)  :: InputC
      INTEGER,ALLOCATABLE :: InputI(:)
      REAL(8),ALLOCATABLE :: InputR(:)

      REAL(8),ALLOCATABLE,TARGET :: SnapUCurr(:,:,:)
      REAL(8),ALLOCATABLE,TARGET :: SnapUWind(:,:,:)
      REAL(8),ALLOCATABLE,TARGET :: SnapVCurr(:,:,:)
      REAL(8),ALLOCATABLE,TARGET :: SnapVWind(:,:,:)
      REAL(8),ALLOCATABLE,TARGET :: SnapWCurr(:,:,:)
      REAL(8)                    :: Times(2)

      CONTAINS



      SUBROUTINE LPT_Read

      USE LPT_Comm_Module
      USE LPT_Data_Module

      IMPLICIT NONE

      INTEGER :: Snap

      REAL(8) :: JunkR

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the LPT_Read routine.")
#endif

      ! Read the user-specified parameters from the input file.
      CALL LPT_Read_Parameters

      ! Read the initial locations of the particles, if necessary.
      IF(INDEX(ParticleInputMethod,"READINITIALCONDITIONS").GT.0)THEN
         CALL LPT_Read_Initial_Conditions
      ELSEIF(INDEX(ParticleInputMethod,"READFROMFILE").GT.0)THEN
         NumParticlesLocal = 0
         WRITE(HeaderParticles,'(A,A)') "Particles from "// &
     &         TRIM(ADJUSTL(ParticleFile))
      ELSE
         NumParticlesLocal = 0
         WRITE(HeaderParticles,'(A)') "Particles from Source"
      ENDIF

      ! Read the unstructured mesh.
      IF(INDEX(MeshFileOrigin,"ADCIRC").GT.0)THEN
         CALL LPT_Read_Unstructured_Mesh
      ELSEIF(INDEX(MeshFileOrigin,"HYCOM").GT.0)THEN
         CALL LPT_Read_HYCOM_Mesh
      ENDIF

      ! These calls will open the files and read the headers,
      ! but not actually read any information.  Some of this info,
      ! such as the number of layers (NumLayers), will be needed
      ! in the next few routines.
      CALL LPT_Read_Velocity_Field("C",Times(1))
      CALL LPT_Read_Velocity_Field("W",Times(1))

#ifdef MPI
      ! Send the number of layers, so the tracking cores
      ! can allocate the arrays.
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = NumLayers
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The number of layers could not be passed to "       &
     &            //"the tracking cores.")
         NumLayers = MessageI(1)
         ! Send the sigma layers if necessary.
         IF(NumLayers.GT.1)THEN
            IF(.NOT.ALLOCATED(Sigma)) ALLOCATE(Sigma(1:NumLayers))
            IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
            ALLOCATE(MessageR(NumLayers))
            MessageR(:) = Sigma(:)
            CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",       &
     &               "The sigma layers could not be passed to the "    &
     &               //"tracking cores.")
            Sigma(:) = MessageR(:)
         ENDIF
      ENDIF
#endif   

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Read routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_Parameters

      USE LPT_Comm_Module
      USE LPT_Data_Module

      IMPLICIT NONE

      CHARACTER(LEN=100)  :: InputFile
      CHARACTER(LEN=100)  :: JunkC
      CHARACTER(LEN=100)  :: JunkC1

      INTEGER             :: IA
      INTEGER             :: IARGC
      INTEGER             :: IERR

      LOGICAL             :: InputFileExists
      LOGICAL             :: InputFileFound

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the LPT_Read_Parameters "      &
     &         //"routine.")
#endif

      ! Obtain the input file name.
      IF(MyRank.EQ.0)THEN
         InputFileFound = .FALSE.
         IF(IARGC().GT.0)THEN
            IA = 0
            DO WHILE (IA.LT.IARGC())
               IA = IA + 1
               CALL GETARG(IA,JunkC)
               IF(JunkC(1:2).EQ."-I")THEN
                  InputFileFound = .TRUE.
                  IA = IA + 1
                  CALL GETARG(IA,JunkC)
                  READ(JunkC,'(A)') InputFile
               ENDIF
            ENDDO
         ENDIF
         IF(.NOT.InputFileFound)THEN
            WRITE(*,'(A,$)') "Enter name of input file: "
            READ(*,'(A)') InputFile
         ENDIF
#ifdef DEBUG
         INQUIRE(FILE=TRIM(InputFile),EXIST=InputFileExists)
         IF(.NOT.InputFileExists)THEN
            CALL LPT_Print(0,"FATAL ERROR","Could not find the input " &
     &               //"file: "//TRIM(InputFile)//".")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(.NOT.UseReaderCore)THEN
         WRITE(MessageC,'(A)') TRIM(InputFile)
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The name of the input file could not be passed to " &
     &            //"the tracker cores.")
         WRITE(InputFile,'(A)') TRIM(MessageC)
      ENDIF
#endif
#if VERBOSE > 0
      CALL LPT_Print(0,"INFO","The input file is "                     &
     &         //TRIM(InputFile)//".")
#endif

      ! Open the input file.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         OPEN(UNIT=11,FILE=TRIM(InputFile),ACTION='READ',              &
     &        DELIM='APOSTROPHE',IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR","The input file "      &
     &            //"could not be opened for reading.")
         ENDIF
#endif
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The input file was opened.")
#endif
      ENDIF

      ! Read the timing information.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=TimingInformation,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR","The timing "          &
     &            //"information could not be read from the input "    &
     &            //"file.")
         ENDIF
#endif
         SimulationLength = SimulationLength/3600.D0
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1))
         MessageR(1) = SimulationLength
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The tracking simulation length could not be "       &
     &            //"passed to the tracking cores.")
         SimulationLength = MessageR(1)
      ENDIF
      IF(UseWriterCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1))
         MessageR(1) = SimulationLength
         CALL LPT_Comm_To_Writer(DummyC,DummyI,MessageR,"R",           &
     &            "The tracking simulation length could not be "       &
     &            //"passed to the writer core.")
         SimulationLength = MessageR(1)
      ENDIF
#endif
#ifdef DEBUG
      IF(SimulationLength.LT.-98.D0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The tracking simulation "     &
     &            //"length must be provided in the input file.")
      ENDIF
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1))
         MessageR(1) = StartingTime
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The starting time of the tracking simulation "      &
     &            //"could not be passed to the tracking cores.")
         StartingTime = MessageR(1)
      ENDIF
      IF(UseWriterCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1))
         MessageR(1) = StartingTime
         CALL LPT_Comm_To_Writer(DummyC,DummyI,MessageR,"R",           &
     &            "The starting time of the tracking simulation "      &
     &            //"could not be passed to the writer core.")
         StartingTime = MessageR(1)
      ENDIF
#endif
#ifdef DEBUG
      IF(StartingTime.LT.-98.D0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The starting time of the "    &
     &            //"tracking simulation must be provided in the "     &
     &            //"input file.")
      ENDIF
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = NumTrackingSnaps
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The number of tracking snaps could not be passed "  &
     &            //"to the tracking cores.")
         NumTrackingSnaps = MessageI(1)
      ENDIF
      IF(UseWriterCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = NumTrackingSnaps
         CALL LPT_Comm_To_Writer(DummyC,MessageI,DummyR,"I",           &
     &            "The number of tracking snaps could not be passed "  &
     &            //"to the writer core.")
         NumTrackingSnaps = MessageI(1)
      ENDIF
#endif
#ifdef DEBUG
      IF(NumTrackingSnaps.LT.-98)THEN
         CALL LPT_Print(0,"FATAL ERROR","The number of tracking snaps "&
     &            //"must be provided in the input file.")
      ENDIF
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1))
         MessageR(1) = MinimumTimeStep
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The minimum time step could not be passed to the "  &
     &            //"tracking cores.")
         MinimumTimeStep = MessageR(1)
      ENDIF
      IF(UseWriterCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1))
         MessageR(1) = MinimumTimeStep
         CALL LPT_Comm_To_Writer(DummyC,DummyI,MessageR,"R",           &
     &            "The minimum time step could not be passed to the "  &
     &            //"writer core.")
         MinimumTimeStep = MessageR(1)
      ENDIF
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = LatticeSearchBins
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The number of lattice search bins could not be "    &
     &            //"passed to the tracking cores.")
         LatticeSearchBins = MessageI(1)
      ENDIF
      IF(UseWriterCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = LatticeSearchBins
         CALL LPT_Comm_To_Writer(DummyC,MessageI,DummyR,"I",           &
     &            "The number of lattice search bins could not be "    &
     &            //"passed to the writer core.")
         LatticeSearchBins = MessageI(1)
      ENDIF
#endif
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         WRITE(JunkC,'(F15.6)') SimulationLength
         CALL LPT_Print(0,"INFO","The tracking simulation length is "  &
     &            //TRIM(ADJUSTL(JunkC))//" hours.")
         WRITE(JunkC,'(F15.6)') StartingTime/3600.D0
         CALL LPT_Print(0,"INFO","The starting time of the tracking "  &
     &            //"is "//TRIM(ADJUSTL(JunkC))//" hours.")
         WRITE(JunkC,'(I24)') NumTrackingSnaps
         CALL LPT_Print(0,"INFO","There are "//TRIM(ADJUSTL(JunkC))    &
     &            //" tracking snaps.")
         WRITE(JunkC,'(F15.6)') MinimumTimeStep
         CALL LPT_Print(0,"INFO","The minimum time step is "           &
     &            //TRIM(ADJUSTL(JunkC))//" hours.")
         WRITE(JunkC,'(I24)') LatticeSearchBins
         CALL LPT_Print(0,"INFO","The lattice search will use "        &
     &            //TRIM(ADJUSTL(JunkC))//" bins.")
      ENDIF
#endif

      ! Read the information about the current velocity field.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=CurrentField,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR","The parameters for "  &
     &               //"the current field could not be read from the " &
     &               //"input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') TRIM(CurrentFile)
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The name of the current velocities file could not " &
     &            //"be passed to the tracking cores.")
         WRITE(CurrentFile,'(A)') TRIM(MessageC)
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') TRIM(CurrentFile)
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The name of the current velocities file could not " &
     &            //"be passed to the writer core.")
         WRITE(CurrentFile,'(A)') TRIM(MessageC)
      ENDIF
#endif
#ifdef DEBUG
      IF(INDEX(CurrentFile,"NULL").GT.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The name of the current "     &
     &            //"velocities file must be provided in the input "   &
     &            //"file.")
      ENDIF
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') TRIM(CurrentFileOrigin)
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The origin of the current velocities file could "   &
     &            //"not be passed to the tracking cores.")
         WRITE(CurrentFileOrigin,'(A)') TRIM(MessageC)
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') TRIM(CurrentFileOrigin)
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The origin of the current velocities file could "   &
     &            //"not be passed to the writer core.")
         WRITE(CurrentFileOrigin,'(A)') TRIM(MessageC)
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(CurrentFileOrigin)
#ifdef DEBUG
      IF(INDEX(CurrentFileOrigin,"ADCIRC").LE.0.AND.                   &
     &   INDEX(CurrentFileOrigin,"HYCOM").LE.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The origin of the current "   &
     &            //"velocities file must be specified as either "     &
     &            //"'ADCIRC' or 'NetCDF' in the input file.")
       ENDIF
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') TRIM(CurrentFileFormat)
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The format of the current velocities file could "   &
     &            //"not be passed to the tracking cores.")
         WRITE(CurrentFileFormat,'(A)') TRIM(MessageC)
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') TRIM(CurrentFileFormat)
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The format of the current velocities file could "   &
     &            //"not be passed to the writer core.")
         WRITE(CurrentFileFormat,'(A)') TRIM(MessageC)
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(CurrentFileFormat)
#ifdef DEBUG
      IF(INDEX(CurrentFileFormat,"NULL").GT.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The format of the current "   &
     &            //"velocities file must be specified in the input "  &
     &            //"file.")
      ENDIF
      IF((INDEX(CurrentFileFormat,"ASCII").LE.0).AND.                  &
     &   (INDEX(CurrentFileFormat,"NETCDF").LE.0))THEN
         CALL LPT_Print(0,"FATAL ERROR","The format of the current "   &
     &            //"velocities file must be specified as either "     &
     &            //"'ASCII' or 'NetCDF' in the input file.")
      ENDIF
      IF(INDEX(CurrentFileOrigin,"HYCOM").GT.0.AND.                    &
     &   INDEX(CurrentFileFormat,"NETCDF").LE.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The current velocities file " &
     &            //"was obtained from HYCOM, so its format must be "  &
     &            //"specified as 'NetCDF' in the input file.")
      ENDIF
#ifndef NETCDF
      IF(INDEX(CurrentFileFormat,"NETCDF").GT.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The current velocities file " &
     &            //"is formatted for NetCDF, but this program has "   &
     &            //"not been compiled with support for NetCDF.")
      ENDIF
#endif
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') TRIM(CurrentDimensions)
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The dimensions of the current velocities could "    &
     &            //"not be passed to the tracking cores.")
         WRITE(CurrentDimensions,'(A)') TRIM(MessageC)
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') TRIM(CurrentDimensions)
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The dimensions of the current velocities could "    &
     &            //"not be passed to the writer core.")
         WRITE(CurrentDimensions,'(A)') TRIM(MessageC)
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(CurrentDimensions)
#ifdef DEBUG
      IF((INDEX(CurrentDimensions,"2D").LE.0).AND.                     &
     &   (INDEX(CurrentDimensions,"3D").LE.0))THEN
         CALL LPT_Print(0,"FATAL ERROR","The current velocity field "  &
     &            //"must be specified in 2D or 3D.")
      ENDIF
      IF(INDEX(CurrentFileOrigin,"HYCOM").GT.0.AND.                    &
     &   INDEX(CurrentDimensions,"3D").LE.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The current velocity field "  &
     &            //"was obtained from HYCOM, so its dimensions must " &
     &            //"specified as 3D in the input file.")
      ENDIF
#endif
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         WRITE(JunkC,'(A)') TRIM(CurrentFile)
         CALL LPT_Print(0,"INFO","The current velocities will be "     &
     &            //"read from "//TRIM(ADJUSTL(JunkC))//".")
         WRITE(JunkC,'(A)') TRIM(CurrentFileOrigin)
         CALL LPT_Print(0,"INFO","The current velocities file was "    &
     &            //"obtained from "//TRIM(ADJUSTL(JunkC))//".")
         WRITE(JunkC,'(A)') TRIM(CurrentFileFormat)
         CALL LPT_Print(0,"INFO","The current velocities file is in "  &
     &            //TRIM(ADJUSTL(JunkC))//" format.")
         WRITE(JunkC,'(A)') TRIM(CurrentDimensions)
         CALL LPT_Print(0,"INFO","The current velocity field is in "   &
     &            //TRIM(ADJUSTL(JunkC))//".")
      ENDIF
#endif

      ! Read the information about the wind velocity field.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=WindField,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"WARNING","The parameters for "      &
     &               //"the wind field could not be read from the "    &
     &               //"input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') TRIM(WindFile)
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The name of the wind velocities file could not be " &
     &            //"passed to the tracking cores.")
         WRITE(WindFile,'(A)') TRIM(MessageC)
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') TRIM(WindFile)
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The name of the wind velocities file could not be " &
     &            //"passed to the writer core.")
         WRITE(WindFile,'(A)') TRIM(MessageC)
      ENDIF
#endif
      IF(INDEX(WindFile,"NULL").LE.0)THEN
#ifdef MPI
         IF(UseReaderCore)THEN
            WRITE(MessageC,'(A)') TRIM(WindFileFormat)
            CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",       &
     &               "The format of the wind velocities file could "   &
     &               //"not be passed to the tracking cores.")
            WRITE(WindFileFormat,'(A)') TRIM(MessageC)
         ENDIF
         IF(UseWriterCore)THEN
            WRITE(MessageC,'(A)') TRIM(WindFileFormat)
            CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",        &
     &               "The format of the wind velocities file could "   &
     &               //"be passed to the writer core.")
            WRITE(WindFileFormat,'(A)') TRIM(MessageC)
         ENDIF
#endif
         CALL LPT_Read_Capitalize_Word(WindFileFormat)
#ifdef DEBUG
         IF(INDEX(WindFileFormat,"NULL").GT.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The format of the wind "   &
     &               //"velocities file must be specified in the "     &
     &               //"input file.")
         ENDIF
         IF((INDEX(WindFileFormat,"ASCII").LE.0).AND.                  &
     &      (INDEX(WindFileFormat,"NETCDF").LE.0))THEN
            CALL LPT_Print(0,"FATAL ERROR","The format of the wind "   &
     &               //"velocities file must be specified as either "  &
     &               //"ASCII or NetCDF.")
         ENDIF
#ifndef NETCDF
         IF(INDEX(WindFileFormat,"NETCDF").GT.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The wind velocities file " &
     &               //"is formatted for NetCDF, but this program has "&
     &               //"not been compiled with support for NetCDF.")
         ENDIF
#endif
#endif
#if VERBOSE > 0
         IF(MyRank.EQ.0)THEN
            WRITE(JunkC,'(A)') TRIM(WindFile)
            CALL LPT_Print(0,"INFO","The wind velocities will be "     &
     &               //"read from "//TRIM(ADJUSTL(JunkC))//".")
            WRITE(JunkC,'(A)') TRIM(WindFileFormat)
            CALL LPT_Print(0,"INFO","The wind velocities file is in "  &
     &               //TRIM(ADJUSTL(JunkC))//" format.")
         ENDIF
#endif
      ENDIF

      ! Read the information about the file containing the initial 
      ! locations of the particles.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=ParticleInput,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR","The information "     &
     &               //"about the particle input could not be read "   &
     &               //"from the input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') ParticleInputMethod
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The method for particle input could not be passed " &
     &            //"to the tracking cores.")
         WRITE(ParticleInputMethod,'(A)') MessageC
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') ParticleInputMethod
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The method for particle input could not be passed " &
     &            //"to the writer core.")
         WRITE(ParticleInputMethod,'(A)') MessageC
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(ParticleInputMethod)
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') ParticleInputCoordinates
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The coordinate system for the particles could "     &
     &            //"not be passed to the tracking cores.")
         WRITE(ParticleInputCoordinates,'(A)') MessageC
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') ParticleInputCoordinates
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The coordinate system for the particles could "     &
     &            //"not be passed to the write core.")
         WRITE(ParticleInputCoordinates,'(A)') MessageC
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(ParticleInputCoordinates)
      IF(INDEX(ParticleInputMethod,"READINITIALCONDITIONS").GT.0.OR.   &
     &   INDEX(ParticleInputMethod,"READFROMFILE").GT.0)THEN
#ifdef MPI
         IF(UseReaderCore)THEN
            WRITE(MessageC,'(A)') ParticleFile
            CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",       &
     &               "The name of the file containing the initial "    &
     &               //"locations of the particles could not be "      &
     &               //"passed to the tracking cores.")
            WRITE(ParticleFile,'(A)') MessageC
         ENDIF
         IF(UseWriterCore)THEN
            WRITE(MessageC,'(A)') ParticleFile
            CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",        &
     &               "The name of the file containing the initial "    &
     &               //"locations of the particles could not be "      &
     &               //"passed to the writer core.")
            WRITE(ParticleFile,'(A)') MessageC
         ENDIF
#endif
#ifdef DEBUG
         IF(INDEX(ParticleFile,"NULL").GT.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The name of the file "     &
     &               //"containing the initial locations of the "      &
     &               //"particles must be specified in the input file.")
         ENDIF
#endif
#ifdef MPI
         IF(UseReaderCore)THEN
            WRITE(MessageC,'(A)') ParticleFileFormat
            CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",       &
     &               "The format of the file containing the initial "  &
     &               //"locations of the particles could not be "      &
     &               //"to the tracking cores.")
            WRITE(ParticleFileFormat,'(A)') MessageC
         ENDIF
         IF(UseWriterCore)THEN
            WRITE(MessageC,'(A)') ParticleFileFormat
            CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",        &
     &               "The format of the file containing the initial "  &
     &               //"locations of the particles could not be "      &
     &               //"passed to the writer core.")
            WRITE(ParticleFileFormat,'(A)') MessageC
         ENDIF
#endif
         CALL LPT_Read_Capitalize_Word(ParticleFileFormat)
#ifdef DEBUG
         IF(INDEX(ParticleFileFormat,"NULL").GT.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The format of the file "   &
     &               //"containing the initial locations of the "      &
     &               //"particles must be specified in the input file.")
         ENDIF
         IF((INDEX(ParticleFileFormat,"ASCII").LE.0).AND.            &
     &      (INDEX(ParticleFileFormat,"NETCDF").LE.0))THEN
            CALL LPT_Print(0,"FATAL ERROR","The format of the file "   &
     &               //"containing the initial conditions of the "     &
     &               //"particles  must be specified as either "       &
     &               //"ASCII or NetCDF.")
         ENDIF
#ifndef NETCDF
         IF(INDEX(ParticleFileFormat,"NETCDF").GT.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The file containing the  " &
     &               //"initial locations of the particles is "        &
     &               //"formatted for NetCDF, but this program has "   &
     &               //"not been compiled with support for NetCDF.")
         ENDIF
#endif
#endif
#if VERBOSE > 0
         IF(MyRank.EQ.0)THEN
            WRITE(JunkC,'(A)') TRIM(ParticleFile)
            CALL LPT_Print(0,"INFO","The initial locations of the "    &
     &               //"particles will be read from "                  &
     &               //TRIM(ADJUSTL(JunkC))//".")
            WRITE(JunkC,'(A)') TRIM(ParticleFileFormat)
            CALL LPT_Print(0,"INFO","The file containing the initial " &
     &               //"locations of the particles is in "             &
     &               //TRIM(ADJUSTL(JunkC))//" format.")
            WRITE(JunkC,'(A)') TRIM(ParticleInputCoordinates)
            CALL LPT_Print(0,"INFO","The particles are given in "      &
     &               //TRIM(ADJUSTL(JunkC))//" coordinates.")
         ENDIF
#endif
      ELSEIF(INDEX(ParticleInputMethod,"SOURCE").GT.0)THEN
#ifdef MPI
         IF(UseReaderCore)THEN
            IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
            ALLOCATE(MessageR(1:3))
            MessageR(1) = ParticleSourceX
            MessageR(2) = ParticleSourceY
            MessageR(3) = ParticleSourceZ
            CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",       &
     &               "The (x,y,z) location of the particle source "    &
     &               //"could not be passed to the tracking cores.")
            ParticleSourceX = MessageR(1)
            ParticleSourceY = MessageR(2)
            ParticleSourceZ = MessageR(3)
         ENDIF
         IF(UseWriterCore)THEN
            IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
            ALLOCATE(MessageR(1:3))
            MessageR(1) = ParticleSourceX
            MessageR(2) = ParticleSourceY
            MessageR(3) = ParticleSourceZ
            CALL LPT_Comm_To_Writer(DummyC,DummyI,MessageR,"R",        &
     &               "The (x,y,z) location of the particle source "    &
     &               //"could not be passed to the writer core.")
            ParticleSourceX = MessageR(1)
            ParticleSourceY = MessageR(2)
            ParticleSourceZ = MessageR(3)
         ENDIF
#endif
#if VERBOSE > 0
         IF(MyRank.EQ.0)THEN
            WRITE(JunkC,'(F15.6)') ParticleSourceX
            CALL LPT_Print(0,"INFO","The particles will originate "    &
     &               //"from a x-location of "                         &
     &               //TRIM(ADJUSTL(JunkC))//".")
            WRITE(JunkC,'(F15.6)') ParticleSourceY
            CALL LPT_Print(0,"INFO","The particles will originate "    &
     &               //"from a y-location of "                         &
     &               //TRIM(ADJUSTL(JunkC))//".")
            WRITE(JunkC,'(F15.6)') ParticleSourceZ
            CALL LPT_Print(0,"INFO","The particles will originate "    &
     &               //"from a z-location of "                         &
     &               //TRIM(ADJUSTL(JunkC))//".")
         ENDIF
#endif
      ELSE
         CALL LPT_Print(0,"FATAL ERROR","The particle input method "   &
     &            //"must be specified as either 'ReadFromFile' or "   &
     &            //"'ReadInitialConditions' or 'Source' in the "      &
     &            //"parameters file.")
      ENDIF

      ! Read the information about the output settings.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=OutputSettings,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"WARNING","The information "         &
     &               //"about the output settings could not be read "  &
     &               //"from the input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') OutputFileName
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The output file name could not be passed to the "   &
     &            //"tracking cores.")
         WRITE(OutputFileName,'(A)') MessageC
         WRITE(MessageC,'(A)') OutputFileFormat
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The output file format could not be passed to the " &
     &            //"tracking cores.")
         WRITE(OutputFileFormat,'(A)') MessageC
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') OutputFileName
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The output file name could not be passed to the "   &
     &            //"writer core.")
         WRITE(OutputFileName,'(A)') MessageC
         WRITE(MessageC,'(A)') OutputFileFormat
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The output file format could not be passed to the " &
     &            //"writer core.")
         WRITE(OutputFileFormat,'(A)') MessageC
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(OutputFileFormat)
#ifdef DEBUG
      IF(INDEX(OutputFileName,"NULL").GT.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The output file name must "   &
     &            //"be specified in the parameters file.")
      ENDIF
#endif
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         CALL LPT_Print(0,"INFO","The output will be written to "      &
     &            //TRIM(ADJUSTL(OutputFileName))//".")
         CALL LPT_Print(0,"INFO","The output will be written in the "  &
     &            //"'"//TRIM(ADJUSTL(OutputFileFormat))//"' format.")
      ENDIF
#endif

      ! Read the information about the unstructured mesh.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=UnstructuredMesh,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR","The information "     &
     &               //"about the file containing the unstructured "   &
     &               //"mesh could not be read from the input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') MeshFile
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The name of the file containing the unstructured "  &
     &            //"mesh could not be passed to the tracking cores.")
         WRITE(MeshFile,'(A)') MessageC
      ENDIF
#endif
#ifdef DEBUG
      IF(INDEX(MeshFile,"NULL").GT.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The name of the file "        &
     &            //"containing the unstructured mesh must be "        &
     &            //"specified in the input file.")
      ENDIF
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') MeshFileOrigin
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The origin of the file containing the unstructured "&
     &            //"mesh could not be passed to the tracking cores.")
         WRITE(MeshFileOrigin,'(A)') MessageC
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(MeshFileOrigin)
#ifdef DEBUG
#ifndef NETCDF
      IF(INDEX(MeshFileOrigin,"HYCOM").GT.0)THEN
         CALL LPT_Print(0,"FATAL ERROR","The mesh file was obtained "  &
     &            //"from HYCOM and is presumed to be in NetCDF "      &
     &            //"format, but this code was not compiled with "     &
     &            //"NetCDF support.")
      ENDIF
#endif
#endif
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') MeshCoordinates
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The coordinate system for the unstructured mesh "   &
     &            //"could not be passed to the tracking cores.")
         WRITE(MeshCoordinates,'(A)') MessageC
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(MeshCoordinates)
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         WRITE(JunkC,'(A)') TRIM(MeshFile)
         CALL LPT_Print(0,"INFO","The unstructured mesh will be read " &
     &            //"from "//TRIM(ADJUSTL(JunkC))//".")
         WRITE(JunkC,'(A)') TRIM(MeshFileOrigin)
         CALL LPT_Print(0,"INFO","The unstructured mesh was obtained " &
     &            //"from "//TRIM(ADJUSTL(JunkC))//".")
         WRITE(JunkC,'(A)') TRIM(MeshCoordinates)
         CALL LPT_Print(0,"INFO","The unstructured mesh is given in "  &
     &            //TRIM(ADJUSTL(JunkC))//" coordinates.")
      ENDIF
#endif

      ! Read the information about the velocity combination.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=VelocityCombination,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR","The information "     &
     &               //"about the velocity combination could not be "  &
     &               //"read from the input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') CombinationMethod
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The method for the velocity combination could not " &
     &            //"be passed to the tracking cores.")
         WRITE(CombinationMethod,'(A)') MessageC
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1))
         MessageR(1) = PercentageCurrent
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The currents percentage to be used in the velocity "&
     &            //"combination could not be passed to the tracking " &
     &            //"cores.")
         PercentageCurrent = MessageR(1)
         MessageR(1) = PercentageWind
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The winds percentage to be used in the velocity "   &
     &            //"combination could not be passed to the tracking " &
     &            //"cores.")
         PercentageWind = MessageR(1)
         MessageR(1) = AngleWind
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The wind adjustment angle to be used in the "       &
     &            //"velocity combination could not be passed to the " &
     &            //"tracking cores.")
         AngleWind = MessageR(1)
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(CombinationMethod)
#ifdef DEBUG
      IF(INDEX(WindFile,"NULL").GT.0.AND.PercentageWind.GT.0.D0)THEN
         CALL LPT_Print(MyRank,"FATAL ERROR","To include a non-zero "  &
     &            //"percentage of the wind velocities, you must "     &
     &            //"specify the information about the wind file.")
      ENDIF
#endif
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         WRITE(JunkC,'(A)') TRIM(CombinationMethod)
         CALL LPT_Print(0,"INFO","The velocities will be combined by " &
     &            //"using the '"//TRIM(ADJUSTL(JunkC))//"' method.")
         WRITE(JunkC,'(F15.6)') PercentageCurrent
         CALL LPT_Print(0,"INFO","The combination will include a "     &
     &            //"current fraction of "//TRIM(ADJUSTL(JunkC))//".")
         WRITE(JunkC,'(F15.6)') PercentageWind
         CALL LPT_Print(0,"INFO","The combination will include a "     &
     &            //"wind fraction of "//TRIM(ADJUSTL(JunkC))//".")
         WRITE(JunkC,'(F15.6)') AngleWind
         CALL LPT_Print(0,"INFO","The combination will include a "     &
     &            //"wind adjustment angle of "//TRIM(ADJUSTL(JunkC))  &
     &            //".")
      ENDIF
#endif

      ! Read the information about the velocity combination.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=Diffusion,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(0,"WARNING","The information "              &
     &               //"about the horizontal diffusion could not be "  &
     &               //"read from the input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') DiffusionMethod
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The method for the diffusion could not be passed "  &
     &            //"to the tracking cores.")
         WRITE(DiffusionMethod,'(A)') MessageC
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(2))
         MessageR(1) = Cx
         MessageR(2) = Cy
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The diffusion scaling coefficients could not be "   &
     &            //"passed to the tracking cores.")
         Cx = MessageR(1)
         Cy = MessageR(2)
         MessageR(1) = Evx
         MessageR(2) = Evy
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The diffusion turbulence coefficients could not "   &
     &            //"be passed to the tracking cores.")
         Evx = MessageR(1)
         Evy = MessageR(2)
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(DiffusionMethod)
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         IF(INDEX(DiffusionMethod,"NULL").LE.0)THEN
            CALL LPT_Print(0,"INFO","The horizontal diffusion will "   &
     &               //"use the "//TRIM(ADJUSTL(DiffusionMethod))      &
     &               //" method.")
            WRITE(JunkC,'(F15.6)') Cx
            WRITE(JunkC1,'(F15.6)') Cy
            CALL LPT_Print(0,"INFO","The horizontal diffusion will "   &
     &               //"use scaling coefficients of "                  &
     &               //TRIM(ADJUSTL(JunkC))//" and "                   &
     &               //TRIM(ADJUSTL(JunkC1))//".")
            WRITE(JunkC,'(F15.6)') Evx
            WRITE(JunkC1,'(F15.6)') Evy
            CALL LPT_Print(0,"INFO","The horizontal diffusion will "   &
     &               //"use turbulence coefficients of "               &
     &               //TRIM(ADJUSTL(JunkC))//" and "                   &
     &               //TRIM(ADJUSTL(JunkC1))//".")
         ELSE
            CALL LPT_Print(0,"INFO","The horizontal diffusion will "   &
     &               //"not be considered.")
         ENDIF
      ENDIF
#endif

      ! Read the information about the water properties.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=WaterProperties,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"WARNING","The information "         &
     &               //"about the water properties could not be read " &
     &               //"from the input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(2))
         MessageR(1) = DensityWater
         MessageR(2) = DynamicViscosityWater
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The water properties could not be passed to the "   &
     &            //"tracking cores.")
         DensityWater = MessageR(1)
         DynamicViscosityWater = MessageR(2)
      ENDIF
#endif
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         WRITE(JunkC,'(F15.6)') DensityWater
         CALL LPT_Print(0,"INFO","The water density is set to "        &
     &            //TRIM(ADJUSTL(JunkC))//" kg/m3.")
         WRITE(JunkC,'(F15.6)') DynamicViscosityWater
         CALL LPT_Print(0,"INFO","The water dynamic viscosity is "     &
     &            //TRIM(ADJUSTL(JunkC))//" Pa-s.")
      ENDIF
#endif

      ! Read the information about the oil physics.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=OilPhysics,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"WARNING","The information "         &
     &               //"about the oil physics could not be read "      &
     &               //"from the input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') BuoyancyMethod
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The method for the oil buoyancy could not be "      &
     &            //"passed to the tracking cores.")
         WRITE(BuoyancyMethod,'(A)') MessageC
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') BuoyancyMethod
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The method for the oil buoyancy could not be "      &
     &            //"passed to the writer core.")
         WRITE(BuoyancyMethod,'(A)') MessageC
      ENDIF
#endif
      CALL LPT_Read_Capitalize_Word(BuoyancyMethod)
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         CALL LPT_Print(0,"INFO","The oil buoyancy will use the "      &
     &            //"method from "//TRIM(ADJUSTL(BuoyancyMethod))//".")
      ENDIF
#endif

#ifdef MPI
      ! The dedicated writer core does not need any other information
      ! from the input parameters file.
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF
#endif

      ! Read the information about the oil properties.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         REWIND(11)
         READ(11,NML=OilProperties,IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"WARNING","The information "         &
     &               //"about the oil properties could not be read "   &
     &               //"from the input file.")
         ENDIF
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(3))
         MessageR(1) = DensityOil
         MessageR(2) = DiameterEffective
         MessageR(3) = InterfacialTension
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The oil properties could not be passed to the "     &
     &            //"tracking cores.")
         DensityOil = MessageR(1)
         DiameterEffective = MessageR(2)
         InterfacialTension = MessageR(3)
      ENDIF
#endif
#if VERBOSE > 0
      IF(MyRank.EQ.0)THEN
         IF(INDEX(BuoyancyMethod,"NULL").LE.0)THEN
            WRITE(JunkC,'(F15.6)') DensityOil
            CALL LPT_Print(0,"INFO","The oil density is set to "       &
     &               //TRIM(ADJUSTL(JunkC))//" kg/m3.")
            WRITE(JunkC,'(F15.6)') DiameterEffective
            CALL LPT_Print(0,"INFO","The oil particle diameters are "  &
     &              //TRIM(ADJUSTL(JunkC))//" m.")
            WRITE(JunkC,'(F15.6)') InterfacialTension
            CALL LPT_Print(0,"INFO","The oil interfacial tension is "  &
     &               //TRIM(ADJUSTL(JunkC))//" N/m.")
         ENDIF
      ENDIF
#endif
 
      ! Close the input file.
      IF(MyRank.EQ.0)THEN
         CLOSE(UNIT=11,STATUS='KEEP',IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The input file could not " &
     &            //"be closed.")
         ENDIF
#endif
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The input file was closed.")
#endif
      ENDIF
      
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the LPT_Read_Parameters "       &
     &               //"routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_Initial_Conditions

      USE LPT_Comm_Module
      USE LPT_Data_Module

      IMPLICIT NONE

      CHARACTER(LEN=100)  :: JunkC
      CHARACTER(LEN=100)  :: JunkC1

      INTEGER             :: IC
      INTEGER             :: IERR
      INTEGER             :: IP
      INTEGER             :: IT
      INTEGER             :: JunkI
      INTEGER,ALLOCATABLE :: ParticlesOnLocalCore(:)

      REAL(8),ALLOCATABLE :: PartDepthGlobal(:)
      REAL(8),ALLOCATABLE :: PartLatGlobal(:)
      REAL(8),ALLOCATABLE :: PartLonGlobal(:)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the "                          &
     &               //"LPT_Read_Initial_Conditions routine.")
#endif

      ! Open the input file.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         OPEN(UNIT=12,FILE=TRIM(ParticleFile),ACTION='READ',         &
     &        IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR","The input file "      &
     &            //"could not be opened for reading.")
         ENDIF
#endif
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The input file was opened.")
#endif
      ENDIF

      ! Read the alpha-numeric header.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         CALL LPT_Read_From_File(12,(/"C"/),                           &
     &            "The header could not be read from the file "        &
     &            //"containing the initial particle locations.")
         WRITE(HeaderParticles,'(A)') TRIM(InputC)
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         WRITE(MessageC,'(A)') TRIM(HeaderParticles)
         CALL LPT_Comm_Distribute(MessageC,DummyI,DummyR,"C",          &
     &            "The header could not be passed to the tracker "     &
     &            //"cores.")
         WRITE(HeaderParticles,'(A)') TRIM(MessageC)
      ENDIF
      IF(UseWriterCore)THEN
         WRITE(MessageC,'(A)') TRIM(HeaderParticles)
         CALL LPT_Comm_To_Writer(MessageC,DummyI,DummyR,"C",           &
     &            "The header could not be passed to the writer core.")
         WRITE(HeaderParticles,'(A)') TRIM(MessageC)
      ENDIF
      ! The dedicated writer core does not need any other information
      ! from the file containing the initial locations of the particles.
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF
#endif

      ! Read the total number of initial particles.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
!        IF(SIZE(InputI).LT.1)THEN
            IF(ALLOCATED(InputI)) DEALLOCATE(InputI)
            ALLOCATE(InputI(1:1))
!        ENDIF
         CALL LPT_Read_From_File(12,(/"I"/),                           &
     &            "The total number of particles could not be read.")
         NumParticlesGlobal = InputI(1)
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = NumParticlesGlobal
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The total number of particles could not be passed " &
     &            //"to the cores.")
         NumParticlesGlobal = MessageI(1)
      ENDIF
#endif
#if VERBOSE > 0
      WRITE(JunkC,'(I24)') NumParticlesGlobal
      CALL LPT_Print(0,"INFO","There is an initial total of "          &
     &         //TRIM(ADJUSTL(JunkC))//" particles.")
#endif

      ! Determine the number of particles to be tracked by each core.
      NumParticlesLocal = INT( NumParticlesGlobal / NumTrackerCores )
      IF(ALLOCATED(ParticlesOnLocalCore))THEN
         DEALLOCATE(ParticlesOnLocalCore)
      ENDIF
      ALLOCATE(ParticlesOnLocalCore(1:NumTrackerCores))
      ParticlesOnLocalCore = NumParticlesLocal
      JunkI = NumParticlesGlobal - SUM(ParticlesOnLocalCore)
      DO IT=1,JunkI
         ParticlesOnLocalCore(IT) = ParticlesOnLocalCore(IT) + 1
      ENDDO
      DO IT=1,NumTrackerCores
         IF(MyRank.EQ.TrackerRanks(IT))THEN
            NumParticlesLocal = ParticlesOnLocalCore(IT)
         ENDIF
      ENDDO
#if VERBOSE > 0
      WRITE(JunkC,'(I24)') NumParticlesLocal
      CALL LPT_Print(0,"INFO","Each tracker core will begin with "     &
     &         //"about "//TRIM(ADJUSTL(JunkC))//" particles.")
#endif

      ! Read the global set of initial particle locations.
      IF(ALLOCATED(PartLonGlobal)) DEALLOCATE(PartLonGlobal)
      ALLOCATE(PartLonGlobal(1:NumParticlesGlobal))
      PartLonGlobal = 0.D0
      IF(ALLOCATED(PartLatGlobal)) DEALLOCATE(PartLatGlobal)
      ALLOCATE(PartLatGlobal(1:NumParticlesGlobal))
      PartLatGlobal = 0.D0
      IF(ALLOCATED(PartDepthGlobal)) DEALLOCATE(PartDepthGlobal)
      ALLOCATE(PartDepthGlobal(1:NumParticlesGlobal))
      PartDepthGlobal = 0.D0
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
#if VERBOSE > 2
         IF(MyRank.EQ.0)THEN
            WRITE(*,'(A,$)') "LPT: INFO: Reading the initial "         &
     &            //"particle locations: +"
         ENDIF
#endif
!        IF(SIZE(InputR).LT.3)THEN
            IF(ALLOCATED(InputR)) DEALLOCATE(InputR)
            ALLOCATE(InputR(1:3))
!        ENDIF
         DO IP=1,NumParticlesGlobal
#if VERBOSE > 3
            WRITE(JunkC,'(I24)') IP
            CALL LPT_Read_From_File(12,(/"R","R","R"/),                &
     &               "The initial location for particle number "       &
     &               //TRIM(ADJUSTL(JunkC))//" could not be read.")
#else
            READ(12,*) InputR(1),InputR(2),InputR(3)
#endif
            PartLonGlobal(IP)   = InputR(1)
            PartLatGlobal(IP)   = InputR(2)
            PartDepthGlobal(IP) = InputR(3)
#if VERBOSE > 2
            IF(MyRank.EQ.0)THEN
               CALL LPT_Progress(IP,NumParticlesGlobal)
            ENDIF
#endif
         ENDDO
      ENDIF
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The initial particle locations were "   &
     &         //"read successfully.")
#endif

      ! Close the file containing the initial particle locations.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         CLOSE(UNIT=12,STATUS='KEEP',IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The file containing the "  &
     &               //"initial particle locations could not "         &
     &               //"be closed.")
         ENDIF
#endif
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The file containing the initial "    &
     &            //"particle locations was closed.")
#endif
      ENDIF

#ifdef MPI
      ! If we used a reader core, then the global information
      ! is known only on that core.  Now send it to the tracking cores.
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1:NumParticlesGlobal))
         MessageR = PartLonGlobal
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The longitudes for the global set of particles "    &
     &            //"could not be passed to the tracker cores.")
         PartLonGlobal = MessageR
         MessageR = PartLatGlobal
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The latitudes for the global set of particles "     &
     &            //"could not be passed to the tracker cores.")
         PartLatGlobal = MessageR
         MessageR = PartDepthGlobal
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The depths for the global set of particles "        &
     &            //"could not be passed to the tracker cores.")
         PartDepthGlobal = MessageR
#if VERBOSE > 0
         CALL LPT_Print(0,"INFO","The initial locations for the "      &
     &            //"global set of particles were sent to the "        &
     &            //"tracker cores.")
#endif
      ENDIF
#endif

      ! Now that the tracker cores know the global information,
      ! they can pull easily the local information.
      DO IT=1,NumTrackerCores
         IF(MyRank.EQ.TrackerRanks(IT))THEN
            IF(ALLOCATED(PartLonLocal)) DEALLOCATE(PartLonLocal)
            ALLOCATE(PartLonLocal(1:NumParticlesLocal))
            IF(ALLOCATED(PartLatLocal)) DEALLOCATE(PartLatLocal)
            ALLOCATE(PartLatLocal(1:NumParticlesLocal))
            IF(ALLOCATED(PartDepthLocal)) DEALLOCATE(PartDepthLocal)
            ALLOCATE(PartDepthLocal(1:NumParticlesLocal))
            IF(ALLOCATED(PartNumberLocal)) DEALLOCATE(PartNumberLocal)
            ALLOCATE(PartNumberLocal(1:NumParticlesLocal))
            IF(INDEX(BuoyancyMethod,"NULL").LE.0)THEN
               IF(ALLOCATED(PartDiameterLocal))THEN
                  DEALLOCATE(PartDiameterLocal)
               ENDIF
               ALLOCATE(PartDiameterLocal(1:NumParticlesLocal))
            ENDIF
            IF(IT-1.EQ.0)THEN
               JunkI = 0
            ELSE
               JunkI = SUM(ParticlesOnLocalCore(1:IT-1))
            ENDIF
            DO IP=1,NumParticlesLocal
               PartLonLocal(IP) = PartLonGlobal(IP+JunkI)
               PartLatLocal(IP) = PartLatGlobal(IP+JunkI)
               PartDepthLocal(IP) = PartDepthGlobal(IP+JunkI)
               PartNumberLocal(IP) = IP+JunkI
               IF(INDEX(BuoyancyMethod,"NULL").LE.0)THEN
                  PartDiameterLocal(IP) = DiameterEffective
               ENDIF
            ENDDO
#if VERBOSE > 1
            WRITE(JunkC,'(I24)') MyRank
            WRITE(JunkC1,'(I24)') NumParticlesLocal
            CALL LPT_Print(MyRank,"INFO","Tracker core "               &
     &               //TRIM(ADJUSTL(JunkC))                            &
     &               //" determined the initial locations for its "    &
     &               //TRIM(ADJUSTL(JunkC1))//" particles.")
#endif
         ENDIF
      ENDDO
      IF(UseReaderCore.AND.(MyRank.EQ.0))THEN
         NumParticlesLocal = 0
      ENDIF

      ! Finally, convert the coordinates if necessary.
      IF(INDEX(ParticleInputCoordinates,"POLAR").GT.0)THEN
         DO IP=1,NumParticlesLocal
            PartLonLocal(IP) = PartLonLocal(IP) * Deg2Rad
            PartLatLocal(IP) = PartLatLocal(IP) * Deg2Rad
            PartLonLocal(IP) = PartLonLocal(IP) * EarthRadius
            PartLatLocal(IP) = PartLatLocal(IP) * EarthRadius
         ENDDO
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The particle locations were "        &
     &            //"converted to Cartesian coordinates.")
#endif
      ENDIF

      ! Clear up the used memory.
      IF(ALLOCATED(ParticlesOnLocalCore))THEN
         DEALLOCATE(ParticlesOnLocalCore)
      ENDIF
      IF(ALLOCATED(PartDepthGlobal)) DEALLOCATE(PartDepthGlobal)
      IF(ALLOCATED(PartLatGlobal))   DEALLOCATE(PartLatGlobal)
      IF(ALLOCATED(PartLonGlobal))   DEALLOCATE(PartLonGlobal)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the "                           &
     &         //"LPT_Read_Initial_Conditions routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_From_Particle_File(Snap)

      USE LPT_Comm_Module
      USE LPT_Data_Module
#ifdef NETCDF
      USE LPT_NetCDF_Module
#endif

      IMPLICIT NONE

      INTEGER,INTENT(IN)  :: Snap

      CHARACTER(LEN=100)  :: JunkC
      CHARACTER(LEN=100)  :: JunkC1

      INTEGER             :: IMOD
      INTEGER             :: IP
      INTEGER             :: IT
      INTEGER             :: JunkI
      INTEGER             :: NumParticlesGlobalTemp
      INTEGER,ALLOCATABLE :: ParticlesOnLocalCore(:)

      LOGICAL,SAVE        :: FirstCall = .TRUE.

      REAL(8),ALLOCATABLE :: PartDepthGlobal(:)
      REAL(8),ALLOCATABLE :: PartDiameterGlobal(:)
      REAL(8),ALLOCATABLE :: PartLatGlobal(:)
      REAL(8),ALLOCATABLE :: PartLonGlobal(:)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the "                          &
     &               //"LPT_Read_From_Particle_File routine.")
#endif

#ifdef MPI
      ! The dedicated writer core does not need any other information
      ! from the file containing the locations of the particles.
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF
#endif

      IF(INDEX(ParticleFileFormat,"ASCII").GT.0)THEN

      ! This section will be added later.
      CONTINUE

      ELSEIF(INDEX(ParticleFileFormat,"NETCDF").GT.0)THEN

      ! Only open the NetCDF file on the first call to this routine.
      ! Otherwise, use the existing identifier to access the file.
      IF(.NOT.FirstCall)THEN

         ! Pull the existing NetCDF identifier.
#ifdef NETCDF
         NC_ID = NC_ID_Particle
#endif

      ELSEIF(FirstCall)THEN

         ! Open the file containing the particle locations.
         IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
#ifdef NETCDF
            CALL LPT_NetCDF_Check( NF90_OPEN(TRIM(ParticleFile),       &
     &               NF90_NOWRITE,NC_ID) )
#endif
#if VERBOSE > 1
            CALL LPT_Print(0,"INFO","The particle file was opened.")
#endif
         ENDIF

      ENDIF

      ! Read the total number of particles.
      NumParticlesGlobalTemp = 0
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
#ifdef NETCDF
!        IF(SIZE(InputI).LT.1)THEN
            IF(ALLOCATED(InputI)) DEALLOCATE(InputI)
            ALLOCATE(InputI(1:1))
!        ENDIF
         CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,                  &
     &            'Particles_per_time_snap',NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,InputI(:),   &
     &            START=(/Snap/),COUNT=(/1/)) )
         NumParticlesGlobalTemp = InputI(1)
#endif
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = NumParticlesGlobalTemp
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The total number of particles could not be passed " &
     &            //"to the cores.")
         NumParticlesGlobalTemp = MessageI(1)
      ENDIF
#endif
#if VERBOSE > 0
      WRITE(JunkC,'(I24)') NumParticlesGlobalTemp
      CALL LPT_Print(0,"INFO","There is now a total of "               &
     &         //TRIM(ADJUSTL(JunkC))//" particles.")
#endif

      ! Read the global set of particle locations.
      IF(ALLOCATED(PartLonGlobal)) DEALLOCATE(PartLonGlobal)
      ALLOCATE(PartLonGlobal(1:NumParticlesGlobalTemp))
      PartLonGlobal = 0.D0
      IF(ALLOCATED(PartLatGlobal)) DEALLOCATE(PartLatGlobal)
      ALLOCATE(PartLatGlobal(1:NumParticlesGlobalTemp))
      PartLatGlobal = 0.D0
      IF(ALLOCATED(PartDepthGlobal)) DEALLOCATE(PartDepthGlobal)
      ALLOCATE(PartDepthGlobal(1:NumParticlesGlobalTemp))
      PartDepthGlobal = 0.D0
      IF(ALLOCATED(PartDiameterGlobal)) DEALLOCATE(PartDiameterGlobal)
      ALLOCATE(PartDiameterGlobal(1:NumParticlesGlobalTemp))
      PartDiameterGlobal = 0.D0
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
#ifdef NETCDF
         CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'X_particle',     &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,             &
     &            PartLonGlobal(:),START=(/1,Snap/),                   &
     &            COUNT=(/NumParticlesGlobalTemp,1/)) )
         CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'Y_particle',     &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,             &
     &            PartLatGlobal(:),START=(/1,Snap/),                   &
     &            COUNT=(/NumParticlesGlobalTemp,1/)) )
         CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'Z_particle',     &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,             &
     &            PartDepthGlobal(:),START=(/1,Snap/),                 &
     &            COUNT=(/NumParticlesGlobalTemp,1/)) )
         IF(INDEX(BuoyancyMethod,"NULL").LE.0)THEN
            CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,               &
     &               'Particle_diam',NC_VAR) )
            CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,          &
     &               PartDiameterGlobal(:),START=(/1,Snap/),           &
     &               COUNT=(/NumParticlesGlobalTemp,1/)) )
         ENDIF
#endif
      ENDIF
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The initial particle locations were "   &
     &         //"read successfully.")
#endif

      ! Close the ASCII/NetCDF conditional.
      ENDIF

#ifdef MPI
      ! If we used a reader core, then the global information
      ! is known only on that core.  Now send it to the tracking cores.
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(1:NumParticlesGlobalTemp))
         MessageR = PartLonGlobal
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The longitudes for the global set of particles "    &
     &            //"could not be passed to the tracker cores.")
         PartLonGlobal = MessageR
         MessageR = PartLatGlobal
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The latitudes for the global set of particles "     &
     &            //"could not be passed to the tracker cores.")
         PartLatGlobal = MessageR
         MessageR = PartDepthGlobal
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The depths for the global set of particles "        &
     &            //"could not be passed to the tracker cores.")
         PartDepthGlobal = MessageR
         IF(INDEX(BuoyancyMethod,"NULL").LE.0)THEN
            MessageR = PartDiameterGlobal
            CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",       &
     &               "The depths for the global set of particles "     &
     &               //"could not be passed to the tracker cores.")
            PartDiameterGlobal = MessageR
         ENDIF
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The initial locations for the "      &
     &            //"global set of particles were sent to the "        &
     &            //"tracker cores.")
#endif
      ENDIF
#endif

      ! If this is the first time snap, then we can use the standard
      ! partition of particles.
      IF(Snap.EQ.1)THEN

      NumParticlesGlobal = NumParticlesGlobalTemp

      ! Determine the number of particles to be tracked by each core.
      NumParticlesLocal = INT( NumParticlesGlobal / NumTrackerCores )
      IF(ALLOCATED(ParticlesOnLocalCore))THEN
         DEALLOCATE(ParticlesOnLocalCore)
      ENDIF
      ALLOCATE(ParticlesOnLocalCore(1:NumTrackerCores))
      ParticlesOnLocalCore = NumParticlesLocal
      JunkI = NumParticlesGlobal - SUM(ParticlesOnLocalCore)
      DO IT=1,JunkI
         ParticlesOnLocalCore(IT) = ParticlesOnLocalCore(IT) + 1
      ENDDO
#if VERBOSE > 0
      WRITE(JunkC,'(I24)') NumParticlesLocal
      CALL LPT_Print(0,"INFO","Each tracker core will begin with "     &
     &         //"about "//TRIM(ADJUSTL(JunkC))//" particles.")
#endif
      NumParticlesLocal = 0
 
      ! Now that the tracker cores know the global information,
      ! they can pull easily the local information.
      DO IT=1,NumTrackerCores
         IF(MyRank.EQ.TrackerRanks(IT))THEN
            IF(IT-1.EQ.0)THEN
               JunkI = 0
            ELSE
               JunkI = SUM(ParticlesOnLocalCore(1:IT-1))
            ENDIF
            DO IP=1,ParticlesOnLocalCore(IT)
               CALL LPT_Oil_Source(PartLonGlobal(IP+JunkI),            &
     &                  PartLatGlobal(IP+JunkI),                       &
     &                  PartDepthGlobal(IP+JunkI),                     &
     &                  PartDiameterGlobal(IP+JunkI),                  &
     &                  IP+JunkI,1)
            ENDDO
#if VERBOSE > 1
            WRITE(JunkC,'(I24)') MyRank
            WRITE(JunkC1,'(I24)') NumParticlesLocal
            CALL LPT_Print(MyRank,"INFO","Tracker core "               &
     &               //TRIM(ADJUSTL(JunkC))                            &
     &               //" determined the initial locations for its "    &
     &               //TRIM(ADJUSTL(JunkC1))//" particles.")
#endif
         ENDIF
      ENDDO

      ! Otherwise, if the global number of particles has not changed
      ! since the last time snap, then we can also use the existing
      ! partition of particles.  But we must not over-write the locations
      ! of particles that we are already tracking.
      ELSEIF(NumParticlesGlobal.EQ.NumParticlesGlobalTemp)THEN

      ! Determine the old breakdown of the number of particles
      ! to be tracked by each core.
      NumParticlesLocal = INT( NumParticlesGlobal / NumTrackerCores )
      IF(ALLOCATED(ParticlesOnLocalCore))THEN
         DEALLOCATE(ParticlesOnLocalCore)
      ENDIF
      ALLOCATE(ParticlesOnLocalCore(1:NumTrackerCores))
      ParticlesOnLocalCore = NumParticlesLocal
      JunkI = NumParticlesGlobal - SUM(ParticlesOnLocalCore)
      DO IT=1,JunkI
         ParticlesOnLocalCore(IT) = ParticlesOnLocalCore(IT) + 1
      ENDDO 
      DO IT=1,NumTrackerCores
         IF(MyRank.EQ.TrackerRanks(IT))THEN
            NumParticlesLocal = ParticlesOnLocalCore(IT)
         ENDIF    
      ENDDO       
#if VERBOSE > 0   
      WRITE(JunkC,'(I24)') NumParticlesLocal
      CALL LPT_Print(0,"INFO","Each tracker core will begin with "     &
     &         //"about "//TRIM(ADJUSTL(JunkC))//" particles.")
#endif               
      IF(UseReaderCore.AND.(MyRank.EQ.0))THEN
         NumParticlesLocal = 0
      ENDIF

      ! Now that the tracker cores know the global information,
      ! they can pull easily the local information.
      DO IT=1,NumTrackerCores
         IF(MyRank.EQ.TrackerRanks(IT))THEN
            IF(IT-1.EQ.0)THEN
               JunkI = 0
            ELSE
               JunkI = SUM(ParticlesOnLocalCore(1:IT-1))
            ENDIF
            DO IP=1,NumParticlesLocal
               IF(PartDomainLocal(IP).EQ.0)THEN
                  PartLonLocal(IP) = PartLonGlobal(IP+JunkI)
                  PartLatLocal(IP) = PartLatGlobal(IP+JunkI)
                  PartDepthLocal(IP) = PartDepthGlobal(IP+JunkI)
                  PartDiameterLocal(IP) = PartDiameterGlobal(IP+JunkI)
                  PartNumberLocal(IP) = IP+JunkI
                  IF(INDEX(ParticleInputCoordinates,"POLAR").GT.0)THEN
                     PartLonLocal(IP) = PartLonLocal(IP) * Deg2Rad
                     PartLatLocal(IP) = PartLatLocal(IP) * Deg2Rad
                     PartLonLocal(IP) = PartLonLocal(IP) * EarthRadius
                     PartLatLocal(IP) = PartLatLocal(IP) * EarthRadius
                  ENDIF
                  CALL LPT_Drog_Find_Element(PartLonLocal(IP),         &
     &                     PartLatLocal(IP),PartDepthLocal(IP),        &
     &                     PartElemJ(IP),PartElemL(IP),                &
     &                     PartDomainLocal(IP))
               ENDIF
            ENDDO
#if VERBOSE > 1
            WRITE(JunkC,'(I24)') MyRank
            WRITE(JunkC1,'(I24)') NumParticlesLocal
            CALL LPT_Print(MyRank,"INFO","Tracker core "               &
     &               //TRIM(ADJUSTL(JunkC))                            &
     &               //" determined the locations for its "            &
     &               //TRIM(ADJUSTL(JunkC1))//" particles.")
#endif
         ENDIF
      ENDDO

      ! Otherwise, if the global number of particles has changed/increased,
      ! then we must use the existing partition for the old set of particles,
      ! and then distribute the new set of particles.
      ELSEIF(NumParticlesGlobal.LT.NumParticlesGlobalTemp)THEN

      ! Determine the old breakdown of the number of particles
      ! to be tracked by each core.
      NumParticlesLocal = INT( NumParticlesGlobal / NumTrackerCores )
      IF(ALLOCATED(ParticlesOnLocalCore))THEN
         DEALLOCATE(ParticlesOnLocalCore)
      ENDIF
      ALLOCATE(ParticlesOnLocalCore(1:NumTrackerCores))
      ParticlesOnLocalCore = NumParticlesLocal
      JunkI = NumParticlesGlobal - SUM(ParticlesOnLocalCore)
      DO IT=1,JunkI
         ParticlesOnLocalCore(IT) = ParticlesOnLocalCore(IT) + 1
      ENDDO
      DO IT=1,NumTrackerCores
         IF(MyRank.EQ.TrackerRanks(IT))THEN
            NumParticlesLocal = ParticlesOnLocalCore(IT)
         ENDIF
      ENDDO
      IF(UseReaderCore.AND.(MyRank.EQ.0))THEN
         NumParticlesLocal = 0
      ENDIF

      ! Start by partitioning the old set of particles, but we must not
      ! over-write the locations of any particles that we are already tracking.
      DO IT=1,NumTrackerCores
         IF(MyRank.EQ.TrackerRanks(IT))THEN
            IF(IT-1.EQ.0)THEN
               JunkI = 0
            ELSE
               JunkI = SUM(ParticlesOnLocalCore(1:IT-1))
            ENDIF
            DO IP=1,NumParticlesLocal
               IF(PartDomainLocal(IP).EQ.0)THEN
                  PartLonLocal(IP) = PartLonGlobal(IP+JunkI)
                  PartLatLocal(IP) = PartLatGlobal(IP+JunkI)
                  PartDepthLocal(IP) = PartDepthGlobal(IP+JunkI)
                  PartDiameterLocal(IP) = PartDiameterGlobal(IP+JunkI)
                  PartNumberLocal(IP) = IP+JunkI
                  IF(INDEX(ParticleInputCoordinates,"POLAR").GT.0)THEN
                     PartLonLocal(IP) = PartLonLocal(IP) * Deg2Rad
                     PartLatLocal(IP) = PartLatLocal(IP) * Deg2Rad
                     PartLonLocal(IP) = PartLonLocal(IP) * EarthRadius
                     PartLatLocal(IP) = PartLatLocal(IP) * EarthRadius
                  ENDIF
                  CALL LPT_Drog_Find_Element(PartLonLocal(IP),         &
     &                     PartLatLocal(IP),PartDepthLocal(IP),        &
     &                     PartElemJ(IP),PartElemL(IP),                &
     &                     PartDomainLocal(IP))
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DO IP=NumParticlesGlobal+1,NumParticlesGlobalTemp
         IMOD = MOD(IP-NumParticlesGlobal,NumTrackerCores)
         IF(IMOD.EQ.0) IMOD = NumTrackerCores
         DO IT=1,NumTrackerCores
            IF(MyRank.EQ.TrackerRanks(IT).AND.IMOD.EQ.IT)THEN
               CALL LPT_Oil_Source(PartLonGlobal(IP),PartLatGlobal(IP),&
     &                  PartDepthGlobal(IP),PartDiameterGlobal(IP),IP,1)
            ENDIF
         ENDDO
      ENDDO
      NumParticlesGlobal = NumParticlesGlobalTemp

      ! Otherwise, if the global number of particles has changed/decreased,
      ! then we're screwed.
      ELSEIF(NumParticlesGlobal.GT.NumParticlesGlobalTemp)THEN

#if VERBOSE > 1
      CALL LPT_Print(MyRank,"FATAL ERROR","The global number of "      &
     &         //"particles has decreased.")
#endif

      ! End the first snap conditional.
      ENDIF

      ! Clear up the used memory.
      IF(ALLOCATED(ParticlesOnLocalCore))THEN
         DEALLOCATE(ParticlesOnLocalCore)
      ENDIF
      IF(ALLOCATED(PartDepthGlobal))    DEALLOCATE(PartDepthGlobal)
      IF(ALLOCATED(PartDiameterGlobal)) DEALLOCATE(PartDiameterGlobal)
      IF(ALLOCATED(PartLatGlobal))      DEALLOCATE(PartLatGlobal)
      IF(ALLOCATED(PartLonGlobal))      DEALLOCATE(PartLonGlobal)

      ! We won't need to open the NetCDF file when we return
      ! to this routine.
      IF(FirstCall)THEN
         FirstCall = .FALSE.
      ENDIF
#ifdef NETCDF
      NC_ID_Particle = NC_ID
#endif

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the "                           &
     &         //"LPT_Read_From_Particle_File routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_Unstructured_Mesh

      USE LPT_Comm_Module
      USE LPT_Data_Module
#ifdef KDTREE
      USE LPT_KDTREE_Module, ONLY: LPT_KDTREE_Initialize
#endif

      IMPLICIT NONE

      CHARACTER(LEN=100) :: JunkC
      CHARACTER(LEN=100) :: JunkC1

      INTEGER            :: IE
      INTEGER            :: IERR
      INTEGER            :: IV

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the "                          &
     &         //"LPT_Read_Unstructured_Mesh routine.")
#endif

#ifdef MPI
      ! The dedicated writer core does not need any information
      ! about the unstructured mesh.
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF
#endif

      ! Open the unstructured mesh file.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         OPEN(UNIT=14,FILE=TRIM(MeshFile),ACTION='READ',IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR","The unstructured "    &
     &            //"mesh file could not be opened for reading.")
         ENDIF
#endif
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The unstructured mesh file was "     &
     &            //"opened.")
#endif
      ENDIF

      ! Read the alpha-numeric header.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         READ(14,'(A)',IOSTAT=IERR) JunkC
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The header could not be "  &
     &               //"read from the unstructured mesh file.")
         ENDIF
#endif
      ENDIF

      ! Read the numbers of elements and vertices.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
!        IF(SIZE(InputI).LT.2)THEN
            IF(ALLOCATED(InputI)) DEALLOCATE(InputI)
            ALLOCATE(InputI(1:2))
!        ENDIF
         CALL LPT_Read_From_File(14,(/"I","I"/),                       &
     &            "The numbers of elements and vertices could not be " &
     &            //"read from the unstructured mesh.")
         NumElems = InputI(1)
         NumVerts = InputI(2)
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = NumElems
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The total number of mesh elements could not be "    &
     &            //"passed to the tracker cores.")
         NumElems = MessageI(1)
         MessageI(1) = NumVerts
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The total number of mesh vertices could not be "    &
     &            //"passed to the tracker cores.")
         NumVerts = MessageI(1)
      ENDIF
#endif
#if VERBOSE > 0
      WRITE(JunkC,'(I24)') NumElems
      WRITE(JunkC1,'(I24)') NumVerts
      CALL LPT_Print(0,"INFO","The unstructured mesh contains "        &
     &               //TRIM(ADJUSTL(JunkC))//" elements and "          &
     &               //TRIM(ADJUSTL(JunkC1))//" vertices.")
#endif

      ! Read the vertex information of longitudes, latitudes and bathy.
      IF(ALLOCATED(MeshLon)) DEALLOCATE(MeshLon)
      ALLOCATE(MeshLon(1:NumVerts))
      IF(ALLOCATED(MeshLat)) DEALLOCATE(MeshLat)
      ALLOCATE(MeshLat(1:NumVerts))
      IF(ALLOCATED(MeshDepth)) DEALLOCATE(MeshDepth)
      ALLOCATE(MeshDepth(1:NumVerts))
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
         WRITE(*,'(A,$)') "LPT: INFO: Reading the vertex table from "  &
     &         //"the mesh file: +"
      ENDIF
#endif
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
!        IF(SIZE(InputI).LT.1)THEN
            IF(ALLOCATED(InputI)) DEALLOCATE(InputI)
            ALLOCATE(InputI(1:1))
!        ENDIF
!        IF(SIZE(InputR).LT.4)THEN
            IF(ALLOCATED(InputR)) DEALLOCATE(InputR)
            ALLOCATE(InputR(1:4))
!        ENDIF
         DO IV=1,NumVerts
#if VERBOSE > 3
            WRITE(JunkC,'(I24)') IV
            CALL LPT_Read_From_File(14,(/"I","R","R","R"/),            &
     &               "The mesh information could not be read for "     &
     &                //"vertex "//TRIM(ADJUSTL(JunkC))//".")
#else
            READ(14,*) InputI(1),InputR(2),InputR(3),InputR(4)
#endif
            MeshLon(IV)   =         InputR(2)
            MeshLat(IV)   =         InputR(3)
            MeshDepth(IV) = -1.D0 * InputR(4)
#if VERBOSE > 2
            IF(MyRank.EQ.0)THEN
               CALL LPT_Progress(IV,NumVerts)
            ENDIF
#endif
         ENDDO
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(NumVerts))
         MessageR = MeshLon
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The mesh longitudes could not be passed to the "    &
     &            //"tracker cores.")
         MeshLon = MessageR
         MessageR = MeshLat
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The mesh latitudes could not be passed to the "     &
     &            //"tracker cores.")
         MeshLat = MessageR
         MessageR = MeshDepth
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The mesh depths could not be passed to the "        &
     &            //"tracker cores.")
         MeshDepth = MessageR
      ENDIF
#endif
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The vertex-based information was read " &
     &         //"from the mesh file.")
#endif

      ! Read the unstructured mesh connectivity table.
      IF(ALLOCATED(MeshConnEV)) DEALLOCATE(MeshConnEV)
      ALLOCATE(MeshConnEV(NumElems,3))
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
         WRITE(*,'(A,$)') "LPT: INFO: Reading the element table from " &
     &         //"the mesh file: +"
      ENDIF
#endif
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
!        IF(SIZE(InputI).LT.5)THEN
            IF(ALLOCATED(InputI)) DEALLOCATE(InputI)
            ALLOCATE(InputI(1:5))
!        ENDIF
         DO IE=1,NumElems
#if VERBOSE > 3
            WRITE(JunkC,'(I24)') IE
            CALL LPT_Read_From_File(14,(/"I","I","I","I","I"/),        &
     &               "The mesh information could not be read for "     &
     &                //"element "//TRIM(ADJUSTL(JunkC))//".")
#else
            READ(14,*) InputI(1),InputI(2),InputI(3),InputI(4),InputI(5)
#endif
            MeshConnEV(IE,1) = InputI(3)
            MeshConnEV(IE,2) = InputI(4)
            MeshConnEV(IE,3) = InputI(5)
#if VERBOSE > 2
            IF(MyRank.EQ.0)THEN
               CALL LPT_Progress(IE,NumElems)
            ENDIF
#endif
         ENDDO
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(NumElems))
         DO IV=1,3
            MessageI(:) = MeshConnEV(:,IV)
            CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",       &
     &               "The connectivity table could not be passed to "  &
     &               //"the tracker cores.")
            MeshConnEV(:,IV) = MessageI(:)
         ENDDO
      ENDIF
#endif
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The element-based information was read "&
     &         //"from the mesh file.")
#endif

      ! Close the unstructured mesh file.
      IF(MyRank.EQ.0)THEN
         CLOSE(UNIT=14,STATUS='KEEP',IOSTAT=IERR)
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(0,"FATAL ERROR","The unstructured mesh "    &
     &               //"file could not be closed.")
         ENDIF
#endif
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The unstructured mesh file was "     &
     &            //"closed.")
#endif
      ENDIF

      ! Finally, convert the coordinates if necessary.
      IF(INDEX(MeshCoordinates,"POLAR").GT.0)THEN
         DO IV=1,NumVerts
            MeshLon(IV) = MeshLon(IV) * Deg2Rad
            MeshLat(IV) = MeshLat(IV) * Deg2Rad
            MeshLon(IV) = MeshLon(IV) * EarthRadius
            MeshLat(IV) = MeshLat(IV) * EarthRadius
         ENDDO
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The unstructured mesh was "          &
     &            //"converted to Cartesian coordinates.")
#endif
      ENDIF

      IF(AmTrackerCore)THEN
         CALL LPT_Drog_MAK_NEINFO
      ENDIF

#ifdef KDTREE
      IF(AmTrackerCore)THEN
         CALL LPT_KDTREE_Initialize
      ENDIF
#endif

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the "                           &
     &         //"LPT_Read_Unstructured_Mesh routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_HYCOM_Mesh

      USE LPT_Comm_Module
      USE LPT_Data_Module
#ifdef KDTREE
      USE LPT_KDTREE_Module, ONLY: LPT_KDTREE_Initialize
#endif
#ifdef NETCDF
      USE LPT_NetCDF_Module
#endif

      IMPLICIT NONE

      CHARACTER(LEN=100) :: JunkC
      CHARACTER(LEN=100) :: JunkC1

      INTEGER :: IDX
      INTEGER :: IE
      INTEGER :: ILAT
      INTEGER :: ILON
      INTEGER :: IV
      INTEGER :: NumDep
      INTEGER :: NumLat
      INTEGER :: NumLon

      REAL(8),ALLOCATABLE :: TempDep(:)
      REAL(8),ALLOCATABLE :: TempLat(:)
      REAL(8),ALLOCATABLE :: TempLon(:)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the "                          &
     &         //"LPT_Read_HYCOM_Mesh routine.")
#endif

#ifdef MPI
      ! The dedicated writer core does not need any information
      ! about the unstructured mesh.
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF
#endif

      ! Open the file containing the HYCOM information.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
#ifdef NETCDF
         CALL LPT_NetCDF_Check( NF90_OPEN(TRIM(MeshFile),              &
     &            NF90_NOWRITE,NC_ID) )
#endif
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The HYCOM mesh file was opened.")
#endif
      ENDIF

      ! Read the dimensions of the HYCOM mesh.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
#ifdef NETCDF
         CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,'Longitude',      &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION(NC_ID,          &
     &            NC_VAR,LEN=NumLon) )
         CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,'Latitude',       &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION(NC_ID,          &
     &            NC_VAR,LEN=NumLat) )
         CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,'Depth',          &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION(NC_ID,          &
     &            NC_VAR,LEN=NumDep) )
#endif
      ENDIF

      ! Determine the numbers of vertices and elements.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         NumVerts = NumLon * NumLat
         NumElems = 2 * (NumLon-1) * (NumLat-1)
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(1))
         MessageI(1) = NumElems
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The total number of mesh elements could not be "    &
     &            //"passed to the tracker cores.")
         NumElems = MessageI(1)
         MessageI(1) = NumVerts
         CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",          &
     &            "The total number of mesh vertices could not be "    &
     &            //"passed to the tracker cores.")
         NumVerts = MessageI(1)
      ENDIF
#endif
#if VERBOSE > 0
      WRITE(JunkC,'(I24)') NumElems
      WRITE(JunkC1,'(I24)') NumVerts
      CALL LPT_Print(0,"INFO","The unstructured mesh contains "        &
     &               //TRIM(ADJUSTL(JunkC))//" elements and "          &
     &               //TRIM(ADJUSTL(JunkC1))//" vertices.")
#endif

      ! Develop the arrays of longitudes, latitudes and depths.
      IF(ALLOCATED(MeshLon)) DEALLOCATE(MeshLon)
      ALLOCATE(MeshLon(1:NumVerts))
      IF(ALLOCATED(MeshLat)) DEALLOCATE(MeshLat)
      ALLOCATE(MeshLat(1:NumVerts))
      IF(ALLOCATED(MeshDepth)) DEALLOCATE(MeshDepth)
      ALLOCATE(MeshDepth(1:NumVerts))
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         IF(ALLOCATED(TempLon)) DEALLOCATE(TempLon)
         ALLOCATE(TempLon(1:NumLon))
#ifdef NETCDF
         CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'Longitude',      &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,TempLon(:),  &
     &            START=(/1/),COUNT=(/NumLon/)) )
#endif
         IF(ALLOCATED(TempLat)) DEALLOCATE(TempLat)
         ALLOCATE(TempLat(1:NumLat))
#ifdef NETCDF
         CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'Latitude',       &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,TempLat(:),  &
     &            START=(/1/),COUNT=(/NumLat/)) )
#endif
         IF(ALLOCATED(TempDep)) DEALLOCATE(TempDep)
         ALLOCATE(TempDep(1:NumDep))
#ifdef NETCDF
         CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'Depth',          &
     &            NC_VAR) )
         CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,TempDep(:),  &
     &            START=(/1/),COUNT=(/NumDep/)) )
#endif
         TempDep(:) = -1.D0 * TempDep(:)
         IV = 0
         DO ILON=1,NumLon
            DO ILAT=1,NumLat
               IV = IV + 1
               MeshLon(IV)   = TempLon(ILON)
               MeshLat(IV)   = TempLat(ILAT)
               MeshDepth(IV) = TempDep(NumDep)
            ENDDO
         ENDDO
         IF(ALLOCATED(TempLon)) DEALLOCATE(TempLon)
         IF(ALLOCATED(TempLat)) DEALLOCATE(TempLat)
         IF(ALLOCATED(TempDep)) DEALLOCATE(TempDep)
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(NumVerts))
         MessageR = MeshLon
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The mesh longitudes could not be passed to the "    &
     &            //"tracker cores.")
         MeshLon = MessageR
         MessageR = MeshLat
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The mesh latitudes could not be passed to the "     &
     &            //"tracker cores.")
         MeshLat = MessageR
         MessageR = MeshDepth
         CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",          &
     &            "The mesh depths could not be passed to the "        &
     &            //"tracker cores.")
         MeshDepth = MessageR
      ENDIF
#endif
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The vertex-based information was read " &
     &         //"from the mesh file.")
#endif

      ! Develop the element-to-vertex connectivity table.
      IF(ALLOCATED(MeshConnEV)) DEALLOCATE(MeshConnEV)
      ALLOCATE(MeshConnEV(NumElems,3))
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
         IE = 0
         DO ILON=1,NumLon-1
            DO ILAT=1,NumLat-1
               ! Lower-left triangle in the rectangle.
               IE = IE + 1
               IDX = (ILON-1)*NumLat + (ILAT  )
               MeshConnEV(IE,1) = IDX
               IDX = (ILON  )*NumLat + (ILAT  )
               MeshConnEV(IE,2) = IDX
               IDX = (ILON-1)*NumLat + (ILAT+1)
               MeshConnEV(IE,3) = IDX
               ! Upper-right triangle in the rectangle.
               IE = IE + 1
               IDX = (ILON  )*NumLat + (ILAT  )
               MeshConnEV(IE,1) = IDX
               IDX = (ILON  )*NumLat + (ILAT+1)
               MeshConnEV(IE,2) = IDX
               IDX = (ILON-1)*NumLat + (ILAT+1)
               MeshConnEV(IE,3) = IDX
            ENDDO
         ENDDO
      ENDIF
#ifdef MPI
      IF(UseReaderCore)THEN
         IF(ALLOCATED(MessageI)) DEALLOCATE(MessageI)
         ALLOCATE(MessageI(NumElems))
         DO IV=1,3
            MessageI(:) = MeshConnEV(:,IV)
            CALL LPT_Comm_Distribute(DummyC,MessageI,DummyR,"I",       &
     &               "The connectivity tables could not be passed to " &
     &               //"the tracker cores.")
            MeshConnEV(:,IV) = MessageI(:)
         ENDDO
      ENDIF
#endif
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The element-based information was read "&
     &         //"from the mesh file.")
#endif

      ! Close the unstructured mesh file.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN
#ifdef NETCDF
         CALL LPT_NetCDF_Check( NF90_CLOSE(NC_ID) )
#endif
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The HYCOM mesh file was closed.")
#endif
      ENDIF

      ! Finally, convert the coordinates if necessary.
      IF(INDEX(MeshCoordinates,"POLAR").GT.0)THEN
         DO IV=1,NumVerts
            MeshLon(IV) = MeshLon(IV) * Deg2Rad
            MeshLat(IV) = MeshLat(IV) * Deg2Rad
            MeshLon(IV) = MeshLon(IV) * EarthRadius
            MeshLat(IV) = MeshLat(IV) * EarthRadius
         ENDDO
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The HYCOM mesh was converted to "    &
     &            //"Cartesian coordinates.")
#endif
      ENDIF

      IF(AmTrackerCore)THEN
         CALL LPT_Drog_MAK_NEINFO
      ENDIF

#ifdef KDTREE
      IF(AmTrackerCore)THEN
         CALL LPT_KDTREE_Initialize
      ENDIF
#endif

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the "                           &
     &         //"LPT_Read_HYCOM_Mesh routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_Update_Velocities(Snap)

      USE LPT_Comm_Module
      USE LPT_Data_Module

      IMPLICIT NONE

      INTEGER,INTENT(IN)       :: Snap

      CHARACTER(LEN=100)       :: JunkC
      CHARACTER(LEN=100)       :: JunkC1
      CHARACTER(LEN=100)       :: JunkC2

      INTEGER                  :: IERR
      INTEGER                  :: IL
      INTEGER                  :: IV

      LOGICAL,SAVE             :: FirstCall = .TRUE.

      REAL(8)                  :: CurrMag
      REAL(8),ALLOCATABLE      :: CurrU(:)
      REAL(8),ALLOCATABLE      :: CurrV(:)
      REAL(8),ALLOCATABLE      :: CurrW(:)
      REAL(8)                  :: SnapTimes(2)
      REAL(8)                  :: SnapUDiff
      REAL(8)                  :: SnapVDiff
      REAL(8)                  :: TimeDiff
      REAL(8)                  :: TimeInc
!     REAL(8),SAVE             :: Times(2)
      REAL(8)                  :: WindU(1)
      REAL(8)                  :: WindV(1)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the "                          &
     &         //"LPT_Read_Update_Velocities routine.")
#endif

#ifdef MPI
      ! The dedicated writer core does not need any information
      ! about the velocity field.
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF
#endif

      ! The goal of this subroutine is to fill arrays
      ! with velocities that are (a) some combination of currents
      ! and winds that are (b) interpolated to the time associated
      ! with the particle tracking step.

      ! Read the currents and winds (if necessary) and interpolate
      ! them to the correct time.
      IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN

         ! CURRENTS

         SnapTimes(:) = Times(:)

         ! Initialize the arrays containing the snaps read from
         ! the input files.
         IF(FirstCall)THEN

            SnapUCurr(:,:,1) = SnapUCurr(:,:,2)
            SnapVCurr(:,:,1) = SnapVCurr(:,:,2)
            SnapWCurr(:,:,1) = SnapWCurr(:,:,2)

            ! Read the second snap.
            CALL LPT_Read_Velocity_Field("C",SnapTimes(2))

         ENDIF

         ! If the time associated with the current particle
         ! tracking step has advanced past the second snap,
         ! then update the second snap.
         DO WHILE((Time1.GE.SnapTimes(2)).AND.                         &
     &            (Time1.LT.StartingTime+(NumTrackingSnaps)*           &
     &                         SimulationLength/NumTrackingSnaps))

            ! Shift the currents backward.
            SnapTimes(1) = SnapTimes(2)
            SnapUCurr(:,:,1) = SnapUCurr(:,:,2)
            SnapVCurr(:,:,1) = SnapVCurr(:,:,2)
            SnapWCurr(:,:,1) = SnapWCurr(:,:,2)

            ! Read the next snap.
            CALL LPT_Read_Velocity_Field("C",SnapTimes(2))

         ENDDO

         ! WINDS

         IF((INDEX(WindFile,"NULL").LE.0).OR.                          &
     &      (INDEX(CombinationMethod,"Original").GT.0.AND.             &
     &       PercentageWind.GT.0.D0))THEN

            SnapTimes(:) = Times(:)

            ! Initialize the arrays containing the snaps read from
            ! the input files.
            IF(FirstCall)THEN

               ! Shift the winds backward.
               SnapUWind(:,:,1) = SnapUWind(:,:,2)
               SnapVWind(:,:,1) = SnapVWind(:,:,2)

               ! Read the second snap.
               CALL LPT_Read_Velocity_Field("W",SnapTimes(2))

            ENDIF

            ! If the time associated with the current particle
            ! tracking step has advanced past the second snap,
            ! then update the second snap.
            DO WHILE((Time1.GE.SnapTimes(2)).AND.                      &
     &               (Time1.LT.StartingTime+(NumTrackingSnaps)*        &
     &                            SimulationLength/NumTrackingSnaps))

               ! Shift the winds backward.
               SnapTimes(1) = SnapTimes(2)
               SnapUWind(:,:,1) = SnapUWind(:,:,2)
               SnapVWind(:,:,1) = SnapVWind(:,:,2)

               ! Read the next snap.
               CALL LPT_Read_Velocity_Field("W",SnapTimes(2))

            ENDDO

         ELSE

            SnapUWind = 0.D0
            SnapVWind = 0.D0

         ENDIF

         ! Develop arrays of current and wind velocities at the time
         ! associated with the current particle tracking step,
         ! by interpolating in time.
#if VERBOSE > 0
         WRITE(JunkC,'(F15.6)') Time1/3600.D0
         CALL LPT_Print(0,"INFO","The particle tracking time is "      &
     &            //TRIM(ADJUSTL(JunkC))//" hours.")
         WRITE(JunkC,'(F15.6)') SnapTimes(1)/3600.D0
         WRITE(JunkC1,'(F15.6)') SnapTimes(2)/3600.D0
         CALL LPT_Print(0,"INFO","The velocities will be "             &
     &            //"interpolated between "//TRIM(ADJUSTL(JunkC))      &
     &            //" and "//TRIM(ADJUSTL(JunkC1))//" hours.")
#endif

         ! Combine the currents and winds into a single velocity
         ! for use in the rest of the code.

         IF(ALLOCATED(CurrU)) DEALLOCATE(CurrU)
         ALLOCATE(CurrU(NumLayers))
         IF(ALLOCATED(CurrV)) DEALLOCATE(CurrV)
         ALLOCATE(CurrV(NumLayers))
         IF(ALLOCATED(CurrW)) DEALLOCATE(CurrW)
         ALLOCATE(CurrW(NumLayers))

         IF(ALLOCATED(VelU)) DEALLOCATE(VelU)
         ALLOCATE(VelU(NumVerts,NumLayers))
         IF(ALLOCATED(VelV)) DEALLOCATE(VelV)
         ALLOCATE(VelV(NumVerts,NumLayers))
         IF(ALLOCATED(VelW)) DEALLOCATE(VelW)
         ALLOCATE(VelW(NumVerts,NumLayers))

         TimeInc  = SnapTimes(2) - SnapTimes(1)
         TimeDiff = SnapTimes(2) - Time1
         Times(:) = SnapTimes(:)

         DO IV=1,NumVerts

            ! Interpolate the currents and winds in time.
            CurrU(:) = SnapUCurr(IV,:,2) - ((TimeDiff*                 &
     &                (SnapUCurr(IV,:,2)-SnapUCurr(IV,:,1)))/TimeInc)
            CurrV(:) = SnapVCurr(IV,:,2) - ((TimeDiff*                 &
     &                (SnapVCurr(IV,:,2)-SnapVCurr(IV,:,1)))/TimeInc)
            CurrW(:) = SnapWCurr(IV,:,2) - ((TimeDiff*                 &
     &                (SnapWCurr(IV,:,2)-SnapWCurr(IV,:,1)))/TimeInc)
            WindU(:) = SnapUWind(IV,:,2) - ((TimeDiff*                 &
     &                (SnapUWind(IV,:,2)-SnapUWind(IV,:,1)))/TimeInc)
            WindV(:) = SnapVWind(IV,:,2) - ((TimeDiff*                 &
     &                (SnapVWind(IV,:,2)-SnapVWind(IV,:,1)))/TimeInc)

            ! Load the 3D currents into the single velocity field.
            VelU(IV,:) = PercentageCurrent * CurrU(:)
            VelV(IV,:) = PercentageCurrent * CurrV(:)
            VelW(IV,:) = PercentageCurrent * CurrW(:)

            IF(INDEX(CombinationMethod,"ORIGINAL").GT.0)THEN

               ! Adjust the horizontal velocities on the top layer
               ! for the contribution from the winds.
               VelU(IV,NumLayers)                                      &
     &                = PercentageCurrent * CurrU(NumLayers)           &
     &                + PercentageWind * WindU(1)                      &
     &                      * COS(AngleWind*Pi/180.D0)                 &
     &                - PercentageWind * WindV(1)                      &
     &                      * SIN(AngleWind*Pi/180.D0)
               VelV(IV,NumLayers)                                      &
     &                = PercentageCurrent * CurrV(NumLayers)           &
     &                + PercentageWind * WindU(1)                      &
     &                      * SIN(AngleWind*Pi/180.D0)                 &
     &                + PercentageWind * WindV(1)                      &
     &                      * COS(AngleWind*Pi/180.D0)

               ! If the current has zero magnitude, then do not allow
               ! the wind alone to move the particles.
               CurrMag = SQRT( CurrU(NumLayers)*CurrU(NumLayers)       &
     &                       + CurrV(NumLayers)*CurrV(NumLayers) )
               IF(CurrMag.LT.0.00001D0)THEN
                  VelU(IV,NumLayers) = 0.D0
                  VelV(IV,NumLayers) = 0.D0
               ENDIF

            ENDIF

         ENDDO

      ENDIF

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The currents and winds were combined "  &
     &         //"into a single velocity field.")
#endif

      ! Send the single velocity field to the tracking cores.
#ifdef MPI
      IF(UseReaderCore)THEN
         ! Send the arrays, layer by layer.
         IF(.NOT.ALLOCATED(VelU)) ALLOCATE(VelU(1:NumVerts,1:NumLayers))
         IF(.NOT.ALLOCATED(VelV)) ALLOCATE(VelV(1:NumVerts,1:NumLayers))
         IF(.NOT.ALLOCATED(VelW)) ALLOCATE(VelW(1:NumVerts,1:NumLayers))
         IF(ALLOCATED(MessageR)) DEALLOCATE(MessageR)
         ALLOCATE(MessageR(NumVerts))
         DO IL=1,NumLayers
            MessageR(:) = VelU(:,IL)
            CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",       &
     &               "The u-component of the velocity field could "    &
     &               //"not be passed to the tracking cores.")
            VelU(:,IL) = MessageR(:)
            MessageR(:) = VelV(:,IL)
            CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",       &
     &               "The v-component of the velocity field could "    &
     &               //"not be passed to the tracking cores.")
            VelV(:,IL) = MessageR(:)
            MessageR(:) = VelW(:,IL)
            CALL LPT_Comm_Distribute(DummyC,DummyI,MessageR,"R",       &
     &               "The w-component of the velocity field could "    &
     &               //"not be passed to the tracking cores.")
            VelW(:,IL) = MessageR(:)
         ENDDO
#if VERBOSE > 1
         CALL LPT_Print(0,"INFO","The velocity field was passed to "   &
     &            //"the tracking cores.")
#endif
      ENDIF
#endif

      IF(FirstCall)THEN
         FirstCall = .FALSE.
      ENDIF

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the "                           &
     &         //"LPT_Read_Update_Velocities routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_Velocity_Field(Field,TimeSnap)

      USE LPT_Comm_Module
      USE LPT_Data_Module
#ifdef NETCDF
      USE LPT_NetCDF_Module
#endif

      IMPLICIT NONE

      CHARACTER(LEN=1),INTENT(IN)  :: Field
      REAL(8),INTENT(OUT)          :: TimeSnap

      CHARACTER(LEN=100)           :: FieldDimensions
      CHARACTER(LEN=100)           :: FieldFile
      CHARACTER(LEN=100)           :: FieldFileError
      CHARACTER(LEN=100)           :: FieldFileFormat
      CHARACTER(LEN=100)           :: FieldFileOrigin
#ifdef NETCDF
      CHARACTER(LEN=100)           :: FieldNameU
      CHARACTER(LEN=100)           :: FieldNameV
      CHARACTER(LEN=100)           :: FieldNameW
#endif
      CHARACTER(LEN=100)           :: JunkC
      CHARACTER(LEN=1),ALLOCATABLE :: VariableTypes(:)

      INTEGER                      :: IERR
      INTEGER                      :: II
      INTEGER                      :: IL
      INTEGER                      :: ILAT
      INTEGER                      :: ILON
      INTEGER                      :: IV
      INTEGER                      :: JunkI
      INTEGER,SAVE                 :: NumLat
      INTEGER,SAVE                 :: NumLon
      INTEGER                      :: NumNonDefaultVerts
      INTEGER                      :: NumSnaps
      INTEGER                      :: UnitNumber

      LOGICAL                      :: FirstCall
      LOGICAL,SAVE                 :: FirstCallCurr = .TRUE.
      LOGICAL,SAVE                 :: FirstCallWind = .TRUE.

      REAL(8)                      :: DefaultValue
      REAL(8),TARGET               :: Dummy(1,1,2)
      REAL(8)                      :: FillValue
      REAL(8)                      :: JunkR(1)
      REAL(8),ALLOCATABLE          :: Temp1(:)
      REAL(8),ALLOCATABLE          :: Temp2(:)
      REAL(8),ALLOCATABLE          :: Temp3(:,:,:)
      REAL(8),ALLOCATABLE          :: Temp4(:,:,:)
      REAL(8),POINTER              :: U(:,:,:)
      REAL(8),POINTER              :: V(:,:,:)
      REAL(8),POINTER              :: W(:,:,:)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the "                          &
     &         //"LPT_Read_Velocity_Field routine.")
#endif

#ifdef MPI
      ! The dedicated writer core does not need any information
      ! about the velocity field.
      IF(UseWriterCore.AND.(MyRank+1.EQ.NumRanks))THEN
         RETURN
      ENDIF
#endif

#ifdef NETCDF
      IF(FirstCallCurr) NC_SnapCurr = 0
      IF(FirstCallWind) NC_SnapWind = 0
#endif

      IF(INDEX(Field,"C").GT.0)THEN
         WRITE(FieldFile,'(A)') TRIM(CurrentFile)
         WRITE(FieldFileOrigin,'(A)') TRIM(CurrentFileOrigin)
         WRITE(FieldFileFormat,'(A)') TRIM(CurrentFileFormat)
         WRITE(FieldDimensions,'(A)') TRIM(CurrentDimensions)
         WRITE(FieldFileError,'(A)') "current"
         FirstCall = FirstCallCurr
         IF(INDEX(FieldDimensions,"2D").GT.0)THEN
            UnitNumber = 64
         ELSEIF(INDEX(FieldDimensions,"3D").GT.0)THEN
            UnitNumber = 45
         ENDIF
#ifdef NETCDF
         IF(INDEX(FieldFileOrigin,"ADCIRC").GT.0)THEN
            IF(INDEX(FieldDimensions,"2D").GT.0)THEN
               WRITE(FieldNameU,'(A)') "u-vel"
               WRITE(FieldNameV,'(A)') "v-vel"
            ELSEIF(INDEX(FieldDimensions,"3D").GT.0)THEN
               WRITE(FieldNameU,'(A)') "u-vel3D"
               WRITE(FieldNameV,'(A)') "v-vel3D"
               WRITE(FieldNameW,'(A)') "w-vel3D"
            ENDIF
         ELSEIF(INDEX(FieldFileOrigin,"HYCOM").GT.0)THEN
            WRITE(FieldNameU,'(A)') "u"
            WRITE(FieldNameV,'(A)') "v"
            WRITE(FieldNameW,'(A)') "w_velocity"
         ENDIF
         NC_ID = NC_ID_Curr
         NC_Snap = NC_SnapCurr
#endif
      ELSEIF(INDEX(Field,"W").GT.0)THEN
         WRITE(FieldFile,'(A)') TRIM(WindFile)
         WRITE(FieldFileOrigin,'(A)') "ADCIRC"
         WRITE(FieldFileFormat,'(A)') TRIM(WindFileFormat)
         WRITE(FieldDimensions,'(A)') "2D"
         WRITE(FieldFileError,'(A)') "wind"
         FirstCall = FirstCallWind
         UnitNumber = 74
#ifdef NETCDF
         WRITE(FieldNameU,'(A)') "windx"
         WRITE(FieldNameV,'(A)') "windy"
         NC_ID = NC_ID_Wind
         NC_Snap = NC_SnapWind
#endif
      ENDIF

      ! On the first call to this routine, we need to open the files
      ! and process the header information.
      IF(FirstCall)THEN

         IF(INDEX(FieldFileOrigin,"ADCIRC").GT.0)THEN

            IF(INDEX(FieldFileFormat,"ASCII").GT.0)THEN

               ! Open the file containing the velocity field.
               IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.           &
     &            (MyRank.EQ.0))THEN
                  OPEN(UNIT=UnitNumber,FILE=TRIM(FieldFile),           &
     &                 ACTION='READ',IOSTAT=IERR)
#ifdef DEBUG
                  IF(IERR.NE.0)THEN
                     CALL LPT_Print(MyRank,"FATAL ERROR","The file "   &
     &                        //"containing the "                      &
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities could not be opened for " &
     &                        //"reading.")
                  ENDIF
#endif
#if VERBOSE > 1
                  CALL LPT_Print(0,"INFO","The file with the "         &
     &                     //TRIM(ADJUSTL(FieldFileError))             &
     &                     //" velocities was opened.")
#endif
               ENDIF

               ! Read the alphanumeric header.
               IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.           &
     &            (MyRank.EQ.0))THEN
                  CALL LPT_Read_From_File(UnitNumber,(/"C"/),          &
     &                     "The header could not be read from the "    &
     &                     //"file containing the "                    &
     &                     //TRIM(ADJUSTL(FieldFileError))             &
     &                     //" velocities.")
                  WRITE(JunkC,'(A)') TRIM(InputC)
               ENDIF

               ! Read the numbers of time snaps and vertices.
               IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.           &
     &            (MyRank.EQ.0))THEN
!                 IF(SIZE(InputI).LT.5)THEN
                     IF(ALLOCATED(InputI)) DEALLOCATE(InputI)
                     ALLOCATE(InputI(1:5))
!                 ENDIF
!                 IF(SIZE(InputR).LT.3)THEN
                     IF(ALLOCATED(InputR)) DEALLOCATE(InputR)
                     ALLOCATE(InputR(1:3))
!                 ENDIF
                  CALL LPT_Read_From_File(UnitNumber,                  &
     &                     (/"I","I","R","I","I"/),                    &
     &                     "The numbers of time snaps and vertices "   &
     &                     //"could not be read from the file "        &
     &                     //"containing the "                         &
     &                     //TRIM(ADJUSTL(FieldFileError))             &
     &                     //" velocities.")
#ifdef DEBUG
                  IF(NumTrackingSnaps.GT.InputI(1))THEN
                     CALL LPT_Print(MyRank,"FATAL ERROR","The file "   &
     &                        //"containing the "                      &
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities does not have enough "    &
     &                        //"time snaps.")
                  ENDIF
                  IF(NumVerts.NE.InputI(2))THEN
                     CALL LPT_Print(MyRank,"FATAL ERROR","The file "   &
     &                        //"containing the "                      &
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities does not have the "       &
     &                        //"correct number of vertices.")
                  ENDIF
#endif
                  NumSnaps = InputI(1)
                  NumVerts = InputI(2)
                  JunkR(1) = InputR(3)
                  JunkI    = InputI(4)
                  IF(INDEX(Field,"C").GT.0)THEN
                     IF(INDEX(FieldDimensions,"2D").GT.0)THEN
                        NumLayers = 1
                     ELSEIF(INDEX(FieldDimensions,"3D").GT.0)THEN
                        NumLayers = InputI(5)
                     ENDIF
                  ENDIF
               ENDIF
#if VERBOSE > 0
               WRITE(JunkC,'(I24)') NumSnaps
               CALL LPT_Print(0,"INFO","The file containing the "      &
     &                  //TRIM(ADJUSTL(FieldFileError))                &
     &                  //" velocities has "//TRIM(ADJUSTL(JunkC))     &
     &                  //" time snaps.")
               IF(INDEX(FieldDimensions,"3D").GT.0)THEN
                  WRITE(JunkC,'(I24)') NumLayers
                  CALL LPT_Print(0,"INFO","The file containing the "   &
     &                     //TRIM(ADJUSTL(FieldFileError))             &
     &                     //" velocities has "//TRIM(ADJUSTL(JunkC))  &
     &                     //" vertical layers.")
               ENDIF
#endif

            ELSEIF(INDEX(FieldFileFormat,"NETCDF").GT.0)THEN
#ifdef NETCDF

               IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.           &
     &            (MyRank.EQ.0))THEN

                  ! Open the file containing the velocity field.
                  CALL LPT_NetCDF_Check( NF90_OPEN(TRIM(FieldFile),    &
     &                     NF90_NOWRITE,NC_ID) )

                  ! Read the number of time snaps.
                  CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,'time',  &
     &                     NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION(NC_ID, &
     &                     NC_VAR,LEN=JunkI) )
#ifdef DEBUG
                  IF(NumTrackingSnaps.GT.JunkI)THEN
                     CALL LPT_Print(MyRank,"FATAL ERROR","The file "   &
     &                        //"containing the "                      &
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities does not have enough "    &
     &                        //"time snaps.")
                  ENDIF
#endif
                  NumSnaps = JunkI
#if VERBOSE > 0
                  WRITE(JunkC,'(I24)') NumSnaps
                  CALL LPT_Print(0,"INFO","The file containing the "   &
     &                     //TRIM(ADJUSTL(FieldFileError))             &
     &                     //" velocities has "//TRIM(ADJUSTL(JunkC))  &
     &                     //" time snaps.")
#endif

                  ! Read the number of vertices.
                  CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,'node',  &
     &                     NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION(NC_ID, &
     &                     NC_VAR,LEN=JunkI) )
#ifdef DEBUG
                  IF(NumVerts.NE.JunkI)THEN
                     CALL LPT_Print(MyRank,"FATAL ERROR","The file "   &
     &                        //"containing the "                      &
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities does not have the "       &
     &                        //"correct number of vertices.")
                  ENDIF
#endif
                  NumVerts = JunkI

                  ! Read the number of vertical layers.
                  IF(INDEX(Field,"C").GT.0)THEN
                     IF(INDEX(FieldDimensions,"2D").GT.0)THEN
                        NumLayers = 1
                     ELSEIF(INDEX(FieldDimensions,"3D").GT.0)THEN
                        CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,   &
     &                           'num_v_nodes',NC_VAR) )
                        CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION( &
     &                           NC_ID,NC_VAR,LEN=NumLayers) )
                     ENDIF
                  ENDIF
#if VERBOSE > 0
                  IF(INDEX(FieldDimensions,"3D").GT.0)THEN
                     WRITE(JunkC,'(I24)') NumLayers
                     CALL LPT_Print(0,"INFO","The file containing the "&
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities has "                     &
     &                        //TRIM(ADJUSTL(JunkC))                   &
     &                        //" vertical layers.")
                  ENDIF
#endif

               ENDIF

#endif
            ENDIF ! ASCII/NetCDF

         ELSEIF(INDEX(FieldFileOrigin,"HYCOM").GT.0)THEN

            IF(INDEX(FieldFileFormat,"NETCDF").GT.0)THEN
#ifdef NETCDF

               IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.           &
     &            (MyRank.EQ.0))THEN

                  ! Open the file containing the velocity field.
                  CALL LPT_NetCDF_Check( NF90_OPEN(TRIM(FieldFile),    &
     &                     NF90_NOWRITE,NC_ID) )

                  ! Read the number of time snaps.
                  CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,'MT',    &
     &                     NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION(NC_ID, &
     &                     NC_VAR,LEN=JunkI) )
#ifdef DEBUG
                  IF(NumTrackingSnaps.GT.JunkI)THEN
                     CALL LPT_Print(MyRank,"FATAL ERROR","The file "   &
     &                        //"containing the "                      &
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities does not have enough "    &
     &                        //"time snaps.")
                  ENDIF
#endif
                  NumSnaps = JunkI
#if VERBOSE > 0
                  WRITE(JunkC,'(I24)') NumSnaps
                  CALL LPT_Print(0,"INFO","The file containing the "   &
     &                     //TRIM(ADJUSTL(FieldFileError))             &
     &                     //" velocities has "//TRIM(ADJUSTL(JunkC))  &
     &                     //" time snaps.")
#endif

                  ! Read the number of vertices.
                  CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,         &
     &                     'Longitude',NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION(NC_ID, &
     &                     NC_VAR,LEN=NumLon) )
                  CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,         &
     &                     'Latitude',NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION(NC_ID, &
     &                     NC_VAR,LEN=NumLat) )
#ifdef DEBUG
                  IF(NumVerts.NE.NumLon*NumLat)THEN
                     CALL LPT_Print(MyRank,"FATAL ERROR","The file "   &
     &                        //"containing the "                      &
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities does not have the "       &
     &                        //"correct number of vertices.")
                  ENDIF
#endif
                  NumVerts = NumLon*NumLat

                  ! Read the number of vertical layers.
                  IF(INDEX(Field,"C").GT.0)THEN
                     IF(INDEX(FieldDimensions,"3D").GT.0)THEN
                        CALL LPT_NetCDF_Check( NF90_INQ_DIMID(NC_ID,   &
     &                           'Depth',NC_VAR) )
                        CALL LPT_NetCDF_Check( NF90_INQUIRE_DIMENSION( &
     &                           NC_ID,NC_VAR,LEN=NumLayers) )
                     ENDIF
                  ENDIF
#if VERBOSE > 0
                  IF(INDEX(FieldDimensions,"3D").GT.0)THEN
                     WRITE(JunkC,'(I24)') NumLayers
                     CALL LPT_Print(0,"INFO","The file containing the "&
     &                        //TRIM(ADJUSTL(FieldFileError))          &
     &                        //" velocities has "                     &
     &                        //TRIM(ADJUSTL(JunkC))                   &
     &                        //" vertical layers.")
                  ENDIF
#endif

               ENDIF

#endif
            ENDIF ! ASCII/NetCDF

         ENDIF ! Origin

      ENDIF ! FirstCall

      ! Now that we know the number of layers,
      ! we can allocate the following arrays.
      IF(.NOT.ALLOCATED(SnapUCurr))THEN
         ALLOCATE(SnapUCurr(1:NumVerts,1:NumLayers,1:2))
      ENDIF
      IF(.NOT.ALLOCATED(SnapVCurr))THEN
         ALLOCATE(SnapVCurr(1:NumVerts,1:NumLayers,1:2))
      ENDIF
      IF(.NOT.ALLOCATED(SnapWCurr))THEN
         ALLOCATE(SnapWCurr(1:NumVerts,1:NumLayers,1:2))
      ENDIF
      IF(.NOT.ALLOCATED(SnapUWind))THEN
         ALLOCATE(SnapUWind(1:NumVerts,1:1,1:2))
      ENDIF
      IF(.NOT.ALLOCATED(SnapVWind))THEN
         ALLOCATE(SnapVWind(1:NumVerts,1:1,1:2))
      ENDIF

      ! For ease of programming, we can use pointers
      ! to refer to the snap arrays in the following section.
      IF(INDEX(Field,"C").GT.0)THEN
         U => SnapUCurr
         V => SnapVCurr
         W => SnapWCurr
      ELSEIF(INDEX(Field,"W").GT.0)THEN
         U => SnapUWind
         V => SnapVWind
         W => Dummy
      ENDIF

      ! Read a snap from the velocity field.
      IF(INDEX(FieldFileOrigin,"ADCIRC").GT.0)THEN

         IF(INDEX(FieldFileFormat,"ASCII").GT.0)THEN

            ! Read the header information for the next time snap.
            IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN

               IF(INDEX(FieldDimensions,"2D").GT.0)THEN

                  CALL LPT_Read_From_File(UnitNumber,(/"C"/),          &
     &                     "The next time snap could not be read "     &
     &                     //"from the file containing the "           &
     &                     //TRIM(ADJUSTL(FieldFileError))             &
     &                     //" velocities.")
                  WRITE(JunkC,'(A)') TRIM(InputC)
                  READ(JunkC,*,END=333,ERR=333) TimeSnap,JunkI,        &
     &                 NumNonDefaultVerts,DefaultValue
                  GOTO 444
 333              CONTINUE
                  READ(JunkC,*) TimeSnap,JunkI
                  NumNonDefaultVerts = NumVerts
                  DefaultValue = 0.D0
 444              CONTINUE
#if VERBOSE > 1
                  WRITE(JunkC,'(I24)') NumNonDefaultVerts
                  CALL LPT_Print(0,"INFO","The next snap contains "    &
     &                     //TRIM(ADJUSTL(JunkC))//" vertices with "   &
     &                     //"non-zero "//TRIM(ADJUSTL(FieldFileError))&
     &                     //"s.")
#endif

               ELSEIF(INDEX(FieldDimensions,"3D").GT.0)THEN

                  IF(ALLOCATED(Sigma)) DEALLOCATE(Sigma)
                  ALLOCATE(Sigma(1:NumLayers))
!                 IF(SIZE(InputI).LT.2)THEN
                     IF(ALLOCATED(InputI)) DEALLOCATE(InputI)
                     ALLOCATE(InputI(1:2))
!                 ENDIF
!                 IF(SIZE(InputR).LT.3*NumLayers+1)THEN
                     IF(ALLOCATED(InputR)) DEALLOCATE(InputR)
                     ALLOCATE(InputR(1:3*NumLayers+1))
!                 ENDIF
#if VERBOSE > 3
                  IF(ALLOCATED(VariableTypes)) DEALLOCATE(VariableTypes)
                  ALLOCATE(VariableTypes(1:3*NumLayers+1))
                  WRITE(VariableTypes(1),'(A)') "R"
                  WRITE(VariableTypes(2),'(A)') "I"
                  DO II=3,3*NumLayers+1
                     WRITE(VariableTypes(II),'(A)') "R"
                  ENDDO
                  CALL LPT_Read_From_File(UnitNumber,VariableTypes,    &
     &                     "The next time snap could not be read "     &
     &                     //"from the file containing the "           &
     &                     //TRIM(ADJUSTL(FieldFileError))             &
     &                     //" velocities.")
#else
                  READ(UnitNumber,*) InputR(1),InputI(2),              &
     &                  (InputR(IL),InputR(IL+1),InputR(IL+2),         &
     &                     IL=3,3*NumLayers-1,3),                      &
     &                  InputR(3*NumLayers),InputR(3*NumLayers+1)
#endif
                  TimeSnap = InputR(1)
                  JunkI    = InputI(2)
                  DO IL=1,NumLayers
                     Sigma(IL) = InputR(3*IL)
                  ENDDO
                  NumNonDefaultVerts = NumVerts
                  DefaultValue = 0.D0

               ENDIF

            ENDIF

            ! Read the next time snap.
#if VERBOSE > 2
            IF(MyRank.EQ.0)THEN
               WRITE(*,'(A,$)') "LPT: INFO: Reading the next snap "    &
     &               //"from the "//TRIM(ADJUSTL(FieldFileError))      &
     &               //" velocity file: +"
            ENDIF
#endif
            IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN

               IF(INDEX(FieldDimensions,"2D").GT.0)THEN

                  U(:,:,2) = DefaultValue
                  V(:,:,2) = DefaultValue
                  W(:,:,2) = 0.D0
!                 IF(SIZE(InputI).LT.1)THEN
                     IF(ALLOCATED(InputI)) DEALLOCATE(InputI)
                     ALLOCATE(InputI(1:1))
!                 ENDIF
!                 IF(SIZE(InputR).LT.3)THEN
                     IF(ALLOCATED(InputR)) DEALLOCATE(InputR)
                     ALLOCATE(InputR(1:3))
!                 ENDIF
                  DO IV=1,NumNonDefaultVerts
#if VERBOSE > 3
                     WRITE(JunkC,'(I24)') IV
                     CALL LPT_Read_From_File(UnitNumber,               &
     &                        (/"I","R","R"/),                         &
     &                        "The "//TRIM(ADJUSTL(FieldFileError))    &
     &                        //" velocities could not be read for "   &
     &                        //"non-zero vertex "                     &
     &                        //TRIM(ADJUSTL(JunkC))//".")
#else
                     READ(UnitNumber,*) InputI(1),InputR(2),InputR(3)
#endif
                     U(InputI(1),1,2) = InputR(2)
                     V(InputI(1),1,2) = InputR(3)
#if VERBOSE > 2
                     IF(MyRank.EQ.0)THEN
                        CALL LPT_Progress(IV,NumNonDefaultVerts)
                     ENDIF
#endif
                  ENDDO

               ELSEIF(INDEX(FieldDimensions,"3D").GT.0)THEN

                  U(:,:,2) = DefaultValue
                  V(:,:,2) = DefaultValue
                  W(:,:,2) = DefaultValue
                  DO IV=1,NumNonDefaultVerts
                     READ(UnitNumber,*) JunkI,(U(JunkI,IL,2),          &
     &                                         V(JunkI,IL,2),          &
     &                                         W(JunkI,IL,2),          &
     &                                         IL=1,NumLayers)
#if VERBOSE > 2
                     IF(MyRank.EQ.0)THEN
                        CALL LPT_Progress(IV,NumNonDefaultVerts)
                     ENDIF
#endif
                  ENDDO

               ENDIF

            ENDIF

         ELSEIF(INDEX(FieldFileFormat,"NETCDF").GT.0)THEN
#ifdef NETCDF

            IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN

               U(:,:,2) = 0.D0
               V(:,:,2) = 0.D0
               W(:,:,2) = 0.D0

               ! Increment the snap for the next go-round.
               NC_Snap = NC_Snap + 1

               ! Read the time associated with this snap.
               CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'time',     &
     &                  NC_VAR) )
               CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,JunkR, &
     &                  START=(/NC_Snap/),COUNT=(/1/)) )
               TimeSnap = JunkR(1)

               IF(INDEX(FieldDimensions,"2D").GT.0)THEN

                  ! Read the x-component.
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,         &
     &                     FieldNameU,NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     U(:,1,2),START=(/1,NC_Snap/),               &
     &                     COUNT=(/NumVerts,1/)) )

                  ! Read the y-component.
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,         &
     &                     FieldNameV,NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     V(:,1,2),START=(/1,NC_Snap/),               &
     &                     COUNT=(/NumVerts,1/)) )

               ELSEIF(INDEX(FieldDimensions,"3D").GT.0)THEN

                  ! Read the sigma layers.
                  IF(ALLOCATED(Sigma)) DEALLOCATE(Sigma)
                  ALLOCATE(Sigma(1:NumLayers))
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'sigma', &
     &                     NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     Sigma(:),START=(/1/),COUNT=(/NumLayers/)) )

                  ! Read the x-component.
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,         &
     &                     FieldNameU,NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     U(:,:,2),START=(/1,1,NC_Snap/),             &
     &                     COUNT=(/NumVerts,NumLayers,1/)) )

                  ! Read the y-component.
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,         &
     &                     FieldNameV,NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     V(:,:,2),START=(/1,1,NC_Snap/),             &
     &                     COUNT=(/NumVerts,NumLayers,1/)) )

                  ! Read the z-component.
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,         &
     &                     FieldNameW,NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     W(:,:,2),START=(/1,1,NC_Snap/),             &
     &                     COUNT=(/NumVerts,NumLayers,1/)) )

               ENDIF

            ENDIF

#endif
         ENDIF

      ELSEIF(INDEX(FieldFileOrigin,"HYCOM").GT.0)THEN

         IF(INDEX(FieldFileFormat,"NETCDF").GT.0)THEN
#ifdef NETCDF

            IF((.NOT.UseReaderCore.AND.AmTrackerCore).OR.(MyRank.EQ.0))THEN

               U(:,:,2) = 0.D0
               V(:,:,2) = 0.D0
               W(:,:,2) = 0.D0

               ! Increment the snap for the next go-round.
               NC_Snap = NC_Snap + 1

               ! Read the time associated with this snap.
               CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'MT',       &
     &                  NC_VAR) )
               CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,JunkR, &
     &                  START=(/NC_Snap/),COUNT=(/1/)) )
               TimeSnap = DBLE(NINT(JunkR(1)*24.D0))*3600.D0

               IF(INDEX(FieldDimensions,"3D").GT.0)THEN

                  ! Read the depths and convert to sigma layers.
                  IF(ALLOCATED(Temp1)) DEALLOCATE(Temp1)
                  ALLOCATE(Temp1(1:NumLayers))
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,'Depth', &
     &                     NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     Temp1(:),START=(/1/),COUNT=(/NumLayers/)) )
                  Temp1(:) = -1.D0 * Temp1(:)
                  IF(ALLOCATED(Temp2)) DEALLOCATE(Temp2)
                  ALLOCATE(Temp2(1:NumLayers))
                  DO IL=1,NumLayers
                     Temp2(IL) = Temp1(NumLayers-IL+1)
                  ENDDO
                  IF(ALLOCATED(Sigma)) DEALLOCATE(Sigma)
                  ALLOCATE(Sigma(1:NumLayers))
                  DO IL=1,NumLayers
                     Sigma(IL) = -1.D0 + 2.D0*(Temp2(IL)-Temp2(1))     &
     &                                /(Temp2(NumLayers)-Temp2(1))
                  ENDDO
                  IF(ALLOCATED(Temp1)) DEALLOCATE(Temp1)
                  IF(ALLOCATED(Temp2)) DEALLOCATE(Temp2)

                  ! Read the x-component.
                  IF(ALLOCATED(Temp3)) DEALLOCATE(Temp3)
                  ALLOCATE(Temp3(1:NumLon,1:NumLat,1:NumLayers))
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,         &
     &                     FieldNameU,NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     Temp3(:,:,:),START=(/1,1,1,NC_Snap/),       &
     &                     COUNT=(/NumLon,NumLat,NumLayers,1/)) )
                  CALL LPT_NetCDF_Check( NF90_GET_ATT(NC_ID,NC_VAR,    &
     &                     "_FillValue",FillValue) )
                  IF(ALLOCATED(Temp4)) DEALLOCATE(Temp4)
                  ALLOCATE(Temp4(1:NumLon,1:NumLat,1:NumLayers))
                  DO IL=1,NumLayers
                     Temp4(:,:,IL) = Temp3(:,:,NumLayers-IL+1)
                  ENDDO
                  IV = 0
                  DO ILON=1,NumLon
                     DO ILAT=1,NumLat
                        IV = IV + 1
                        DO IL=1,NumLayers
                           IF(Temp4(ILON,ILAT,IL).EQ.FillValue)THEN
                              Temp4(ILON,ILAT,IL) = 0.D0
                           ENDIF
                           U(IV,IL,2) = Temp4(ILON,ILAT,IL)
                        ENDDO
                     ENDDO
                  ENDDO
                  IF(ALLOCATED(Temp3)) DEALLOCATE(Temp3)
                  IF(ALLOCATED(Temp4)) DEALLOCATE(Temp4)

                  ! Read the y-component.
                  IF(ALLOCATED(Temp3)) DEALLOCATE(Temp3)
                  ALLOCATE(Temp3(1:NumLon,1:NumLat,1:NumLayers))
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,         &
     &                     FieldNameV,NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     Temp3(:,:,:),START=(/1,1,1,NC_Snap/),       &
     &                     COUNT=(/NumLon,NumLat,NumLayers,1/)) )
                  CALL LPT_NetCDF_Check( NF90_GET_ATT(NC_ID,NC_VAR,    &
     &                     "_FillValue",FillValue) )
                  IF(ALLOCATED(Temp4)) DEALLOCATE(Temp4)
                  ALLOCATE(Temp4(1:NumLon,1:NumLat,1:NumLayers))
                  DO IL=1,NumLayers
                     Temp4(:,:,IL) = Temp3(:,:,NumLayers-IL+1)
                  ENDDO
                  IV = 0
                  DO ILON=1,NumLon
                     DO ILAT=1,NumLat
                        IV = IV + 1
                        DO IL=1,NumLayers
                           IF(Temp4(ILON,ILAT,IL).EQ.FillValue)THEN
                              Temp4(ILON,ILAT,IL) = 0.D0
                           ENDIF
                           V(IV,IL,2) = Temp4(ILON,ILAT,IL)
                        ENDDO
                     ENDDO
                  ENDDO
                  IF(ALLOCATED(Temp3)) DEALLOCATE(Temp3)
                  IF(ALLOCATED(Temp4)) DEALLOCATE(Temp4)
      
                  ! Read the z-component.
                  IF(ALLOCATED(Temp3)) DEALLOCATE(Temp3)
                  ALLOCATE(Temp3(1:NumLon,1:NumLat,1:NumLayers))
                  CALL LPT_NetCDF_Check( NF90_INQ_VARID(NC_ID,         &
     &                     FieldNameW,NC_VAR) )
                  CALL LPT_NetCDF_Check( NF90_GET_VAR(NC_ID,NC_VAR,    &
     &                     Temp3(:,:,:),START=(/1,1,1,NC_Snap/),       &
     &                     COUNT=(/NumLon,NumLat,NumLayers,1/)) )
                  CALL LPT_NetCDF_Check( NF90_GET_ATT(NC_ID,NC_VAR,    &
     &                     "_FillValue",FillValue) )
                  IF(ALLOCATED(Temp4)) DEALLOCATE(Temp4)
                  ALLOCATE(Temp4(1:NumLon,1:NumLat,1:NumLayers))
                  DO IL=1,NumLayers
                     Temp4(:,:,IL) = Temp3(:,:,NumLayers-IL+1)
                  ENDDO
                  IV = 0
                  DO ILON=1,NumLon
                     DO ILAT=1,NumLat
                        IV = IV + 1
                        DO IL=1,NumLayers
                           IF(Temp4(ILON,ILAT,IL).EQ.FillValue)THEN
                              Temp4(ILON,ILAT,IL) = 0.D0
                           ENDIF
                           W(IV,IL,2) = Temp4(ILON,ILAT,IL)
                        ENDDO
                     ENDDO
                  ENDDO
                  IF(ALLOCATED(Temp3)) DEALLOCATE(Temp3)
                  IF(ALLOCATED(Temp4)) DEALLOCATE(Temp4)

               ENDIF

            ENDIF

#endif
         ENDIF

      ENDIF
#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","The next snap of "                      &
     &         //TRIM(ADJUSTL(FieldFileError))                         &
     &         //" velocities was read successfully.")
#endif

      IF(INDEX(Field,"C").GT.0)THEN
         FirstCallCurr = .FALSE.
         NULLIFY(U)
         NULLIFY(V)
         NULLIFY(W)
#ifdef NETCDF
         NC_ID_Curr = NC_ID
         NC_SnapCurr = NC_Snap
#endif
      ELSEIF(INDEX(Field,"W").GT.0)THEN
         FirstCallWind = .FALSE.
         NULLIFY(U)
         NULLIFY(V)
         NULLIFY(W)
#ifdef NETCDF
         NC_ID_Wind = NC_ID
         NC_SnapWind = NC_Snap
#endif
      ENDIF

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the "                           &
     &         //"LPT_Read_Velocity_Field routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_From_File(UnitNumber,VariableTypes,          &
     &               ErrorMessage)

      USE LPT_Comm_Module

      IMPLICIT NONE

      CHARACTER(*),INTENT(IN) :: ErrorMessage
      CHARACTER(*),INTENT(IN) :: VariableTypes(:)

      INTEGER,INTENT(IN)      :: UnitNumber

      CHARACTER(LEN=500)      :: JunkC
      CHARACTER(LEN=500)      :: JunkC1

      INTEGER                 :: IERR
      INTEGER                 :: IV

      READ(UnitNumber,'(A)',IOSTAT=IERR) JunkC
#ifdef DEBUG
      IF(IERR.NE.0)THEN
         CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
      ENDIF
#endif

      DO IV=1,SIZE(VariableTypes)
         IF(INDEX(VariableTypes(IV),"C").EQ.1)THEN
            READ(JunkC,'(A)',IOSTAT=IERR) InputC
         ELSEIF(INDEX(VariableTypes(IV),"I").EQ.1)THEN
            READ(JunkC,*,IOSTAT=IERR) InputI(IV)
         ELSEIF(INDEX(VariableTypes(IV),"R").EQ.1)THEN
            READ(JunkC,*,IOSTAT=IERR) InputR(IV)
         ENDIF
#ifdef DEBUG
         IF(IERR.NE.0)THEN
            CALL LPT_Print(MyRank,"FATAL ERROR",ErrorMessage)
         ENDIF
#endif
         IF(IV.LT.SIZE(VariableTypes))THEN
            ! Remove any leading blanks.
            WRITE(JunkC,'(A)') TRIM(ADJUSTL(JunkC))
            ! Move up the remaining entries.
            WRITE(JunkC1,'(A)') JunkC(INDEX(JunkC," "):LEN_TRIM(JunkC))
            WRITE(JunkC,'(A)') TRIM(ADJUSTL(JunkC1))
         ENDIF
      ENDDO

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Read_Capitalize_Word(Word)

      IMPLICIT NONE

      CHARACTER(100),INTENT(INOUT) :: Word

      INTEGER                      :: I

      DO I=1,LEN_TRIM(Word)

         IF(Word(I:I).EQ."a") Word(I:I) = "A"
         IF(Word(I:I).EQ."b") Word(I:I) = "B"
         IF(Word(I:I).EQ."c") Word(I:I) = "C"
         IF(Word(I:I).EQ."d") Word(I:I) = "D"
         IF(Word(I:I).EQ."e") Word(I:I) = "E"
         IF(Word(I:I).EQ."f") Word(I:I) = "F"
         IF(Word(I:I).EQ."g") Word(I:I) = "G"
         IF(Word(I:I).EQ."h") Word(I:I) = "H"
         IF(Word(I:I).EQ."i") Word(I:I) = "I"
         IF(Word(I:I).EQ."j") Word(I:I) = "J"
         IF(Word(I:I).EQ."k") Word(I:I) = "K"
         IF(Word(I:I).EQ."l") Word(I:I) = "L"
         IF(Word(I:I).EQ."m") Word(I:I) = "M"
         IF(Word(I:I).EQ."n") Word(I:I) = "N"
         IF(Word(I:I).EQ."o") Word(I:I) = "O"
         IF(Word(I:I).EQ."p") Word(I:I) = "P"
         IF(Word(I:I).EQ."q") Word(I:I) = "Q"
         IF(Word(I:I).EQ."r") Word(I:I) = "R"
         IF(Word(I:I).EQ."s") Word(I:I) = "S"
         IF(Word(I:I).EQ."t") Word(I:I) = "T"
         IF(Word(I:I).EQ."u") Word(I:I) = "U"
         IF(Word(I:I).EQ."v") Word(I:I) = "V"
         IF(Word(I:I).EQ."w") Word(I:I) = "W"
         IF(Word(I:I).EQ."x") Word(I:I) = "X"
         IF(Word(I:I).EQ."y") Word(I:I) = "Y"
         IF(Word(I:I).EQ."z") Word(I:I) = "Z"

      ENDDO

      END SUBROUTINE



      END MODULE


