      MODULE LPT_Data_Module

      IMPLICIT NONE

      REAL(8),PARAMETER   :: EarthRadius = 6.3675D6
      REAL(8),PARAMETER   :: Gravity = 9.805D0

      CHARACTER(LEN=100)  :: HeaderParticles

      INTEGER,ALLOCATABLE :: BoundarySegmentElement(:)
      INTEGER,ALLOCATABLE :: BoundarySegmentFollowing(:)
      INTEGER,ALLOCATABLE :: BoundarySegmentPreceding(:)
      INTEGER,ALLOCATABLE :: BoundarySegmentVerts(:,:)
      INTEGER,ALLOCATABLE :: MapOpenMP(:)
      INTEGER,ALLOCATABLE :: MeshConnEE(:,:)
      INTEGER,ALLOCATABLE :: MeshConnEV(:,:)
      INTEGER,ALLOCATABLE :: MeshConnVE(:,:)
      INTEGER             :: NumBoundarySegments
      INTEGER             :: NumElems
      INTEGER             :: NumLayers
      INTEGER             :: NumParticlesGlobal
      INTEGER             :: NumParticlesLocal
      INTEGER             :: NumVerts = 0
      INTEGER             :: OutputTimeStep
      INTEGER,ALLOCATABLE :: PartBoundary(:)
      INTEGER,ALLOCATABLE :: PartBoundaryCode(:)
      INTEGER,ALLOCATABLE :: PartElemJ(:)
      INTEGER,ALLOCATABLE :: PartElemL(:)
      INTEGER,ALLOCATABLE :: PartDomainLocal(:)
      INTEGER,ALLOCATABLE :: PartNumberLocal(:)

      REAL(8)             :: Deg2Rad
      REAL(8),ALLOCATABLE :: ElementArea_AR(:)
      REAL(8),ALLOCATABLE :: ElementArea_T(:,:)
      REAL(8),ALLOCATABLE :: MeshDepth(:)
      REAL(8),ALLOCATABLE :: MeshLat(:)
      REAL(8),ALLOCATABLE :: MeshLon(:)
      REAL(8),ALLOCATABLE :: PartDepthLocal(:)
      REAL(8),ALLOCATABLE :: PartDiameterLocal(:)
      REAL(8),ALLOCATABLE :: PartLatLocal(:)
      REAL(8),ALLOCATABLE :: PartLonLocal(:)
      REAL(8)             :: Pi
      REAL(8)             :: RequiredAccuracy = 1.D0
      REAL(8),ALLOCATABLE :: ShapeFunc_A(:,:)
      REAL(8),ALLOCATABLE :: ShapeFunc_A0(:,:)
      REAL(8),ALLOCATABLE :: ShapeFunc_B(:,:)
      REAL(8),ALLOCATABLE :: Sigma(:)
      REAL(8)             :: Time1
      REAL(8)             :: Time2
      REAL(8),ALLOCATABLE :: TimeStepSizes(:)
      REAL(8),ALLOCATABLE :: VelU(:,:)
      REAL(8),ALLOCATABLE :: VelV(:,:)
      REAL(8),ALLOCATABLE :: VelW(:,:)

      ! ESSENTIAL INPUT
      ! If these parameters are not specified in the input file,
      ! then the program will print a fatal error message and stop.
      ! Thus there is no need to set default values for these parameters,
      ! because their values must be specified by the user.

      INTEGER :: NumTrackingSnaps = -99
      REAL(8) :: SimulationLength = -99999999.D0
      REAL(8) :: StartingTime = -99.D0
      REAL(8) :: MinimumTimeStep = 0.01D0
      INTEGER :: LatticeSearchBins = 1000
      NAMELIST /TimingInformation/ NumTrackingSnaps,                   &
     &                             SimulationLength,                   &
     &                             StartingTime,                       &
     &                             MinimumTimeStep,                    &
     &                             LatticeSearchBins

      CHARACTER(LEN=100) :: CurrentFile = "NULL"
      CHARACTER(LEN=100) :: CurrentFileOrigin = "ADCIRC"
      CHARACTER(LEN=100) :: CurrentFileFormat = "NULL"
      CHARACTER(LEN=100) :: CurrentDimensions = "NULL"
      NAMELIST /CurrentField/ CurrentFile,                             &
     &                        CurrentFileOrigin,                       &
     &                        CurrentFileFormat,                       &
     &                        CurrentDimensions

      CHARACTER(LEN=100) :: WindFile = "NULL"
      CHARACTER(LEN=100) :: WindFileFormat = "NULL"
      NAMELIST /WindField/ WindFile,                                   &
     &                     WindFileFormat

      CHARACTER(LEN=100) :: ParticleInputMethod = "NULL"
      CHARACTER(LEN=100) :: ParticleInputCoordinates = "Polar"
      CHARACTER(LEN=100) :: ParticleFile = "NULL"
      CHARACTER(LEN=100) :: ParticleFileFormat = "NULL"
      REAL(8) :: ParticleSourceX
      REAL(8) :: ParticleSourceY
      REAL(8) :: ParticleSourceZ
      NAMELIST /ParticleInput/ ParticleInputMethod,                    &
                               ParticleInputCoordinates,               &
     &                         ParticleFile,                           &
     &                         ParticleFileFormat,                     &
     &                         ParticleSourceX,                        &
     &                         ParticleSourceY,                        &
     &                         ParticleSourceZ

      CHARACTER(LEN=100) :: MeshFile = "NULL"
      CHARACTER(LEN=100) :: MeshFileOrigin = "ADCIRC"
      CHARACTER(LEN=100) :: MeshCoordinates = "Polar"
      NAMELIST /UnstructuredMesh/ MeshFile,                            &
     &                            MeshFileOrigin,                      &
     &                            MeshCoordinates


      CHARACTER(LEN=100) :: CombinationMethod = "Original"
      REAL(8) :: PercentageCurrent = 1.D0
      REAL(8) :: PercentageWind = 0.D0
      REAL(8) :: AngleWind = 0.07D0
      NAMELIST /VelocityCombination/ CombinationMethod,                &
     &                               PercentageCurrent,                &
     &                               PercentageWind,                   &
     &                               AngleWind

      CHARACTER(LEN=100) :: DiffusionMethod = "NULL"
      REAL(8) :: Cx  = 12.D0
      REAL(8) :: Cy  = 12.D0
      REAL(8) :: Evx = 10.D0
      REAL(8) :: Evy = 10.D0
      NAMELIST /Diffusion/ DiffusionMethod, Cx, Cy, Evx, Evy

      REAL(8) :: DensityWater = 998.2071D0 ! km/m3 at 20C
      REAL(8) :: DynamicViscosityWater = 0.001002D0 ! Pa-s at 20C
      NAMELIST /WaterProperties/ DensityWater,                         &
     &                           DynamicViscosityWater

      REAL(8) :: DensityOil = 858.D0 ! kg/m3 from Socolofsky et al. (2011)
      REAL(8) :: DiameterEffective = 0.00005D0 ! m 
      REAL(8) :: InterfacialTension = 0.023D0 ! N/m from Belore et al. (2011)
      NAMELIST /OilProperties/ DensityOil,                             &
     &                         DiameterEffective,                      &
     &                         InterfacialTension

      CHARACTER(LEN=100) :: BuoyancyMethod = "NULL"
      NAMELIST /OilPhysics/ BuoyancyMethod

      CHARACTER(LEN=100) :: OutputFileName = "NULL"
      CHARACTER(LEN=100) :: OutputFileFormat = "ASCII"
      NAMELIST /OutputSettings/ OutputFileName,                        &
     &                          OutputFileFormat

      END MODULE



      MODULE LPT_Data_Lattice_Table

      USE LPT_Data_Module, ONLY: NDIV => LatticeSearchBins

      IMPLICIT NONE

      INTEGER,ALLOCATABLE :: NE_PIECE(:,:)
      INTEGER,ALLOCATABLE :: NE_PIECE_INDEX(:,:)
      INTEGER,ALLOCATABLE :: NE_PIECE_LIST(:)

      REAL(8)             :: DX(1:2)
      REAL(8)             :: XMIN(1:2)

      END MODULE LPT_Data_Lattice_Table

