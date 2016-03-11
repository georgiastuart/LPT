      SUBROUTINE LPT_Oil_Buoyancy(DiameterEffective,VelocityTerminal)

      USE LPT_Data_Module, ONLY: BuoyancyMethod,                       &
     &                           DensityOil,                           &
     &                           DensityWater,                         &
     &                           DynamicViscosityWater,                &
     &                           Gravity,                              &
     &                           InterfacialTension

      IMPLICIT NONE

      REAL(8),INTENT(IN)  :: DiameterEffective
      REAL(8),INTENT(OUT) :: VelocityTerminal

      REAL(8) :: DensityDifference
      REAL(8) :: DiameterCritical

      REAL(8) :: A1
      REAL(8) :: A2
      REAL(8) :: B1
      REAL(8) :: B2
      REAL(8) :: EO
      REAL(8) :: H
      REAL(8) :: J
      REAL(8) :: M
      REAL(8) :: ND
      REAL(8) :: R
      REAL(8) :: W
      REAL(8) :: X1
      REAL(8) :: X2
      REAL(8) :: Y1
      REAL(8) :: Y2

      DensityDifference = ABS(DensityOil - DensityWater)

      ! Two-Equation Approach
      IF(INDEX(BuoyancyMethod,"CLIFT1978").GT.0)THEN

         ! Compute the critical diameter.
         DiameterCritical = (9.52D0                                    &
     &      *(DynamicViscosityWater**(2.D0/3.D0)))                     &
     &      /((Gravity*DensityWater*DensityDifference)**(1.D0/3.D0))

         ! Small Reynolds numbers.
         IF(DiameterEffective.LE.DiameterCritical)THEN

            VelocityTerminal = (Gravity*(DiameterEffective**2.D0)      &
     &           *DensityDifference)/(18.D0*DynamicViscosityWater)

         ! Larger Reynolds numbers.
         ELSEIF(DiameterEffective.GT.DiameterCritical)THEN

            VelocityTerminal = SQRT((8.D0*Gravity*DiameterEffective    &
     &            *DensityDifference)/(3.D0*DensityWater))

         ENDIF

      ! Integrated Approach
      ELSEIF(INDEX(BuoyancyMethod,"ZHENG2000").GT.0)THEN

         ! Spherical shape, small size range
         IF(DiameterEffective.LE.0.001D0)THEN

            ND = 4.D0*DensityWater*DensityDifference*Gravity           &
     &        *(DiameterEffective**3.D0)/(3.D0                         &
     &        *(DynamicViscosityWater**2.D0))
            W = LOG10(ND)
            IF(ND.LE.73.D0)THEN
               R = ND/24.D0 - (1.7569D0*10.D0**-4.D0)*(ND**2.D0)       &
     &            + (6.9252D0*10.D0**-7.D0)*(ND**3.D0)                 &
     &            - (2.3027D0*10.D0**-10.D0)*(ND**4.D0)
            ELSEIF(ND.LE.580.D0)THEN
               R = -1.7095D0 + 1.33438D0*W - 0.11591D0*(W**2.D0)
               R = 10.D0**R
            ELSE
               R = -1.81391D0 + 1.34671D0*W - 0.12427D0*(W**2.D0)      &
     &            + 0.006344D0*(W**3.D0)
               R = 10.D0**R
            ENDIF
            VelocityTerminal = (R*DynamicViscosityWater)/(DensityWater &
     &         *DiameterEffective)

         ELSEIF(DiameterEffective.GT.0.001D0)THEN

            ! Compute X1, Y1
            X1 = LOG10(DiameterEffective)
            M = (Gravity*(DynamicViscosityWater**4.D0)                 &
     &         *DensityDifference)/((DensityWater**2.D0)               &
     &         *(InterfacialTension**3.D0))
            H = 59.3D0
            J = 0.94D0*(H**0.757D0)
            VelocityTerminal = DynamicViscosityWater/(DensityWater     &
     &         *DiameterEffective)*(M**-0.149D0)*(J-0.857D0)
            Y1 = LOG10(VelocityTerminal)

            ! Compute X2, Y2
            X2 = LOG10(0.015D0)
            EO = (Gravity*DensityDifference*(0.015D0**2.D0))           &
     &          /InterfacialTension
            M = (Gravity*(DynamicViscosityWater**4.D0)                 &
     &         *DensityDifference)/((DensityWater**2.D0)               &
     &         *(InterfacialTension**3.D0))
            H = (4.D0/3.D0)*EO*(M**-0.149D0)*(1.D0**-0.14D0)
            IF(H.LE.59.3D0)THEN
               J = 0.94D0*(H**0.757D0)
            ELSE
               J = 3.42D0*(H**0.441D0)
            ENDIF
            VelocityTerminal = DynamicViscosityWater/(DensityWater     &
     &         *0.015D0)*(M**-0.149D0)*(J-0.857D0)
            Y2 = LOG10(VelocityTerminal)

            ! Compute A1, B1, A2, B2
            A1 = 0.5D0
            B1 = LOG10(0.711D0*SQRT(Gravity*DensityDifference          &
     &          /DensityWater))
            A2 = (Y2-Y1)/(X2-X1)
            B2 = Y1 - A2*X1

            ! Compute the critical diameter.
            DiameterCritical = 10.D0**((B2-B1)/(A1-A1))

            ! Ellipsoidal shape, intermediate size range.
            IF(DiameterEffective.LE.DiameterCritical)THEN

               EO = (Gravity*DensityDifference                         &
     &            *(DiameterEffective**2.D0))/InterfacialTension
               M = (Gravity*(DynamicViscosityWater**4.D0)              &
     &            *DensityDifference)/((DensityWater**2.D0)            &
     &            *(InterfacialTension**3.D0))
               H = (4.D0/3.D0)*EO*(M**-0.149D0)*(1.D0**-0.14D0)
               IF(H.LE.59.3D0)THEN
                  J = 0.94D0*(H**0.757D0)
               ELSE
                  J = 3.42D0*(H**0.441D0)
               ENDIF
               VelocityTerminal = DynamicViscosityWater/(DensityWater  &
     &            *DiameterEffective)*(M**-0.149D0)*(J-0.857D0)

            ! Spherical-cap share, large size range.
            ELSEIF(DiameterEffective.GT.DiameterCritical)THEN

               VelocityTerminal = 0.711D0*SQRT(Gravity                 &
     &            *DiameterEffective*DensityDifference/DensityWater)

            ENDIF

         ENDIF

      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE



      SUBROUTINE LPT_Oil_Source(XNEW,YNEW,ZNEW,DNEW,NNEW,SourceType)

      USE LPT_Data_Lattice_Table
      USE LPT_Data_Module, ONLY: Deg2Rad,                              &
     &                           EarthRadius,                          &
     &                           NumParticlesGlobal,                   &
     &                           NumParticlesLocal,                    &
     &                           NumTrackingSnaps,                     &
     &                           PartBoundary,                         &
     &                           PartBoundaryCode,                     &
     &                           PartElemJ,                            &
     &                           PartElemL,                            &
     &                           PartDepthLocal,                       &
     &                           ParticleInputCoordinates,             &
     &                           PartDiameterLocal,                    &
     &                           PartDomainLocal,                      &
     &                           PartNumberLocal,                      &
     &                           PartLatLocal,                         &
     &                           PartLonLocal,                         &
     &                           SimulationLength,                     &
     &                           TimeStepSizes

      IMPLICIT NONE

      INTEGER,INTENT(IN)  :: NNEW
      INTEGER,INTENT(IN)  :: SourceType

      REAL(8),INTENT(IN)  :: DNEW
      REAL(8),INTENT(IN)  :: XNEW
      REAL(8),INTENT(IN)  :: YNEW
      REAL(8),INTENT(IN)  :: ZNEW

      INTEGER             :: ArrayLength
      INTEGER             :: ICHECK
      INTEGER             :: JJDR
      INTEGER             :: LLDR
      INTEGER,ALLOCATABLE :: TempI(:)

      REAL(8),ALLOCATABLE :: TempR(:)
      REAL(8)             :: XSTART
      REAL(8)             :: YSTART
      REAL(8)             :: ZSTART

      ! Increase the number of particles.
      NumParticlesLocal = NumParticlesLocal + 1
      IF(SourceType.EQ.0)THEN
         NumParticlesGlobal = NumParticlesGlobal + 1
      ENDIF

      ! Expand the existing integer arrays.
      IF(ALLOCATED(PartElemJ))THEN
         ArrayLength = SIZE(PartElemJ)
         ALLOCATE(TempI(1:ArrayLength+1))
         TempI(1:ArrayLength) = PartElemJ(1:ArrayLength)
         DEALLOCATE(PartElemJ)
         ALLOCATE(PartElemJ(1:ArrayLength+1))
         PartElemJ(1:ArrayLength) = TempI(1:ArrayLength)
         DEALLOCATE(TempI)
      ELSE
         ALLOCATE(PartElemJ(1:1))
      ENDIF
      IF(ALLOCATED(PartElemL))THEN
         ArrayLength = SIZE(PartElemL)
         ALLOCATE(TempI(1:ArrayLength+1))
         TempI(1:ArrayLength) = PartElemL(1:ArrayLength)
         DEALLOCATE(PartElemL)
         ALLOCATE(PartElemL(1:ArrayLength+1))
         PartElemL(1:ArrayLength) = TempI(1:ArrayLength)
         DEALLOCATE(TempI)
      ELSE
         ALLOCATE(PartElemL(1:1))
      ENDIF
      IF(ALLOCATED(PartBoundary))THEN
         ArrayLength = SIZE(PartBoundary)
         ALLOCATE(TempI(1:ArrayLength+1))
         TempI(1:ArrayLength) = PartBoundary(1:ArrayLength)
         DEALLOCATE(PartBoundary)
         ALLOCATE(PartBoundary(1:ArrayLength+1))
         PartBoundary(1:ArrayLength) = TempI(1:ArrayLength)
         DEALLOCATE(TempI)
      ELSE
         ALLOCATE(PartBoundary(1:1))
      ENDIF
      IF(ALLOCATED(PartBoundaryCode))THEN
         ArrayLength = SIZE(PartBoundaryCode)
         ALLOCATE(TempI(1:ArrayLength+1))
         TempI(1:ArrayLength) = PartBoundaryCode(1:ArrayLength)
         DEALLOCATE(PartBoundaryCode)
         ALLOCATE(PartBoundaryCode(1:ArrayLength+1))
         PartBoundaryCode(1:ArrayLength) = TempI(1:ArrayLength)
         DEALLOCATE(TempI)
      ELSE
         ALLOCATE(PartBoundaryCode(1:1))
      ENDIF
      IF(ALLOCATED(PartNumberLocal))THEN
         ArrayLength = SIZE(PartNumberLocal)
         ALLOCATE(TempI(1:ArrayLength+1))
         TempI(1:ArrayLength) = PartNumberLocal(1:ArrayLength)
         DEALLOCATE(PartNumberLocal)
         ALLOCATE(PartNumberLocal(1:ArrayLength+1))
         PartNumberLocal(1:ArrayLength) = TempI(1:ArrayLength)
         DEALLOCATE(TempI)
      ELSE
         ALLOCATE(PartNumberLocal(1:1))
      ENDIF
      IF(ALLOCATED(PartDomainLocal))THEN
         ArrayLength = SIZE(PartDomainLocal)
         ALLOCATE(TempI(1:ArrayLength+1))
         TempI(1:ArrayLength) = PartDomainLocal(1:ArrayLength)
         DEALLOCATE(PartDomainLocal)
         ALLOCATE(PartDomainLocal(1:ArrayLength+1))
         PartDomainLocal(1:ArrayLength) = TempI(1:ArrayLength)
         DEALLOCATE(TempI)
      ELSE
         ALLOCATE(PartDomainLocal(1:1))
      ENDIF

      ! Expand the existing real arrays.
      IF(ALLOCATED(PartLonLocal))THEN
         ArrayLength = SIZE(PartLonLocal)
         ALLOCATE(TempR(1:ArrayLength+1))
         TempR(1:ArrayLength) = PartLonLocal(1:ArrayLength)
         DEALLOCATE(PartLonLocal)
         ALLOCATE(PartLonLocal(1:ArrayLength+1))
         PartLonLocal(1:ArrayLength) = TempR(1:ArrayLength)
         DEALLOCATE(TempR)
      ELSE
         ALLOCATE(PartLonLocal(1:1))
      ENDIF
      IF(ALLOCATED(PartLatLocal))THEN
         ArrayLength = SIZE(PartLatLocal)
         ALLOCATE(TempR(1:ArrayLength+1))
         TempR(1:ArrayLength) = PartLatLocal(1:ArrayLength)
         DEALLOCATE(PartLatLocal)
         ALLOCATE(PartLatLocal(1:ArrayLength+1))
         PartLatLocal(1:ArrayLength) = TempR(1:ArrayLength)
         DEALLOCATE(TempR)
      ELSE
         ALLOCATE(PartLatLocal(1:1))
      ENDIF
      IF(ALLOCATED(PartDepthLocal))THEN
         ArrayLength = SIZE(PartDepthLocal)
         ALLOCATE(TempR(1:ArrayLength+1))
         TempR(1:ArrayLength) = PartDepthLocal(1:ArrayLength)
         DEALLOCATE(PartDepthLocal)
         ALLOCATE(PartDepthLocal(1:ArrayLength+1))
         PartDepthLocal(1:ArrayLength) = TempR(1:ArrayLength)
         DEALLOCATE(TempR)
      ELSE
         ALLOCATE(PartDepthLocal(1:1))
      ENDIF
      IF(ALLOCATED(PartDiameterLocal))THEN
         ArrayLength = SIZE(PartDiameterLocal)
         ALLOCATE(TempR(1:ArrayLength+1))
         TempR(1:ArrayLength) = PartDiameterLocal(1:ArrayLength)
         DEALLOCATE(PartDiameterLocal)
         ALLOCATE(PartDiameterLocal(1:ArrayLength+1))
         PartDiameterLocal(1:ArrayLength) = TempR(1:ArrayLength)
         DEALLOCATE(TempR)
      ELSE
         ALLOCATE(PartDiameterLocal(1:1))
      ENDIF
      IF(ALLOCATED(TimeStepSizes))THEN
         ArrayLength = SIZE(TimeStepSizes)
         ALLOCATE(TempR(1:ArrayLength+1))
         TempR(1:ArrayLength) = TimeStepSizes(1:ArrayLength)
         DEALLOCATE(TimeStepSizes)
         ALLOCATE(TimeStepSizes(1:ArrayLength+1))
         TimeStepSizes(1:ArrayLength) = TempR(1:ArrayLength)
         DEALLOCATE(TempR)
      ELSE
         ALLOCATE(TimeStepSizes(1:1))
      ENDIF

      ! Add the new particle.
      XSTART = XNEW
      YSTART = YNEW
      ZSTART = ZNEW

      ! Convert the coordinates if necessary.
      IF(INDEX(ParticleInputCoordinates,"POLAR").GT.0)THEN
         XSTART = XSTART * Deg2Rad
         YSTART = YSTART * Deg2Rad
         XSTART = XSTART * EarthRadius
         YSTART = YSTART * EarthRadius
      ENDIF

      ! Determine the horizontal and vertical elements
      ! containing the new particle.
      ICHECK = 0
      CALL LPT_Drog_Find_Element(XSTART,YSTART,ZSTART,JJDR,LLDR,ICHECK)

      ! Assign the new variables.
      PartLonLocal(NumParticlesLocal) = XSTART
      PartLatLocal(NumParticlesLocal) = YSTART
      PartDepthLocal(NumParticlesLocal) = ZSTART
      PartDiameterLocal(NumParticlesLocal) = DNEW
      PartNumberLocal(NumParticlesLocal) = NNEW
      PartDomainLocal(NumParticlesLocal) = ICHECK
      PartElemJ(NumParticlesLocal) = JJDR
      PartElemL(NumParticlesLocal) = LLDR
      TimeStepSizes(NumParticlesLocal) = DABS(SimulationLength         &
     &             /FLOAT(NumTrackingSnaps)/10.d0)
      PartBoundary(NumParticlesLocal) = 0
      PartBoundaryCode(NumParticlesLocal) = 0

      ! We're done here.
      RETURN

      END SUBROUTINE

