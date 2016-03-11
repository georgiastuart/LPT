! Casey 120815: Added the KDTREE2 search algorithm.
!
!     The implementation is VERY similar to the implementation
!     within ADCIRC, provided by Chris Massey.  The search tree
!     is initialized after the mesh is read, and then the search
!     routine can be called from anywhere in the LPT code.
!
!     The KDTREE2 search algorithm requires the kdtree2.f90
!     source file, and it can be triggered at compilation
!     with the -DKDTREE pre-compiler flag.
!
!     It was found that the the KDTREE2 search algorithm
!     is MUCH SLOWER than the lattice search algorithm that was
!     implemented by Seizo Tanaka.  Conceptually, this makes sense,
!     because the lattice search algorithm should only be operating
!     on a small subset of the mesh elements.  No matter how fast
!     the KDTREE2 search algorithm is, it is still operating
!     on the full set of the mesh elements.
!
!     For this reason, the KDTREE2 search algorithm was NOT
!     implemented fully in the LPT code.  It is only available
!     as a search option to determine the starting element
!     for each of the particles.
!
!     THUS THE KDTREE2 SEARCH ALGORITHM IS NOT RECOMMENDED.



      MODULE LPT_KDTREE_Module

#ifdef KDTREE
      USE kdtree2_module
#endif

      IMPLICIT NONE

#ifdef KDTREE
      REAL(8),ALLOCATABLE :: ElemCenter(:,:)
      REAL(8),ALLOCATABLE :: ElemRadii(:)

      TYPE(KDTREE2),POINTER            :: KD_Tree
      TYPE(KDTREE2_RESULT),ALLOCATABLE :: KD_Result(:)
#endif



#ifdef KDTREE
      CONTAINS



      SUBROUTINE LPT_KDTREE_Initialize

      USE LPT_Comm_Module, ONLY: MyRank
      USE LPT_Data_Module, ONLY: MeshConnEV,                           &
     &                           MeshLat,                              &
     &                           MeshLon,                              &
     &                           NumElems

      IMPLICIT NONE

      INTEGER :: IE
      INTEGER :: IV
      INTEGER :: Vert(2)

      REAL(8) :: EdgeLength(3)
      REAL(8) :: ElemArea
      REAL(8) :: Lat(3)
      REAL(8) :: Lon(3)

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Entering the "                          &
     &         //"LPT_KDTREE_Initialize routine.")
#endif

      IF(.NOT.ALLOCATED(ElemCenter))THEN
         ALLOCATE(ElemCenter(1:2,1:NumElems))
      ENDIF
      IF(.NOT.ALLOCATED(ElemRadii))THEN
         ALLOCATE(ElemRadii(1:NumElems))
      ENDIF
#if VERBOSE > 2
      IF(MyRank.EQ.0)THEN
         WRITE(*,'(A,$)') "LPT: INFO: Computing the element centers "  &
     &         //"and radii: +"
      ENDIF
#endif
      DO IE=1,NumElems
         DO IV=1,3
            Lon(IV)  = MeshLon(MeshConnEV(IE,IV))
            Lat(IV)  = MeshLat(MeshConnEV(IE,IV))
         ENDDO
         ElemCenter(1,IE) = (Lon(1)+Lon(2)+Lon(3))/3.D0
         ElemCenter(2,IE) = (Lat(1)+Lat(2)+Lat(3))/3.D0
         DO IV=1,3
            Vert(1) = IV
            Vert(2) = IV+1
            IF(Vert(2).GT.3) Vert(2) = 1
            EdgeLength(IV) = SQRT( (Lon(Vert(2))-Lon(Vert(1)))**2      &
     &                           + (Lat(Vert(2))-Lat(Vert(1)))**2 )
         ENDDO
         ElemArea = (Lon(1)-Lon(3)) * (Lat(2)-Lat(1))                  &
     &            + (Lon(3)-Lon(2)) * (Lat(1)-Lat(3))
         ElemRadii(IE) = 1.5D0                                         &
     &            * (EdgeLength(1)*EdgeLength(2)*EdgeLength(3))        &
     &            / (2.D0 * ElemArea)
#if VERBOSE > 2
         IF(MyRank.EQ.0)THEN
            CALL LPT_Progress(IE,NumElems)
         ENDIF
#endif
      ENDDO

      KD_Tree => kdtree2_create(ElemCenter,rearrange=.TRUE.,           &
     &                          sort=.TRUE.)

      IF(.NOT.ALLOCATED(KD_Result))THEN
         ALLOCATE(KD_Result(1:NumElems))
      ENDIF

#if VERBOSE > 1
      CALL LPT_Print(0,"INFO","Exiting the "                           &
     &         //"LPT_KDTREE_Initialize routine.")
#endif

      ! We're done here.
      RETURN

      END SUBROUTINE
#endif



#ifdef KDTREE
      SUBROUTINE LPT_KDTREE_Search(Lon,Lat,Element,FoundFlag)

      USE LPT_Data_Module, ONLY: MeshConnEV,                           &
     &                           MeshLat,                              &
     &                           MeshLon,                              &
     &                           NumElems

      IMPLICIT NONE

      INTEGER,INTENT(OUT) :: Element
      INTEGER,INTENT(OUT) :: FoundFlag

      REAL(8),INTENT(IN) :: Lat
      REAL(8),INTENT(IN) :: Lon

      INTEGER :: ClosestElem
      INTEGER :: ElemNumber
      INTEGER :: IE
      INTEGER :: IV
      INTEGER :: SearchElems
      INTEGER :: Vert(3)

      LOGICAL :: ElemFound

      REAL(8) :: A1
      REAL(8) :: A2
      REAL(8) :: A3
      REAL(8) :: AA
      REAL(8) :: AE
      REAL(8) :: AREASK
      REAL(8) :: Dist
      REAL(8) :: ElemMin(2)
      REAL(8) :: VertLat(3)
      REAL(8) :: VertLon(3)
      REAL(8) :: X1
      REAL(8) :: X2
      REAL(8) :: X3
      REAL(8) :: Y1
      REAL(8) :: Y2
      REAL(8) :: Y3

      REAL(8),PARAMETER :: Tolerance = 1.0D-5

      ElemFound = .FALSE.
!     SearchElems = NumElems
      SearchElems = 20

      CALL kdtree2_n_nearest(tp=KD_Tree,qv=(/Lon,Lat/),                &
     &                  nn=NumElems,results=KD_Result)

      IE = 1
      ClosestElem = KD_Result(IE)%idx

      ElemMin = MINVAL(SQRT(KD_Result(1:SearchElems)%dis)              &
     &                      - ElemRadii(KD_Result(1:SearchElems)%idx) )

      IF(ElemMin(1).LE.0.0D0)THEN
         DO WHILE((.NOT.ElemFound).AND.(IE.LE.SearchElems))
            ElemNumber = KD_Result(IE)%idx
            Dist = SQRT(KD_Result(IE)%dis)
            IF(Dist-ElemRadii(ElemNumber).LE.0.0D0)THEN
               DO IV=1,3
                  Vert(IV) = MeshConnEV(ElemNumber,IV)
                  VertLon(IV) = MeshLon(Vert(IV))
                  VertLat(IV) = MeshLat(Vert(IV))
               ENDDO
               X1=VertLon(1)
               X2=VertLon(2)
               X3=VertLon(3)
               Y1=VertLat(1)
               Y2=VertLat(2)
               Y3=VertLat(3)
               A1=(Lon-X3)*(Y2-Y3)+(X2-X3)*(Y3-Lat)
               A2=(Lon-X1)*(Y3-Y1)-(Lat-Y1)*(X3-X1)
               A3=(Lat-Y1)*(X2-X1)-(Lon-X1)*(Y2-Y1)
               AA=ABS(A1)+ABS(A2)+ABS(A3)
               AREASK=X2*Y3+X1*Y2+X3*Y1-Y1*X2-Y2*X3-Y3*X1
               AE=ABS(AA-AREASK)/AREASK
               IF (AE.LT.Tolerance)THEN
                  ElemFound = .TRUE.
                  ClosestElem = ElemNumber
                  Element = ClosestElem
                  FoundFlag = 1
               ELSE
                  IE = IE + 1
               ENDIF
            ELSE
               IE = IE + 1
            ENDIF
         ENDDO
      ENDIF

      IF(.NOT.ElemFound)THEN
         Element = 0
         FoundFlag = 0
      ENDIF

      ! We're done here.
      RETURN

      END SUBROUTINE
#endif



      END MODULE
