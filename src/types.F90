module types
use const, only: dp
implicit none
type :: tbc
   integer                  :: conType
   !< >0: BLOCKNUMMER
   !< BOUNDARY CONDITION, see CONSTS (BC_)

!!!!!!!! BLOCKCONNECTION
   integer                  :: conCPU
   !< CPU auf dem Sicher der Block befindet
   integer                  :: conFace
   !< Face of neighbouring BLOCK
   integer                  :: permutation
   !< Permutation of Block Connection

!!!!!!!! WALL
   integer                  :: WallType
   !< WallTYpe see CONST *_Wall
   
   real ( kind = dp )       :: WallTemp
   real ( kind = dp )       :: WallHeatFlux

!!!!!!!! INFLOW
   
   integer                  :: InflowType
   !< Art der Einströmung
   
end type tbc

type :: tblock
   integer                  :: nPkt(3)
   
   integer                  :: nCell(3)
   
   integer                  :: nOut(3)
   
   real(kind=dp),allocatable :: Q(:,:,:,:)
   !< Cell Centered Variables (I,J,K,VAR-INDEX)

   real(kind=dp),allocatable :: Res(:,:,:,:)
   !< Cell Residual Variables (I,J,K,VAR-INDEX)

   real(kind=dp),allocatable :: Flux(:,:,:,:,:)
   !< Cell Centered Variables (I,J,K,VAR-INDEX,FLUX_DIR)

   real(kind=dp),allocatable :: xyz(:,:,:,:)
   !< GITTERPUNKTE -POSITION (I,J,K,COORD)

   real(kind=dp),allocatable :: dt(:,:,:)
   !< Cell-Timestep (I,J,K)


   real(kind=dp),allocatable :: schwerpunkt(:,:,:,:)
   
   real(kind=dp),allocatable :: metric1(:,:,:,:,:)
   !< CELL Centered METRIKEN (I,J,K,X/Y/Z,I/J/K)
   !< dx / di @(i,j,k)

   real(kind=dp),allocatable :: metric2(:,:,:,:,:)
   !< CELL Centered Inversed METRIKEN (I,J,K,I/J/K,X/Y/Z)
   !< di / dx @(i,j,k)
   
   real(kind=dp),allocatable :: Jac (:,:,:)
   
   real(kind=dp),allocatable :: JacI (:,:,:)
   
   real(kind=dp),allocatable :: Edge_Len (:,:,:,:)
   !< Edge Length (I,J,K, EDGE_DIR)
   !< EDGE_DIR: 1 ist in i-Richtung (West&Ost), 2 ist in j-Richtung (Süd&Nord),...
   !<  i+1 West ist Ost, j+1 Süd ist Nord, k+1 Back ist Front

   real(kind=dp),allocatable :: Edge_Vec (:,:,:,:,:)
   !< Normalenvektor der Edge ( COORD,I,J,K, EDGE_DIR)
   !< Richtung ist für Ost&Süd Vektoren richtig, West&Nord muss mit -1 multipliziert werden.
   !< EDGE_DIR: 1 ist in i-Richtung (West&Ost), 2 ist in j-Richtung (Süd&Nord),...
   !< COORD gibt X,y,z wert an
   real(kind=dp),allocatable :: Area (:,:,:)

   type(tbc), allocatable    :: face(:)
   !< GIBT BLOCKVERBDINUNG DES BLOCKES AN: (FACE)
   !< FACE: 1 = W, 2 = E, 3 = S, 4 = N, 5 = B, 6 = F

end type tblock
end module types
