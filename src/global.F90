module global
use const, only: dp
use types
implicit none
public 
save
type(tblock),allocatable            :: block(:)

integer                             :: Dimen
!< Dimension des Gitters
integer                             :: nFaces
!< Anzahl der Rechenblock-Randflächen
integer                             :: nBlock
!< Anzahl der Rechenblöcke
end module global
