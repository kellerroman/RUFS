module global
use const
use types
implicit none
public 
save
type(tblock),allocatable            :: block(:)

character( len = io_len_VarName), allocatable   :: VarName(:)

integer                             :: Dimen
!< Dimension des Gitters
integer                             :: nFaces
!< Anzahl der Rechenblock-Randflächen
integer                             :: nCorners
!< Anzahl der Rechenblock-Randecken
integer                             :: nBlock
!< Anzahl der Rechenblöcke

logical :: res_out
!< Is set true if Residual is outputted to screen/file
!< is set in iter_control
logical :: sol_out
!< is True if in this iteration the solution is written to the output file
!< is set in iter_control
logical :: stop_iter
!< is True if the iteration is stopped
!< is set in iter_control
logical :: write_sol_header
!< is true is header must be written, in case of restart can be false to append to current file
!< is set in input_control

integer                             :: iteration
!< is the current iteration

end module global
