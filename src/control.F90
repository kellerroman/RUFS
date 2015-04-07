module control
implicit none
save
public
!!! CONFIG VARIABLES FROM INPUT FILE
integer :: space_Order

character ( len = 100 )             :: file_git_in = "git.bin"

character ( len = 100 )             :: file_bc_in  = "bc.bin"

!!! CALCULATED VARIABLES 
integer :: n_BC_Cells

integer :: nVar


end module control