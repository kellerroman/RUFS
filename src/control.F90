module control
use const
implicit none
save
public
!!! CONFIG VARIABLES FROM INPUT FILE
integer :: space_Order

character ( len = len_str_filename ):: file_git_in          = "git.bin"
!<
character ( len = len_str_filename ):: file_bc_in           = "bc.bin"

character ( len = len_str_filename ):: file_sol_out         = "sol.bin"

integer                             :: control_num_iteration
integer                             :: control_sol_out
integer                             :: control_res_out

!!! CALCULATED VARIABLES 
integer :: n_BC_Cells

integer :: nVar
!< Anzahl der Variablen im Variablenvektor

integer :: sol_out_nVar
!< Anzahl der Variablen in der Ausgabedatei

end module control
