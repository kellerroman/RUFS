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
integer                             :: control_bc_cells_out
!< defines if bc_cells are outputtet to sol file
!< 0: not outputted
!< 1: outputted

!!! CALCULATED VARIABLES 
integer :: n_BC_Cells

integer :: nVar
!< Anzahl der Variablen im Variablenvektor Q

integer :: sol_out_nVar
!< Anzahl der Variablen in der Ausgabedatei

end module control
