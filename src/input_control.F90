subroutine input_control()
use control
use global
implicit none

! this is read in input_control
n_BC_Cells = 2

nVar = 4

control_num_iteration = 10
control_sol_out = 5
control_res_out = 5
iteration = 0

write_sol_header = .true.

sol_out_nVar = 7
allocate (VarName(sol_out_nVar))
VarName(1) = VarName_Jac
VarName(2) = VarName_SwpX
VarName(3) = VarName_SwpY
VarName(4) = VarName_M2XI
VarName(5) = VarName_M2XJ
VarName(6) = VarName_M2YI
VarName(7) = VarName_M2YJ

end subroutine input_control
