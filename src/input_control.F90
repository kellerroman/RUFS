subroutine input_control()
use control
use global
implicit none

integer :: t

! this is read in input_control
n_BC_Cells = 1

nVar = 4

control_num_iteration = 10000
control_sol_out = 1000
control_res_out = 500
iteration = 0
time = 0.0E0_dp
const_dt = .true.
control_CFL = 0.1E0_dp


write_sol_header = .true.

control_bc_cells_out = 0

sol_out_nVar = 6
t = 1
allocate (VarName(sol_out_nVar))
VarName(t) = VarName_Rho; t = t + 1
VarName(t) = VarName_SpU; t = t + 1
VarName(t) = VarName_SpV; t = t + 1
VarName(t) = VarName_Ene; t = t + 1
!VarName(t) = VarName_SwpX; t = t + 1
!VarName(t) = VarName_SwpY; t = t + 1
VarName(t) = VarName_JacI; t = t + 1
VarName(t) = VarName_Jac; t = t + 1
!VarName(t) = VarName_M1XI; t = t + 1
!VarName(t) = VarName_M1XJ; t = t + 1
!VarName(t) = VarName_M1YI; t = t + 1
!VarName(t) = VarName_M1YJ; t = t + 1
!VarName(t) = VarName_M2XI; t = t + 1
!VarName(t) = VarName_M2XJ; t = t + 1
!VarName(t) = VarName_M2YI; t = t + 1
!VarName(t) = VarName_M2YJ; t = t + 1

end subroutine input_control
