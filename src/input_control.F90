subroutine input_control()
use control
use global
implicit none
integer, parameter :: fu = 111
integer, parameter :: cfg_version = 2
integer :: t, version
logical :: fexists

iteration = 0
time = 0.0E0_dp
sol_out_nVar = 0

write_sol_header = .true.
control_bc_cells_out = 0

inquire(file=trim(file_cfg_in),exist=fexists)

if(.not. fexists) then
   call error("CONFIG DATEI konnte nicht gefunden werden: "//TRIM(file_cfg_in),__FILE__,__LINE__)
end if

open(fu,file=trim(file_cfg_in),form="UNFORMATTED",access="STREAM",status="OLD")

read(fu) version

if (version /= cfg_version) then
   call error("CONFIG FILE VERSION IS WRONG",__FILE__,__LINE__)
end if

read(fu)   control_num_iteration,control_sol_out,control_res_out &
          ,control_equation,control_dimension, space_disc, space_order &
          ,control_dt_method, control_riemann_solver,control_CFL &
          ,control_timestep

call wr("CONFIG INFORMATION",2)
write(*,'(A40," = ",I0)') "Number of Iterations", control_num_iteration
write(*,'(A40," = ",I0)') "Output Solution", control_sol_out
write(*,'(A40," = ",I0)') "Output Residual", control_res_out
write(*,'(A40," = ",A)') "Riemann Solver", string_Riemann_Solver(control_riemann_solver)
write(*,'(A40," = ",A)') "Time Integration Theme",string_dt_solver(control_dt_method)
write(*,'(A40," = ",ES10.4)') "Timestep",control_timestep

if ( control_equation == 1 .and. control_dimension <= 2 ) then
   nVar = 4
   sol_out_nVar = sol_out_nVar + 4
end if

if (space_order < 3) then
   n_BC_Cells = 1
else
   n_BC_Cells = 2
end if

t = 1
allocate ( VarName ( sol_out_nVar ) )
if ( control_equation == 1 .and. control_dimension <= 2 ) then
   VarName(t) = VarName_Rho; t = t + 1
   VarName(t) = VarName_SpU; t = t + 1
   VarName(t) = VarName_SpV; t = t + 1
   VarName(t) = VarName_Ene; t = t + 1
end if
!VarName(t) = VarName_EdLeW; t = t + 1
!VarName(t) = VarName_EdLeE; t = t + 1
!VarName(t) = VarName_EdAnW; t = t + 1
!VarName(t) = VarName_EdAnE; t = t + 1
!VarName(t) = VarName_EdAnS; t = t + 1
!VarName(t) = VarName_EdAnN; t = t + 1
!VarName(t) = VarName_SwpX; t = t + 1
!VarName(t) = VarName_SwpY; t = t + 1
!VarName(t) = VarName_Area; t = t + 1
!VarName(t) = VarName_JacI; t = t + 1
!VarName(t) = VarName_Jac; t = t + 1
!VarName(t) = VarName_M1XI; t = t + 1
!VarName(t) = VarName_M1XJ; t = t + 1
!VarName(t) = VarName_M1YI; t = t + 1
!VarName(t) = VarName_M1YJ; t = t + 1
!VarName(t) = VarName_M2XI; t = t + 1
!VarName(t) = VarName_M2XJ; t = t + 1
!VarName(t) = VarName_M2YI; t = t + 1
!VarName(t) = VarName_M2YJ; t = t + 1

!write(*,*) control_num_iteration,control_sol_out,control_res_out &
!          ,equation,soldim, space_disc, space_order, control_CFL
end subroutine input_control
