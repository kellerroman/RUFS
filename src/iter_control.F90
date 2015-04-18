subroutine iter_control()
use control
use global
implicit none
   sol_out = .false.
   res_out = .false.
   iteration = iteration + 1
   if (iteration >= control_num_iteration) then
      stop_iter = .true.
   end if
   if (mod(iteration,control_sol_out) == 0 .or. stop_iter) then
      sol_out = .true.
   end if
   if (mod(iteration,control_res_out) == 0 .or. iteration <= 50) then
      res_out = .true.
   end if

end subroutine iter_control
