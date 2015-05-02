program solver
use const
use global
implicit none
call wr(achar(27)//"[31mSolver V0.0.1 "//__DATE__//" "//__TIME__//achar(27)//"[0m",1)

call input()
call init()

!!!! MAIN LOOP

call wr(achar(27)//"[31mStart Calculation"//achar(27)//"[0m",2)
stop_iter = .false.

main_loop: do while (.not. stop_iter)
   call iter_control !! SETTING VARIABLES FOR EACH ITERATION

   call set_boundary()
   call calc_fluxes()
!   call calc_residual()

   call write_residual()
!   call calc_timestep()
!   call update_sol()

   call write_solution()
end do main_loop



call end_program(EXIT_WITHOUT_ERROR)

end program
