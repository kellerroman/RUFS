program solver
use const
use global
implicit none
call wr("Solver V0.0.1 "//__DATE__//" "//__TIME__,1)

call input()
call init()

!!!! MAIN LOOP

call wr("Start Calculation",2)

stop_iter = .false.

main_loop: do while (.not. stop_iter)
   call iter_control !! SETTING VARIABLES FOR EACH ITERATION

!   call set_boundary()
!   call calc_fluxes()
!   call calc_residual()

   call write_residual()

!   call update_sol()

   call write_solution()
end do main_loop



call end_program(EXIT_WITHOUT_ERROR)

end program
