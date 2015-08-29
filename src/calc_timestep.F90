subroutine calc_timestep
   use global
   use types
   use const
   use control
   implicit none

   integer :: b,i,j,k

   real(kind=dp) :: dt_min

   real(kind=dp) :: u, v,rho,x,y,dt

   do b = 1, nBlock
      if (control_dt_method == 1) then
         block(b) % dt = control_timestep
         time  = time + control_timestep
      else
         do k= 1, block(b) % nCell(3)
            do j= 1, block(b) % nCell(2)
               do i= 1, block(b) % nCell(1)
                  rho = block(b) % Q(i,j,k,1)
                  u = block(b) % Q(i,j,k,2) / rho
                  v = block(b) % Q(i,j,k,3) / rho
                  x = block(b) % len_dt(1,i,j,k)
                  y = block(b) % len_dt(2,i,j,k)
                  dt = min(x/u,y/v) * control_CFL
                  block(b) % dt(i,j,k) = dt
               end do
            end do
         end do
      end if
   end do

   if (control_dt_method == 3) then
      dt_min = 1E10_dp
      !!!MIN dt aller blöcke finden
      do b = 1, nBlock
         dt_min = min(dt_min,minval(block(b) % dt))
      end do
      !!! hier muss komuniziert werden im falle von MPI
      do b = 1, nBlock
         block(b) % dt = dt_min
      end do
      time  = time + dt_min
   end if

end subroutine calc_timestep
