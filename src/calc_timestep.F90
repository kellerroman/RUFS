subroutine calc_timestep
   use global
   use types
   use const
   use control
   implicit none

   integer :: b,i,j,k

   real(kind=dp) :: dt_min

   real(kind=dp) :: us, vs, ce, sq, p, rho, dt

   do b = 1, nBlock
      if (const_dt) then
         block(b) % dt = control_timestep
      else
         do k= 1, block(b) % nCell(3)
            do j= 1, block(b) % nCell(2)
               do i= 1, block(b) % nCell(1)
                  rho = block(b)%Q(i,j,k,1)
                  us = abs(block(b)%metric2(i,j,k,1,1) * block(b)%Q(i,j,k,2)           &
                     +     block(b)%metric2(i,j,k,1,2) * block(b)%Q(i,j,k,3))          &
                     /     rho
                  vs = abs(block(b)%metric2(i,j,k,2,1) * block(b)%Q(i,j,k,2)           &
                     +     block(b)%metric2(i,j,k,2,2) * block(b)%Q(i,j,k,3))          &
                     /     rho
                  p  = (0.4E0_dp) * (block(b)%Q(i,j,k,4)                               &
                     -  0.5E0_dp  *  rho * (block(b)%Q(i,j,k,2) * block(b)%Q(i,j,k,2)  &
                     + block(b)%Q(i,j,k,3) * block(b)%Q(i,j,k,3)) )
                  ce = sqrt( 1.4E0_dp * p / rho )
                  block(b) % dt(i,j,k) = control_CFL / ( us + vs + ce &
                     * sqrt( block(b)%metric2(i,j,k,1,1) * block(b)%metric2(i,j,k,1,1) &
                           + block(b)%metric2(i,j,k,1,2) * block(b)%metric2(i,j,k,1,2) &
                           + block(b)%metric2(i,j,k,2,1) * block(b)%metric2(i,j,k,2,1) &
                           + block(b)%metric2(i,j,k,2,2) * block(b)%metric2(i,j,k,2,2) &
                           + 2.d0 * abs( block(b)%metric2(i,j,k,1,1) * block(b)%metric2(i,j,k,2,1) &
                                       + block(b)%metric2(i,j,k,1,2) * block(b)%metric2(i,j,k,2,2) &
                                       ) &
                           ) &
                     )
               end do
            end do
         end do
      end if
   end do

   if (const_dt) then
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
