subroutine calc_residual
   use global
   use types
   use control
   use const
   implicit none

   integer :: b,i,j,k
   minRes = 1E10_dp
   maxRes = -1E10_dp
   avgRes = 0.0E0_dp
   do b = 1,nBlock
      do k = 1, block(b) % nCell(3)
         do j = 1, block(b) % nCell(2)
            do i = 1, block(b) % nCell(1)
               block(b) % Res(i,j,k,:) = + block(b) % Flux(i+1,j,k,:,1) - block(b) % Flux(i,j,k,:,1) &
                                         + block(b) % Flux(i,j+1,k,:,2) - block(b) % Flux(i,j,k,:,2)
!               if (j == 1) then
!                  write(*,*) i,block(b) % Res(i,j,k,:)
!               end if
            end do
         end do
      end do
      minRes = min(minRes,minval(abs(block(b) % Res(:,:,:,1) )))
      maxRes = max(minRes,maxval(abs(block(b) % Res(:,:,:,1) )))
      avgRes = avgRes + sum(abs(block(b) % Res(:,:,:,1) ))
   end do
   avgRes = avgRes / nCell
end subroutine calc_residual
