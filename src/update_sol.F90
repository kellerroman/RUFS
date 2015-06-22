subroutine update_sol
   use global
   use const
   use control
   implicit none
   integer :: b,i,j,k

   do b = 1,nBlock
      do k = 1, block(b) % nCell(3)
         do j = 1, block(b) % nCell(2)
            do i = 1, block(b) % nCell(1)
               block(b) % Q(i,j,k,:) =  block(b) % Q    (i,j,k,:) &
                                     -  block(b) % dt   (i,j,k)   &
                                     *  block(b) % Area (i,j,k)   &
                                     *  block(b) % Res  (i,j,k,:)
            end do
         end do
      end do
   end do
end subroutine update_sol
