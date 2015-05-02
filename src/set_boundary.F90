subroutine set_boundary()
   use global
   implicit none
   integer :: b,f,var
   integer :: i,j,k
   integer :: i1,j1,k1
   integer :: idir,jdir,kdir
   integer :: idir_loop,jdir_loop,kdir_loop
   integer :: i_start,j_start,k_start
   integer :: i_end,j_end,k_end

   do b = 1, nBlock
      do f = 1, nFaces
         call loop_indices_for_bc_faces(f,b                          &
                                       ,i_start,i_end,idir,idir_loop &
                                       ,j_start,j_end,jdir,jdir_loop &
                                       ,k_start,k_end,kdir,kdir_loop )

         !Boundary, extrapolate Values
!         if (block(b) % face(f) % conType == BC_WALL) then
         if (block(b) % face(f) % conType < 0 ) then
            do k = k_start,k_end, kdir_loop
               if (kdir /= 0) then
                  k1 = k_start - (k - k_start - kdir)
               else
                  k1 = k
               end if
               do j = j_start,j_end,jdir_loop
                  if (jdir /= 0) then
                     j1 = j_start - (j - j_start - jdir)
                  else
                     j1 = j
                  end if
                  do i = i_start,i_end,idir_loop
                     if (idir /= 0) then
                        i1 = i_start - (i - i_start - idir)
                     else
                        i1 = i
                     end if
                     block(b) % Q (i,j,k,1) = block(b) % Q (i1,j1,k1,1)
                     block(b) % Q (i,j,k,2:3) = - block(b) % Q (i1,j1,k1,2:3)
                     block(b) % Q (i,j,k,4) = block(b) % Q (i1,j1,k1,4)




                  end do
               end do
            end do
         !Angrenzender Block, Werte von anderem Block
         else !conType < 0
            call error ("Angrenzende Bloecke noch nicht implementiert",__FILE__,__LINE__)
         end if !conType < 0

      end do !f = 1,nFaces
   end do !b = 1,nBlock

end subroutine set_boundary
