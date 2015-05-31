subroutine set_boundary()
   use global
   use control
   implicit none
   integer :: b,f,var
   integer :: i,j,k
   integer :: i1,j1,k1
   integer :: i2,j2,k2, edge
   integer :: idir,jdir,kdir
   integer :: idir_loop,jdir_loop,kdir_loop
   integer :: i_start,j_start,k_start
   integer :: i_end,j_end,k_end
   real(kind = dp) :: un,u,v


   do b = 1, nBlock
      do f = 1, nFaces
         call loop_indices_for_bc_faces(f,b                          &
                                       ,i_start,i_end,idir,idir_loop &
                                       ,j_start,j_end,jdir,jdir_loop &
                                       ,k_start,k_end,kdir,kdir_loop )

         !Boundary, extrapolate Values
         if (block(b) % face(f) % conType == BC_WALL) then
!         if (block(b) % face(f) % conType < 0 ) then
!!!!! EULER EQUATION
            if (control_equation == 1) then
                do k = k_start,k_end, kdir_loop
                   if (kdir /= 0) then
                      k1 = k_start - (k - k_start - kdir)
                      edge = 3
                      if ( f == 5) then
                        k2 = 1
                      else
                        k2 = block(b) % nPkt(3)
                      end if
                   else
                      k1 = k
                      k2 = k
                   end if
                   do j = j_start,j_end,jdir_loop
                      if (jdir /= 0) then
                         j1 = j_start - (j - j_start - jdir)
                         edge = 2
                         if (f == 3)  then
                            j2 = 1
                         else
                            j2 = block(b) % nPkt(2)
                        end if
                      else
                         j1 = j
                         j2 = j
                      end if
                      do i = i_start,i_end,idir_loop
                         if (idir /= 0) then
                            i1 = i_start - (i - i_start - idir)
                            edge = 1
                            if ( f == 1) then
                                i2 = 1
                             else
                                i2 = block(b) % nPkt(1)
                            end if
                         else
                            i1 = i
                            i2 = i
                         end if
                         block(b) % Q (i,j,k,1) = block(b) % Q (i1,j1,k1,1)
                         u = block(b) % Q (i1,j1,k1,2) ! / block(b) % Q (i,j,k,1)
                         v = block(b) % Q (i1,j1,k1,3) ! / block(b) % Q (i,j,k,1)
                         un = + u * block(b) % edge_vec(1,i2,j2,k2,edge) &
                                  + v * block(b) % edge_vec(2,i2,j2,k2,edge)
                         block(b) % Q (i,j,k,2) = block(b) % Q (i1,j1,k1,2) &
                                        - 2.0E0_dp * un * block(b) % edge_vec(1,i2,j2,k2,edge)

                         block(b) % Q (i,j,k,3) = block(b) % Q (i1,j1,k1,3) &
                                        - 2.0E0_dp * un * block(b) % edge_vec(2,i2,j2,k2,edge)

                         block(b) % Q (i,j,k,4) = block(b) % Q (i1,j1,k1,4)

!                         if ( i == 20 .and. j == 0) then
!
!                         write(*,*) i,j,i1,j1,i2,j2,edge
!                         write(*,*) block(b) % edge_vec(:,i2,j2,k2,edge),un
!                         write(*,*) block(b) % Q (i,j,k,:),sqrt(block(b) % Q (i,j,k,2)**2+block(b) % Q (i,j,k,3)**2) !&
!                         !                 / block(b) % Q (i,j,k,1)
!                         write(*,*) block(b) % Q (i1,j1,k1,:),sqrt(block(b) % Q (i1,j1,k1,2)**2+block(b) % Q (i1,j1,k1,3)**2) !&
!                         !                 / block(b) % Q (i,j,k,1)
!!                         write(*,*) 0.4*(block(b) % Q (i,j,k,4)-0.5*block(b) % Q (i,j,k,1) * &
!!                                 (block(b) % Q (i,j,k,2)**2+block(b) % Q (i,j,k,3)**2)) &
!!                                 ,0.4*(block(b) % Q (i1,j1,k1,4)-0.5*block(b) % Q (i1,j1,k1,1) * &
!!                                 (block(b) % Q (i1,j1,k1,2)**2+block(b) % Q (i1,j1,k1,3)**2))
!
!                         end if
!                         write(*,*) i,un,-un*block(b) % edge_vec(1,i2,j2,k2,edge) &
!                                              ,-un*block(b) % edge_vec(2,i2,j2,k2,edge) &
!                                              ,block(b) % Q (i,j,k,2:3)&
!                                              ,sqrt(block(b) % Q (i,j,k,2)*block(b) % Q (i,j,k,2)&
!                                              + block(b) % Q (i,j,k,3)*block(b) % Q (i,j,k,3))

                      end do
                   end do
                end do

 !!!!! NAVIE STOKES EQUATION
            else if ( control_equation == 2) then
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
            end if !control_equation == 1
         else if (block(b) % face(f) % conType == BC_INFLOW) then

            !!überschall inflow bleibt constant
         else if (block(b) % face(f) % conType == BC_OUTFLOW) then

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
                     block(b) % Q (i,j,k,2:3) = block(b) % Q (i1,j1,k1,2:3)
                     block(b) % Q (i,j,k,4) = block(b) % Q (i1,j1,k1,4)
                  end do
               end do
            end do
         else if (block(b) % face(f) % conType == BC_SYMMETRY) then
!         if (block(b) % face(f) % conType < 0 ) then
            do k = k_start,k_end, kdir_loop
                                  if (kdir /= 0) then
                      k1 = k_start - (k - k_start - kdir)
                      edge = 3
                      if ( f == 5) then
                        k2 = 1
                      else
                        k2 = block(b) % nPkt(3)
                      end if
                   else
                      k1 = k
                      k2 = k
                   end if
                   do j = j_start,j_end,jdir_loop
                      if (jdir /= 0) then
                         j1 = j_start - (j - j_start - jdir)
                         edge = 2
                         if (f == 3)  then
                            j2 = 1
                         else
                            j2 = block(b) % nPkt(2)
                        end if
                      else
                         j1 = j
                         j2 = j
                      end if
                      do i = i_start,i_end,idir_loop
                         if (idir /= 0) then
                            i1 = i_start - (i - i_start - idir)
                            edge = 1
                            if ( f == 1) then
                                i2 = 1
                             else
                                i2 = block(b) % nPkt(1)
                            end if
                         else
                            i1 = i
                            i2 = i
                         end if
                     block(b) % Q (i,j,k,1) = block(b) % Q (i1,j1,k1,1)
                     un = + block(b) % Q (i1,j1,k1,2) * block(b) % edge_vec(1,i2,j2,k2,edge) &
                              + block(b) % Q (i1,j1,k1,3) * block(b) % edge_vec(2,i2,j2,k2,edge)

                     block(b) % Q (i,j,k,2) = block(b) % Q (i1,j1,k1,2) &
                                        + 2.0E0_dp * un * block(b) % edge_vec(1,i2,j2,k2,edge)
                     block(b) % Q (i,j,k,3) = block(b) % Q (i1,j1,k1,3) &
                                        + 2.0E0_dp * un * block(b) % edge_vec(2,i2,j2,k2,edge)
                     block(b) % Q (i,j,k,4) = block(b) % Q (i1,j1,k1,4)

                  end do
               end do
            end do
         else if (block(b) % face(f) % conType < 0 ) then
            write(*,*) "BLOCK:",b,BC_TYPE(-block(b) % face(f) % conType)//"@"//FACES(f)
!            call error ("Randtype noch nicht impementiert: "&
!            //BC_TYPE(-block(b) % face(f) % conType)//"@"//FACES(f),__FILE__,__LINE__)

         !Angrenzender Block, Werte von anderem Block
         else !conType < 0
            call error ("Angrenzende Bloecke noch nicht implementiert",__FILE__,__LINE__)
         end if !conType < 0

      end do !f = 1,nFaces
   end do !b = 1,nBlock

end subroutine set_boundary
