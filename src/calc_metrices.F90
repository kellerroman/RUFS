subroutine calc_metrices()
   use global
   use const
   use control
   implicit none
   integer :: i,j,k,b,f,d1,d2
   integer :: idir,jdir,kdir
   integer :: idir_loop,jdir_loop,kdir_loop
   integer :: i_start,j_start,k_start
   integer :: i_end,j_end,k_end
   real(kind=dp) :: rt1,rt2
   integer :: ndir(3), nlv(3)

   do b = 1,nBlock
      do k = 1,block(b) % nCell(3)
         do j = 1,block(b) % nCell(2)
            do i = 1,block(b) % nCell(1)
               block(b) % schwerpunkt(i,j,k,:) =  ( block(b) % xyz(i  ,j  ,k  ,:) &
                                                  + block(b) % xyz(i+1,j  ,k  ,:) &
                                                  + block(b) % xyz(i  ,j+1,k  ,:) &
                                                  + block(b) % xyz(i+1,j+1,k  ,:) &
                                                  ) * 0.25E0_dp

               block(b) % Area(i,j,k) = 2.0E0_dp / abs (                                                 &
                                      + ( block(b) % xyz(i  ,j  ,k  ,2) - block(b) % xyz(i+1,j+1,k  ,2)) &
                                      * ( block(b) % xyz(i  ,j+1,k  ,1) - block(b) % xyz(i+1,j  ,k  ,1)) &
                                      + ( block(b) % xyz(i+1,j  ,k  ,2) - block(b) % xyz(i  ,j+1,k  ,2)) &
                                      * ( block(b) % xyz(i  ,j  ,k  ,1) - block(b) % xyz(i+1,j+1,k  ,1)) )

            end do
         end do
      end do

      do k = 1,block(b) % nCell(3)
         do j = 1,block(b) % nCell(2)
            do i = 1,block(b) % nPkt(1)
               block(b) % Edge_Len(i,j,k,1) = sqrt ( &
                                              ( block(b) % xyz(i  ,j+1,k  ,1) - block(b) % xyz(i  ,j  ,k  ,1)) &
                                            * ( block(b) % xyz(i  ,j+1,k  ,1) - block(b) % xyz(i  ,j  ,k  ,1)) &
                                            + ( block(b) % xyz(i  ,j+1,k  ,2) - block(b) % xyz(i  ,j  ,k  ,2)) &
                                            * ( block(b) % xyz(i  ,j+1,k  ,2) - block(b) % xyz(i  ,j  ,k  ,2)) )
               ! Normalenvector ist dx = y2-y1
               block(b) % Edge_Vec(1,i,j,k,1) = ( block(b) % xyz(i  ,j+1,k  ,2) - block(b) % xyz(i  ,j  ,k  ,2)) &
                                              / block(b) % Edge_Len(i,j,k,1)
               block(b) % Edge_Vec(2,i,j,k,1) = - ( block(b) % xyz(i  ,j+1,k  ,1) - block(b) % xyz(i  ,j  ,k  ,1))&
                                              / block(b) % Edge_Len(i,j,k,1)
!               if (i == 24) then
!                  write(*,*) j,block(b) % Edge_Vec(:,i,j,k,1),block(b) % xyz(i  ,j+1,k  ,2) , block(b) % xyz(i  ,j  ,k  ,2)&
!                  ,block(b) % Edge_Len(i,j,k,1)
!                  end if
            end do
         end do
      end do
      do k = 1,block(b) % nCell(3)
         do j = 1,block(b) % nPkt(2)
            do i = 1,block(b) % nCell(1)
               block(b) % Edge_Len(i,j,k,2) = sqrt ( &
                                              ( block(b) % xyz(i+1,j  ,k  ,1) - block(b) % xyz(i  ,j  ,k  ,1)) &
                                            * ( block(b) % xyz(i+1,j  ,k  ,1) - block(b) % xyz(i  ,j  ,k  ,1)) &
                                            + ( block(b) % xyz(i+1,j  ,k  ,2) - block(b) % xyz(i  ,j  ,k  ,2)) &
                                            * ( block(b) % xyz(i+1,j  ,k  ,2) - block(b) % xyz(i  ,j  ,k  ,2)) )
               ! Normalenvector ist dx = y2-y1
               block(b) % Edge_Vec(1,i,j,k,2) = - ( block(b) % xyz(i+1,j  ,k  ,2) - block(b) % xyz(i  ,j  ,k  ,2)) &
                                              / block(b) % Edge_Len(i,j,k,2)
               block(b) % Edge_Vec(2,i,j,k,2) = ( block(b) % xyz(i+1,j  ,k  ,1) - block(b) % xyz(i  ,j  ,k  ,1))&
                                              / block(b) % Edge_Len(i,j,k,2)

            end do
         end do
      end do
      do d1 = 1,dimen
         do k = 1,block(b) % nCell(3)
            do j = 1,block(b) % nCell(2)
               do i = 1,block(b) % nCell(1)

                  rt1 =       max( block(b) % xyz(i  ,j  ,k,d1) &
                                 , block(b) % xyz(i+1,j  ,k,d1) &
                                 , block(b) % xyz(i  ,j+1,k,d1) &
                                 , block(b) % xyz(i+1,j+1,k,d1) )
                  rt2 =       min( block(b) % xyz(i  ,j  ,k,d1) &
                                 , block(b) % xyz(i+1,j  ,k,d1) &
                                 , block(b) % xyz(i  ,j+1,k,d1) &
                                 , block(b) % xyz(i+1,j+1,k,d1) )
                  block(b) % Len_dt(d1,i,j,k) = rt1-rt2
               end do
            end do
         end do
      end do
!      do k = 1,block(b) % nCell(3)
!         do j = 1,block(b) % nCell(2)
!            do i = 1,block(b) % nCell(1)
!               write(*,'(3(I4,1X),15(F8.5,1X))') i,j,k &
!               ,block(b) % schwerpunkt(i,j,k,:) &
!               ,1.0E0_dp / block(b) % Area(i,j,k) &
!               ,block(b) % Edge_Len(i,j,k,:) &
!               ,block(b) % Edge_Vec(:,i,j,k,1)&
!               ,block(b) % Edge_Vec(:,i+1,j,k,1)&
!               ,block(b) % Edge_Vec(:,i,j,k,2)&
!               ,block(b) % Edge_Vec(:,i,j+1,k,2)
!            end do
!!         stop
!         end do
!      end do
!      stop

!      i =50
!      j= 2
!      k = 1
!      write(*,*) block(b) % edge_len(i-1,j,k,1) ,block(b) % edge_len(i,j,k,1) &
!      ,block(b) % edge_len(i,j-1,k,2) ,block(b) % edge_len(i,j,k,2)
!      write(*,*) block(b) % edge_Vec(:,i-1,j,k,1) ,block(b) % edge_Vec(:,i,j,k,1) &
!      ,block(b) % edge_Vec(:,i,j-1,k,2) ,block(b) % edge_Vec(:,i,j,k,2)
!      stop

      do f = 1,nFaces

         call loop_indices_for_bc_faces(f,b                          &
                                       ,i_start,i_end,idir,idir_loop &
                                       ,j_start,j_end,jdir,jdir_loop &
                                       ,k_start,k_end,kdir,kdir_loop )

         !Boundary, extrapolate Schwerpunkt
         if (block(b) % face(f) % conType < 0) then
            do k = k_start,k_end, kdir_loop
               do j = j_start,j_end,jdir_loop
                  do i = i_start,i_end,idir_loop

                     block(b) % schwerpunkt(i,j,k,:) = 2.0D0 &
                                                     * block(b) % schwerpunkt(i+idir,j+jdir,k+kdir,:) &
                                                     - block(b) % schwerpunkt(i+idir*2,j+jdir*2,k+kdir*2,:)
                  end do
               end do
            end do
         !Angrenzender Block, Schwerpunkt von anderem Block
         else !conType < 0
            call error ("Angrenzende Bloecke noch nicht implementiert",__FILE__,__LINE__)
         end if !conType < 0
      end do !f = 1,nFaces


      !!! ECKPUNKTE
      if (Dimen > 2) call error("3D NOCH NICHT IMPLEMENTIERT ECKPUNKTE",__FILE__,__LINE__)
!      if (n_BC_Cells > 1) call error ("n_BC_CELLS >1 für ECKPUNKTE NICHT IMPLEMENTIERT",__FILE__,__LINE__)
!
!      k = 1
!      do f = 1, nCorners
!         if (f == 1 .or. f == 4 .or. f == 5 .or. f == 7) then
!            i_start  = 0
!            i_end    = -n_BC_Cells + 1
!            idir = 1
!         else
!            i_start  = block(b)%nCell(1)+1 !n_BC_Cells
!            i_end    = block(b)%nCell(1) + n_BC_Cells
!            idir = -1
!         end if
!
!         if (f == 1 .or. f == 2 .or. f == 5 .or. f == 6) then
!            j_start  = 0
!            j_end    =
!         else
!            j = block(b)%nCell(2)+1 !n_BC_Cells
!         end if
!
!         block(b) % schwerpunkt(i,j,k,:) = 2.0D0 &
!                                         * block(b) % schwerpunkt(i+idir  ,j,k,:) &
!                                         - block(b) % schwerpunkt(i+idir*2,j,k,:)
!
!      end do !f = 1, nCorners
   end do !b = 1,nBlock

end subroutine calc_metrices
