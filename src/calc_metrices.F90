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
               block(b) % Edge_Vec(1,i,j,k,2) = ( block(b) % xyz(i+1,j  ,k  ,2) - block(b) % xyz(i  ,j  ,k  ,2)) &
                                              / block(b) % Edge_Len(i,j,k,2)
               block(b) % Edge_Vec(2,i,j,k,2) = - ( block(b) % xyz(i+1,j  ,k  ,1) - block(b) % xyz(i  ,j  ,k  ,1))&
                                              / block(b) % Edge_Len(i,j,k,2)

            end do
         end do
      end do



      
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
      if (n_BC_Cells > 1) call error ("n_BC_CELLS >1 für ECKPUNKTE NICHT IMPLEMENTIERT",__FILE__,__LINE__)

      k = 1
      do f = 1, nCorners
         if (f == 1 .or. f == 4 .or. f == 5 .or. f == 7) then
            i = 0
            idir = 1
         else
            i = block(b)%nCell(1)+1 !n_BC_Cells
            idir = -1
         end if

         if (f == 1 .or. f == 2 .or. f == 5 .or. f == 6) then
            j = 0
         else
            j = block(b)%nCell(2)+1 !n_BC_Cells
         end if
         block(b) % schwerpunkt(i,j,k,:) = 2.0D0 &
                                         * block(b) % schwerpunkt(i+idir  ,j,k,:) &
                                         - block(b) % schwerpunkt(i+idir*2,j,k,:)

      end do !f = 1, nCorners




      !do k = 1-n_BC_Cells, block(b)%nCell(3)+n_BC_Cells
       k = 1
         nlv(3) = k
         do j = 1-n_BC_Cells, block(b)%nCell(2)+n_BC_Cells
            nlv(2) = j
            do i = 1-n_BC_Cells, block(b)%nCell(1)+n_BC_Cells
               nlv(1) = i
               do d1 = 1, Dimen ! physikalische Dimensionen x,y,z
                  do d2 = 1,Dimen  ! Gitter dimensionen i,j,k eta xi phi
                     ndir = 0
                     ndir(d2) = 1
                     if (nlv(d2) == 1-n_BC_Cells) then

                        block(b) % metric1 (i,j,k,d1,d2) = 0.5E0_dp * ( &
                         - 3.0E0_dp * block(b) % schwerpunkt(i          ,j          ,k          ,d1) &
                         + 4.0E0_dp * block(b) % schwerpunkt(i+ndir(1)  ,j+ndir(2)  ,k+ndir(3)  ,d1) &
                         -            block(b) % schwerpunkt(i+ndir(1)*2,j+ndir(2)*2,k+ndir(3)*2,d1) )

                     else if (nlv(d2) == block(b)%nCell(d2)+n_BC_Cells) then
                        block(b) % metric1 (i,j,k,d1,d2) = 0.5E0_dp * ( &
                         + 3.0E0_dp * block(b) % schwerpunkt(i          ,j          ,k          ,d1) &
                         - 4.0E0_dp * block(b) % schwerpunkt(i-ndir(1)  ,j-ndir(2)  ,k-ndir(3)  ,d1) &
                         +            block(b) % schwerpunkt(i-ndir(1)*2,j-ndir(2)*2,k-ndir(3)*2,d1) )
                     else
                        block(b) % metric1 (i,j,k,d1,d2) = 0.5E0_dp * ( &
                           block(b) % schwerpunkt(i+ndir(1),j+ndir(2),k+ndir(3),d1) &
                         - block(b) % schwerpunkt(i-ndir(1),j-ndir(2),k-ndir(3),d1) )
                     end if
                  end do !d2 = 1,Dimen
               end do !d1 = 1, Dimen

               block(b) % JacI(i,j,k) = block(b) % metric1(i,j,k,1,1) &
                                      * block(b) % metric1(i,j,k,2,2) &
                                      - block(b) % metric1(i,j,k,1,2) &
                                      * block(b) % metric1(i,j,k,2,1)

               block(b) % Jac(i,j,k) = 1.0E0_dp / block(b) % JacI(i,j,k)

               do d1 = 1, Dimen
                  do d2 = 1,Dimen
                     block(b) % metric2(i,j,k,d1,d2) = block(b) % Jac(i,j,k) &
                                 * block(b) % metric1(i,j,k,3-d1,3-d2)
                     if (d2 /= d1) &
                        block(b) % metric2(i,j,k,d1,d2) = - block(b) % metric2(i,j,k,d1,d2)
                  end do
               end do


            end do !i = 1-n_BC_Cells, block(b)%nCell(1)+n_BC_Cells
         end do !j = 1-n_BC_Cells, block(b)%nCell(2)+n_BC_Cells
   !   end do !k = 1-n_BC_Cells, block(b)%nCell(3)+n_BC_Cells
   end do !b = 1,nBlock

end subroutine calc_metrices
