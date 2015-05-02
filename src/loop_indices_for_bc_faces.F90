subroutine loop_indices_for_bc_faces(f,b                          &
                                    ,i_start,i_end,idir,idir_loop &
                                    ,j_start,j_end,jdir,jdir_loop &
                                    ,k_start,k_end,kdir,kdir_loop )
   use control, only : n_BC_Cells
   use global
   implicit none
   integer, intent(in)  :: f
   integer, intent(in)  :: b
   integer, intent(out) :: idir
   integer, intent(out) :: jdir
   integer, intent(out) :: kdir
   integer, intent(out) :: idir_loop
   integer, intent(out) :: jdir_loop
   integer, intent(out) :: kdir_loop
   integer, intent(out) :: i_start
   integer, intent(out) :: j_start
   integer, intent(out) :: k_start
   integer, intent(out) :: i_end
   integer, intent(out) :: j_end
   integer, intent(out) :: k_end

   if (f <= 4) then
      k_start = 1
      k_end = block(b) % nCell(3)
      kdir_loop = 1
      kdir = 0
   else if (f == 5) then
      k_start = 0
      k_end = k_start - n_BC_Cells + 1
      kdir = 1
      kdir_loop = -1
   else if (f == 6) then
      k_start = block(b) % nCell(3) + 1
      k_end = k_start + n_BC_Cells -1
      kdir = -1
      kdir_loop = 1
   end if

   if (f <= 2) then
      j_start = 1
      j_end = block(b) % nCell(2)
      jdir_loop = 1
      jdir = 0
   else if (f == 3) then
      j_start = 0
      j_end = j_start - n_BC_Cells + 1
      jdir = 1
      jdir_loop = -1
   else if (f == 4) then
      j_start = block(b) % nCell(2) + 1
      j_end = j_start + n_BC_Cells -1
      jdir = -1
      jdir_loop = 1
   end if

   if (f >= 3) then
      i_start = 1
      i_end = block(b) % nCell(1)
      idir_loop = 1
      idir = 0
   else if (f == 1) then
      i_start = 0
      i_end = i_start - n_BC_Cells + 1
      idir = 1
      idir_loop = -1
   else if (f == 2) then
      i_start = block(b) % nCell(1) + 1
      i_end = i_start + n_BC_Cells - 1
      idir = -1
      idir_loop = 1
   end if

end subroutine loop_indices_for_bc_faces
