subroutine allocate_vars()
use control
use global
implicit none
integer :: b
do b = 1,nBlock
   allocate ( block(b) % schwerpunkt(1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells &
                                    ,Dimen))

   allocate ( block(b) % metric1    (1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells &
                                    ,Dimen,Dimen))

   allocate ( block(b) % metric2    (1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells &
                                    ,Dimen,Dimen))

   allocate ( block(b) % Jac        (1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells ))

   allocate ( block(b) % JacI       (1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells ))

   allocate ( block(b) % Edge_Len   (1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells &
                                    ,Dimen))

   allocate ( block(b) % Edge_Vec   (1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells &
                                    ,Dimen,Dimen))

   allocate ( block(b) % Area       (1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells ))

   allocate ( block(b) % Q          (1-n_BC_Cells : block(b)%nCell(1)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(2)+n_BC_Cells &
                                    ,1-n_BC_Cells : block(b)%nCell(3)+n_BC_Cells &
                                    ,nVar))

   allocate ( block(b) % Res        (1            : block(b)%nCell(1)             &
                                    ,1            : block(b)%nCell(2)             &
                                    ,1            : block(b)%nCell(3)             &
                                    ,nVar))

   allocate ( block(b) % Flux       (1            : block(b)%nPkt(1)             &
                                    ,1            : block(b)%nPkt(2)             &
                                    ,1            : block(b)%nPkt(3)             &
                                    ,nVar,Dimen))

   allocate ( block(b) % dt         (1            : block(b)%nCell(1)             &
                                    ,1            : block(b)%nCell(2)             &
                                    ,1            : block(b)%nCell(3)             ))



   block(b) % Q(1:block(b)%nCell(1),1:block(b)%nCell(2),1:block(b)%nCell(3),1) = 1.0E0_dp
   block(b) % Q(1:block(b)%nCell(1),1:block(b)%nCell(2),1:block(b)%nCell(3),2) = 0.0E0_dp
   block(b) % Q(1:block(b)%nCell(1),1:block(b)%nCell(2),1:block(b)%nCell(3),3) = 0.0E0_dp
   block(b) % Q(1:block(b)%nCell(1),1:block(b)%nCell(2),1:block(b)%nCell(3),4) = 1.0E0_dp / 0.4E0_dp

   block(b) % Q(block(b)%nCell(1)/2+1 : block(b)%nCell(1),1:block(b)%nCell(2),1:block(b)%nCell(3),1) = 0.125E0_dp
   block(b) % Q(block(b)%nCell(1)/2+1 : block(b)%nCell(1),1:block(b)%nCell(2),1:block(b)%nCell(3),4) = 0.1E0_dp / 0.4E0_dp
end do

end subroutine allocate_vars
