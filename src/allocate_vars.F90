subroutine allocate_vars()
use control
use global
implicit none
integer :: b,j
real(kind=dp) :: p,gamma,rho,u,v,gm1
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

   allocate ( block(b) % Edge_Len   (1            : block(b)%nPkt (1)            &
                                    ,1            : block(b)%nPkt (2)            &
                                    ,1            : block(b)%nPkt (3)            &
                                    ,Dimen))

   allocate ( block(b) % Edge_Vec   (Dimen                                       &
                                    ,1            : block(b)%nPkt (1)            &
                                    ,1            : block(b)%nPkt (2)            &
                                    ,1            : block(b)%nPkt (3)            &
                                    ,Dimen))

   allocate ( block(b) % Area       (1            : block(b)%nCell(1)            &
                                    ,1            : block(b)%nCell(2)            &
                                    ,1            : block(b)%nCell(3)            ))

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

   allocate ( block(b) % Q_rec (1 : block(b)%nPkt(1) &
                                                   ,1 : block(b)%nPkt(2) &
                                                   ,1 : block(b)%nPkt(3) &
                                                   ,nVar,Dimen,2))


   block(b) % Q(:,:,:,1) = 1.0E0_dp
   block(b) % Q(:,:,:,2) = 0.0E0_dp
   block(b) % Q(:,:,:,3) = 0.0E0_dp
   block(b) % Q(:,:,:,4) = 1.0E0_dp / 0.4E0_dp
if ( block(b)%nCell(1) > block(b)%nCell(2)) then
!!!!!! x - TUBE
   block(b) % Q(block(b)%nCell(1)/2+1 : ,:,:,1) = 0.125E0_dp
   block(b) % Q(block(b)%nCell(1)/2+1 : ,:,:,4) = 0.1E0_dp / 0.4E0_dp
else if ( block(b)%nCell(1) < block(b)%nCell(2)) then
!!!!!!! y - TUBE
   block(b) % Q(:, block(b)%nCell(2)/2+1:, :, 1) = 0.125E0_dp
   block(b) % Q(:, block(b)%nCell(2)/2+1:, :, 4) = 0.1E0_dp / 0.4E0_dp
else
!!!!!!! dia - TUBE
    do j = 1, block(b)%nPkt(2)
       block(b) % Q(block(b)%nPkt(1) - j : block(b)%nPkt(1),j,:,1) = 0.125E0_dp
       block(b) % Q(block(b)%nPkt(1) - j : block(b)%nPkt(1),j,:,4) = 0.1E0_dp / 0.4E0_dp
   end do
end if

   Gamma  = 1.4E0_dp
   gm1 = gamma - 1.0E0_dp
   rho = 4.0E0_dp
   u = 1.0E0_dp
   v = 0.0E0_dp
   p = 1/ gamma
   block(b) % Q(:,:,:,1) = rho
   block(b) % Q(:,:,:,2) = u * rho
   block(b) % Q(:,:,:,3) = v * rho
   block(b) % Q(:,:,:,4) = p / gm1 + 0.5E0_dp*rho*(u*u+v*v)

   !(1.E0_dp/(1.4E0_dp*0.4E0_dp)+2.E0_dp)

!   block(b) % Q(10,10,1,2) = 2.E0_dp
!   block(b) % Q(10,10,1, 4) = (1.E0_dp/(1.4E0_dp*0.4E0_dp)+4.E0_dp)
end do

end subroutine allocate_vars
