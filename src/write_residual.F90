subroutine write_residual()
use control
use global
implicit none
if (res_out) then
   write(*,*) iteration
end if
end subroutine write_residual
