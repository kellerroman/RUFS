subroutine write_residual()
use control
use global
implicit none

if (res_out) then
   write(*,'(1X,I7,4(1X,ES10.3))') iteration,time,minRes,maxRes,avgRes
end if
end subroutine write_residual
