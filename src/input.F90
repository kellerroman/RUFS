subroutine input()
use global
use const
use control
implicit none
integer :: b,f


call input_control()

call read_grid()

call read_bc()

call wr("DATA IN",2)
IF(Dimen == 3) THEN
  WRITE(*,'("3D SIMULATION")')
ELSE IF (Dimen == 2) THEN
  WRITE(*,'("2D SIMULATION")')
else
  WRITE(*,'("1D SIMULATION")')
END IF
write(*,'("BLOCKS ON FILE:         ",I0)') nBlock
write(*,'(4(A5,1X),6(A2,1X))') "BLOCK","NCI","NCJ","NCK",FACES(1:nFaces)
do b = 1, nBlock
   write(*,'(4(I5,1X))',ADVANCE="NO") b, block(b) % nCell(:)
   do f = 1 , nFaces
      if (block(b) % face(f) % conType > 0) THEN
         write(*,'(I2,1X)',ADVANCE="NO") block(b) % face(f) % conType
      else if (block(b) % face(f) % conType == 0) THEN
         write(*,'(A2,1X)',ADVANCE="NO") "-"
      else
         write(*,'(A2,1X)',ADVANCE="NO") BC_TYPE(-block(b) % face(f) % conType)
      end if
   end do
   write(*,*)
end do

end subroutine input
