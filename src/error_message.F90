subroutine error(text,file,line)
use const
implicit none
character ( len = *) ,intent(in) :: text
character ( len = *) ,intent(in) :: file
integer              ,intent(in) :: line

write(*,'(A,1X,A,1X,A,1X,I0)') "ERROR IN",file,"@",line
write(*,'(A)') text
call end_program(EXIT_WITH_ERROR)
end subroutine error
