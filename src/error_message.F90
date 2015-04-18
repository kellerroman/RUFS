subroutine error(text,file,line)
use const
implicit none
character ( len = *) ,intent(in) :: text
character ( len = *) ,intent(in) :: file
integer              ,intent(in) :: line

write(*,'(2(/100("!")))')
write(*,'(10("!"),2X,A,1X,A,1X,A,1X,I0)') "ERROR IN",file(8:),"@",line
write(*,'(10("!"),2X,A)') text
write(*,'(100("!")/100("!"))')
call end_program(EXIT_WITH_ERROR)
end subroutine error
