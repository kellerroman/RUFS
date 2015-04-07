subroutine end_program(ExitCode)
use const
implicit none
integer, intent(in) :: ExitCode

!!!! DEALOCATE ARRAYS IF ALLOCATED


if (ExitCode == EXIT_WITHOUT_ERROR) then
   call wr("Program Ended Regularily",1)
else if (ExitCode == EXIT_WITH_ERROR) then
   call wr("ERROR",1)
   stop 1
end if


end subroutine end_program