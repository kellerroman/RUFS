program gridgen
implicit none

integer, parameter :: imax = 51
integer, parameter :: jmax = 51
integer, parameter :: fu_git = 99
integer, parameter :: Version = 1000
INTEGER, PARAMETER :: ioout = 10

integer, parameter :: Dimen = 2
integer, parameter :: nBlock = 1

integer, parameter :: bc(4) = (/-2,-1,-3,-3/)
real(kind= 8), parameter :: r1 = 1.0D0
real(kind= 8), parameter :: r2 = 2.0D0
real(kind=8), parameter :: a2d = 3.1415927D0 / 180.0D0
real(kind= 8), parameter :: a1 = a2d * -90.0D0
real(kind= 8), parameter :: a2 = a2d *    90.0D0
real(kind=8) :: xyz (imax,jmax,Dimen)
real(kind = 8) :: radius, angle
integer :: i,j,d

write(*,*) "SIMPLE GRID GEN vor Reentry capsula ( circle)"

do i = 1,imax
   radius = r1 + (r2-r1)/(imax-1) * (i-1)
   do j = 1,jmax
      angle = a1 + (a2-a1)/(jmax-1) * (j-1)
      xyz(i,j,1) = - radius * cos ( angle)

      xyz(i,j,2) = radius * sin(angle)
   end do
end do
open (fu_git,file="git.bin",form="UNFORMATTED",access="STREAM",status="replace")

open (ioout,file="bc.bin",form="unformatted",access="stream",status="replace")

write (fu_git) Version,Dimen,nBlock
write (ioout) Version,Dimen,nBlock

i = 3*4
write(ioout) i

write (fu_git) imax,jmax
write (ioout)  bc
write (fu_git) (((xyz(i,j,d),d= 1,Dimen) &
                          ,i= 1, imax) &
                          ,j= 1, jmax) 

close (fu_git)
close (ioout)
write(*,*) "done"

end program
