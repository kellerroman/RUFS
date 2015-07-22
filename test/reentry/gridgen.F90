program gridgen
implicit none

integer, parameter :: imax = 71
integer, parameter :: jmax = 51
integer, parameter :: fu_git = 99
integer, parameter :: Version = 1000
INTEGER, PARAMETER :: ioout = 10

integer, parameter :: Dimen = 2
integer, parameter :: nBlock = 1

integer, parameter :: bc(4) = (/-2,-1,-3,-3/)
real(kind= 8), parameter :: r1 = 0.5D0
real(kind= 8), parameter :: r2 = 1.5D0
real(kind= 8), parameter :: r3 = 3.5D0
real(kind=8), parameter :: a2d = 3.1415927D0 / 180.0D0
real(kind= 8), parameter :: a1 = a2d * (- 90.0D0)
real(kind= 8), parameter :: a2 = a2d *    90.0D0
real(kind=8) :: xyz (imax,jmax,Dimen)
real(kind = 8) :: radius, angle,radius2
integer :: i,j,d

write(*,*) "SIMPLE GRID GEN vor Reentry capsula ( circle)"

do i = 1,imax

   do j = 1,jmax
      angle = a1 + (a2-a1)/(jmax-1) * (j-1)
      radius = r2 + (r3-r2) * (1.0D0-(a2-abs(angle))/a2)**2
      radius2 =  radius + (r1-radius) * (i-1)/(imax-1)
!      write(*,*) 1.0D0-(a2-abs(angle))/a2
!      write(*,*) angle,( (abs(angle)-a2) / a2 +3.0D0) / 2.0D0
      xyz(i,j,1) = - radius2 * cos ( angle)

      xyz(i,j,2) = radius2 * sin(angle)
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
