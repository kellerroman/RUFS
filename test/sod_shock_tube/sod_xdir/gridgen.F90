program gridgen
implicit none

integer, parameter :: imax = 101
integer, parameter :: jmax = 3
integer, parameter :: fu_git = 99
integer, parameter :: Version = 1000
INTEGER, PARAMETER :: ioout = 10

integer, parameter :: Dimen = 2
integer, parameter :: nBlock = 1

integer, parameter :: bc(4) = (/-2,-3,-4,-4/)

real(kind=8) :: xyz (imax,jmax,Dimen)

integer :: i,j,d

write(*,*) "SIMPLE GRID GEN"

do i = 1,imax
   do j = 1,jmax
      xyz(i,j,1) = 1.0D0/dble(imax-1) * dble(i-1)

      xyz(i,j,2) = 1.0D0/dble(jmax-1) * dble(j-1)+5.0D-1/dble(imax-1) * dble(i-1)
!      write(*,*) i,j,xyz(i,j,1),xyz(i,j,2)
   end do
end do
do d = 1,2
   write(*,'(A3)',advance = "no") "J\I"
   do i = 1,imax
      write(*,'(I9,1X)',advance= "no") I
   end do
   write(*,*)
   do j = jmax,1,-1
      write(*,'(I3)',advance="no") j
      do i = 1,imax
         write(*,'(1PF10.5)',advance= "no") xyz(i,j,d)
      end do
      write(*,*)
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
