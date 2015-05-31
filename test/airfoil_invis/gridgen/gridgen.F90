program gridgen
implicit none

integer, parameter :: imax = 71
integer, parameter :: jmax = 48
integer, parameter :: fu_git = 99
integer, parameter :: Version = 1000
INTEGER, PARAMETER :: ioout = 10

integer, parameter :: Dimen = 2
integer, parameter :: nBlock = 1

integer, parameter :: bc(4) = (/-2,-3,-1,-4/)

real(kind=8) :: xyz (imax,jmax,2)

integer :: i,j,d

write(*,*) "AIRFOIL converter"

open(unit=2,file="fort.1.fine")
do i = 1,imax
    do j = 1, jmax
        read(2,1000) xyz(i,j,:)
    end do
    read(2,2000)
end do
close(2)

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

1000 format(1x,e11.4,5x,e11.4)
2000 format (1x)
end program
