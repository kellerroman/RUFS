subroutine read_bc()
use global
use control
implicit none
integer, parameter :: fu_bc = 99
integer, parameter :: Version = 1000
integer ::  fileVersion, fileNB,fileDimen

logical :: fexists

integer :: b,f, temp

inquire(file=trim(file_bc_in), exist=fexists)
if (.not. fexists) then
  call error("Boundary Connection File nicht gefunden: "//trim(file_bc_in),__FILE__,__LINE__)
end if

open(fu_bc, file=trim(file_bc_in),form="UNFORMATTED",access="STREAM",status="OLD")
read(fu_bc) fileVersion,fileDimen,fileNB
 
if (fileVersion /= Version) then
   call error("Boundary Connection File Version stimmt nicht",__FILE__,__LINE__)
end if

if (fileNB /= nBlock) then
   call error("Blockanzahl stimmt nicht",__FILE__,__LINE__)
end if

if (fileDimen /= Dimen) then
   call error("Block Connection Dimensionen stimmen nicht",__FILE__,__LINE__)
end if

! reading start byte to read for each processor
! UNNECESSARY FOR NONPARALLEL
read(fu_bc) temp

do b = 1, nBlock
   allocate ( block(b) % face(Dimen * 2) )
end do

read(fu_bc) ((block(b)% face(f) % conType, f = 1, Dimen * 2),b = 1, nBlock)

close(fu_bc)
end subroutine read_bc
