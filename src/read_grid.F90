subroutine read_grid()
!use const, only: dp
use global
use control
implicit none
integer, parameter :: fu_git = 99
integer, parameter :: Version = 1000
integer :: fileVersion
integer :: b, i ,j ,k ,d
logical :: fexists
inquire(file=trim(file_git_in),exist=fexists)
if(.not. fexists) then
  call error("GITTER INPUT DATEI konnte nicht gefunden werden: "//TRIM(file_git_in),__FILE__,__LINE__)
end if

open(fu_git,file=trim(file_git_in),form="UNFORMATTED",access="STREAM",status="OLD")

read(fu_git) fileVersion,Dimen,nBlock

nFaces = Dimen * 2
nCorners = 2**Dimen

if (fileVersion /= Version) then
call error("Gitterdatei Version stimmt nicht",__FILE__,__LINE__)
end if

allocate(block(nBlock))

do b = 1,nBlock
   block(b) % nPkt = 1
   block(b) % nCell = 1
   block(b) % nOut = 1
   read(fu_git) block(b) % nPkt(1:Dimen)
    
   block(b) % nCell(1:Dimen) = block(b) % nPkt(1:Dimen)-1

   allocate (block(b) % xyz(1-n_BC_Cells : block(b)%nPkt(1)+n_BC_Cells &
                           ,1-n_BC_Cells : block(b)%nPkt(2)+n_BC_Cells &
                           ,1-n_BC_Cells : block(b)%nPkt(3)+n_BC_Cells &
                           ,Dimen))

end do !B= 1,NB  SCHLEIFE ÜBER BLÖCKE GIT

do b = 1, nBlock
   read(fu_git) ((((block(b) % xyz(i,j,k,d),d=1,Dimen) &
                                           ,i= 1,block(b) % nPkt(1)) &
                                           ,j= 1,block(b) % nPkt(2)) &
                                           ,k= 1,block(b) % nPkt(3))
end do

close ( fu_git)

end subroutine read_grid
