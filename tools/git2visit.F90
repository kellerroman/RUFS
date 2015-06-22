program git2visit

   use const
   use types
   implicit none

include "silo_f9x.inc"
   integer, parameter :: fu = 99
   integer, parameter :: fg = 101
   integer, parameter :: git_Version = 1000
   integer :: arg_count
   character(len = 100) :: arg
   character(len = 100) :: git_file
   character(len = 100) :: out_file
   integer :: sol_file_version,Dimen,nVar,nBlock,Precision,len_VarName
   integer :: rg_Version,rg_Dimen,rg_nBlock
   type(tblock), allocatable :: block(:)

   integer :: status
   integer :: nCell

   integer :: nBCC_out(3,2)
   !< Gibt die Ausgeschriebenen Randzellen an
   integer :: b,var,i,j,k,d

!!!!! SILO
   integer :: dbfile, ierr, err, optlistid

   character(len = 8), allocatable ::  blocknames(:)
   integer, allocatable ::             lblocknames(:)
   integer, allocatable ::             meshtypes(:)

   git_file = "git.bin"
   out_file = "git.silo"

   call wr("CONVERT GIT 2 VisIt",1)

   open( unit = fg , file = trim(git_file)                             &
       , form = "unformatted", access = "stream", status = "old")

   read (fg) rg_Version,Dimen,nBlock

   if (rg_Version /= git_Version) then
      write(*,*) "ERROR: GIT Version wrong!",rg_Version,git_Version
      STOP 1
   end if

   if (Dimen == 1) then
      write(*,'(10X,A)') "1D Calculation"
   else if (Dimen == 2) then
      write(*,'(10X,A)') "2D Calculation"
   else if (Dimen == 3) then
      write(*,'(10X,A)') "3D Calculation"
   else
      write(*,*) "ERROR: Dimen not Supported",Dimen
      STOP
   end if

   allocate(block(nBlock))
   allocate (blocknames  (nBlock))
   allocate (lblocknames (nBlock))
   allocate (meshtypes   (nBlock))
   lblocknames = 9
   meshtypes = DB_QUAD_CURV
   do b = 1, nBlock
      write(blocknames(b),'("block",I3.3)') b
      lblocknames(b) = len_trim(blocknames(b))
   end do

   call wr("BLOCKS",3)
   write(*,'(A9,1X,A3,1X,3(A4,1X))') "","#", "I","J","K"
   nCell = 0
   do b = 1,nBlock
      block(b) % nCell = 1
      block(b) % nPkt = 1
      read(fg) block(b) % nPkt(1:Dimen)
      block(b) % nCell(1:Dimen) = block(b) % nPkt(1:Dimen) -1
      write(*,'(10X,I3,1X,3(I4,1X))') b,block(b) % nCell(1:Dimen)


      nCell = nCell + product(block(b) % nCell)
      allocate ( block(b) % xyz(block(b) % nPkt(1),block(b) % nPkt(2),block(b) % nPkt(3),Dimen))
   end do
   write(*,'(10X,A,1X,I0)') "Total Number of Cells:" ,nCell

   do b = 1,nBlock
      read(fg) ((((block(b) % xyz(i,j,k,d),d=1,Dimen) &
                                           ,i= 1,block(b) % nPkt(1)) &
                                           ,j= 1,block(b) % nPkt(2)) &
                                           ,k= 1,block(b) % nPkt(3))
   end do

      ierr = dbcreate(trim(out_file), len_trim(out_file)   &
                     , DB_CLOBBER, DB_LOCAL                &
                     ,"Comment about the data", 22, DB_HDF5, dbfile)
      if(dbfile.eq.-1) then
         write(*,*) "Could not create Silo file"
         STOP 1
      endif


      do b = 1, nBlock
         write(*,*) len_trim(blocknames(b))
         err = dbputqm (dbfile, trim(blocknames(b)), len_trim(blocknames(b)) &
                       , "xc", 2,"yc", 2, "zc", 2 &
                       ,block(b) % xyz(:,:,:,1) &
                       ,block(b) % xyz(:,:,:,2) &
                       ,DB_F77NULL, block(b) % nPKT(1:Dimen), Dimen&
                       ,DB_DOUBLE, DB_NONCOLLINEAR, DB_F77NULL, ierr)
      end do
      ! Set the maximum string length to 20 since that's how long our strings are
      err = dbset2dstrlen(8)
      err = dbputmmesh(dbfile, "quadmesh", 8, nBlock, blocknames, lblocknames, meshtypes, DB_F77NULL, ierr)
      ierr = dbclose(dbfile)

   close(fg)
   write(*,*)
   call wr("CONVERT GIT done",1)
end program git2visit


