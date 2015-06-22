program sol2visit


   use const
   use types
   implicit none

include "silo_f9x.inc"
   integer, parameter :: fu = 99
   integer, parameter :: fg = 101
   integer, parameter :: git_Version = 1000
   integer :: arg_count
   character(len = 100) :: arg
   character(len = 100) :: sol_file
   character(len = 100) :: git_file
   character(len = 100) :: out_file
   character(len = io_marker_len) :: marker
   character(len =io_len_VarName), allocatable :: VarName(:)
   integer :: sol_file_version,Dimen,nVar,nBlock,Precision,len_VarName
   integer :: rg_Version,rg_Dimen,rg_nBlock
   type(tblock), allocatable :: block(:)

   character(len=*), parameter :: marker_error = '(/2(100("!")/),10X,"ERROR@",I0,": ",A," not found!:""",A,""""/2(100("!")/))'

   integer :: iteration
   integer :: num_ausgabe
   !< Anzahl der Ausgabepunkt ider der Datei
   integer :: status
   integer :: nCell
   integer :: num_Pkt(3)
   integer :: nBCC_out(3,2)
   !< Gibt die Ausgeschriebenen Randzellen an
   integer :: b,var,i,j,k,d

!!!!! SILO
   integer :: dbfile, ierr, err, optlistid

   character(len = 8), allocatable ::  blocknames(:)
!   integer, allocatable ::             lblocknames(:)
!   integer, allocatable ::             meshtypes(:)

   sol_file = "sol.bin"
   git_file = "git.bin"
   out_file = "sol.silo"

   call wr("CONVERT SOLUTION 2 VisIt",1)

   open( unit = fg , file = trim(git_file)                             &
       , form = "unformatted", access = "stream", status = "old")

   read (fg) rg_Version,rg_Dimen,rg_nBlock

   if (rg_Version /= git_Version) then
      write(*,*) "ERROR: GIT Version wrong!",rg_Version,git_Version
      STOP 1
   end if

   arg_count=command_argument_count()
   open( unit = fu , file = trim(sol_file)                             &
       , form = "unformatted", access = "stream", status = "old")

   read(fu) marker
   if (marker /= io_marker_file_start) then
      write(*,marker_error) __LINE__,"IO_MARKER_FILE_START",marker
      STOP 1
   end if
   read(fu) marker
   if (marker /= io_marker_header_start) then
      write(*,marker_error) __LINE__,"IO_MARKER_HEADER_START",marker
      STOP 1
   end if

   read(fu) sol_file_version,Dimen,nVar,nBlock,Precision,len_VarName

   if (sol_file_version /= io_sol_version) then
      write(*,*) "ERROR: File Version wrong!",sol_file_version,io_sol_version
      STOP 1
   end if

   if (Precision /= dp ) then
      write(*,*) "ERROR: Real-Precision wrong!",Precision,dp
      STOP 1
   end if

   if (len_VarName /= io_len_VarName) then
      write(*,*) "ERROR: VarName-Length wrong!",len_VarName,io_len_VarName
      STOP 1
   end if
   call wr("HEADER INFO",2)

   if (Dimen /= rg_Dimen) then
      write(*,*) "Grid and Sol Dimension is different",rg_Dimen,Dimen
      stop 1
   end if

   if (rg_nBlock /= nBlock) then
      write(*,*) "Grid and Sol Number of Blocks is different",rg_nBlock,nBlock
      stop 1
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
   allocate(VarName(nVar))
   call wr("BLOCKS",3)
   write(*,'(A9,1X,A3,1X,3(A4,1X))') "","#", "I","J","K"
   nCell = 0
   do b = 1,nBlock
      block(b) % nCell = 1
      block(b) % nPkt = 1
      read(fu) block(b) % nCell(1:Dimen)
      read(fg) block(b) % nPkt(1:Dimen)
      write(*,'(10X,I3,1X,3(I4,1X))') b,block(b) % nCell(1:Dimen)


      nCell = nCell + product(block(b) % nCell)
   end do
   write(*,'(10X,A,1X,I0)') "Total Number of Cells:" ,nCell
   nBCC_out = 0
   do b = 1,nBlock
      read(fu) nBCC_out(1:Dimen,:)
      allocate ( block(b) % Q(block(b) % nCell(1) &
                                                  ,block(b) % nCell(2) &
                                                  ,block(b) % nCell(3), nVar))
      allocate ( block(b) % xyz( block(b) % nPkt(1) + nBCC_out(1,1) + nBCC_out(1,2) &
                                                     , block(b) % nPkt(2) + nBCC_out(2,1) + nBCC_out(2,2) &
                                                     , block(b) % nPkt(3) + nBCC_out(3,1) + nBCC_out(3,2), Dimen))
      do var = 1,Dimen
         if (block(b) % nCell(var)-sum(nBCC_out(var,:)) /= block(b) % nPkt(var)-1 ) then
            write(*,'(A,1X,I0,1X,A,I0,1X,I0)') "Dimensions of Block",b,"are different in Git and Sol:" &
                                             , block(b) % nPkt(var)-1,block(b) % nCell(var)
            stop 1
         end if
      end do

      read(fg) ((((block(b) % xyz(i,j,k,d),d=1,Dimen) &
                                           ,i = 1+nBCC_out(1,1), block(b) % nPkt(1)+nBCC_out(1,1) ) &
                                           ,j = 1+nBCC_out(2,1), block(b) % nPkt(2)+nBCC_out(2,1) ) &
                                           ,k= 1+nBCC_out(3,1), block(b) % nPkt(3)+nBCC_out(3,1) )
   end do
   do b  = 1, nBlock
      do i = 1+nBCC_out(1,1), block(b) % nPkt(1)+nBCC_out(1,1)
         do j = nBCC_out(2,1),1,-1
            block(b) % xyz(i,j,1,:) = 2.0D0 * block(b) % xyz(i,j+1,1,:) - block(b) % xyz(i,j+2,1,:)
!            write(*,*) block(b) % xyz(i,j+2,1,:),block(b) % xyz(i,j+1,1,:)
!            write(*,*) i,j,block(b) % xyz(i,j,1,:)
         end do
         do j = block(b) % nPkt(2)+nBCC_out(2,1)+1 , block(b) % nPkt(2)+nBCC_out(2,1)+nBCC_out(2,2)
            block(b) % xyz(i,j,1,:) = 2.0D0 * block(b) % xyz(i,j -1,1,:) - block(b) % xyz(i,j-2,1,:)
!            write(*,*) i,j,block(b) % xyz(i,j,1,:)
         end do
      end do
      do i = nBCC_out(1,1), 1, -1
         do j =  1, block(b) % nPkt(2)+ nBCC_out(2,1) + nBCC_out(2,2)
            block(b) % xyz(i,j,1,:) = 2 * block(b) % xyz(i+1,j,1,:) - block(b) % xyz(i+2,j,1,:)
!            write(*,*) i,j,block(b) % xyz(i,j,1,:)
         end do
      end do
      do i = block(b) % nPkt(1)+1+nBCC_out(1,1) , block(b) % nPkt(1)+nBCC_out(1,1)+nBCC_out(1,2)
         do j =  1, block(b) % nPkt(2)+ nBCC_out(2,1) + nBCC_out(2,2)
            block(b) % xyz(i,j,1,:) = 2 * block(b) % xyz(i -1,j,1,:) - block(b) % xyz(i-2,j,1,:)
!            write(*,*) i,j,block(b) % xyz(i,j,1,:)
         end do
      end do
   end do

   read(fu) VarName

   call wr("Variable on File",3)
   write(*,'(12X)',advance= "no")
   do var = 1,nVar
      if (var>1) write(*,'(", ")',advance="NO")
      write(*,'(A)',advance="NO") trim(VarName(var))
   end do
   write(*,*)
   read(fu) marker
   if (marker /= io_marker_header_end) then
      write(*,marker_error) __LINE__,"IO_MARKER_HEADER_END",marker
      STOP 1
   end if
call wr("Iterations",2)
!write(*,'(10X)',advance= "no")
num_ausgabe = 0
do
!   read(fu,iostat=status,pos = mypos) marker

   read(fu,iostat=status) marker

   if (status /= 0 .or. marker == io_marker_file_end) then
      if (status /= 0) then
         Write(*,*) "end of file"
      end if
      exit
   else
      num_ausgabe = num_ausgabe + 1
      if (marker /= io_marker_solution_start) then
         write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_START",marker
         STOP 1
      end if
      read(fu) marker
      if (marker /= io_marker_solution_header_start) then
         write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_HEADER_START",marker
         STOP 1
      end if

      read(fu) iteration
!      write(*,'(I0,1X)',advance = "no") iteration

      read(fu) marker
      if (marker /= io_marker_solution_header_end) then
         write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_HEADER_END",marker
         STOP 1
      end if
      write(out_file,'("sol",I6.6,".silo")') iteration
      ierr = dbcreate(trim(out_file), len_trim(out_file)   &
                     , DB_CLOBBER, DB_LOCAL                &
                     ,"Comment about the data", 22, DB_HDF5, dbfile)
      if(dbfile.eq.-1) then
         write(*,*) "Could not create Silo file"
         STOP 1
         endif

      err = dbmkoptlist(1, optlistid)
      err = dbaddiopt(optlistid, DBOPT_CYCLE, iteration)
      err = dbfreeoptlist(optlistid)

      do b = 1, nBlock
         do i = 1, Dimen
            num_pkt(i) = block(b) % nPkt(i) + nBCC_out(i,1) + nBCC_out(i,2)
         end do
         err = dbputqm (dbfile, "quadmesh", 8 &
                       , "xc", 2,"yc", 2, "zc", 2 &
                       ,block(b) % xyz(:,:,:,1) &
                       ,block(b) % xyz(:,:,:,2) &
                       ,DB_F77NULL, num_pkt(1:Dimen), Dimen &
                       ,DB_DOUBLE, DB_NONCOLLINEAR, DB_F77NULL, ierr)
         read(fu) marker
         if (marker /= io_marker_solution_block_start) then
            write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_BLOCK_START",marker
            STOP 1
         end if
         do var = 1, nVar
            read(fu) marker
            if (marker /= io_marker_solution_var_start) then
               write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_VAR_START",marker
               STOP 1
            end if

            read(fu) block(b) % Q(:,:,:,var)

            !!!!! Zeichne SCHWERPUNKT GITTER
            if (var >= 2 ) then
               if (VarName(var-1) == VarName_SwpX .and. &
                   VarName(var)       == VarName_SwpY .and. &
                   Dimen == 2) then
!                   do i = 1,Dimen
!                     cells(i) = block(b) % nCell(i) + nBCC_out(i,1) + nBCC_out(i,2)
!                  end do
                  err = dbputqm (dbfile, "schwerpunkte", 12 &
                                , "xc", 2,"yc", 2, "zc", 2  &
                                ,block(b) % Q(:,:,:,var-1)    &
                                ,block(b) % Q(:,:,:,var)    &
                                ,DB_F77NULL, block(b) % nCell(1:Dimen), Dimen  &
                                ,DB_DOUBLE, DB_NONCOLLINEAR, DB_F77NULL, ierr)
               end if
            end if
!            do i = 1,Dimen
!               cells(i) = block(b) % nCell(i) + nBCC_out(i,1) + nBCC_out(i,2)
!            end do
            err = dbputqv1(dbfile, trim(VarName(var)), len_trim(VarName(var)) &
                          ,"quadmesh", 8, block(b) % Q(:,:,:,var), block(b) % nCell(1:Dimen), Dimen &
                          , DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, ierr)

            read(fu) marker
            if (marker /= io_marker_solution_var_end) then
               write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_VAR_END",marker
               STOP 1
            end if

         end do
         read(fu) marker
         if (marker /= io_marker_solution_block_end) then
            write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_BLOCK_END",marker
            STOP 1
         end if
      end do

      read(fu) marker
      if (status /= 0 .or. marker /= io_marker_solution_end) then
         write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_END",marker
         STOP 1
      end if
      ierr = dbclose(dbfile)
   end if

end do

   close(fu)
   close(fg)
   write(*,*)
   call wr("CONVERT SOLUTION done",1)
end program sol2visit
