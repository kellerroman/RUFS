program extract_1D
   use const
   use types
   implicit none
   integer, parameter :: fu = 99
   integer, parameter :: fo = 101

   integer :: arg_count
   character(len = 100) :: arg
   character(len = 100) :: sol_file
   character(len = 100) :: git_file
   character(len = 100) :: out_file
   character(len = io_marker_len) :: marker
   character(len =io_len_VarName), allocatable :: VarName(:)
   integer :: sol_file_version,Dimen,nVar,nBlock,Precision,len_VarName
   type(tblock), allocatable :: block(:)

   character(len=*), parameter :: marker_error = '(/2(100("!")/),10X,"ERROR@",I0,": ",A," not found!:""",A,""""/2(100("!")/))'

   integer :: iteration
   integer :: num_ausgabe
   !< Anzahl der Ausgabepunkt ider der Datei
   integer :: status
   integer :: nCell

   integer :: nBCC_out(3,2)
   !< Gibt die Ausgeschriebenen Randzellen an
   integer :: b,var,i,j

   integer :: var_dir
   !< Welche Richtung soll variert werden

   integer :: istart,iend
   integer :: jstart,jend

   var_dir = 1
   sol_file = "sol.bin"
   git_file = "git.bin"
   out_file = "extract_1dXXXXXX.dat"
   write(*,*)
   call wr("Extract 1D",1)
    i = 1
    DO
        CALL get_command_argument(i, arg)
        IF (LEN_TRIM(arg) == 0) EXIT
        if (trim(arg) == "x") then
            var_dir = 1
        else if (trim(arg) == "y") then
            var_dir = 2
        else if (trim(arg) == "xy") then
            var_dir = 3
        end if
        i = i+1
    END DO

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
      read(fu) block(b) % nCell(1:Dimen)

      write(*,'(10X,I3,1X,3(I4,1X))') b,block(b) % nCell(1:Dimen)
      nCell = nCell + product(block(b) % nCell)
      allocate ( block(b) % Q(block(b) % nCell(1),block(b) % nCell(2),block(b) % nCell(3),nVar))
   end do
   write(*,'(10X,A,1X,I0)') "Total Number of Cells:" ,nCell

   do b = 1,nBlock
      read(fu) nBCC_out(1:Dimen,:)
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
      write(*,'(I0,1X)',advance = "no") iteration


      read(fu) marker
      if (marker /= io_marker_solution_header_end) then
         write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_HEADER_END",marker
         STOP 1
      end if

      do b = 1, nBlock
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
            read(fu) marker
            if (marker /= io_marker_solution_var_end) then
               write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_VAR_END",marker
               STOP 1
            end if

         end do !var = 1, nVar
         read(fu) marker
         if (marker /= io_marker_solution_block_end) then
            write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_BLOCK_END",marker
            STOP 1
         end if
         write(out_file,'("extract_1d_",I6.6,".dat")') iteration
         open(unit = fo , file = trim(out_file))

         write(fo,'(A)',advance="no") "I,X"
         do var = 1, nVar
            write(fo,'(",",A)',advance="no") trim(varname(var))
         end do
        write(arg,'("(I3,",I0,"("","",F20.12))")') nVar+1
        write(fo,*)
        if (var_dir == 1) then
            j = block(b) %  nCell(2)/2
            do i = 1,block(b) %  nCell(1)
               write(fo,arg) i-nBCC_out(1,1) &
               ,DBLE((i-nBCC_out(1,1)))/DBLE(block(b) %  nCell(1)-nBCC_out(1,1))&
               ,(block(b) % Q(i,j,1,var),var=1,nVar)
            end do
         else if (var_dir == 2) then
            i = block(b) %  nCell(1)/2
            do j = 1,block(b) %  nCell(2)
               write(fo,arg) j-nBCC_out(2,1)&
                ,DBLE((j-nBCC_out(2,1)))/DBLE(block(b) %  nCell(2)-nBCC_out(2,1))&
                ,(block(b) % Q(i,j,1,var),var=1,nVar)
            end do
         else if (var_dir == 3) then
            do j = 1,min(block(b) %  nCell(1),block(b) %  nCell(2))
               i = j
               write(fo,arg) j-nBCC_out(2,1) &
               ,DBLE(j-nBCC_out(2,1))/DBLE(block(b) %  nCell(2)-nBCC_out(2,1)) &
               ,(block(b) % Q(i,j,1,var),var=1,nVar)
            end do
         end if
         close (fo)
      end do !b = 1, nBlock

      read(fu) marker
      if (status /= 0 .or. marker /= io_marker_solution_end) then
         write(*,marker_error) __LINE__,"IO_MARKER_SOLUTION_END",marker
         STOP 1
      end if

   end if

end do

   close(fu)
   write(*,*)
   call wr("SHOW SOLUTION done",1)
end program extract_1D


