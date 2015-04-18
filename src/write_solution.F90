subroutine write_solution()
use const
use control
use global

implicit none
integer , parameter :: fu = 99

integer              :: b,var,i,j,k


if (sol_out) then
   call wr("Writing Solution to File",3)

   if (write_sol_header) then
      write_sol_header = .false.
      open( unit = fu , file = trim(file_sol_out)                             &
          , form = "unformatted", access = "stream", status = "replace"       )
      write(fu) io_marker_file_start
      write(fu) io_marker_header_start
      write(fu) io_sol_version
      write(fu) Dimen
      write(fu) sol_out_nVar
      write(fu) nBlock
      write(fu) dp !SIZE OF SOLUTION VARIABLES
      write(fu) io_len_VarName !LEN OF VARIABLE NAMES
      do b = 1, nBlock
         write(fu) block(b) % nCell(1:Dimen)
      end do
      write(fu) VarName
      write(fu) io_marker_header_end

   else
      open( unit = fu , file = trim(file_sol_out)                             &
          , form = "unformatted", access = "stream", status = "old"           &
          , position = "append"                                               )

   end if
   write(fu) io_marker_solution_start
   write(fu) io_marker_solution_header_start
   write(fu) iteration
   write(fu) io_marker_solution_header_end
   do b = 1, nBlock
      write(fu) io_marker_solution_block_start
      do var = 1, sol_out_nVar
         write(fu) io_marker_solution_var_start
         select case(VarName(var))
            case(VarName_Jac)
               write(fu) (((block(b) % Jac(i,j,k),i = 1,block(b) % nCell(1) ) &
                                                 ,j = 1,block(b) % nCell(2) ) &
                                                 ,k = 1,block(b) % nCell(3) )
            case(VarName_SwpX)
               write(fu) (((block(b) % schwerpunkt(i,j,k,1),i = 1,block(b) % nCell(1) ) &
                                                 ,j = 1,block(b) % nCell(2) ) &
                                                 ,k = 1,block(b) % nCell(3) )
            case(VarName_SwpY)
               write(fu) (((block(b) % schwerpunkt(i,j,k,2),i = 1,block(b) % nCell(1) ) &
                                                 ,j = 1,block(b) % nCell(2) ) &
                                                 ,k = 1,block(b) % nCell(3) )
            case(VarName_M2XI)
               write(fu) (((block(b) % metric2(i,j,k,1,1),i = 1,block(b) % nCell(1) ) &
                                                 ,j = 1,block(b) % nCell(2) ) &
                                                 ,k = 1,block(b) % nCell(3) )
            case(VarName_M2XJ)
               write(fu) (((block(b) % metric2(i,j,k,1,2),i = 1,block(b) % nCell(1) ) &
                                                 ,j = 1,block(b) % nCell(2) ) &
                                                 ,k = 1,block(b) % nCell(3) )
            case(VarName_M2YI)
               write(fu) (((block(b) % metric2(i,j,k,2,1),i = 1,block(b) % nCell(1) ) &
                                                 ,j = 1,block(b) % nCell(2) ) &
                                                 ,k = 1,block(b) % nCell(3) )
            case(VarName_M2YJ)
               write(fu) (((block(b) % metric2(i,j,k,2,2),i = 1,block(b) % nCell(1) ) &
                                                 ,j = 1,block(b) % nCell(2) ) &
                                                 ,k = 1,block(b) % nCell(3) )

            case default
               call error("VAR TO WRITE UNKNOWN: "//trim(varname(var)),__FILE__,__LINE__)
         end select
         write(fu) io_marker_solution_var_end
      end do
      write(fu) io_marker_solution_block_end
   end do
   write(fu) io_marker_solution_end

   if (iteration == control_num_iteration) then ! Try to set the end_marker at the last iteration
      write(fu) io_marker_file_end
   end if
   close (fu)
end if
end subroutine write_solution
