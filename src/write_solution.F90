subroutine write_solution()
use const
use control
use global

implicit none
integer , parameter :: fu = 99

integer              :: b,var,i,j,k

integer              :: nBCC_out(3,2)
!< ANZAHL DER RANDZELLEN DIE MIT AUSGEGEBEN WERDEN SOLLEN
!< (Gitter Richtung (i,j,k); Anfang(1)/Ende(2)


nBCC_out = 0

if  (control_bc_cells_out == 1) then
   ! Output all boundary values
   nBCC_out(1:Dimen,:) = n_BC_Cells
end if


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
         do var = 1,Dimen
            write(fu) block(b) % nCell(var)+sum(nBCC_out(var,:))
         end do
      end do
      do b = 1,nBlock
         write(fu) nBCC_out(1:Dimen,:)
      end do
      if (sol_out_nVar > 0) then
         write(fu) VarName(1:sol_out_nVar)
      end if
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
            case(VarName_Rho)
               write(fu) (((block(b) % Q(i,j,k,1)           ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_SpU)
               write(fu) (((block(b) % Q(i,j,k,2)           ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_SpV)
               write(fu) (((block(b) % Q(i,j,k,3)           ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_Ene)
               write(fu) (((block(b) % Q(i,j,k,4)           ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_Pre)
               write(fu) ((((0.4E0_dp*(block(b) % Q(i,j,k,4) -0.5E0_dp &
                                                            * (block(b) % Q(i,j,k,2)*block(b) % Q(i,j,k,2)&
                                                            + block(b) % Q(i,j,k,3)*block(b) % Q(i,j,k,3))&
                                                            / block(b) % Q(i,j,k,1) ) )&
                                                             ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_SwpX)
               write(fu) (((block(b) % schwerpunkt(i,j,k,1) ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_SwpY)
               write(fu) (((block(b) % schwerpunkt(i,j,k,2) ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )

            case(VarName_Area)
               write(fu) (((block(b) % area(i,j,k)          ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )

            case(VarName_EdLeW)
               write(fu) (((block(b) % Edge_Len(i  ,j,k,1)  ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_EdLeE)
               write(fu) (((block(b) % Edge_Len(i+1,j,k,1)  ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )

            case(VarName_EdLeS)
               write(fu) (((block(b) % Edge_Len(i,j  ,k,2)  ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_EdLeN)
               write(fu) (((block(b) % Edge_Len(i,j+1,k,2)  ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )

            case(VarName_EdAnW)
               write(fu) (((Normvec_Angle(-block(b) % Edge_Vec(i,j,k,1,:)) &
                                                            ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_EdAnE)
               write(fu) (((Normvec_Angle(block(b) % Edge_Vec(i+1,j,k,1,:)) &
                                                            ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )

            case(VarName_EdAnS)
               write(fu) (((Normvec_Angle(block(b) % Edge_Vec(i,j,k,2,:)) &
                                                            ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case(VarName_EdAnN)
               write(fu) (((Normvec_Angle(-block(b) % Edge_Vec(i,j+1,k,2,:)) &
                                                            ,i = 1-nBCC_out(1,1),block(b) % nCell(1)+nBCC_out(1,2) ) &
                                                            ,j = 1-nBCC_out(2,1),block(b) % nCell(2)+nBCC_out(2,2) ) &
                                                            ,k = 1-nBCC_out(3,1),block(b) % nCell(3)+nBCC_out(3,2) )
            case default
               call error("VAR TO WRITE UNKNOWN: "//trim(varname(var)),__FILE__,__LINE__)
         end select
         write(fu) io_marker_solution_var_end
      end do
      write(fu) io_marker_solution_block_end
   end do
   write(fu) io_marker_solution_end

   if (iteration >= control_num_iteration) then ! Try to set the end_marker at the last iteration
      write(fu) io_marker_file_end
   end if
   close (fu)
end if

contains

function Normvec_Angle(vec) result( winkel)

implicit none
real(kind=dp), intent(in) :: vec(2)
real(kind=dp) :: winkel
real(kind=dp), parameter :: a2d = 180E0_dp / 3.1415927E0_dp
!! NORMIERUNG:
! ( 1, 0) =   0°
! ( 0, 1) =  90°
! (-1, 0) = 180°
! (-1,-1) = 270°
winkel = acos(vec(1)) * a2d
if (vec(2) < 0 ) then
   winkel = 360E0_dp - winkel
end if
end function

end subroutine write_solution
