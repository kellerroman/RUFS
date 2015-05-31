module const
implicit none
public
save
integer, parameter :: dp = KIND(1.0D+0) ! DOUBLE PRECISION
integer, parameter :: ip = KIND(1)      ! INTEGER PRECISION
integer, parameter :: len_str_filename = 100

integer, parameter ::  WEST_FACE          = 1
integer, parameter ::  EAST_FACE          = 2
integer, parameter :: SOUTH_FACE          = 3
integer, parameter :: NORTH_FACE          = 4
integer, parameter ::  BACK_FACE          = 5
integer, parameter :: FRONT_FACE          = 6

integer, parameter ::    BC_WALL          = -1
integer, parameter ::    BC_INFLOW        = -2
integer, parameter ::    BC_OUTFLOW       = -3
integer, parameter ::    BC_SYMMETRY     = -4

integer, parameter ::    isothermal_Wall  = 1
integer, parameter ::     adiabatic_Wall  = 2
integer, parameter ::      heatflux_Wall  = 3
integer, parameter ::var_isothermal_Wall  = 4
integer, parameter ::  var_heatflux_Wall  = 5

integer, parameter :: EXIT_WITHOUT_ERROR  = 0
integer, parameter :: EXIT_WITH_ERROR     = 1


integer, parameter :: inflow_super_p_t_u  = 1
!< Inflow Supersonic mit räum. konstantem
!< Druck, Temperatur und Geschwindigkeit(Vektor)
integer, parameter :: inflow_sub_t_u      = 2
!< Subsonishcer Inflow mit räuml. konstanter Temp
!< und Geschwindigkeit (Vektor)
integer, parameter :: inflow_sub_t_up     = 3
!< Subsonishcer Inflow mit räuml. konstanter Temp
!< und variabler Geschwindigkeit (Vektor)
integer, parameter :: inflow_sub_mf_t_u   = 4
integer, parameter :: inflow_sub_mf_t_up  = 5

character ( len = 1 ), parameter :: FACES(6) = (/"W","E","S","N","B","F"/)
character ( len = 2 ), parameter :: BC_TYPE(4) = (/"Wa","If","Of","Sy"/)
!<              -1: WALL
!<              -2: INFLOW
!<              -3: OUTFLOW
!<              -4: SYMMETRIE
character ( len = 2 ), parameter :: marker_prefix = "$%"
!character ( len = 4 ), parameter :: marker_suffix = "$"//CHAR(13)//CHAR(11)//CHAR (0) !!! MAKES THE SOL_FILE LITTLE READABLE IN ASCII (CAT)
character ( len = 1 ), parameter :: marker_suffix = "$"

character ( len = 4 ), parameter :: string_Riemann_Solver(3) = (/"RHLL","Roe ","AUSM"/)
character ( len = 9 ), parameter :: string_dt_solver(3) = (/"given dt ","given cfl","calc  dt "/)
integer, parameter :: io_marker_len = len(marker_prefix)+1+len(marker_suffix)

character ( len = io_marker_len ), parameter :: io_marker_file_start           = marker_prefix//"f"//marker_suffix
character ( len = io_marker_len ), parameter :: io_marker_header_start         = marker_prefix//"h"//marker_suffix
character ( len = io_marker_len ), parameter :: io_marker_header_end           = marker_prefix//"H"//marker_suffix

character ( len = io_marker_len ), parameter :: io_marker_solution_start       = marker_prefix//"s"//marker_suffix

character ( len = io_marker_len ), parameter :: io_marker_solution_header_start     = marker_prefix//"i"//marker_suffix
character ( len = io_marker_len ), parameter :: io_marker_solution_header_end       = marker_prefix//"I"//marker_suffix

character ( len = io_marker_len ), parameter :: io_marker_solution_block_start     = marker_prefix//"b"//marker_suffix
character ( len = io_marker_len ), parameter :: io_marker_solution_block_end       = marker_prefix//"B"//marker_suffix

character ( len = io_marker_len ), parameter :: io_marker_solution_var_start     = marker_prefix//"v"//marker_suffix
character ( len = io_marker_len ), parameter :: io_marker_solution_var_end       = marker_prefix//"V"//marker_suffix

character ( len = io_marker_len ), parameter :: io_marker_solution_end         = marker_prefix//"S"//marker_suffix


character ( len = io_marker_len ), parameter :: io_marker_file_end             = marker_prefix//"F"//marker_suffix

integer, parameter               :: io_sol_version = 1

integer, parameter               :: io_len_VarName = 20

character ( len = io_len_Varname ), parameter :: VarName_Rho        = "Dichte"
character ( len = io_len_Varname ), parameter :: VarName_SpU        = "Geschw_U"
character ( len = io_len_Varname ), parameter :: VarName_SpV        = "Geschw_V"
character ( len = io_len_Varname ), parameter :: VarName_Ene        = "Energie"
character ( len = io_len_Varname ), parameter :: VarName_SwpX       = "Schwerpunkt_X"
character ( len = io_len_Varname ), parameter :: VarName_SwpY       = "Schwerpunkt_Y"
character ( len = io_len_Varname ), parameter :: VarName_SwpZ       = "Schwerpunkt_Z"
character ( len = io_len_Varname ), parameter :: VarName_Area       = "Flaeche_Volumen"
character ( len = io_len_Varname ), parameter :: VarName_EdLeN      = "Edge_Length_North"
character ( len = io_len_Varname ), parameter :: VarName_EdLeE      = "Edge_Length_East"
character ( len = io_len_Varname ), parameter :: VarName_EdLeS      = "Edge_Length_South"
character ( len = io_len_Varname ), parameter :: VarName_EdLeW      = "Edge_Length_West"
character ( len = io_len_Varname ), parameter :: VarName_EdAnN      = "Edge_Angle_North"
character ( len = io_len_Varname ), parameter :: VarName_EdAnE      = "Edge_Angle_East"
character ( len = io_len_Varname ), parameter :: VarName_EdAnS      = "Edge_Angle_South"
character ( len = io_len_Varname ), parameter :: VarName_EdAnW      = "Edge_Angle_West"


character ( len = io_len_Varname ), parameter :: VarName_Jac        = "Jacobi"
character ( len = io_len_Varname ), parameter :: VarName_JacI       = "JacobiInv"
character ( len = io_len_Varname ), parameter :: VarName_M1XI       = "M1_x_i"
character ( len = io_len_Varname ), parameter :: VarName_M1XJ       = "M1_x_j"
character ( len = io_len_Varname ), parameter :: VarName_M1YI       = "M1_y_i"
character ( len = io_len_Varname ), parameter :: VarName_M1YJ       = "M1_y_j"
character ( len = io_len_Varname ), parameter :: VarName_M2XI       = "M2_x_i"
character ( len = io_len_Varname ), parameter :: VarName_M2XJ       = "M2_x_j"
character ( len = io_len_Varname ), parameter :: VarName_M2YI       = "M2_y_i"
character ( len = io_len_Varname ), parameter :: VarName_M2YJ       = "M2_y_j"
end module const
