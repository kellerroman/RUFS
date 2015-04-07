module const
implicit none
public
save
integer, parameter :: dp = KIND(1.0D+0)

integer, parameter ::  WEST_FACE          = 1
integer, parameter ::  EAST_FACE          = 2
integer, parameter :: SOUTH_FACE          = 3
integer, parameter :: NORTH_FACE          = 4
integer, parameter ::  BACK_FACE          = 5
integer, parameter :: FRONT_FACE          = 6

integer, parameter ::    BC_WALL          = -1
integer, parameter ::    BC_INFLOW        = -2
integer, parameter ::    BC_OUTFLOW       = -3
integer, parameter ::    BC_SYMMERTIE     = -4

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

end module const
