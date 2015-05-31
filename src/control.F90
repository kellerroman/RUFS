module control
use const
implicit none
save
public
!!! CONFIG VARIABLES FROM INPUT FILE
integer :: space_Order,space_disc

character ( len = len_str_filename ):: file_cfg_in          = "config.bin"
!<
character ( len = len_str_filename ):: file_git_in          = "git.bin"
!<
character ( len = len_str_filename ):: file_bc_in           = "bc.bin"
!<
character ( len = len_str_filename ):: file_sol_out         = "sol.bin"
!<
character ( len = len_str_filename ):: file_ani_out         = "ani.bin"
!<

integer                             :: control_num_iteration
integer                             :: control_sol_out
integer                             :: control_res_out
integer                             :: control_bc_cells_out
!< defines if bc_cells are outputtet to sol file
!< 0: not outputted
!< 1: outputted
integer                             :: control_dimension
!< Gibt die Art der Simulation an
!< 1 : 2D Simulation
!< 2 : 2D Rotationssymmetrisch e Simulation
!< 3 : 3D Simulation
integer                             :: control_riemann_solver
!< Gibt die Riemann Löser an
!< 1 : RHLL
!< 2 : Roe
!< 3 : AUSM
real(kind=dp)                       :: control_CFL
real(kind=dp)                       :: control_timestep
integer                             :: control_equation
!< Gibt das Physikalische Modell an:
!< 1: Euler - Gleichung
!< 2: Navier-Stokes
!< [ 3: Wärmeleitung ]

integer :: control_dt_method
!< Gibt an wie der Zeitschritt berechnet wird:
!< 1: Es wird ein vorgegebener Zeitschritt verwendet
!< 2: Es wird mittels einer CFL Nummer ein Lokaler Zeitschritt in jeder Zelle berechnet
!< 3: Es wird mittels einer CFL Nummer ein Lokaler Zeitschritt in jeder Zelle berechnet,
!<     und alle Zellen verwenden den kleinsten Zeitschritt aller Zellen

!!! CALCULATED VARIABLES 
integer :: n_BC_Cells

integer :: nVar
!< Anzahl der Variablen im Variablenvektor Q

integer :: sol_out_nVar
!< Anzahl der Variablen in der Ausgabedatei

end module control
