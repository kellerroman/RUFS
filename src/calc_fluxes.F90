subroutine calc_fluxes()
   use global
   use const
   use control
   implicit none
   real(kind = dp),parameter :: dh = 1.5E0_dp
   real(kind = dp),parameter :: eh = 0.5E0_dp
   real(kind = dp),parameter :: epsi = 1.0E-10_dp
   real(kind = dp) :: met(2,2)
   real(kind = dp) :: Jac
   real(kind = dp) :: sp(2,2)
   real(kind = dp) :: UL(nVar),UR(nVar)
   real(kind = dp) :: flux(nVar)
   real(kind = dp) :: R(nVar)
   integer :: b,i,j,k,d1,d2

   if (space_order == 2 .and. space_disc == 0) then
      do b = 1,nBlock
         do k = 1, block(b) % nCell(3)
            do j = 1, block(b) % nCell(2)
               do i = 1, block(b) % nPkt(1)
                  block(b) % Q_rec(i,j,k,:,1,1) = dh * block(b) % Q(i-1,j,k,:) -eh *  block(b) % Q(i-2,j,k,:)
                  block(b) % Q_rec(i,j,k,:,1,2) = dh * block(b) % Q(i     ,j,k,:) -eh *  block(b) % Q(i+1,j,k,:)
                  end do
            end do
         end do
         do k = 1, block(b) % nCell(3)
            do j = 1, block(b) % nPkt(2)
               do i = 1, block(b) % nCell(1)
                  block(b) % Q_rec(i,j,k,:,2,1) = dh * block(b) % Q(i,j-1,k,:) -eh *  block(b) % Q(i,j-2,k,:)
                  block(b) % Q_rec(i,j,k,:,2,2) = dh * block(b) % Q(i,j     ,k,:) -eh *  block(b) % Q(i,j+1,k,:)
               end do
            end do
         end do
      end do
   else if (space_order == 2 .and. space_disc == 1) then ! MUSCL 2nd ORDER
      do b = 1,nBlock
         do k = 1, block(b) % nCell(3)
            do j = 1, block(b) % nCell(2)
               do i = 1, block(b) % nPkt(1)
                  R = (block(b) % Q(i+1,j,k,:) - block(b) % Q(i,j,k,:))  &
                      / (block(b) % Q(i,j,k,:) - block(b) % Q(i-1,j,k,:)+epsi)
                  call minmod(R)
                  block(b) % Q_rec(i,j,k,:,1,1) = block(b) % Q(i-1,j,k,:) + eh * R &
                                                                      * (block(b) % Q(i-1,j,k,:) - block(b) % Q(i-2,j,k,:))
                  R = ( block(b) % Q(i+1,j,k,:) - block(b) % Q(i,j,k,:) )  &
                      / ( block(b) % Q(i+2,j,k,:) - block(b) % Q(i+1,j,k,:) + epsi )
                  call minmod(R)
                  block(b) % Q_rec(i,j,k,:,1,2) = block(b) % Q(i,j,k,:) - eh * R &
                                                                      * (block(b) % Q(i+1,j,k,:) - block(b) % Q(i,j,k,:))
                  end do
            end do
         end do
         do k = 1, block(b) % nCell(3)
            do j = 1, block(b) % nPkt(2)
               do i = 1, block(b) % nCell(1)

                  R = (block(b) % Q(i,j+1,k,:) - block(b) % Q(i,j,k,:))  &
                      / (block(b) % Q(i,j,k,:) - block(b) % Q(i,j-1,k,:)+epsi)
                  call minmod(R)
                  block(b) % Q_rec(i,j,k,:,2,1) = block(b) % Q(i,j-1,k,:) + eh * R &
                                                                      * (block(b) % Q(i,j-1,k,:) - block(b) % Q(i,j-2,k,:))

                  R = ( block(b) % Q(i,j+1,k,:) - block(b) % Q(i,j,k,:) )  &
                      / ( block(b) % Q(i,j+2,k,:) - block(b) % Q(i,j+1,k,:) + epsi )
                  call minmod(R)
                  block(b) % Q_rec(i,j,k,:,2,2) = block(b) % Q(i,j,k,:) - eh * R &
                                                                      * (block(b) % Q(i,j+1,k,:) - block(b) % Q(i,j,k,:))
!                  block(b) % Q_rec(i,j,k,:,2,1) = dh * block(b) % Q(i,j-1,k,:) -eh *  block(b) % Q(i,j-2,k,:)
!                  block(b) % Q_rec(i,j,k,:,2,2) = dh * block(b) % Q(i,j     ,k,:) -eh *  block(b) % Q(i,j+1,k,:)
               end do
            end do
         end do
      end do
   else ! 1st order
      do b = 1,nBlock
         do k = 1, block(b) % nCell(3)
            do j = 1, block(b) % nCell(2)
               do i = 1, block(b) % nPkt(1)
                  block(b) % Q_rec(i,j,k,:,1,1) = block(b) % Q(i-1,j,k,:)
                  block(b) % Q_rec(i,j,k,:,1,2) = block(b) % Q(i,j,k,:)
               end do
            end do
         end do
         do k = 1, block(b) % nCell(3)
            do j = 1, block(b) % nPkt(2)
               do i = 1, block(b) % nCell(1)
                  block(b) % Q_rec(i,j,k,:,2,1) = block(b) % Q(i,j-1,k,:)
                  block(b) % Q_rec(i,j,k,:,2,2) = block(b) % Q(i,j,k,:)
               end do
            end do
         end do
      end do
   end if




   do b = 1,nBlock

      do k = 1, block(b) % nCell(3)
         do j = 1, block(b) % nCell(2)
            do i = 1, block(b) % nPkt(1)
               UL = block(b) % Q_rec(i,j,k,:,1,1)
               UR = block(b) % Q_rec(i,j,k,:,1,2)
               if ( control_riemann_solver == 1) then
                  block(b) % Flux(i,j,k,:,1) = Rotated_RHLL &
                                                                  (UL, UR, block(b) % Edge_Vec(:,i,j,k,1) )
               else if ( control_riemann_solver == 2) then
                  block(b) % Flux(i,j,k,:,1) = Roe &
                                                                  (UL, UR, block(b) % Edge_Vec(:,i,j,k,1) )
               end if
               block(b) % Flux(i,j,k,:,1) = block(b) % Flux(i,j,k,:,1) * block(b) % Edge_Len(i,j,k,1)
!               if (i == 24) then
!                  write(*,*) j,block(b) % Edge_Vec(:,i,j,k,1),block(b) % Edge_Vec(:,i,j,k,2)
!                  end if
            end do
         end do
      end do

      do k = 1, block(b) % nCell(3)
         do j = 1, block(b) % nPkt(2)
            do i = 1, block(b) % nCell(1)
               UL = block(b) % Q_rec(i,j,k,:,2,1)
               UR = block(b) % Q_rec(i,j,k,:,2,2)

               if ( control_riemann_solver == 1) then
                  block(b) % Flux(i,j,k,:,2) = Rotated_RHLL &
                                                                  (UL, UR, block(b) % Edge_Vec(:,i,j,k,2) )

               else if ( control_riemann_solver == 2) then
                  block(b) % Flux(i,j,k,:,2) = Roe &
                                                                  (UL, UR, block(b) % Edge_Vec(:,i,j,k,2) )

               end if

               block(b) % Flux(i,j,k,:,2) =block(b) % Flux(i,j,k,:,2) * block(b) % Edge_Len(i,j,k,2)
!               if (j == 24) then
!                  write(*,*) i,block(b) % Edge_Vec(:,i,j,k,1),block(b) % Edge_Vec(:,i,j,k,2)
!               end if
            end do
         end do
      end do

   end do

contains

function AUSM(uL,uR,vec)
use const, only : dp
implicit none
 real(kind = dp) :: uL(4), uR(4)    !  Input: conservative variables rho*[1, u, v, E]
 real(kind = dp) :: vec(2)             !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real(kind = dp) :: AUSM (4) ! Output: AUSM  flux function.

end function AUSM

!****************************************************************
!* -- Roe's Flux Function ---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!****************************************************************
 function Roe(uL, uR,vec)
use const, only : dp
implicit none
 real(kind = dp) :: uL(4), uR(4) !  Input: conservative variables rho*[1, u, v, E]
 real(kind = dp) :: vec(2)             !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real(kind = dp) :: Roe(4)       ! Output: Roe flux function (upwind)
!Local constants
 real(kind = dp) :: nx, ny       !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real(kind = dp) :: gamma                          ! Ratio of specific heat.
 real(kind = dp) :: zero, fifth, half, one, two    ! Numbers
!Local variables
 real(kind = dp) :: tx, ty       ! Tangent vector (perpendicular to the face normal)
 real(kind = dp) :: vxL, vxR, vyL, vyR             ! Velocity components.
 real(kind = dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
 real(kind = dp) :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
 real(kind = dp) :: aL, aR, HL, HR                 ! Speeds of sound.
 real(kind = dp) :: RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
 real(kind = dp) :: drho,dvx,dvy,dvn,dvt,dpr,dV(4)  ! Wave strenghs
 real(kind = dp) :: ws(4),dws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
 real(kind = dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 integer :: i, j
     nx = vec(1)
     ny = vec(2)
!Constants.
     gamma = 1.4
     zero = 0.0
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0

!Tangent vector (Do you like it? Actually, Roe flux can be implemented
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  tx = -ny
  ty = nx
!Primitive and other variables.
!  Left state

If ( ul(1) == 0.0E0_dp .or. ur(1) == 0.0E0_dp) then
   write(*,*) "Roe",ul,ur
end if

    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
     vnL = vxL*nx+vyL*ny
     vtL = vxL*tx+vyL*ty
      pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
     vnR = vxR*nx+vyR*ny
     vtR = vxR*tx+vyR*ty
      pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR
!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
    vx = (vxL+RT*vxR)/(one+RT)
    vy = (vyL+RT*vyR)/(one+RT)
     H = ( HL+RT* HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
    vn = vx*nx+vy*ny
    vt = vx*tx+vy*ty
!Wave Strengths
   drho = rhoR - rhoL
     dpr =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dpr - rho*a*dvn )/(two*a*a)
  dV(2) = rho*dvt/a
  dV(3) =  drho - dpr/(a*a)
  dV(4) = (dpr + rho*a*dvn )/(two*a*a)

!Wave Speed
  ws(1) = abs(vn-a)
  ws(2) = abs(vn)
  ws(3) = abs(vn)
  ws(4) = abs(vn+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one
  Rv(2,1) = vx - a*nx
  Rv(3,1) = vy - a*ny
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = vt*a

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx
  Rv(3,4) = vy + a*ny
  Rv(4,4) =  H + vn*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*vnL
  fL(2) = rhoL*vnL * vxL + pL*nx
  fL(3) = rhoL*vnL * vyL + pL*ny
  fL(4) = rhoL*vnL *  HL

  fR(1) = rhoR*vnR
  fR(2) = rhoR*vnR * vxR + pR*nx
  fR(3) = rhoR*vnR * vyR + pR*ny
  fR(4) = rhoR*vnR *  HR

  Roe = half * (fL + fR - diss)

 end function Roe

!****************************************************************
!* -- Rotated-Roe-HLL Flux Function ---
!*
!* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
!* Resolving, Rotated-Hybrid Riemann Solvers,
!* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
!*
!* The most robust Riemann solver known to the author (in terms of nonlinear
!* instability such as carbuncle).
!*
!* NB: At a boundary face, need to switch to a geometric normal vector:
!*               (nx2,ny2)=(nx, ny), (nx1,ny1)=(-ny,nx).
!*     This is not implemented here. It requires information on whether
!*     the geometric normal, (nx,ny), is on a boundary face or not.
!*     It shouldn't be difficult for you to implement it.
!*
!* Katate Masatsuka, February 2010. http://www.cfdbooks.com
!****************************************************************
 function Rotated_RHLL(uL, uR, vec)
 use const, only : dp
 implicit none
 real(kind = dp) :: uL(4), uR(4)    !  Input: conservative variables rho*[1, u, v, E]
 real(kind = dp) :: vec(2)             !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real(kind = dp) :: Rotated_RHLL(4) ! Output: Rotated_RHLL flux function.

!Local constants
 real(kind = dp) :: gamma                          ! Ratio of specific heat.
 real(kind = dp) :: zero, fifth, half, one, two    ! Numbers
 real(kind = dp) :: eps                            !
!Local variables
 real(kind = dp) :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
 real(kind = dp) :: tx, ty                         ! Tangent vector (taken as n1)
 real(kind = dp) :: alpha1, alpha2                 ! Projections of the new normals
 real(kind = dp) :: vxL, vxR, vyL, vyR             ! Velocity components.
 real(kind = dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
 real(kind = dp) :: vnL, vnR, vtL, vtR             ! Normal and tagent velocities
 real(kind = dp) :: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
 real(kind = dp) :: RT,rho,vx,vy,H,a               ! Roe-averages
 real(kind = dp) :: vn, vt                         ! Normal and tagent velocities(Roe-average)
 real(kind = dp) :: drho,dvx,dvy,dvn,dvt,dp1,dV(4)  ! Wave strenghs
 real(kind = dp) :: abs_dq                         ! Magnitude of the velocity difference
 real(kind = dp) :: abs_ws(4),ws(4),dws(4), Rv(4,4)! Wave speeds and right-eigevectors
 real(kind = dp) :: SRp,SLm                        ! Wave speeds for the HLL part
 real(kind = dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 real(kind = dp) :: temp
 real(kind = dp) :: nx, ny
 integer :: i, j

!Constants.
     nx = vec(1)
     ny = vec(2)
     gamma = 1.4E0_dp
      zero = 0.0E0_dp
     fifth = 0.2E0_dp
      half = 0.5E0_dp
       one = 1.0E0_dp
       two = 2.0E0_dp
       eps = 1.0e-5_dp ! 1.0e-12 in the original paper (double precision)

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
      pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
      pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

     vnL = vxL*nx + vyL*ny
     vnR = vxR*nx + vyR*ny

!Compute the flux.
   fL(1) = rhoL*vnL
   fL(2) = rhoL*vnL * vxL + pL*nx
   fL(3) = rhoL*vnL * vyL + pL*ny
   fL(4) = rhoL*vnL *  HL

   fR(1) = rhoR*vnR
   fR(2) = rhoR*vnR * vxR + pR*nx
   fR(3) = rhoR*vnR * vyR + pR*ny
   fR(4) = rhoR*vnR *  HR

!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(NB: n1 and n2 may need to be frozen at some point during
!     a steady calculation to fully make it converge. For time-accurate
!     calculation, this is fine.)
! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (vxR-vxL)**2+(vyR-vyL)**2 )
  if ( abs_dq > eps) then
       nx1 = (vxR-vxL)/abs_dq
       ny1 = (vyR-vyL)/abs_dq
  else
    nx1 = -ny
    ny1 =  nx
  endif

       nx2  = nx
       ny2 = ny
       nx1 = -ny
       ny1 = nx

    alpha1 = nx * nx1 + ny * ny1
!   To make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
    alpha1 = temp * alpha1

! Take n2 as perpendicular to n1.
       nx2 = -ny1
       ny2 =  nx1

       nx2  = nx
       ny2 = ny
       nx1 = -ny
       ny1 = nx

    alpha2 = nx * nx2 + ny * ny2
!   To make alpha2 always positive.
      temp = sign(one,alpha2)
       nx2 = temp * nx2
       ny2 = temp * ny2
    alpha2 = temp * alpha2

!Now we are going to compute the Roe flux with n2 as the normal
!and n1 as the tagent vector, with modified wave speeds (5.12)

!Compute the Roe Averages
     RT = sqrt(rhoR/rhoL)
    rho = RT*rhoL
     vx = (vxL+RT*vxR)/(one+RT)
     vy = (vyL+RT*vyR)/(one+RT)
      H = ( HL+RT* HR)/(one+RT)
      a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
     vn = vx*nx2+vy*ny2
     vt = vx*nx1+vy*ny1

!Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = vxL*nx2 + vyL*ny2
    vnR = vxR*nx2 + vyR*ny2
    vtL = vxL*nx1 + vyL*ny1
    vtR = vxR*nx1 + vyR*ny1

   drho = rhoR - rhoL
     dp1 =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dp1 - rho*a*dvn )/(two*a*a)
  dV(2) =  rho*dvt/a
  dV(3) =  drho - dp1/(a*a)
  dV(4) = (dp1 + rho*a*dvn )/(two*a*a)

!Wave Speeds for Roe flux part.
    ws(1) = vn-a
    ws(2) = vn
    ws(3) = vn
    ws(4) = vn+a
  abs_ws  = abs(ws)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
!only for the nonlinear fields.
  dws(1) = fifth
   if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
  dws(4) = fifth
   if (abs_ws(4)<dws(4)) abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))

!HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
   SRp = max( zero, vtR + aR, vt + a)
   SLm = min( zero, vtL - aL, vt - a)

!Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
   ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm)

!Right Eigenvectors: with n2 as normal and n1 as tangent.
  tx = nx1
  ty = ny1

  Rv(1,1) = one
  Rv(2,1) = vx - a*nx2
  Rv(3,1) = vy - a*ny2
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = a*vt

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx2
  Rv(3,4) = vy + a*ny2
  Rv(4,4) =  H + vn*a

!Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the Rotated-RHLL flux.
  Rotated_RHLL = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

 end function Rotated_RHLL



end subroutine calc_fluxes

subroutine minmod ( a )
   use const, only: dp
   use control, only : nVar
   implicit none
   real(kind = dp), intent(inout) ::a(nVar)
   integer :: n
   do n = 1, nVar
         a(n) = max(0.0E0_dp,min(1.0E0_dp,a(n)))
   end do
end subroutine minmod
