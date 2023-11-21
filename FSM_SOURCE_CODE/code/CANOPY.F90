!-----------------------------------------------------------------------
! Mass balance of canopy snow
!-----------------------------------------------------------------------
subroutine CANOPY(Eveg,unload,intcpt,Sbveg)

use MODCONF, only: CRO_ON

use MODTILE, only: tthresh 

use CONSTANTS, only: &
  eps,               &! Ratio of molecular weights of water and dry air
  Rair,              &! Gas constant for air (J/K/kg)
  Tm,                &! Melting point (K)
  vkman               ! Von Karman constant

use DRIVING, only: &
  dt,                &! Timestep (s)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Ta,                &! Air temperature (K)
  Sf,                &! Snowfall rate (kg/m2/s)
  zU,                &! Wind measurement height (m)
  Ua                  ! Wind speed (m/s)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  cveg,              &! Vegetation turbulent transfer coefficient
  rchd,              &! Ratio of displacement height to canopy height
  rchz,              &! Ratio of roughness length to canopy height
  tcnc,              &! Canopy unloading time scale for cold snow (s)
  tcnm,              &! Canopy unloading time scale for melting snow (s)
  psf,               &! Scaling factor for solid precipitation (within forest stand, at min fveg)
  psr                 ! Range of solid precipitation (within forest stand, spread min-max CC)
  
use PARAMMAPS, only: &
  VAI,               &! Vegetation area index
  scap,              &! Canopy snow capacity (kg/m^2)
  z0sf                ! Snow-free surface roughness length (m)

use STATE_VARIABLES, only: &
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tveg                ! Vegetation temperature (K)

use LANDUSE, only: &
  fveg,              &! Canopy cover fraction
  pmultf,            &! Precip multiplier applied to open-area snowfall 
  hcan,              &! Canopy height (m)
  tilefrac            ! Grid cell tile fraction

implicit none

real, intent(in) :: &
  Eveg(Nx,Ny)         ! Moisture flux from vegetation (kg/m^2/s)

real, intent(out) :: &
  intcpt(Nx,Ny),     &! Canopy interception (kg/m^2)
  Sbveg(Nx,Ny),      &! Sublimation from the vegetation (kg/m^2)
  unload(Nx,Ny)       ! Snow mass unloaded from canopy (kg/m^2)

real :: &
  Evegs,             &! Canopy snow sublimation rate (kg/m^2/s)
  Qs,                &! Saturation humidity
  dh,                &! Displacement height (m)
  rho,               &! Air density (kg/m^3)
  zU1,               &! Wind measurement height with offset (m)
  z0,                &! Roughness length for momentum (m)
  tunl                ! Canopy snow unloading timescale (s)

integer :: & 
  i,j                 ! Grid coordinates

do j = 1, Ny
do i = 1, Nx 

  unload(i,j) = 0
  intcpt(i,j) = 0
  Sbveg(i,j)  = 0

  if (tilefrac(i,j) < tthresh) goto 1 ! exclude points outside tile of interest
  
  if (fveg(i,j) > epsilon(fveg(i,j))) then
    
    ! rescale precipitation to correct back precip multiplier applied to open area 
    Sf(i,j) = pmultf(i,j)*Sf(i,j)
    ! interception
    intcpt(i,j) = (scap(i,j) - Sveg(i,j))*(1 - exp(-fveg(i,j)*Sf(i,j)*dt/scap(i,j)))
    Sveg(i,j) = Sveg(i,j) + intcpt(i,j)
    Sf(i,j) = Sf(i,j) - intcpt(i,j)/dt 
    Sf(i,j) = (psf - psr*fveg(i,j))*Sf(i,j) ! including preferential deposition in canopy gaps; might have to be revisited to ensure mass conservation, potentially integrate with pmultf

    ! sublimation
    Evegs = 0
    if (Sveg(i,j) > epsilon(Sveg(i,j)) .or. Tveg(i,j) < Tm) then 
      if (.NOT. CRO_ON) then  
        Evegs = Eveg(i,j)
      else 
        call QSAT(Ps(i,j),Ta(i,j),Qs)
        !Evegs = 0.002/3600*Ua(i,j)*Ps(i,j)/eps*(Qs-Qa(i,j))                  ! Lundquist formulation
        !Evegs = 0.002/3600*2.5/VAI(i,j)*Ua(i,j)*Ps(i,j)/eps*(Qs-Qa(i,j))     ! Lundquist formulation, attempting to 
        
        z0 = ((rchz*hcan(i,j))**fveg(i,j)) * (z0sf(i,j)**(1 - fveg(i,j)))
        dh = rchd * hcan(i,j) 
        zU1 = zU + hcan(i,j)
        rho = Ps(i,j) / (Rair*Ta(i,j))
        Evegs = rho*sqrt((vkman / log((zU1 - dh)/z0))*Ua(i,j))*VAI(i,j)/cveg*(Qs-Qa(i,j))
      end if
    end if
    Sveg(i,j) = Sveg(i,j) - Evegs*dt
    Sbveg(i,j) = Evegs*dt
    if (Sveg(i,j) < 0) Sbveg(i,j) =  Sbveg(i,j) + Sveg(i,j)
    Sveg(i,j) = max(Sveg(i,j), 0.)

    ! unloading
    tunl = tcnc
    if (Tveg(i,j) >= Tm) tunl = tcnm
    tunl = max(tunl, dt)
    unload(i,j) = Sveg(i,j)*dt/tunl
    Sveg(i,j) = Sveg(i,j) - unload(i,j)
    
    ! ensure that vegetation capacity is not exceeded / unload excess snow 
    if (Sveg(i,j) > scap(i,j)) then
      unload(i,j) = unload(i,j) + (Sveg(i,j)-scap(i,j))
      Sveg(i,j) = scap(i,j)
    end if 
  end if
 
  1 continue ! exclude points
  
end do
end do
 
end subroutine CANOPY
