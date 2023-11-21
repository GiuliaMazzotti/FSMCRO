!-----------------------------------------------------------------------
! Read point driving data
!-----------------------------------------------------------------------
subroutine DRIVE(EoR)

use MODCONF, only: CANMOD,CRO_ON,MET1D_ON,TVT_ON

use MODPERT, only: Z0PERT,WCPERT,FSPERT,ALPERT,SLPERT

use CONSTANTS, only: &
  eps,               &! Ratio of molecular weights of water and dry air
  e0,                &! Saturation vapour pressure at Tm (Pa)
  Tm                  ! Melting point (K)

!! ATTENTION: Sf and Rf are initially hourly precipitation (from OSHD matlab scripts), converted here to a precipitation rate
use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour,              &! Hour of day
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Sf,                &! Snowfall rate (kg/m2/s)
  Sf24h,             &! Snowfall 24hr (kg/m2)
  Ta,                &! Air temperature (K)
  Tv,                &! Time-varying transmissivity for direct SWR (-)
  RH,                &! Relative humidity (%)
  Ua,                &! Wind speed (m/s)
  z0P,               &! z0 perturbations
  wcP,               &! liquid water capacity perturbations
  fsP,               &! fresh snow density perturbations
  alP,               &! albedo perturbations
  slP                 ! settling perturbations

use GRID, only: &
  Nx,Ny               ! Grid dimensions
  
use LANDUSE, only: &
  vfhp                ! Hemispherical sky-view fraction including canopy

implicit none

logical, intent(out) :: &
  EoR                 ! End-of-run flag

real*4 :: &
  es(Nx,Ny),      &! Saturation vapour pressure (Pa)
  Tc(Nx,Ny),      &! Temperature (C)
  tmpmet           ! Temporary variable for meteo data

integer :: i,j,where,eastatus 

! FSM driving data
inquire(unit=800, pos=where)
read(800,pos=where,IOSTAT=eastatus) year
inquire(unit=801, pos=where)
read(801,pos=where,IOSTAT=eastatus) month
inquire(unit=802, pos=where)
read(802,pos=where,IOSTAT=eastatus) day
inquire(unit=803, pos=where)
read(803,pos=where,IOSTAT=eastatus) hour
inquire(unit=804, pos=where)

if(MET1D_ON) then
  read(804,pos=where,IOSTAT=eastatus) tmpmet
  Sdir(:,:) = tmpmet
  inquire(unit=805, pos=where)
  read(805,pos=where,IOSTAT=eastatus) tmpmet
  Sdif(:,:) = tmpmet
  inquire(unit=806, pos=where)
  read(806,pos=where,IOSTAT=eastatus) tmpmet 
  LW(:,:) = tmpmet
  inquire(unit=807, pos=where)
  read(807,pos=where,IOSTAT=eastatus) tmpmet 
  Sf(:,:) = tmpmet
  inquire(unit=808, pos=where)
  read(808,pos=where,IOSTAT=eastatus) tmpmet 
  Rf(:,:) = tmpmet
  inquire(unit=809, pos=where)
  read(809,pos=where,IOSTAT=eastatus) tmpmet 
  Ta(:,:) = tmpmet
  inquire(unit=810, pos=where)
  read(810,pos=where,IOSTAT=eastatus) tmpmet 
  RH(:,:) = tmpmet
  inquire(unit=811, pos=where)
  read(811,pos=where,IOSTAT=eastatus) tmpmet 
  Ua(:,:) = tmpmet
  inquire(unit=812, pos=where)
  read(812,pos=where,IOSTAT=eastatus) tmpmet 
  Ps(:,:) = tmpmet
  if (.NOT. CRO_ON) then
    inquire(unit=813, pos=where)
    read(813,pos=where,IOSTAT=eastatus) tmpmet 
    Sf24h(:,:) = tmpmet
  endif
else ! original implementation
  read(804,pos=where,IOSTAT=eastatus) ((Sdir(i,j),j=1,Ny),i=1,Nx)
  inquire(unit=805, pos=where)
  read(805,pos=where,IOSTAT=eastatus) ((Sdif(i,j),j=1,Ny),i=1,Nx)
  inquire(unit=806, pos=where)
  read(806,pos=where,IOSTAT=eastatus) ((LW(i,j),j=1,Ny),i=1,Nx)
  inquire(unit=807, pos=where)
  read(807,pos=where,IOSTAT=eastatus) ((Sf(i,j),j=1,Ny),i=1,Nx)
  inquire(unit=808, pos=where)
  read(808,pos=where,IOSTAT=eastatus) ((Rf(i,j),j=1,Ny),i=1,Nx)
  inquire(unit=809, pos=where)
  read(809,pos=where,IOSTAT=eastatus) ((Ta(i,j),j=1,Ny),i=1,Nx)
  inquire(unit=810, pos=where)
  read(810,pos=where,IOSTAT=eastatus) ((RH(i,j),j=1,Ny),i=1,Nx)
  inquire(unit=811, pos=where)
  read(811,pos=where,IOSTAT=eastatus) ((Ua(i,j),j=1,Ny),i=1,Nx)
  inquire(unit=812, pos=where)
  read(812,pos=where,IOSTAT=eastatus) ((Ps(i,j),j=1,Ny),i=1,Nx)
  if (.NOT. CRO_ON) then
    inquire(unit=813, pos=where)
    read(813,pos=where,IOSTAT=eastatus) ((Sf24h(i,j),j=1,Ny),i=1,Nx)
  endif
endif

if (CANMOD == 1) then
  if (TVT_ON) then
    inquire(unit=814, pos=where)
    read(814,pos=where,IOSTAT=eastatus) ((Tv(i,j),j=1,Ny),i=1,Nx)
  else 
    Tv(:,:) = vfhp
  endif
endif

if (Z0PERT) then
  inquire(unit=815, pos=where)
  read(815,pos=where,IOSTAT=eastatus) ((z0P(i,j),j=1,Ny),i=1,Nx)
endif

if (WCPERT) then
  inquire(unit=816, pos=where)
  read(816,pos=where,IOSTAT=eastatus) ((wcP(i,j),j=1,Ny),i=1,Nx)
endif

if (FSPERT) then
  inquire(unit=817, pos=where)
  read(817,pos=where,IOSTAT=eastatus) ((fsP(i,j),j=1,Ny),i=1,Nx)
endif

if (ALPERT) then
  inquire(unit=818, pos=where)
  read(818,pos=where,IOSTAT=eastatus) ((alP(i,j),j=1,Ny),i=1,Nx)
endif

if (SLPERT) then
  inquire(unit=819, pos=where)
  read(819,pos=where,IOSTAT=eastatus) ((slP(i,j),j=1,Ny),i=1,Nx)
endif

! Check for end of driving data file, otherwise perform unit conversions
if (eastatus < 0) then
  EoR = .true. ! End of file
else 

  ! Meteo data adaptations 
  Ua = max(Ua, 0.1)

  do i = 1,Nx
  do j = 1,Ny
    Sf(i,j) = Sf(i,j)/dt
    Rf(i,j) = Rf(i,j)/dt
    Tc(i,j) = Ta(i,j) - Tm
    es(i,j) = e0*exp(17.5043*Tc(i,j)/(241.3 + Tc(i,j)))
    Qa(i,j) = (RH(i,j)/100)*eps*es(i,j)/Ps(i,j)
  end do
  end do
  
endif

write(*,*) 'timestep: ', year, '-', month, '-', day, '-', hour, ' ', EoR

return

end subroutine DRIVE
