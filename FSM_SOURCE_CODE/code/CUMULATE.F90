!-----------------------------------------------------------------------
! Cumulate fluxes
!-----------------------------------------------------------------------
subroutine CUMULATE(Roff, meltflux_out,Sbsrf,Sdirt,Sdift,LWt,asrf_out,Melt, &
                    Esrf,Eveg,Gsoil,Hsrf,intcpt,KH,KHa,Khg,KHv,KWg,KWv,  &
                    LE,LEsrf,LWsci,LWveg,Rnet,Rsrf,Sbveg,H,Swsci,SWsrf,  &
                    SWveg,Usc,unload)
                    
use MODE_WRITE, only: WRITE_2D

use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour                ! Hour of the day
  
use MODCONF, only: FOR_HN, CRO_ON

use MODOUTPUT, only: &
  WRITE_DIAG_TABLE,  &
  WRITE_STATE_TABLE, &
  OUT3D_ON

use GRID, only: &
  Nsmax,         &! Maximum number of snow layers
  Nsoil,         &! Number of soil layers
  Nx,Ny           ! Grid dimensions

use STATE_VARIABLES, only: &
  albs,              &! Snow albedo
  fsnow,             &! Snow cover fraction
  Qcan,              &! Canopy air space humidity
  Ds,                &! Snow layer thicknesses (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  snowdepthmin,      &! min Snow at time step of swemin (m)
  snowdepthmax,      &! max Snow at time step of swemax (m)
  swemin,            &! min swe during season (mm)
  swemax,            &! max swe during season (mm)  
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsrf,              &! Surface skin temperature (K)
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  Tveg                ! Vegetation temperature (K)

USE CRO_DIAG, only: &
  ZP_SNOWTYPEMEPRA, &
  ZP_SNOWSSA

implicit none

real, intent(in) :: &
  asrf_out(Nx,Ny),   &! Surface albedo
  intcpt(Nx,Ny),     &! Canopy interception (kg/m^2)  
  Roff(Nx,Ny),       &! Total runoff (kg/m^2)
  meltflux_out(Nx,Ny),  &! Runoff from snowmelt at base of snow (kg/m^2)
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)$
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHg(Nx,Ny),        &! Eddy diffusivity for heat from the ground (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny),        &! Eddy diffusivity for water from vegetation (m/s)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  LWsci(Nx,Ny),      &! Subcanopy incoming SWR (W/m^2)
  LWt(Nx,Ny),        &! Incoming longwave radiation corrected for subgrid topography (W/m^2)
  LWveg(Nx,Ny),      &! Net longwave radiation absorbed by vegetation (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Rsrf(Nx,Ny),       &! Net radiation absorbed by the surface (W/m^2
  Sbsrf(Nx,Ny),      &! Sublimation from the snow surface (kg/m^2)
  Sbveg(Nx,Ny),      &! Sublimation from the vegetation (kg/m^2)
  Sdirt(Nx,Ny),      &! Incoming direct beam radiation corrected for subgrid topography (W/m^2)
  Sdift(Nx,Ny),      &! Incoming diffuse radiation corrected for subgrid topography (W/m^2)
  SWsci(Nx,Ny),      &! Subcanopy incoming SWR (W/m^2)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  unload(Nx,Ny),     &! Snow mass unloaded from canopy (kg/m^2)
  Usc(Nx,Ny)          ! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)               

real :: &
  Sliq_out(Nx,Ny),   &!
  snowdepth(Nx,Ny),  &! Snow depth (m)
  SWE(Nx,Ny) ,   &         ! Snow water equivalent (kg/m^2)
  tmpvar(Nx,Ny)    ! Tmp 2D var to store individual layers of 3D var

integer :: i,j,k,where,ii

! BC just in case, these sums should be performed only until Nsnow.
do j = 1,Ny
  do i = 1,Nx
    Sliq_out(i,j) = sum(Sliq(:,i,j))
    snowdepth(i,j) = sum(Ds(:,i,j)) * fsnow(i,j)
    SWE(i,j) = sum(Sice(:,i,j)) + sum(Sliq(:,i,j))
  end do
end do

inquire(unit=1401, pos=where)
write(1401,pos=where) year
inquire(unit=1402, pos=where)
write(1402,pos=where) month
inquire(unit=1403, pos=where)
write(1403,pos=where) day
inquire(unit=1404, pos=where)
write(1404,pos=where) hour

ii = 1
if (write_diag_table(ii)) call WRITE_2D(Roff,1404 + ii)        ! Total runoff, snow and bare soil (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(snowdepth, 1404 + ii)  ! Snow depth (m)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(SWE, 1404 + ii)        ! Snow water equivalent (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sliq_out, 1404 + ii)   ! Sliq_out(i,j) = Sliq(1,i,j) + Sliq(2,i,j) + Sliq(3,i,j)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sdirt,1404 + ii)       ! Incoming direct beam radiation corrected for subgrid topography (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sdift, 1404 + ii)      ! Incoming diffuse radiation corrected for subgrid topography (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LWt, 1404 + ii)        ! Incoming longwave radiation corrected for subgrid topography (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(meltflux_out, 1404 + ii)  ! Runoff from snowmelt at base of snow (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sbsrf, 1404 + ii)      ! Snow sublimation rate (kg/m^2/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(asrf_out, 1404 + ii)   ! Surface albedo
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Melt, 1404 + ii)       ! Surface melt rate (kg/m^2/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Esrf, 1404 + ii)       ! Moisture flux from the surface (kg/m^2/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Eveg, 1404 + ii)       ! Moisture flux from vegetation (kg/m^2/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Gsoil, 1404 + ii)      ! Heat flux into soil (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Hsrf, 1404 + ii)       ! Sensible heat flux from the surface (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(intcpt, 1404 + ii)     ! Canopy interception (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KH, 1404 + ii)         ! Eddy diffusivity for heat to the atmosphere (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KHa, 1404 + ii)        ! Eddy diffusivity for heat from the canopy air space (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KHg, 1404 + ii)        ! Eddy diffusivity for heat from the ground (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KHv, 1404 + ii)        ! Eddy diffusivity for heat from vegetation (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KWg, 1404 + ii)        ! Eddy diffusivity for water from the ground (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(KWv, 1404 + ii)        ! Eddy diffusivity for water from vegetation (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LE, 1404 + ii)         ! Latent heat flux to the atmosphere (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LEsrf, 1404 + ii)      ! Latent heat flux from the surface (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LWsci, 1404 + ii)      ! Subcanopy incoming longwave radiation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(LWveg, 1404 + ii)      ! Net longwave radiation absorbed by vegetation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Rnet, 1404 + ii)       ! Net radiation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Rsrf, 1404 + ii)       ! Net radiation at surface (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Sbveg, 1404 + ii)      ! Sublimation from vegetation (kg/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(H, 1404 + ii)          ! Sensible heat flux to the atmosphere (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(SWsci, 1404 + ii)      ! Subcanopy incoming shortwave radiation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(SWsrf, 1404 + ii)      ! Net SW radiation absorbed by the surface (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(SWveg, 1404 + ii)      ! Net SW radiation absorbed by vegetation (W/m^2)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(Usc, 1404 + ii)        ! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)
ii = ii + 1
if (write_diag_table(ii)) call WRITE_2D(unload, 1404 + ii)     ! Snow mass unloaded from canopy (kg/m^2)


! If necessary, write state vars into results files:
ii = 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(albs,1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(fsnow, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Tsrf, 1500 + ii)
ii = ii + 1
! if (CANMOD==0) then BC this is commented on purpose. DO NOT introduce this if.
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(snowdepthmin, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(snowdepthmax, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(swemin, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(swemax, 1500 + ii)
ii = ii + 1
! else BC this is commented on purpose. DO NOT introduce this if.
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Qcan, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Sveg, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Tcan, 1500 + ii)
ii = ii + 1
if (WRITE_STATE_TABLE(ii)) call WRITE_2D(Tveg, 1500 + ii)
ii = ii + 1
! endif BC this is commented on purpose. DO NOT introduce this if.

! write 01h-07h-13h-19h states (UTC+1) for the subsequent HN model (<-> 00UTC,06UTC,12UTC,18UTC)
if (FOR_HN) then
  if (hour > 0.5 .and. hour < 1.5) then
    write(1227) ((Ds(1,i,j),j=1,Ny),i=1,Nx)
    write(1228) ((Tsrf(i,j),j=1,Ny),i=1,Nx)
    write(1229) ((Tsnow(1,i,j),j=1,Ny),i=1,Nx)
    write(1230) ((Tsoil(1,i,j),j=1,Ny),i=1,Nx)

    close(1227)
    close(1228)
    close(1229)
    close(1230)
  endif
  if (hour > 6.5 .and. hour < 7.5) then
    write(1231) ((Ds(1,i,j),j=1,Ny),i=1,Nx)
    write(1232) ((Tsrf(i,j),j=1,Ny),i=1,Nx)
    write(1233) ((Tsnow(1,i,j),j=1,Ny),i=1,Nx)
    write(1234) ((Tsoil(1,i,j),j=1,Ny),i=1,Nx)

    close(1231)
    close(1232)
    close(1233)
    close(1234)
  endif
  if (hour > 12.5 .and. hour < 13.5) then
    write(1235) ((Ds(1,i,j),j=1,Ny),i=1,Nx)
    write(1236) ((Tsrf(i,j),j=1,Ny),i=1,Nx)
    write(1237) ((Tsnow(1,i,j),j=1,Ny),i=1,Nx)
    write(1238) ((Tsoil(1,i,j),j=1,Ny),i=1,Nx)

    close(1235)
    close(1236)
    close(1237)
    close(1238)
  endif
  if (hour > 18.5 .and. hour < 19.5) then
    write(1239) ((Ds(1,i,j),j=1,Ny),i=1,Nx)
    write(1240) ((Tsrf(i,j),j=1,Ny),i=1,Nx)
    write(1241) ((Tsnow(1,i,j),j=1,Ny),i=1,Nx)
    write(1242) ((Tsoil(1,i,j),j=1,Ny),i=1,Nx)

    close(1239)
    close(1240)
    close(1241)
    close(1242)
  endif
endif

if (OUT3D_ON) then
  
  do k = 1,Nsmax 
    tmpvar = Ds(k,:,:)
    call WRITE_2D(tmpvar, 1600 + k)
    tmpvar = Sice(k,:,:)
    call WRITE_2D(tmpvar, 1700 + k)
    tmpvar = Sliq(k,:,:)
    call WRITE_2D(tmpvar, 1800 + k)
    tmpvar = Tsnow(k,:,:)
    call WRITE_2D(tmpvar, 1900 + k)
    if (CRO_ON) then
      tmpvar(1,:) = ZP_SNOWTYPEMEPRA(:,k)
      call WRITE_2D(tmpvar, 2000 + k)
      tmpvar(1,:) = ZP_SNOWSSA(:,k)
      call WRITE_2D(tmpvar, 2100 + k)
    end if 
  end do

  do k = 1,Nsoil
    tmpvar = Tsoil(k,:,:)
    call WRITE_2D(tmpvar, 2000-k)
  end do

end if

end subroutine CUMULATE
