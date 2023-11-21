!-----------------------------------------------------------------------
! flag presence or absence of snow 
! create the snow-on mask
! consider creating this mask earlier on and skipping all other subroutines
! if MEB is switched off
! note that this is only separated from CRO_COUP because the size of the internal
! variables in CRO_COUP depends on ISIZE_SNOW, so this needs to be assigned in advance
!-----------------------------------------------------------------------
subroutine CRO_CHECKIN(crosnowmask, ISIZE_SNOW, NMASK)

use CRO_CSTS, only : &
  EXSNOWDMIN,         &
  EXRHOSMAX_ES 

use CRO_GRID, only : &
    KSIZE1

use CRO_STATE_VARS, only : &
  CROSNOW_ON

use DRIVING, only: &
  dt, &
  Rf, &
  Ta, &
  Sf                  

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Nsnow               ! Number of snow layers

implicit none

integer, intent(out) :: &
  crosnowmask(Nx,Ny),   &   ! Array that keeps track of points with snow - to enter crocus routine or not
  ISIZE_SNOW, &! Number of points with snow at timestep 
  NMASK(KSIZE1) ! Index correspondence array 
real ::  &
  zsnowfall


integer :: &
i,j                 ! Point counters


! initialize mask the crocus way 
ISIZE_SNOW = 0
NMASK(:) = 0

do j = 1, Ny
do i = 1, Nx
  zsnowfall = Sf(i,j) * dt/EXRHOSMAX_ES 
  if (Nsnow(i,j) > 0 .OR. zsnowfall >= EXSNOWDMIN) then 
    crosnowmask(i,j) = 1
    ISIZE_SNOW = ISIZE_SNOW + 1
    NMASK(ISIZE_SNOW) = j
    ! if snowpack existing: set snow flag to true, no update to state vars needed
    ! leave state vars as they are because it means crocus was running in previous 
    ! timestep 
  else
    crosnowmask(i,j) = 0
    ! no crocus snow 
  end if 
end do
end do 

end subroutine CRO_CHECKIN