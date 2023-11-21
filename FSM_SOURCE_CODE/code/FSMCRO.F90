!-----------------------------------------------------------------------
! FSMCRO: Combination of FSM2 and Crocus snow models
! Code structure based on Flexible Snow Model (FSM version 2.0)
!
! Original implementation> 
! Richard Essery
! School of GeoSciences
! University of Edinburgh
!
! FSM2 developments
! FSM modifications by OSHD@SLF (Adam Winstral, Jan Magnusson, Nora Helbig, 
!   Louis Queno, Bertrand Cluzet, Giulia Mazzotti)
! FSMCRO implementation and Crocus coupling by Giulia Mazzotti, SLF/CEN-MeteoFrance
!   relies on EXT_CROCUS code by Matthieu Lafaysse, Vincent Vionnet, Rafife Nehili, CEN
!
! This program replaces the program in the original version (FSM2.F90 subroutine)
!-----------------------------------------------------------------------
program FSMCRO

use DIAGNOSTICS, only: &
  Nave                ! Number of timesteps in average outputs

implicit none

logical :: EoR        ! End-of-run flag

integer :: i          ! Timestep counter

call SETUP

call CRO_SETUP

! Loop over timesteps
EoR = .false.
do
  do i = 1, Nave
    call DRIVE(EoR)
    if (EoR) goto 1
    call PHYSICS
  end do
end do

1 continue

call DUMP

end program FSMCRO


