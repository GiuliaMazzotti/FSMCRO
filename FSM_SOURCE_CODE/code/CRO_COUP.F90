SUBROUTINE CRO_COUP(unload,ksoil,SWsrf,Rsrf,Hsrf,LEsrf,G, &
    Gsoil, meltflux_out, Roff, Sbsrf,SWsci,LWsci,Usc, &
    crosnowmask,ISIZE_SNOW,NMASK)   
! Crocus-FSM coupler, modified from PROG.F90 as provided in EXT_CROCUS
! State variables that are updated by SNOW.F90 in FSM2 are also updated here. 

! Import CROCUS modules
USE MODD_TYPE_DATE_SURF
USE MODI_SNOWCRO        
USE MODI_SNOWCRO_DIAG     

! Import CRO modules  
USE CRO_CSTS
USE CRO_SETTINGS 
USE CRO_DRIVE 
USE CRO_GRID
USE CRO_PARAMS
USE CRO_SURF
USE CRO_STATE_VARS 
USE CRO_DUMMIES
USE CRO_DIAG

! Import FSM modules 
use MODTILE, only: TILE
USE CONSTANTS, only: sb, pi, Tm
USE GRID, only: Nx,Ny,Nsoil
USE DRIVING, only: year,month,day,hour,Ps,Sf,Rf,Ta,Sdir,Sdif,Qa,Ua,Lw,dt
USE LANDUSE, only: vfhp
USE STATE_VARIABLES, only: albs,Ds,fsnow,Nsnow,Sliq,Sice,Tsnow,Tsoil,Tcan,Qcan
USE MODD_SNOW_PAR,  ONLY : X_RI_MAX

IMPLICIT NONE

TYPE(DATE_TIME)  :: TPTIME     ! current date and time

! Function arguments = FSM variable names 
real, intent(in) :: &
    unload(Nx,Ny),      &
    ksoil(Nsoil,Nx,Ny), &
    Usc(Nx,Ny),         &
    G(Nx,Ny)     

real, intent(inout) :: &
    SWsrf(Nx,Ny),       &
    Rsrf(Nx,Ny),        &
    Hsrf(Nx,Ny),        &
    LEsrf(Nx,Ny),       &
    LWsci(Nx,Ny),       &
    SWsci(Nx,Ny)

real, intent(out) :: &
    Gsoil(Nx,Ny),        & ! Heat flux into soil (W/m^2)
    Roff(Nx,Ny),         & ! Total runoff (kg/m^2)
    meltflux_out(Nx,Ny), & ! Runoff from snowmelt at base of snow (kg/m^2)
    Sbsrf(Nx,Ny)           ! Sublimation from the snow surface (kg/m^2)

! Function arguments: snow mask
integer, intent(in) :: &
    ISIZE_SNOW, &! Number of points with snow at timestep 
    NMASK(KSIZE1) ! Index correspondence array 

integer, intent(inout) :: &
    crosnowmask(Nx,Ny) 
  
! Internal variables = EXTCROCUS variable names
! include everything that requires coupling & update at every time step
! additional to arguments of SNOWCRO, arguments of radiation transfer function and dummies
real(8), DIMENSION(KSIZE1) :: & 
    ZP_EXNS,        & ! Exner function at surface
    ZP_EXNA,        & ! Exner function at lowest atmos level
	ZP_HSNOW,       & ! sensible heat flux from snow (W/m2)
    ZP_LEL3L,       & ! liquid water evaporation flux from snow (always 0 W/m2 in Crocus)
    ZP_LES3L,       & ! sublimation heat flux from snow (W/m2)
    ZP_PEW_B_COEF,  & ! B-wind coefficient (m/s)
    ZP_PET_B_COEF,  & ! B-air temperature coefficient
    ZP_PEQ_B_COEF,  & ! B-air specific humidity coefficient
    ZP_RHOA,        & ! Air density
    ZP_RNSNOW,      & ! net radiative flux from snow (W/m2)
    ZP_SOILCOND,    & ! soil thermal conductivity [W/(m K)] 
    ZP_SWNETSNOW,   & ! net shortwave radiation entering top of snowpack(W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
    ZP_LWNETSNOW,   & ! net longwave radiation entering top of snowpack(W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
    ZP_SWNETSNOWS,  & ! net shortwave radiation in uppermost layer of snowpack (W m-2) Imposed if MEB=T, diagnosed if MEB=F; Used for surface energy budget diagnostics
    ZP_ZENITH,      & ! solar zenith angle
    dsout,          & ! dummy for DZ coupling & masking 
    sweout            ! dummy for SWE coupling & masking 
    
! Internal variables = packed variables 
real(8), DIMENSION(ISIZE_SNOW) :: & 
    X_PEW_A_COEF,X_PET_A_COEF,X_PEQ_A_COEF,X_PEW_B_COEF, &
    X_PET_B_COEF,X_PEQ_B_COEF,X_SNOWALB,X_PS,X_SRSNOW,X_RRSNOW, &
    X_PSN3L,X_TA,X_SW_RAD,X_QA,X_VMOD,X_LW_RAD,X_RHOA,X_UREF, &
    X_EXNS,X_EXNA,X_DIRCOSZW,X_ZREF,X_Z0NAT,X_Z0EFF,X_Z0HNAT, &
    X_ALB,X_SOILCOND,X_THRUFAL,X_GRNDFLUX,X_EVAPCOR,X_GFLXCOR, &
    X_SWNETSNOW,X_SWNETSNOWS,X_LWNETSNOW,X_RNSNOW,X_HSNOW, &
    X_GFLUXSNOW,X_HPSNOW,X_LES3L,X_LEL3L,X_EVAP,X_SNDRIFT,X_RI, &
    X_EMISNOW,X_CDSNOW,X_USTARSNOW,X_CHSNOW,X_SNOWHMASS,X_QS, &
    X_VEGTYPE,X_ZENITH,X_ANGL_ILLUM,X_LAT,X_LON,X_GSFCSNOW,  &
    X_SNOWMAK,X_SNDPT_12H,X_SNDPT_1DY,X_SNDPT_3DY,X_SNDPT_5DY,  &
    X_SNDPT_7DY,X_SNSWE_1DY,X_SNSWE_3DY,X_SNSWE_5DY,X_SNSWE_7DY, & 
    X_SNRAM_SONDE,X_SN_WETTHCKN,X_SN_REFRZNTHCKN,X_DEP_HIG,  &
    X_DEP_MOD,X_DEP_SUP,X_DEP_TOT,X_DEP_HUM,X_ACC_LEV,X_NAT_LEV, &
    X_PRO_SUP_TYP,X_PRO_INF_TYP,X_AVA_TYP

real(8), DIMENSION(ISIZE_SNOW,KSIZE2) :: &     
    X_SNOWSWE, X_SNOWRHO, X_SNOWHEAT, X_SNOWDIAMOPT, X_SNOWSPHERI, &
    X_SNOWHIST, X_SNOWAGE, X_SNOWLIQ, X_SNOWTEMP, X_SNOWDZ, & 
    X_SNOWDEND,X_SNOWSPHER,X_SNOWSIZE,X_SNOWSSA, X_SNOWTYPEMEPRA,    &
    X_SNOWRAM,X_SNOWSHEAR,X_ACC_RAT,X_NAT_RAT

real(8), DIMENSION(ISIZE_SNOW,NIMPUR) :: &   
    X_IMPDRY,X_IMPWET

real(8), DIMENSION(ISIZE_SNOW,KSIZE2,NIMPUR) :: &   
    X_SNOWIMPUR,X_SNOWIMP_CONC

real(8), DIMENSION(ISIZE_SNOW,KSIZE3) :: & 
    X_TG, X_D_G

real(8), DIMENSION(ISIZE_SNOW,KSIZE4) :: &    
    X_BLOWSNW

real(8), DIMENSION(ISIZE_SNOW,2) :: &    
    X_DIR_SW,X_SCA_SW,X_SPEC_ALB,X_DIFF_RATIO,X_SPEC_TOT

! Loop counters
integer :: &
    JI, JJ,          &  ! Mask indices conversions
    j,               &  ! Point counters
    ik,              &  ! Third dim counter
    k                   ! Level counter

! VARIABLE COUPLING AND UPDATES - BEFORE CALLING CROCUS
! Contains all variables that need to be updated at every timestep, prior to 
! being passed to Crocus (intent IN, intent INOUT in SNOWCRO)

! CRO_DRIVE
TPTIME%TDATE%YEAR  = year
TPTIME%TDATE%MONTH = month
TPTIME%TDATE%DAY   = day
TPTIME%TIME         = 1*3600. !hour*3600.
do j = 1, KSIZE1
    ZP_PS(j)     = Ps(1,j)
    ZP_SRSNOW(j) = Sf(1,j)                  ! check slope correction  
    ZP_RRSNOW(j) = Rf(1,j)   
    if (TILE == 'open') then                ! IN NEXT VERSION: instead of this, call pre-processor if tile == forest
        ZP_TA(j)     = Ta(1,j)
        ZP_SW_RAD(j) = Sdir(1,j) + Sdif(1,j) ! check slope / terrain correction 
        ZP_QA(j)     = Qa(1,j)
        ZP_VMOD(j)   = Ua(1,j)
        ZP_LW_RAD(j) = Lw(1,j)               ! potentially change to LWt            ! check slope correction   
    else ! TILE == 'forest'
        ! Note that if MEB is switched on, sub-canopy meteo not really used because fluxes are prescribed - hence whether the sub-canopy meteo
        ! or the driving data is used here should not make a difference
        ZP_SW_RAD(j) = SWsci(1,j)           
        ZP_TA(j)     = Ta(1,j)       
        ZP_QA(j)     = Qa(1,j)
        ZP_VMOD(j)   = min(Usc(1,j), 0.1)  ! min(0.4*Ua(1,j),0.1) 
        ZP_LW_RAD(j) = LW(1,j)*vfhp(1,j)+sb*(ZP_TA(j))**4*(1-vfhp(1,j)) ! no distinction between terrain and canopy elements with this approach             
        if (Ta(1,j) >= Tm) then
            ZP_RRSNOW(j) = ZP_RRSNOW(j) + unload(1,j)/dt
        else
            ZP_SRSNOW(j) = ZP_SRSNOW(j) + unload(1,j)/dt
        end if 
    end if 
    ! check if updating zT and zU is needed?
end do 

! Compute zenith angle using stripped-down version of the Surfex function
CALL SUNPOS(year, month, day, TPTIME%TIME, ZP_LON, ZP_LAT, ZP_ZENITH)

! Atmospheric coupling functions 
if (OMEB) then ! all of this is not needed
    ZP_PEW_B_COEF = EXUNDEF
    ZP_PET_B_COEF = EXUNDEF
    ZP_PEQ_B_COEF = EXUNDEF
    ZP_EXNS       = EXUNDEF
    ZP_EXNA       = EXUNDEF
    ZP_RHOA       = EXUNDEF
else 
    ZP_PEW_B_COEF = ZP_VMOD
    ZP_PET_B_COEF = ZP_TA/(ZP_PS/EXP00)**(EXRD/EXCPD)
    ZP_PEQ_B_COEF = ZP_QA  
    ZP_EXNS       = (ZP_PS/EXP00)**(EXRD/EXCPD)
    ZP_EXNA       = ZP_EXNS
    ZP_RHOA       = ZP_PS/(EXRD*ZP_TA*(1.+((EXRV/EXRD)-1.)*ZP_QA)+EXG*ZP_ZREF)
end if 

! Fluxes coupling -> NOT YET OPERATIONAL
if (OMEB) then     
    ! call crocus radiation routines to compute PTAU_N
    do j = 1,KSIZE1
        ZP_SWNETSNOW(j)  = SWsrf(1,j)
        ZP_SWNETSNOWS(j) = SWsrf(1,j)!*(1-PTAU_N(j,1))
        ZP_HSNOW(j)      = Hsrf(1,j)
        ZP_LES3L(j)      = LEsrf(1,j)
        ZP_LEL3L(j)      = 0.
        ZP_RNSNOW(j)     = Rsrf(1,j)
        ZP_LWNETSNOW(j)  = Rsrf(1,j) - SWsrf(1,j) 
    end do 
else ! if MEB switched off, initialization with 0 within SNOWCRO
    ZP_GFLUXSNOW = 0.
    ZP_HSNOW     = 0.
    ZP_LEL3L     = 0.
    ZP_LES3L     = 0.
    ZP_RNSNOW    = 0.
    ZP_SWNETSNOW =  0.
    ZP_SWNETSNOWS = 0. 
    ZP_LWNETSNOW  = 0.
end if

! State variables coupling (internal)
do j = 1,KSIZE1
    ZP_SOILCOND(j) = ksoil(1,1,j)
end do

! Intent in, has been changed between SNOWCRO and now
do j = 1, KSIZE1
    ZP_TG(j,1) = Tsoil(1,1,j)
end do

! PACK VARIABLES 
! Currently, this change of array size is applied to all variables for the sake of simplicity 
! (and to keep a coding that is somewhat consistent with SURFEX / Crocus)
! This can be sorted much more easily and eleganty in the longer term 
! In Crocus, some variables are only added after the Crocus step to be passed on
! further; in a code cleaning effort, I should really figure out what is needed where and\
! remove obsolete and duplicate variables
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    X_PEW_A_COEF(JJ) = ZP_PEW_A_COEF(JI) 
    X_PET_A_COEF(JJ) = ZP_PET_A_COEF(JI) 
    X_PEQ_A_COEF(JJ) = ZP_PEQ_A_COEF(JI) 
    X_PEW_B_COEF(JJ) = ZP_PEW_B_COEF(JI) 
    X_PET_B_COEF(JJ) = ZP_PET_B_COEF(JI) 
    X_PEQ_B_COEF(JJ) = ZP_PEQ_B_COEF(JI) 
    X_SNOWALB(JJ)    = ZP_SNOWALB(JI) 
    X_PS(JJ)         = ZP_PS(JI)   
    X_SRSNOW(JJ)     = ZP_SRSNOW(JI)          
    X_RRSNOW(JJ)     = ZP_RRSNOW(JI)  
    X_PSN3L(JJ)      = ZP_PSN3L(JI)    
    X_TA(JJ)         = ZP_TA(JI)     
    X_SW_RAD(JJ)     = ZP_SW_RAD(JI)  
    X_QA(JJ)         = ZP_QA(JI) 
    X_VMOD(JJ)       = ZP_VMOD(JI)  
    X_LW_RAD(JJ)     = ZP_LW_RAD(JI)       
    X_RHOA(JJ)       = ZP_RHOA(JI)  
    X_UREF(JJ)       = ZP_UREF(JI) 
    X_EXNS(JJ)       = ZP_EXNS(JI) 
    X_EXNA(JJ)       = ZP_EXNA(JI) 
    X_DIRCOSZW(JJ)   = ZP_DIRCOSZW(JI)   
    X_ZREF(JJ)       = ZP_ZREF(JI) 
    X_Z0NAT(JJ)      = ZP_Z0NAT(JI) 
    X_Z0EFF(JJ)      = ZP_Z0EFF(JI) 
    X_Z0HNAT(JJ)     = ZP_Z0HNAT(JI) 
    X_ALB(JJ)        = ZP_ALB(JI) 
    X_SOILCOND(JJ)   = ZP_SOILCOND(JI) 
    X_THRUFAL(JJ)    = ZP_THRUFAL(JI)  
    X_GRNDFLUX(JJ)   = ZP_GRNDFLUX(JI) 
    X_EVAPCOR(JJ)    = ZP_EVAPCOR(JI)   
    X_GFLXCOR(JJ)    = ZP_GFLXCOR(JI)      
    X_SWNETSNOW(JJ)  = ZP_SWNETSNOW(JI)  
    X_SWNETSNOWS(JJ) = ZP_SWNETSNOWS(JI) 
    X_LWNETSNOW(JJ)  = ZP_LWNETSNOW(JI) 
    X_RNSNOW(JJ)     = ZP_RNSNOW(JI)    
    X_HSNOW(JJ)      = ZP_HSNOW(JI)     
    X_GFLUXSNOW(JJ)  = ZP_GFLUXSNOW(JI)     
    X_HPSNOW(JJ)     = ZP_HPSNOW(JI)        
    X_LES3L(JJ)      = ZP_LES3L(JI)    
    X_LEL3L(JJ)      = ZP_LEL3L(JI)       
    X_EVAP(JJ)       = ZP_EVAP(JI)             
    X_SNDRIFT(JJ)    = ZP_SNDRIFT(JI)     
    X_RI(JJ)         = ZP_RI(JI)     
    X_EMISNOW(JJ)    = ZP_EMISNOW(JI)   
    X_CDSNOW(JJ)     = ZP_CDSNOW(JI)   
    X_USTARSNOW(JJ)  = ZP_USTARSNOW(JI)  
    X_CHSNOW(JJ)     = ZP_CHSNOW(JI)   
    X_SNOWHMASS(JJ)  = ZP_SNOWHMASS(JI)     
    X_QS(JJ)         = ZP_QS(JI)           
    X_VEGTYPE(JJ)    = ZP_VEGTYPE(JI)       
    X_ZENITH(JJ)     = ZP_ZENITH(JI) 
    X_ANGL_ILLUM(JJ) = ZP_ANGL_ILLUM(JI) 
    X_LAT(JJ)        = ZP_LAT(JI) 
    X_LON(JJ)        = ZP_LON(JI) 
    X_GSFCSNOW(JJ)   = ZP_GSFCSNOW(JI) 
    X_SNOWMAK(JJ)    = ZP_SNOWMAK(JI) 
    X_SNDPT_12H(JJ)  = ZP_SNDPT_12H(JI) 
    X_SNDPT_1DY(JJ)  = ZP_SNDPT_1DY(JI) 
    X_SNDPT_3DY(JJ)  = ZP_SNDPT_3DY(JI) 
    X_SNDPT_5DY(JJ)  = ZP_SNDPT_5DY(JI) 
    X_SNDPT_7DY(JJ)  = ZP_SNDPT_7DY(JI) 
    X_SNSWE_1DY(JJ)  = ZP_SNSWE_1DY(JI) 
    X_SNSWE_3DY(JJ)  = ZP_SNSWE_3DY(JI) 
    X_SNSWE_5DY(JJ)  = ZP_SNSWE_5DY(JI) 
    X_SNSWE_7DY(JJ)  = ZP_SNSWE_7DY(JI) 
    X_SNRAM_SONDE(JJ) = ZP_SNRAM_SONDE(JI) 
    X_SN_WETTHCKN(JJ) =  ZP_SN_WETTHCKN(JI) 
    X_SN_REFRZNTHCKN(JJ) = ZP_SN_REFRZNTHCKN(JI)
    X_DEP_HIG(JJ)    = ZP_DEP_HIG(JI) 
    X_DEP_MOD(JJ)    = ZP_DEP_MOD(JI) 
    X_DEP_SUP(JJ)    = ZP_DEP_SUP(JI) 
    X_DEP_TOT(JJ)    = ZP_DEP_TOT(JI) 
    X_DEP_HUM(JJ)    = ZP_DEP_HUM(JI) 
    X_ACC_LEV(JJ)    = ZP_ACC_LEV(JI) 
    X_NAT_LEV(JJ)    = ZP_NAT_LEV(JI) 
    X_PRO_SUP_TYP(JJ) = ZP_PRO_SUP_TYP(JI) 
    X_PRO_INF_TYP(JJ) = ZP_PRO_INF_TYP(JI) 
    X_AVA_TYP(JJ)    = ZP_AVA_TYP(JI)  
end do

do k = 1,KSIZE2
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    X_SNOWSWE(JJ,k)   = ZP_SNOWSWE(JI,k)
    X_SNOWRHO(JJ,k)   = ZP_SNOWRHO(JI,k)
    X_SNOWHEAT(JJ,k)  = ZP_SNOWHEAT(JI,k)
    X_SNOWDIAMOPT(JJ,k) = ZP_SNOWDIAMOPT(JI,k)
    X_SNOWSPHERI(JJ,k) = ZP_SNOWSPHERI(JI,k) 
    X_SNOWHIST(JJ,k)  = ZP_SNOWHIST(JI,k) 
    X_SNOWAGE(JJ,k)   = ZP_SNOWAGE(JI,k) 
    X_SNOWLIQ(JJ,k)   = ZP_SNOWLIQ(JI,k)
    X_SNOWTEMP(JJ,k)  = ZP_SNOWTEMP(JI,k)
    X_SNOWDZ(JJ,k)    = ZP_SNOWDZ(JI,k) 
    X_SNOWDEND(JJ,k)  = ZP_SNOWDEND(JI,k) 
    X_SNOWSPHER(JJ,k) = ZP_SNOWSPHER(JI,k) 
    X_SNOWSIZE(JJ,k)  = ZP_SNOWSIZE(JI,k) 
    X_SNOWSSA(JJ,k)   = ZP_SNOWSSA(JI,k) 
    X_SNOWTYPEMEPRA(JJ,k) = ZP_SNOWTYPEMEPRA(JI,k) 
    X_SNOWRAM(JJ,k)   = ZP_SNOWRAM(JI,k) 
    X_SNOWSHEAR(JJ,k) = ZP_SNOWSHEAR(JI,k) 
    X_ACC_RAT(JJ,k)   = ZP_ACC_RAT(JI,k) 
    X_NAT_RAT(JJ,k)   = ZP_NAT_RAT(JI,k) 
end do 
end do

do ik = 1,NIMPUR
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    X_IMPDRY(JJ,ik)   = ZP_IMPDRY(JI,ik)
    X_IMPWET(JJ,ik)   = ZP_IMPWET(JI,ik) 
    do k = 1,KSIZE2
        X_SNOWIMPUR(JJ,k,ik) = ZP_SNOWIMPUR(JI,k,ik) 
        X_SNOWIMP_CONC(JJ,k,ik) = ZP_SNOWIMP_CONC(JI,k,ik)
    end do
end do 
end do 

do ik = 1,KSIZE3
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    X_TG(JJ,ik)  = ZP_TG(JI,ik)
    X_D_G(JJ,ik) = ZP_D_G(JI,ik)
end do
end do

do ik = 1,KSIZE4
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    X_BLOWSNW(JJ,ik) = ZP_BLOWSNW(JI,ik) 
end do 
end do

do ik = 1,2
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    X_DIR_SW(JJ,ik)     = ZP_DIR_SW(JI,ik)
    X_SCA_SW(JJ,ik)     = ZP_SCA_SW(JI,ik)
    X_SPEC_ALB(JJ,ik)   = ZP_SPEC_ALB(JI,ik)
    X_DIFF_RATIO(JJ,ik) = ZP_DIFF_RATIO(JI,ik)
    X_SPEC_TOT(JJ,ik)   = ZP_SPEC_TOT(JI,ik)
end do 
end do

write(*,*) CSNOWMETAMO, CSNOWFALL, CSNOWCOND
! SNOW MODEL CALL ONLY WHERE THERE IS SNOW  
CALL SNOWCRO(CSNOWRES, TPTIME, OMEB, LGLACIER, HIMPLICIT_WIND,                &
                X_PEW_A_COEF, X_PEW_B_COEF, X_PET_A_COEF, X_PEQ_A_COEF,   &
                X_PET_B_COEF, X_PEQ_B_COEF, X_SNOWSWE, X_SNOWRHO,         &
                X_SNOWHEAT, X_SNOWALB, X_SNOWDIAMOPT, X_SNOWSPHERI,          &
                X_SNOWHIST, X_SNOWAGE, X_SNOWIMPUR, PTSTEP, X_PS,          &
                X_SRSNOW, X_RRSNOW, X_PSN3L, X_TA, X_TG(:,1), X_SW_RAD, &
                X_QA,X_VMOD, X_LW_RAD, X_RHOA, X_UREF, X_EXNS, X_EXNA, &
                X_DIRCOSZW, X_ZREF, X_Z0NAT, X_Z0EFF, X_Z0HNAT,          &
                X_ALB, X_SOILCOND, X_D_G(:,1), X_SNOWLIQ, X_SNOWTEMP,    &
                X_SNOWDZ, X_THRUFAL, X_GRNDFLUX, X_EVAPCOR, X_GFLXCOR,   &
                X_SWNETSNOW, X_SWNETSNOWS, X_LWNETSNOW, X_RNSNOW,         &
                X_HSNOW, X_GFLUXSNOW, X_HPSNOW, X_LES3L, X_LEL3L,        &
                X_EVAP, X_SNDRIFT, X_RI, X_EMISNOW, X_CDSNOW,             &
                X_USTARSNOW, X_CHSNOW, X_SNOWHMASS, X_QS, X_VEGTYPE,     &
                X_ZENITH, X_ANGL_ILLUM, X_LAT, X_LON, X_BLOWSNW,         &
                CSNOWDRIFT, LSNOWDRIFT_SUBLIM, LSNOW_ABS_ZENITH,              &
                CSNOWMETAMO,CSNOWRAD, LATMORAD, X_DIR_SW, X_SCA_SW,         &
                X_SPEC_ALB, X_DIFF_RATIO, X_SPEC_TOT, X_GSFCSNOW,           &
                X_IMPWET, X_IMPDRY, CSNOWFALL, CSNOWCOND, CSNOWHOLD,         &
                CSNOWCOMP, CSNOWZREF, X_SNOWMAK, LSNOWCOMPACT_BOOL,          &
                LSNOWMAK_BOOL, LSNOWTILLER, LSELF_PROD, LSNOWMAK_PROP)


! Purely diagnostic, this routine can be called only at output time steps
CALL SNOWCRO_DIAG(CSNOWHOLD, CSNOWMETAMO, X_SNOWDZ, X_SNOWSWE, &
                X_SNOWRHO, X_SNOWDIAMOPT, X_SNOWSPHERI, X_SNOWAGE,  &
                X_SNOWHIST, X_SNOWTEMP, X_SNOWLIQ, X_DIRCOSZW, X_SNOWIMPUR, X_SNOWDEND, &
                X_SNOWSPHER, X_SNOWSIZE, X_SNOWSSA, X_SNOWTYPEMEPRA, X_SNOWRAM, X_SNOWSHEAR, &
                X_ACC_RAT, X_NAT_RAT, X_SNDPT_12H, X_SNDPT_1DY, X_SNDPT_3DY, X_SNDPT_5DY, &
                X_SNDPT_7DY, X_SNSWE_1DY, X_SNSWE_3DY, X_SNSWE_5DY, X_SNSWE_7DY, &
                X_SNRAM_SONDE, X_SN_WETTHCKN, X_SN_REFRZNTHCKN, X_SNOWIMP_CONC,  &
                X_DEP_HIG, X_DEP_MOD, X_DEP_SUP, X_DEP_TOT, X_DEP_HUM, &
                X_ACC_LEV, X_NAT_LEV, X_PRO_SUP_TYP, X_PRO_INF_TYP, X_AVA_TYP)

! UNPACK VARIABLES -> only those that are really outputs
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    ZP_SNOWALB(JI)   = X_SNOWALB(JJ) 
    ZP_THRUFAL(JI)   = X_THRUFAL(JJ)    
    ZP_GRNDFLUX(JI)  = X_GRNDFLUX(JJ)
    ZP_GFLXCOR(JI)   = X_GFLXCOR(JJ) 
    ZP_GFLUXSNOW(JI) = X_GFLUXSNOW(JJ)       
    ZP_EVAP(JI)      = X_EVAP(JJ)                   
    ZP_SNDRIFT(JI)   = X_SNDRIFT(JJ)        
    ZP_EMISNOW(JI)   = X_EMISNOW(JJ)       
    ZP_QS(JI)        = X_QS(JJ)                
    ZP_GSFCSNOW(JI)  = X_GSFCSNOW(JJ)  
    ! additional diagnostics
    ZP_HSNOW(JI)     = X_HSNOW(JJ)  
    ZP_LES3L(JI)     = X_LES3L(JJ)
    ZP_RNSNOW(JI)    = X_RNSNOW(JJ)
end do

do k = 1,KSIZE2
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    ZP_SNOWSWE(JI,k)   = X_SNOWSWE(JJ,k)    
    ZP_SNOWRHO(JI,k)   = X_SNOWRHO(JJ,k)   
    ZP_SNOWHEAT(JI,k)  = X_SNOWHEAT(JJ,k)  
    ZP_SNOWDIAMOPT(JI,k) = X_SNOWDIAMOPT(JJ,k) 
    ZP_SNOWSPHERI(JI,k) = X_SNOWSPHERI(JJ,k) 
    ZP_SNOWHIST(JI,k)  = X_SNOWHIST(JJ,k)  
    ZP_SNOWAGE(JI,k)   = X_SNOWAGE(JJ,k)   
    ZP_SNOWLIQ(JI,k)   = X_SNOWLIQ(JJ,k)   
    ZP_SNOWTEMP(JI,k)  = X_SNOWTEMP(JJ,k)  
    ZP_SNOWDZ(JI,k)    = X_SNOWDZ(JJ,k)    

    ZP_SNOWSSA(JI,k)    = X_SNOWSSA(JJ,k)
    ZP_SNOWTYPEMEPRA(JI,k) =  X_SNOWTYPEMEPRA(JJ,k)
end do 
end do

do ik = 1,NIMPUR
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    do k = 1,KSIZE2
        ZP_SNOWIMPUR(JI,k,ik) = X_SNOWIMPUR(JJ,k,ik) 
        ZP_SNOWIMP_CONC(JI,k,ik) = X_SNOWIMP_CONC(JJ,k,ik) 
    end do
end do 
end do 

! VARIABLE COUPLING AND UPDATES - AFTER CROCUS UPDATES 
! In case of vanishing snowpack, remove spurious snow amounts 
do j = 1, KSIZE1
    dsout(j) = 0.
    sweout(j) = 0.
    do k = 1, KSIZE2
        dsout(j) = dsout(j) + ZP_SNOWDZ(j,k)
        sweout(j) = sweout(j) + ZP_SNOWSWE(j,k)
    end do
    if (dsout(j) > 0. .AND. dsout(j) < 1.1*EXSNOWDMIN) then 
        dsout = 0.
        ZP_THRUFAL(j) = ZP_THRUFAL(j) + sweout(j)/PTSTEP
        ! ZP_GRNDFLUX(j) !!!! check these corrections with matthieu
        ! ZP_EVAP(j)
        ZP_SNOWSWE(j,:) = 0
    end if 
end do

! Update FSM snow state variables 
fsnow = 0
Nsnow = 0
do j = 1, KSIZE1
    albs(1,j) = ZP_SNOWALB(j) 
    if (dsout(j) > 0.) then 
        fsnow(1,j) = 1
        do k = 1,KSIZE2
            if (ZP_SNOWSWE(j,k) > 0.) then
                Nsnow(1,j) = Nsnow(1,j) + 1
            end if 
        end do
    end if 
end do

do k = 1, KSIZE2
    do j = 1, KSIZE1
        Ds(k,1,j)    = ZP_SNOWDZ(j,k)
        Sliq(k,1,j)  = ZP_SNOWLIQ(j,k) 
        Sice(k,1,j)  = ZP_SNOWSWE(j,k) - ZP_SNOWLIQ(j,k) 
        if (Sice(k,1,j) > 0.) then
            Tsnow(k,1,j) = ZP_SNOWTEMP(j,k) 
        else 
            Tsnow(k,1,j) = Tm
        end if 
    end do
end do

! Other coupled FSM variables 
do j = 1, KSIZE1
    Gsoil(1,j)        = ZP_GRNDFLUX(j) ! ground flux already includes radiation
    meltflux_out(1,j) = (ZP_THRUFAL(j) + ZP_EVAPCOR(j))*dt
    Roff(1,j)         = meltflux_out(1,j)*dt
    Sbsrf(1,j)        = ZP_EVAP(j)
end do

! Reinitialize Crocus variable where there is no (more) snow
do j = 1,KSIZE1
    if (dsout(j) < epsilon(dsout(j))) then 
        ZP_SNOWSWE(j,:) = 0
        ZP_SNOWDZ(j,:) = 0
        ZP_SNOWHEAT(j,:) = EXUNDEF
        ZP_SNOWRHO(j,:) = EXUNDEF
        ZP_SNOWAGE(j,:) = EXUNDEF
        ZP_SNOWDIAMOPT(j,:) = EXUNDEF
        ZP_SNOWSPHERI(j,:) = EXUNDEF
        ZP_SNOWHIST(j,:) = EXUNDEF
        ZP_SNOWLIQ(j,:)= EXUNDEF
        ZP_SNOWTEMP(j,:) = EXUNDEF
     !   ZP_SNOWALB(j) = EXUNDEF
        ZP_SNOWIMPUR(j,:,:) = EXUNDEF
        Ds(:,1,j)    = 0.
        Sliq(:,1,j)  = 0.
        Sice(:,1,j)  = 0.
        ZP_SNOWTYPEMEPRA(j,:) = EXUNDEF
        ZP_SNOWSSA(j,:) = EXUNDEF
    end if 
end do

! update other states and fluxes that are going to be output as diagnostics
do JJ = 1,ISIZE_SNOW
    JI = NMASK(JJ)
    Hsrf(1,JI)  = ZP_HSNOW(JI)    
    LEsrf(1,JI) = ZP_LES3L(JI)
    Rsrf(1,JI)  = ZP_RNSNOW(JI)    
end do 

do j = 1, KSIZE1
    SWsci(1,j)  = ZP_SW_RAD(j)
    LWsci(1,j)  = ZP_LW_RAD(j)
end do 

! update crosnowmask for use in SNOW.F90
do j = 1,KSIZE1
    if (Nsnow(1,j) > 0) then 
        crosnowmask(1,j) = 1
    else
        crosnowmask(1,j) = 0
    end if 
end do 


END SUBROUTINE CRO_COUP
