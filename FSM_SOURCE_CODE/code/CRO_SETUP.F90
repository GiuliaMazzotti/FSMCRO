SUBROUTINE CRO_SETUP

! Import CROCUS modules (from CRO_MODUL)
use CRO_CSTS
USE CRO_SETTINGS
USE CRO_DRIVE
USE CRO_GRID
USE CRO_PARAMS
USE CRO_SURF
USE CRO_STATE_VARS
USE CRO_DUMMIES
USE CRO_DIAG
USE CRO_NLSTPARAMS

! Import FSM modules (from MODULES): parameters that need coupling but are not changed in each time step 
USE MODCONF, only:          CRO_ON
USE MODTILE, only:          TILE
USE CONSTANTS, only:        undef, iundef
USE GRID, only:             Ny, Nsmax, Nsoil, Dzsoil
USE PARAMMAPS, only:        alb0
USE LANDUSE, only:          lat, lon
USE DRIVING, only:          dt, zT, zU
use SOILPARAMS, only:       hcap_soil
USE STATE_VARIABLES, only:  Ds, Sliq, Sice, Tsnow, Tsoil, albs

USE MODD_SNOW_PAR,  only : X_RI_MAX, XVAGING_NOGLACIER

IMPLICIT NONE

integer :: &
  j,               &  ! Point counters
  k                   ! Level counter

character(len=200) :: nlst_file

character(len=20) :: XCSNOWFALL, XCSNOWRES, XCSNOWCOND, XCSNOWHOLD, &
                    XCSNOWCOMP, XCSNOWRAD, XCSNOWMETAMO

real:: XX_RI_MAX, XXCVHEATF, XXVAGING_NOGLACIER

namelist  /nam_crooptions/ XCSNOWFALL, XCSNOWRES, XCSNOWCOND, XCSNOWHOLD, &
                          XCSNOWCOMP, XCSNOWRAD, XCSNOWMETAMO! Namelist defining ESCROC options 
namelist  /nam_croparams/ XX_RI_MAX, XXCVHEATF, XXVAGING_NOGLACIER ! Namelist defining ESCROC parameters 

if (CRO_ON) then 

  call getarg(1, nlst_file)
  open(5000, file = nlst_file)

  ! Assign values to modules 

  ! CRO_SETTINGS 
  ! Options are set to operational version 
  XCSNOWRES       = 'RIL'   !'RIL', 'M98', Deterministic Crocus: RI1, operational: RI2
  XCSNOWFALL      = 'V12'   !'V12','S02','A76','P75','NZE'
  XCSNOWCOND      = 'Y81'   !'Y81', 'I02','C11'
  XCSNOWHOLD      = 'B92'   !'B92', 'SPK','O04','B02'
  XCSNOWCOMP      = 'B92'   !'B92'
  XCSNOWMETAMO    = 'B21'   !'B92','C13','F06','T07','S-C','S-F', new: 'B21'; Deterministic Crocus: C13
  XCSNOWRAD       = 'B92'   !'B92','B93','T17'; Deterministic Crocus: B60, 
  CSNOWZREF      = 'CST'   !'CST','VAR'
  CSNOWDRIFT     = 'VI13'  !'NONE','DFLT','VI13','GA01'
  HIMPLICIT_WIND = 'OLD'   ! Direct wind implicitation
  OMEB               = .FALSE.
  LGLACIER           = .FALSE.
  if (TILE == 'glacier') LGLACIER = .TRUE.
  LSNOWDRIFT_SUBLIM  = .FALSE.
  LSNOW_ABS_ZENITH   = .FALSE.
  LATMORAD           = .FALSE.
  LSNOWCOMPACT_BOOL  = .FALSE.
  LSNOWMAK_BOOL      = .FALSE.
  LSNOWTILLER        = .FALSE.
  LSELF_PROD         = .FALSE.
  LSNOWMAK_PROP      = .FALSE.

  ! check for options input from namelist 
  read(5000, nam_crooptions)
  CSNOWFALL  = XCSNOWFALL
  CSNOWRES   = XCSNOWRES
  CSNOWCOND  = XCSNOWCOND
  CSNOWHOLD  = XCSNOWHOLD
  CSNOWCOMP  = XCSNOWCOMP
  CSNOWRAD   = XCSNOWRAD
  CSNOWMETAMO  = XCSNOWMETAMO

  ! CRO_GRID
  KSIZE1 = Ny  
  KSIZE2 = Nsmax 
  KSIZE3 = Nsoil 
  KSIZE4 = 1
  NIMPUR = 0

  ! CRO_DRIVE 
  PTSTEP = dt  
  allocate(ZP_ZREF(KSIZE1)) 
  allocate(ZP_UREF(KSIZE1)) 
  allocate(ZP_PS(KSIZE1)) 
  allocate(ZP_SRSNOW(KSIZE1)) 
  allocate(ZP_RRSNOW(KSIZE1)) 
  allocate(ZP_SW_RAD(KSIZE1)) 
  allocate(ZP_LW_RAD(KSIZE1)) 
  allocate(ZP_TA(KSIZE1)) 
  allocate(ZP_QA(KSIZE1)) 
  allocate(ZP_VMOD(KSIZE1)) 
  ZP_ZREF(:) = zT !1.5 SURFEX default
  ZP_UREF(:) = zU !2.  SURFEX default

  ! CRO_PARAMS
  allocate(ZP_ALB(KSIZE1)) 
  allocate(ZP_Z0NAT(KSIZE1)) 
  allocate(ZP_Z0EFF(KSIZE1)) 
  allocate(ZP_Z0HNAT(KSIZE1)) 
  allocate(ZP_PEW_A_COEF(KSIZE1)) 
  allocate(ZP_PET_A_COEF(KSIZE1)) 
  allocate(ZP_PEQ_A_COEF(KSIZE1)) 
  ZP_ALB(:)        = alb0(1,1)
  ZP_Z0NAT(:)      = 0.001    ! XZ0SN from NAM_SURF_CSTSm pure snow surface, potentially couple to FSM value later
  ZP_Z0HNAT(:)     = 0.0001   ! XZ0HSN
  ZP_Z0EFF(:)      = ZP_Z0NAT 
  ZP_PEW_A_COEF(:) = 0.       ! = 0 in case of explicit coupling / offline 
  ZP_PET_A_COEF(:) = 0.       ! = 0 in case of explicit coupling / offline 
  ZP_PEQ_A_COEF(:) = 0.       ! = 0 in case of explicit coupling / offline 

  ! CRO_SURF
  allocate(ZP_D_G(KSIZE1,KSIZE3)) 
  allocate(ZP_LAT(KSIZE1)) 
  allocate(ZP_LON(KSIZE1)) 
  allocate(ZP_DIRCOSZW(KSIZE1)) 
  do k = 1, Nsoil
    ZP_D_G(:,k) = Dzsoil(k)
  end do
  ZP_LAT(:)       = lat(1,:)
  ZP_LON(:)       = lon(1,:)
  ZP_DIRCOSZW(:)  = 1.         !'slope',ACOS(PDIRCOSZW(JJ))*(180./XPI),"deg" - check coupling with slopemu

  ! CRO_STATE_VARS
  allocate(ZP_SNOWSWE(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWRHO(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWHEAT(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWDIAMOPT(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWSPHERI(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWHIST(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWAGE(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWIMPUR(KSIZE1,KSIZE2,NIMPUR)) 
  allocate(ZP_SNOWDZ(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWTEMP(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWLIQ(KSIZE1,KSIZE2)) 
  allocate(ZP_SNOWALB(KSIZE1))
  allocate(ZP_TG(KSIZE1,KSIZE3))  
  allocate(CROSNOW_ON(KSIZE1))
  ! use values from NAM_PREP_ISBA_SNOW where available
  ZP_SNOWSWE(:,:)     = 0.    
  ZP_SNOWRHO(:,:)     = 300.
  ZP_SNOWHEAT(:,:)    = 0.
  ZP_SNOWDIAMOPT(:,:)   = 0. 
  ZP_SNOWSPHERI(:,:)   = 0.
  ZP_SNOWHIST(:,:)    = 0.
  ZP_SNOWAGE(:,:)     = 0.
  ZP_SNOWIMPUR(:,:,:) = 0.
  ZP_SNOWDZ(:,:)      = 0.
  ZP_SNOWTEMP(:,:)    = 273.16
  ZP_SNOWLIQ(:,:)     = 0.
  ZP_SNOWALB(:)       = 0.5
  CROSNOW_ON          = .TRUE. ! For now this is a dummy, have not implemented simultaneous treatment of multiple points zet

  ! CRO_DUMMIES
  ! Need to exist in FSMCRO but in fact SNOWCRO internal variables, unaffected by the coupling
  ! Some variables really set to dummies because relevant crocus option not used 
  allocate(ZP_ANGL_ILLUM(KSIZE1))
  allocate(ZP_BLOWSNW(KSIZE1,KSIZE4))
  allocate(ZP_CDSNOW(KSIZE1))
  allocate(ZP_CHSNOW(KSIZE1))
  allocate(ZP_DIFF_RATIO(KSIZE1,2)) 
  allocate(ZP_DIR_SW(KSIZE1,2)) 
  allocate(ZP_EMISNOW(KSIZE1))
  allocate(ZP_EVAP(KSIZE1))
  allocate(ZP_EVAPCOR(KSIZE1))
  allocate(ZP_GFLXCOR(KSIZE1))
  allocate(ZP_GFLUXSNOW(KSIZE1))
  allocate(ZP_GRNDFLUX(KSIZE1))
  allocate(ZP_GSFCSNOW(KSIZE1))
  allocate(ZP_HPSNOW(KSIZE1))
  allocate(ZP_IMPDRY(KSIZE1,NIMPUR)) 
  allocate(ZP_IMPWET(KSIZE1,NIMPUR)) 
  allocate(ZP_PSN3L(KSIZE1)) 
  allocate(ZP_QS(KSIZE1))          
  allocate(ZP_RI(KSIZE1))   
  allocate(ZP_SCA_SW(KSIZE1,2)) 
  allocate(ZP_SNDRIFT(KSIZE1))       
  allocate(ZP_SNOWHMASS(KSIZE1)) 
  allocate(ZP_SNOWMAK(KSIZE1)) 
  allocate(ZP_SPEC_ALB(KSIZE1,2))
  allocate(ZP_SPEC_TOT(KSIZE1,2))    
  allocate(ZP_THRUFAL(KSIZE1))  
  allocate(ZP_VEGTYPE(KSIZE1)) 
  allocate(ZP_USTARSNOW(KSIZE1)) 

  ! initializations as in snow3L_isba.F90
  ZP_ANGL_ILLUM   = EXUNDEF
  ZP_BLOWSNW      = 0.
  ZP_CDSNOW       = 0.
  ZP_CHSNOW       = 0.
  ZP_DIFF_RATIO   = EXUNDEF   ! in snow3L_isba undef
  ZP_DIR_SW       = 0.
  ZP_EMISNOW      = EXEMISSN  ! from snow3L_isba
  ZP_EVAP         = 0.        ! init 0 within SNOWCRO
  ZP_EVAPCOR      = 0.
  ZP_GFLXCOR      = 0.
  ZP_GFLUXSNOW    = 0.
  ZP_GRNDFLUX     = 0.        ! init 0 within SNOWCRO if MEB off, but fine for MEB too
  ZP_GSFCSNOW     = 0.
  ZP_HPSNOW       = 0.
  ZP_IMPDRY       = 0.
  ZP_IMPWET       = 0.
  ZP_PSN3L        = 1.
  ZP_QS           = EXUNDEF   ! in snow3L_isba undef, init 0 in SNOWCRO
  ZP_RI           = EXUNDEF   ! in snow3L_isba undef
  ZP_SCA_SW       = 1. 
  ZP_SNDRIFT      = 0.
  ZP_SNOWHMASS    = 0.        ! init 0 within SNOWCRO
  ZP_SNOWMAK      = 0.
  ZP_SPEC_ALB     = EXUNDEF   ! in snow3L_isba undef
  ZP_SPEC_TOT     = EXUNDEF   ! in snow3L_isba undef
  ZP_THRUFAL      = 0.        ! init 0 within SNOWCRO
  ZP_VEGTYPE      = 0.
  ZP_USTARSNOW    = 0.

  ! CRO_DIAG
  allocate(ZP_SNOWDEND(KSIZE1,KSIZE2))
  allocate(ZP_SNOWSPHER(KSIZE1,KSIZE2))
  allocate(ZP_SNOWSIZE(KSIZE1,KSIZE2))
  allocate(ZP_SNOWSSA(KSIZE1,KSIZE2))
  allocate(ZP_SNOWTYPEMEPRA(KSIZE1,KSIZE2))
  allocate(ZP_SNOWRAM(KSIZE1,KSIZE2))
  allocate(ZP_SNOWSHEAR(KSIZE1,KSIZE2))
  allocate(ZP_ACC_RAT(KSIZE1,KSIZE2))
  allocate(ZP_NAT_RAT(KSIZE1,KSIZE2))
  allocate(ZP_SNDPT_12H(KSIZE1)) 
  allocate(ZP_SNDPT_1DY(KSIZE1)) 
  allocate(ZP_SNDPT_3DY(KSIZE1)) 
  allocate(ZP_SNDPT_5DY(KSIZE1)) 
  allocate(ZP_SNDPT_7DY(KSIZE1)) 
  allocate(ZP_SNSWE_1DY(KSIZE1)) 
  allocate(ZP_SNSWE_3DY(KSIZE1)) 
  allocate(ZP_SNSWE_5DY(KSIZE1)) 
  allocate(ZP_SNSWE_7DY(KSIZE1)) 
  allocate(ZP_SNRAM_SONDE(KSIZE1)) 
  allocate(ZP_SN_WETTHCKN(KSIZE1)) 
  allocate(ZP_SN_REFRZNTHCKN(KSIZE1)) 
  allocate(ZP_DEP_HIG(KSIZE1)) 
  allocate(ZP_DEP_MOD(KSIZE1)) 
  allocate(ZP_DEP_SUP(KSIZE1)) 
  allocate(ZP_DEP_TOT(KSIZE1)) 
  allocate(ZP_DEP_HUM(KSIZE1)) 
  allocate(ZP_ACC_LEV(KSIZE1)) 
  allocate(ZP_NAT_LEV(KSIZE1)) 
  allocate(ZP_PRO_SUP_TYP(KSIZE1)) 
  allocate(ZP_PRO_INF_TYP(KSIZE1)) 
  allocate(ZP_AVA_TYP(KSIZE1)) 
  allocate(ZP_SNOWIMP_CONC(KSIZE1,KSIZE2,NIMPUR))

  ZP_SNOWDEND = EXUNDEF
  ZP_SNOWSPHER= EXUNDEF
  ZP_SNOWSIZE= EXUNDEF
  ZP_SNOWSSA= EXUNDEF
  ZP_SNOWTYPEMEPRA= EXUNDEF
  ZP_SNOWRAM= EXUNDEF
  ZP_SNOWSHEAR= EXUNDEF
  ZP_ACC_RAT= EXUNDEF
  ZP_NAT_RAT= EXUNDEF
  ZP_SNDPT_12H= EXUNDEF
  ZP_SNDPT_1DY= EXUNDEF
  ZP_SNDPT_3DY= EXUNDEF
  ZP_SNDPT_5DY= EXUNDEF
  ZP_SNDPT_7DY= EXUNDEF
  ZP_SNSWE_1DY= EXUNDEF
  ZP_SNSWE_3DY= EXUNDEF
  ZP_SNSWE_5DY= EXUNDEF
  ZP_SNSWE_7DY= EXUNDEF
  ZP_SNRAM_SONDE= EXUNDEF
  ZP_SN_WETTHCKN= EXUNDEF
  ZP_SN_REFRZNTHCKN= EXUNDEF
  ZP_DEP_HIG= EXUNDEF
  ZP_DEP_MOD= EXUNDEF
  ZP_DEP_SUP= EXUNDEF
  ZP_DEP_TOT= EXUNDEF
  ZP_DEP_HUM= EXUNDEF
  ZP_ACC_LEV= EXUNDEF
  ZP_NAT_LEV= EXUNDEF
  ZP_PRO_SUP_TYP= EXUNDEF
  ZP_PRO_INF_TYP= EXUNDEF
  ZP_AVA_TYP= EXUNDEF  
  ZP_SNOWIMP_CONC = EXUNDEF

  ! Initialization of Surfex constants
  CALL INI_CSTS

  ! Overwrite constants that are provided in the namelist (i.e. those dependent on the ensemble member)
  XX_RI_MAX = 0.026 ! RI2 option
  XXCVHEATF = 0.6 ! CV30000 option
  XXVAGING_NOGLACIER = 60 ! B60 option

  read(5000, nam_croparams)
  X_RI_MAX = XX_RI_MAX
  X_CVHEATF = XXCVHEATF 
  XVAGING_NOGLACIER = XXVAGING_NOGLACIER

  close(5000)

  ! apply the correction to the soil heat capacity used by the FSM soil scheme and initialized in SETUP
  hcap_soil(:,:) = X_CVHEATF/0.6*hcap_soil(:,:)  ! no change from the FSM behaviour for default X_CVHEATF value
  
end if 

END SUBROUTINE CRO_SETUP
