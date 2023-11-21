MODULE CRO_MODE_RAD
! subroutines taken from crocus to calculate radiation transfer through first snowpack layer
! GM, slight adaptations, no straightforward version control intended.

CONTAINS
!---------------------------------------------------------------------------------------

SUBROUTINE SNOWCROALB_SPECTRAL_BANDS_MEB(PVEGTYPE,PSNOWALB,PSNOWRHO,PSNOWAGE,  &
                                    PSNOWDIAMOPT,PSNOWSPHERI,PSNOWIMPUR, PPS,  &
                                    PSW_RAD,PPSN,PSNOWDZ,PZENITH,PSOILALB,     &
                                    PSNOWALBVIS,PSNOWALBNIR,PSNOWALBFIR,       &
                                    PTAU_N,HSNOWMETAMO,HSNOWRAD,OATMORAD)
! Split Total snow albedo into N-spectral bands. NOTE currently MEB only uses 2 bands of the 3 possible.
! Adapted GM for use in FSM_CRO (removed T17 options)

USE CRO_CSTS,         ONLY : EXUNDEF
USE CRO_RADTRANS,     ONLY : XSW_WGHT_VIS,XSW_WGHT_NIR
USE MODD_SNOW_METAMO, ONLY : XSNOWDZMIN
USE MODD_SNOW_PAR,    ONLY : NSPEC_BAND_SNOW
!
IMPLICIT NONE
!*      0.1    declarations of arguments
REAL, DIMENSION(:),   INTENT(IN)  :: PVEGTYPE, PSNOWALB      
REAL, DIMENSION(:,:), INTENT(IN)  :: PSNOWRHO, PSNOWDZ, PSNOWDIAMOPT, PSNOWSPHERI  
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSNOWIMPUR
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH, PSOILALB, PPSN      
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWAGE      
REAL, DIMENSION(:),   INTENT(IN)    :: PPS, PSW_RAD        
CHARACTER(3),         INTENT(IN)    :: HSNOWMETAMO, HSNOWRAD
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWALBVIS, PSNOWALBNIR, PSNOWALBFIR  
REAL, DIMENSION(:,:), INTENT(OUT)   :: PTAU_N        
LOGICAL, INTENT(IN)                 :: OATMORAD  

!*      0.2    declarations of local variables
INTEGER                             :: JJ, JI, INI, INLVLS
INTEGER, DIMENSION(SIZE(PZENITH))   :: INLVLS_USE
REAL, DIMENSION(SIZE(PPS))          :: ZWORK, ZWORKA, ZAGE
REAL, DIMENSION(SIZE(PPS))          :: ZPROJLAT, ZDSGRAIN, ZBETA1, ZBETA2, ZBETA3, &
    ZOPTICALPATH1, ZOPTICALPATH2, ZOPTICALPATH3
REAL, DIMENSION(SIZE(PPS))          :: ZPERMSNOWFRAC
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSNOWDZ,ZSNOWRHO,ZSNOWDIAMOPT,ZSNOWSPHERI,ZSNOWAGE
REAL, DIMENSION(SIZE(PPS),NSPEC_BAND_SNOW)       :: ZSPECTRALALBEDO ! spectral albedo 1=VIS, 2=NIR, 3=UV
REAL:: ZDEFAULTG1,ZDEFAULTG2                  !For now these values are constant
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWG0 ! asymmetry parameter of snow grains at nr=1.3 and at non absorbing wavelengths (no unit) (npoints,nlayer)
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWY0 ! Value of y of snow grains at nr=1.3 (no unit
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWW0 ! Value of W of snow grains at nr=1.3 (no unit)
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWB0 ! absorption enhancement parameter of snow grains at nr=1.3 and at non absorbing wavelengths 
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZRADSINK
REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZRADXS,ZZENITH,ZSW_RAD,ZSNOWALB

INI    = SIZE(PSNOWDZ,1)
INLVLS = SIZE(PSNOWDZ,2)
ZDEFAULTG1=0.5
ZDEFAULTG2=0.5  
INLVLS_USE(:) = 0
DO JJ = 1,SIZE(PSNOWDZ(:,:),2)
    DO JI = 1,SIZE(PSNOWDZ(:,:),1)
        IF ( PSNOWDZ(JI,JJ)>0. ) THEN 
            INLVLS_USE(JI) = JJ
            ZSNOWRHO(JI,JJ)=PSNOWRHO(JI,JJ)
            ZSNOWAGE(JI,JJ)=PSNOWAGE(JI,JJ)
            ZSNOWDIAMOPT(JI,JJ)=PSNOWDIAMOPT(JI,JJ)
            ZSNOWSPHERI(JI,JJ)=PSNOWSPHERI(JI,JJ)
        ELSE
            ! default values to avoid numerical problems in case of no snow
            ZSNOWRHO(JI,JJ)=400.
            ZSNOWAGE(JI,JJ)=10.
            ZSNOWDIAMOPT(JI,JJ)=ZDEFAULTG1
            ZSNOWSPHERI(JI,JJ)=ZDEFAULTG2
        ENDIF
    ENDDO  !  end loop snow layers
ENDDO    ! end loop grid points

IF ( (HSNOWRAD=="T17")) THEN
! deleted GM because not used yet
ELSE
! 1) Spectral albedo
! ------------------
ZWORK(:)         = 0.0
ZWORKA(:)        = PSNOWALB(:)
ZPERMSNOWFRAC(:) = PVEGTYPE(:)

CALL SNOWCROALB(.FALSE.,                                             &
ZWORKA,ZSPECTRALALBEDO,PSNOWDZ(:,1),ZSNOWRHO(:,1:2),       &
ZPERMSNOWFRAC,ZSNOWDIAMOPT(:,1),ZSNOWSPHERI(:,1),               &
ZSNOWAGE(:,1),ZSNOWDIAMOPT(:,2),ZSNOWSPHERI(:,2),ZSNOWAGE(:,2), &
PPS, PZENITH, INLVLS_USE, HSNOWMETAMO) 

! Since we only consider VIS and NIR bands for soil and veg in MEB currently:
! (also note, here PSNOWALB doesn't evolve...we just diagnose spectral components).
WHERE(PSNOWALB(:)/=EXUNDEF)
PSNOWALBVIS(:) = ZSPECTRALALBEDO(:,1)
! We diagnose NIR albedo such that total albedo is conserved
! (using just 2 spectral bands in MEB)
PSNOWALBNIR(:) = (PSNOWALB(:) - XSW_WGHT_VIS*PSNOWALBVIS(:))/XSW_WGHT_NIR
! currently NOT used by MEB
PSNOWALBFIR(:) = EXUNDEF                                     
! For the surface layer absorbtion computation:
ZSPECTRALALBEDO(:,1) = PSNOWALBVIS(:)
ZSPECTRALALBEDO(:,2) = PSNOWALBNIR(:)
ZSPECTRALALBEDO(:,3) = PSNOWALBFIR(:)
ELSEWHERE
PSNOWALBVIS(:) = EXUNDEF
PSNOWALBNIR(:) = EXUNDEF
PSNOWALBFIR(:) = EXUNDEF
END WHERE
!
! Snow optical grain diameter (no age dependency over polar regions):
! 2) SW absorption in uppermost snow layer 
! ----------------------------------------
! For now, consider just 2 bands with MEB, so renormalize:
ZSPECTRALALBEDO(:,1) = ZSPECTRALALBEDO(:,1)
ZSPECTRALALBEDO(:,2) = (PSNOWALB(:) - XSW_WGHT_VIS*ZSPECTRALALBEDO(:,1))/XSW_WGHT_NIR
! Adjust thickness to be as in snow computations:
DO JJ=1,INLVLS
DO JI=1,INI
ZSNOWDZ(JI,JJ) = PSNOWDZ(JI,JJ)/MAX(1.E-4,PPSN(JI))
ENDDO
ENDDO

CALL SNOWCRORADTRANS(XSNOWDZMIN, ZSPECTRALALBEDO, ZSNOWDZ, ZSNOWRHO, &
                ZPERMSNOWFRAC, PZENITH,  ZSNOWAGE,         &
                ZSNOWDIAMOPT, ZSNOWSPHERI, HSNOWMETAMO,       &
                INLVLS_USE, PTAU_N)

! Note that because we force a snow thickness to compute tramission, 
! a bogus value ( < 0) can be computed despite the non-existence of snow.
! To check/prevent any problems, make a simple check:
PTAU_N(:,:) = MAX(0., PTAU_N(:,:))
END IF
!
END SUBROUTINE SNOWCROALB_SPECTRAL_BANDS_MEB
!===============================================================================

SUBROUTINE SNOWCRORADTRANS(PSNOWDZMIN, PSPECTRALALBEDO, PSNOWDZ, PSNOWRHO,  &
                        PPERMSNOWFRAC, PZENITH,  PSNOWAGE, PSNOWDIAMOPT,    &
                        PSNOWSPHERI, HSNOWMETAMO, KNLVLS_USE, PRADTRANS)
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave (solar) radiation
!     through the snowpack (using a form of Beer's Law: exponential
!     decay of radiation with increasing snow depth).

USE CRO_CSTS, ONLY : EXUNDEF
USE MODD_SNOW_PAR, ONLY : XVSPEC1,XVSPEC2,XVSPEC3,XVBETA1,XVBETA2, &
                          XVBETA4,XVBETA3,XVBETA5, XMINCOSZEN
USE CRO_RADTRANS,  ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
!
IMPLICIT NONE
!*      0.1    declarations of arguments
REAL(8),                 INTENT(IN)    :: PSNOWDZMIN
REAL, DIMENSION(:),   INTENT(IN)    :: PPERMSNOWFRAC
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWDZ, PSNOWAGE, PSNOWDIAMOPT, PSNOWSPHERI
REAL, DIMENSION(:,:), INTENT(IN)    :: PSPECTRALALBEDO
CHARACTER(3),         INTENT(IN)    :: HSNOWMETAMO
INTEGER,DIMENSION(:), INTENT(IN)    :: KNLVLS_USE
REAL, DIMENSION(:,:), INTENT(OUT)   :: PRADTRANS
!*      0.2    declarations of local variables
INTEGER                              :: JJ, JI
INTEGER                              :: INI
INTEGER                              :: INLVLS
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZRADTOT, ZPROJLAT, ZCOSZEN
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZOPTICALPATH1, ZOPTICALPATH2, ZOPTICALPATH3
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDSGRAIN, ZCOEF, ZSNOWDZ, ZAGE
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZBETA1, ZBETA2, ZBETA3, ZWORK
! 0. Initialization:
INI    = SIZE(PSNOWDZ(:,:),1)
INLVLS = SIZE(PSNOWDZ(:,:),2)
! 1. Vanishingly thin snowpack check:
!    For vanishingly thin snowpacks, much of the radiation can pass through snowpack into underlying soil, making
!    a large (albeit temporary) thermal gradient: by imposing a minimum thickness, this increases the radiation absorbtion
!    for vanishingly thin snowpacks.
ZSNOWDZ(:,:) = MAX(PSNOWDZMIN, PSNOWDZ(:,:))
! 2. Extinction of net shortwave radiation
!   Fn of snow depth and density (Loth and Graf 1993: SNOWCVEXT => from Bohren and Barkstrom 1974
!   SNOWAGRAIN and SNOWBGRAIN=> from Jordan 1976)
!   Coefficient for taking into account the increase of path length of rays in snow due to zenithal angle
ZCOSZEN(:)=MAX(XMINCOSZEN,COS(PZENITH(:)))
! Compensate partly the fact that albedo formulation does not account for zenithal angle (polar or glacier regions)
ZPROJLAT(:)=(1.0-PPERMSNOWFRAC(:))+PPERMSNOWFRAC(:)/ZCOSZEN(:)
ZDSGRAIN = PSNOWDIAMOPT
! Extinction coefficient from Brun et al. (1989):
ZWORK(:,:)=SQRT(ZDSGRAIN(:,:))
ZBETA1(:,:)=MAX(XVBETA1*PSNOWRHO(:,:)/ZWORK(:,:),XVBETA2)
ZBETA2(:,:)=MAX(XVBETA3*PSNOWRHO(:,:)/ZWORK(:,:),XVBETA4)
ZBETA3(:,:)=XVBETA5
ZOPTICALPATH1(:) = 0.0
ZOPTICALPATH2(:) = 0.0
ZOPTICALPATH3(:) = 0.0
DO JJ=1,INLVLS
DO JI=1,INI
    IF (JJ<=KNLVLS_USE(JI)) THEN
        ZOPTICALPATH1(JI) = ZOPTICALPATH1(JI) + ZBETA1(JI,JJ)*ZSNOWDZ(JI,JJ)
        ZOPTICALPATH2(JI) = ZOPTICALPATH2(JI) + ZBETA2(JI,JJ)*ZSNOWDZ(JI,JJ)
        ZCOEF (JI,JJ) = XSW_WGHT_VIS*(1.0-PSPECTRALALBEDO(JI,1))*EXP(-ZOPTICALPATH1(JI)*ZPROJLAT(JI)) &
            + XSW_WGHT_NIR*(1.0-PSPECTRALALBEDO(JI,2))*EXP(-ZOPTICALPATH2(JI)*ZPROJLAT(JI)) 
    ELSE
        ZCOEF(JI,JJ)=0.
    ENDIF
ENDDO
ENDDO
! 3. Radiation trans at base of each layer
! NOTE, at level=0, rad = Swnet*(1-alb)  so radcoef(0)=1 implicitly
PRADTRANS(:,:)  = ZCOEF(:,:)
END SUBROUTINE SNOWCRORADTRANS
!===============================================================================

SUBROUTINE SNOWCROALB(OGLACIER,                            &
    PALBEDOSC,PSPECTRALALBEDO,PSNOWDZ,          &
    PSNOWRHO,PPERMSNOWFRAC,                     &
    PSNOWDIAMOPT_TOP,PSNOWSPHERI_TOP,PSNOWAGE_TOP, &
    PSNOWDIAMOPT_BOT,PSNOWSPHERI_BOT,PSNOWAGE_BOT, &
    PPS, PZENITH, KNLVLS_USE ,HSNOWMETAMO       ) 
    !!   DESC: see snowcro, copied from MEB & slightly modified by GM to remove unused modules
USE MODD_SNOW_PAR, ONLY : XANSMAX, XANSMIN,XAGLAMIN, XAGLAMAX, &
        XVRPRE1,XVRPRE2,XVAGING_NOGLACIER,   &
        XVAGING_GLACIER, XVSPEC1,XVSPEC2,    &
        XVSPEC3, XVW1,XVW2,XVD1,XVD2
!
IMPLICIT NONE
!*      0.1    declarations of arguments
LOGICAL, INTENT(IN)               :: OGLACIER    ! True = Over permanent snow and ice, initialise WGI=WSAT,Hsnow>=10m and allow 0.8<SNOALB<0.85 False = No specific treatment
REAL, DIMENSION(:), INTENT(IN)    :: PSNOWDZ,PPERMSNOWFRAC
REAL,DIMENSION(:,:), INTENT(IN)   :: PSNOWRHO ! For now only the 2 first layers are required
REAL, DIMENSION(:), INTENT(INOUT) :: PALBEDOSC
REAL, DIMENSION(:,:), INTENT(OUT) :: PSPECTRALALBEDO   ! Albedo in the different spectral bands
REAL, DIMENSION(:), INTENT(IN)    :: PSNOWDIAMOPT_TOP,PSNOWSPHERI_TOP,PSNOWAGE_TOP, &
                PSNOWDIAMOPT_BOT,PSNOWSPHERI_BOT,PSNOWAGE_BOT, PPS 
INTEGER, DIMENSION(:), INTENT(IN) :: KNLVLS_USE                
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH ! solar zenith angle for future use
CHARACTER(3),INTENT(IN)           :: HSNOWMETAMO ! metamorphism scheme
!*      0.2    declarations of local variables
REAL, DIMENSION(3,SIZE(PSNOWRHO,1)) :: ZALB_TOP, ZALB_BOT
REAL, DIMENSION(SIZE(PSNOWRHO,1))   :: ZANSMIN, ZANSMAX, ZMIN, ZMAX
REAL, DIMENSION(SIZE(PSNOWRHO,1))   :: ZFAC_TOP, ZFAC_BOT
REAL, DIMENSION(SIZE(PALBEDOSC))  :: ZVAGE1   
INTEGER         :: JJ   ! looping indexes
!-------------------------------------------------------------------------------
! Initialize:
IF ( OGLACIER ) THEN
    ZANSMIN(:) = XAGLAMIN * PPERMSNOWFRAC(:) + XANSMIN * (1.0-PPERMSNOWFRAC(:))
    ZANSMAX(:) = XAGLAMAX * PPERMSNOWFRAC(:) + XANSMAX * (1.0-PPERMSNOWFRAC(:))  
    ZVAGE1(:)  = XVAGING_GLACIER * PPERMSNOWFRAC(:) + XVAGING_NOGLACIER * (1.0-PPERMSNOWFRAC(:)) 
ELSE
    ZANSMIN(:) = XANSMIN
    ZANSMAX(:) = XANSMAX   
    ZVAGE1(:)  = XVAGING_NOGLACIER
ENDIF

! coherence control
! to remove when initialization routines will be updated
IF ( MINVAL(PSNOWAGE_BOT)<0. ) THEN
    print*, ('FATAL ERROR in SNOWCROALB during prep for SNOWCRO: Snow layer age inconsistent')
END IF
!
DO JJ=1, SIZE(PALBEDOSC)
    IF ( KNLVLS_USE(JJ)==0 ) THEN
        ! case with no snow on the ground 
        PALBEDOSC(JJ) = ZANSMIN(JJ)
    ELSE
    !
        CALL GET_ALB(JJ,PSNOWRHO(JJ,1),PPS(JJ),ZVAGE1(JJ),PSNOWDIAMOPT_TOP(JJ),&
        PSNOWSPHERI_TOP(JJ),PSNOWAGE_TOP(JJ),ZALB_TOP(:,JJ),HSNOWMETAMO)
        !                  
        IF ( KNLVLS_USE(JJ)>=2 ) THEN !modif ML
            ! second surface layer when it exists 
            CALL GET_ALB(JJ,PSNOWRHO(JJ,2),PPS(JJ),ZVAGE1(JJ),PSNOWDIAMOPT_BOT(JJ),&
            PSNOWSPHERI_BOT(JJ),MIN(365.,PSNOWAGE_BOT(JJ)),ZALB_BOT(:,JJ),HSNOWMETAMO)
        ELSE
            ! when it does not exist, the second surface layer gets top layer albedo   
            ZALB_BOT(:,JJ) = ZALB_TOP(:,JJ)
        ENDIF
        ! 
        ! computation of spectral albedo over 3 bands taking into account the respective
        ! depths of top layers 
        ZMIN(JJ) = MIN( 1., PSNOWDZ(JJ)/XVD1 )
        ZMAX(JJ) = MAX( 0., (PSNOWDZ(JJ)-XVD1)/XVD2 )
        ZFAC_TOP(JJ) = XVW1 * ZMIN(JJ) + XVW2 * MIN( 1., ZMAX(JJ) )
        ZFAC_BOT(JJ) = XVW1 * ( 1. - ZMIN(JJ) ) + XVW2 * ( 1. - MIN( 1., ZMAX(JJ) ) )
        PSPECTRALALBEDO(JJ,1) = ZFAC_TOP(JJ) * ZALB_TOP(1,JJ) + ZFAC_BOT(JJ) * ZALB_BOT(1,JJ) 
        PSPECTRALALBEDO(JJ,2) = ZFAC_TOP(JJ) * ZALB_TOP(2,JJ) + ZFAC_BOT(JJ) * ZALB_BOT(2,JJ)
        PSPECTRALALBEDO(JJ,3) = ZFAC_TOP(JJ) * ZALB_TOP(3,JJ) + ZFAC_BOT(JJ) * ZALB_BOT(3,JJ)
        !
        ! arbitrarily specified spectral distribution  
        ! to be changed when solar radiation distribution is an input variable 
        PALBEDOSC(JJ) = XVSPEC1 * PSPECTRALALBEDO(JJ,1) + &
        XVSPEC2 * PSPECTRALALBEDO(JJ,2) + &
        XVSPEC3 * PSPECTRALALBEDO(JJ,3) 
    !    
    ENDIF ! end case with snow on the ground
ENDDO ! end loop grid points
END SUBROUTINE SNOWCROALB
!===============================================================================

SUBROUTINE GET_ALB(KJ,PSNOWRHO_IN,PPS_IN,PVAGE1,PSNOWDIAMOPT,PSNOWSPHERI,PSNOWAGE,PALB,&
    HSNOWMETAMO)

USE MODD_SNOW_PAR, ONLY : XALBICE1, XALBICE2, XALBICE3,   &
           XRHOTHRESHOLD_ICE,              &
           XVALB2, XVALB3, XVALB4, XVALB5, &
           XVALB6, XVALB7, XVALB8, XVALB9, &
           XVALB10, XVALB11, XVDIOP1,      &
           XVRPRE1, XVRPRE2, XVPRES1

IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KJ
REAL, INTENT(IN) :: PSNOWRHO_IN, PPS_IN
REAL, INTENT(IN) :: PVAGE1
REAL, INTENT(IN) :: PSNOWDIAMOPT, PSNOWSPHERI, PSNOWAGE
REAL, DIMENSION(3), INTENT(OUT) :: PALB
CHARACTER(3),INTENT(IN)::HSNOWMETAMO
REAL :: ZDIAM, ZDIAM_SQRT
!
IF ( PSNOWRHO_IN<XRHOTHRESHOLD_ICE ) THEN
    ! Normal case (snow)
    ZDIAM = PSNOWDIAMOPT
    ZDIAM_SQRT = SQRT(ZDIAM)
    PALB(1) = MIN( XVALB2 - XVALB3*ZDIAM_SQRT, XVALB4 )
    PALB(2) = MAX( 0.3, XVALB5 - XVALB6*ZDIAM_SQRT )
    ZDIAM   = MIN( ZDIAM, XVDIOP1 )
    ZDIAM_SQRT = SQRT(ZDIAM)  
    PALB(3) = MAX( 0., XVALB7*ZDIAM - XVALB8*ZDIAM_SQRT + XVALB9 ) 
    PALB(1) = MAX( XVALB11, PALB(1) - MIN( MAX(PPS_IN/XVPRES1,XVRPRE1), XVRPRE2 ) * &
        XVALB10 * PSNOWAGE / PVAGE1 )
ELSE
    ! Prescribed spectral albedo for surface ice
    PALB(1) = XALBICE1
    PALB(2) = XALBICE2
    PALB(3) = XALBICE3
ENDIF
END SUBROUTINE GET_ALB
!---------------------------------------------------------------------------------------

END MODULE CRO_MODE_RAD