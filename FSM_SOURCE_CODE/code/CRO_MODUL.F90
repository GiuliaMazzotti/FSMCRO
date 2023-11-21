! Additional modules specific to the FSMCRO coupling

!-----------------------------------------------------------------------
! Physical constants
! some constants used in crocus that are initialized outside of ini_csts
! add 'E' as prefix to mark this is the externaliyed version. there might
! be some redundancy, I didnt check exhaustively
!-----------------------------------------------------------------------
MODULE CRO_CSTS 
    real(8), parameter :: &
        EXDAY        = 86400.,                  & ! Day duration 
        EXPI         = 2.*ASIN(1.),             & ! Pi
        EXAVOGADRO   = 6.0221367E+23,           & ! Avogadro number
        EXBOLTZ      = 1.380658E-23,            & ! Boltzman constant 
        EXG          = 9.80665,                 & ! Gravity constant
        EXMD         = 28.9644E-3,              & ! Molar mass of dry air 
        EXMV         = 18.0153E-3,              & ! Molar mass of vapor
        EXP00        = 1.E5,                    & ! Reference pressure
        EXRD         = EXAVOGADRO*EXBOLTZ/EXMD, & ! Gas constant for dry air
        EXRV         = EXAVOGADRO*EXBOLTZ/EXMV, & ! Gas constant for vapor
        EXCPD        = 7.*EXRD/2.,              & ! Cpd (dry air)         
        EXSNOWDMIN   = 0.000001,                & ! modd_snow_par: ISBA-ES Minimum total snow depth for thermal calculations. Used to prevent numerical problems as snow becomes vanishingly thin. 
        EXRHOSMAX_ES = 750.,                    & ! maximum values of the density of snow for ISBA-ES snow option(kg m-3)
        EXEMISSN     = 0.99,                    & ! Snow emissivity as set in ini_surf_csts.F90
        EXUNDEF      = 1.E+20                   ! Value given to undefined variables
    END MODULE CRO_CSTS

!-----------------------------------------------------------------------
! Model configuration
!-----------------------------------------------------------------------
MODULE CRO_SETTINGS
    CHARACTER(3) ::     &
        CSNOWCOMP,          & ! Compaction option: 'B92', 'S14' Schleef 2014, 'T11' Teufelsbauer 2011
        CSNOWCOND,          & ! Thermal conductivity scheme; HSNOWCOND=Y81 default Crocus from Yen et al. 1981; I02 ISBA_ES snow conductivity param (Boone et al. 2002); C11 Calonne et al. 2011 
        CSNOWFALL,          & ! Falling snow scheme; HSNOWFALL=V12 Vionnet et al. 2012 from Brun et al. 1989; A76 Anderson et al. 1976; S02 Lehning el al. 2002; P75 Pahaut 1975; NZE Constant 200 kg/m3 
        CSNOWHOLD,          & ! Liquid water content scheme; B92 default Crocus from Brun et al. 1992 or Vionnet et al. 2012; B02 ISBA_ES  parametrization (Boone et al. 2002); O04 CLM parametrization (Oleson et al 2004); 02 SNOWPACK aprametrization (Lehning et al 2002)
        CSNOWMETAMO,        & ! Metamorphism scheme; HSNOWMETAMO=C13 Carmagnola et al 2014; T07 Taillandier et al 2007; F06 Flanner et al 2006
        CSNOWRAD,           & ! Radiative transfer scheme; HSNOWRAD=B92 Brun et al 1992; =T17 (Tuzet et al. 2017) (Libois et al. 2013) TARTES with impurities content scheme
        CSNOWRES,           & ! Turbulent exchange option from ISBA-SNOW3L
        CSNOWZREF             ! Reference height; HSNOWZREF='CST' constant reference height from the snow surface; 'VAR' variable reference height from the snow surface (i.e. constant from the ground)                                                                
    CHARACTER(4) ::     &
        CSNOWDRIFT,         & ! Snowdrift scheme (mechanical transformation of snow grains and compaction of falling snow properties)       
        HIMPLICIT_WIND        ! Wind implicitation option
    LOGICAL ::          &
        LATMORAD,           & ! Activate atmotartes scheme
        LGLACIER,           & ! Simulation over glacier/permanent snow or not   
        LSELF_PROD,         & ! GM??? 
        LSNOW_ABS_ZENITH,   & ! Activate absorption of solar radiation for polar regions (not physical but better performance)
        LSNOWCOMPACT_BOOL,  & ! GM???
        LSNOWDRIFT_SUBLIM,  & ! Activate sublimation during drift
        LSNOWMAK_BOOL,      & ! Activate snowmaking option			
        LSNOWMAK_PROP,      & ! Specified properties of artificial snow
        LSNOWTILLER,        & ! Activate tiller
        OMEB                  ! Canopy activated or not (affects use of subroutines)
END MODULE CRO_SETTINGS

!-----------------------------------------------------------------------
! Meteorological driving variables
!-----------------------------------------------------------------------
MODULE CRO_DRIVE
    REAL(8) ::   PTSTEP     ! time step of the integration (REAL in snowcro)
    REAL(8), allocatable :: & 
        ZP_ZREF(:),    & ! Reference height of the first atmospheric level
        ZP_UREF(:),    & ! Wind reference height
        ZP_PS(:),      & ! surface pressure
        ZP_SRSNOW(:),  & ! rain rate [kg/(m2 s)]
        ZP_RRSNOW(:),  & ! snow rate (SWE) [kg/(m2 s)]
        ZP_SW_RAD(:),  & ! incoming solar radiation (W/m2)
        ZP_LW_RAD(:),  & ! atmospheric infrared radiation (W/m2)
        ZP_TA(:),      & ! atmospheric temperature at level za (K)
        ZP_QA(:),      & ! atmospheric specific humidity at level za
        ZP_VMOD(:)       ! modulus of the wind parallel to the orography (m/s)
END MODULE CRO_DRIVE

!-----------------------------------------------------------------------
! Grid parameters
!-----------------------------------------------------------------------
module CRO_GRID
    INTEGER ::      & 
        KSIZE1,         & ! number of modelled points
        KSIZE2,         & ! maximum number of snow layers
        KSIZE3,         & ! number of ground layers
        KSIZE4,         & ! number of blowing snow variables
        NIMPUR            ! number of impurities types
end module CRO_GRID

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
MODULE CRO_PARAMS
    REAL(8), allocatable :: &
        ZP_ALB(:),          & ! soil/vegetation albedo    
        ZP_Z0NAT(:),        & ! PZ0 = grid box average roughness length
        ZP_Z0EFF(:),        & ! PZ0EFF = roughness length for momentum                            
        ZP_Z0HNAT(:)! PZ0H = grid box average roughness length for heat   
    REAL(8), allocatable :: &
        ZP_PEW_A_COEF(:),   & ! A wind coefficient (m2s/kg)
        ZP_PET_A_COEF(:),   & ! A air temperature coefficient
        ZP_PEQ_A_COEF(:)      ! A air specific humidity coefficient  
END MODULE CRO_PARAMS

!-----------------------------------------------------------------------
! Spatial surface characteristics 
!-----------------------------------------------------------------------
MODULE CRO_SURF 
! Stuff with intent in
    REAL(8), allocatable :: &
        ZP_D_G(:,:),    & ! Assumed first soil layer thickness (m) Used to calculate ground/snow heat flux
        ZP_LAT(:),      & ! LAT/LON after packing
        ZP_LON(:),      & !  
        ZP_DIRCOSZW(:)    ! Cosinus of the angle between the normal to the surface and the vertical     
END MODULE CRO_SURF

!-----------------------------------------------------------------------
! Model state variables used in Crocus
!-----------------------------------------------------------------------
MODULE CRO_STATE_VARS
    REAL(8), allocatable ::  &
        ZP_SNOWSWE(:,:),        & ! Snow layer(s) Water Equivalent (SWE:kg m-2)  !! note: this variable is used to track if theres snow or not
        ZP_SNOWRHO(:,:),        & ! Snow layer(s) averaged density (kg/m3)
        ZP_SNOWHEAT(:,:),       & ! Snow layer(s) heat content (J/m2)
        ZP_SNOWDIAMOPT(:,:),      & ! PSNOWDIAMOPT = Snow layers grain feature 1
        ZP_SNOWSPHERI(:,:),      & ! PSNOWSPHERI = Snow layer grain feature 2
        ZP_SNOWHIST(:,:),       & ! Snow layer grain historical parameter (only for non dendritic snow)
        ZP_SNOWAGE(:,:),        & ! Snow grain age 
        ZP_SNOWIMPUR(:,:,:),    & ! Snow impurity content (g) (LOCATION,LAYER,NIMPUR)) Impur type :1/BC 2/Dust
        ZP_SNOWDZ (:,:),        & ! Snow layer(s) thickness (m)
        ZP_SNOWTEMP(:,:),       & ! Snow layer(s) temperature (m)
        ZP_SNOWLIQ(:,:),        & ! Snow layer(s) liquid water content (m)
        ZP_SNOWALB(:),          & ! Prognostic surface snow albedo (does not include anything but the actual snow cover) dim: KSIZE1
        ZP_TG(:,:)                ! Surface soil temperature (effective temperature the of layer lying below snow); dims: KSIZE1,KSIZE3
    LOGICAL, allocatable ::  &
        CROSNOW_ON(:)             ! Flag use of Crocus snow for individual point (either start of a snowpack or existing snowpack when CRO_ON, equivalent to mask used in snow3L_isba)
END MODULE CRO_STATE_VARS

MODULE CRO_DUMMIES ! SNOWCRO arguments that can be initialized with a dummy value, no proper coupling
    REAL(8), allocatable ::    &
        ZP_ANGL_ILLUM(:),       & ! Angle entre le soleil et la normal au sol et le soleil (=zenith sans pente au sol) utilisé dans TARTES
        ZP_BLOWSNW(:,:),        & ! Properties of deposited blowing snow (from Sytron or Meso-NH/Crocus); 1 : Deposition flux (kg/m2/s); 2 : Density of deposited snow (kg/m3); 3 : SGRA1 dims: KSIZE1,KSIZE4                        
        ZP_CDSNOW(:),           & ! drag coefficient for momentum over snow
        ZP_CHSNOW(:),           & ! drag coefficient for heat over snow
        ZP_DIFF_RATIO(:,:),     & ! F.T diffuse to total irradiance ratio
        ZP_DIR_SW(:,:),         & ! F.T  direct spectral irradiance (W/m2/um) / used when Tartes is activated / (npoints,nbands)
        ZP_EMISNOW(:),          & ! snow surface emissivity
        ZP_EVAP(:),             & ! total evaporative flux (kg/m2/s)
        ZP_EVAPCOR(:),          & ! evaporation/sublimation correction term: extract any evaporation exceeding the actual snow cover (as snow vanishes) and apply it as a surface soil water sink. [kg/(m2 s)]
        ZP_GFLUXSNOW(:),        & ! net heat flux from snow (W/m2)
        ZP_GFLXCOR(:),          & ! flux correction to underlying soil for vanishing snowpack (to put any energy excess from snow to soil) (W/m2)
        ZP_GRNDFLUX(:),         & ! soil/snow interface heat flux (W/m2)
        ZP_GSFCSNOW(:),         & ! heat flux between the surface and sub-surface snow layers (2nd layer) (W/m2)) - used within crocus, but not output further (SNOWFLUX in org code)
        ZP_HPSNOW(:),           & ! heat release from rainfall (W/m2)
        ZP_IMPDRY(:,:),         & ! Dry deposit coefficient from Forcing File(g/m²/s)
        ZP_IMPWET(:,:),         & ! Wet deposit coefficient from Forcing File(g/m²/s)
        ZP_PSN3L(:),            & ! 'snow fraction', treated as meteo input for debug purposes in crocus. Different use if MEB "on". In this case, it is only used for Tg update since only this variable has a sub-grid relevance.
        ZP_QS(:),               & ! surface humidity
        ZP_RI(:),               & ! Richardson number
        ZP_SCA_SW(:,:),         & ! F.T diffuse spectral irradiance (W/m2/um) / used when Tartes is activated
        ZP_SNDRIFT(:),          & ! blowing snow sublimation (kg/m2/s)    
        ZP_SNOWHMASS(:),        & ! heat content change due to mass changes in snowpack (J/m2): for budget calculations only.
        ZP_SNOWMAK(:),          & ! Snowmaking thickness (m)
        ZP_SPEC_ALB(:,:),       & ! F.T spectral albedo
        ZP_SPEC_TOT(:,:),       & ! total incoming spectral irradiance (npoints,nbands)
        ZP_THRUFAL(:),          & ! rate that liquid water leaves snow pack: paritioned into soil infiltration/runoff by ISBA [kg/(m2 s)] 
        ZP_VEGTYPE(:),          & ! PPERMSNOWFRAC  = fraction of permanet snow/ice ???GM why such a change of name???
        ZP_USTARSNOW(:)           ! friction velocity over snow (m/s) 
END MODULE CRO_DUMMIES 

MODULE CRO_DIAG ! output of snowcro_diag
    REAL(8), allocatable ::    & 
        ZP_SNOWDEND(:,:),       & ! dendricity (-)
        ZP_SNOWSPHER(:,:),       & ! sphericity (-)
        ZP_SNOWSIZE(:,:),       & ! grain size (m)
        ZP_SNOWSSA(:,:),       & ! specific surface area (m2/kg)
        ZP_SNOWTYPEMEPRA(:,:),       & ! snow type (-)
        ZP_SNOWRAM(:,:),       &  ! ram penetration strength (kgf = 9.81 N)
        ZP_SNOWSHEAR(:,:),       &  ! shear strength (kgf/dm2 = 0.981 kPa)
        ZP_ACC_RAT(:,:),       &  ! accidental ratio shear strength/stress
        ZP_NAT_RAT(:,:),       & ! natural ratio shear strength/stress
        ZP_SNDPT_12H(:),       & ! slope-perpendicular thickness of snow with age <= 12h  (m) (projection in reproj_diag_isban)
        ZP_SNDPT_1DY(:),       &  ! slope-perpendicular thickness of snow with age <= 1 day  (m) (projection in reproj_diag_isban)
        ZP_SNDPT_3DY(:),       & ! slope-perpendicular thickness of snow with age <= 3 days (m) (projection in reproj_diag_isban)
        ZP_SNDPT_5DY(:),       &  ! slope-perpendicular thickness of snow with age <= 5 days (m) (projection in reproj_diag_isban)
        ZP_SNDPT_7DY(:),       &  ! slope-perpendicular thickness of snow with age <= 7 days (m) (projection in reproj_diag_isban)
        ZP_SNSWE_1DY(:),       &  ! swe with age <= 1 day  (kg m-2)
        ZP_SNSWE_3DY(:),       &   ! swe with age <= 3 days (kg m-2)
        ZP_SNSWE_5DY(:),       & ! swe with age <= 5 days (kg m-2)
        ZP_SNSWE_7DY(:),       &  ! swe with age <= 7 days (kg m-2)
        ZP_SNRAM_SONDE(:),       &  ! ramsonde top penetration slope-perpendicular thickness (m) (projection in reproj_diag_isban)
        ZP_SN_WETTHCKN(:),       & ! top continous wet snow slope-perpendicular thickness (m) (projection in reproj_diag_isban)
        ZP_SN_REFRZNTHCKN(:),       & ! top continous refrozen snow slope-perpendicular thickness (m) (projection in reproj_diag_isban)
        ZP_DEP_HIG(:),       &  ! vertical depth of high instability (m)
        ZP_DEP_MOD(:),       &  ! vertical depth of moderate instability (m)
        ZP_DEP_SUP(:),       & ! vertical depth of superior profile bottom (m)
        ZP_DEP_TOT(:),       &     ! total vertical snow depth (m)
        ZP_DEP_HUM(:),       &      ! vertical thickness of the uppest continuous block of humid snow in the sup profile
        ZP_ACC_LEV(:),       &             ! accidental risk level (0-4)
        ZP_NAT_LEV(:),       &  ! natural risk level (0-6)
        ZP_PRO_SUP_TYP(:),       &  ! type of superior profile (0, 4, 5, 6)
        ZP_PRO_INF_TYP(:),       & ! type of inferior profile (0, 1, 6)
        ZP_AVA_TYP(:),       &! type of avalanche (0-6)   
        ZP_SNOWIMP_CONC(:,:,:)   ! concentration of snow impurities (? to be checked...)
END MODULE CRO_DIAG

MODULE CRO_NLSTPARAMS ! Input parameters, i.e. parameters that are editable by user through the namelist but specific to the coupling
    REAL(8)  :: &
        X_CVHEATF         ! factor to modulate the soil heat capacity, aimed at mimicking the ESCRO option
END MODULE CRO_NLSTPARAMS













 











