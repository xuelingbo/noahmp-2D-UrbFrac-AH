module ResistanceCanopyStomataMedlynMod

!!! Compute canopy stomatal resistance and foliage photosynthesis based on Medlyn et al. (2011) scheme
!!! gs = g0 + 1.6 * (1 + g1/sqrt(Ds)) * An / cs, with Ds the vapor pressure deficit at leaf surface [kPa]
!!! Photosynthesis is the same as the Ball-Berry scheme; only the stomatal conductance closure differs

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine ResistanceCanopyStomataMedlyn(noahmp, IndexShade)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: STOMATA
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Medlyn closure: adapted from Noah-MP-ML2 LeafPhotosynthesisMod (after CLM-ml v0, G. Bonan)
! -------------------------------------------------------------------------

    implicit none

    integer          , intent(in   ) :: IndexShade            ! index for sunlit/shaded (0=sunlit;1=shaded)
    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: IndIter               ! iteration index
    integer, parameter               :: NumIter = 3           ! number of iterations
    real(kind=kind_noahmp)           :: RadPhotoActAbsTmp     ! temporary absorbed par for leaves [W/m2]
    real(kind=kind_noahmp)           :: ResistanceStomataTmp  ! temporary leaf stomatal resistance [s/m]
    real(kind=kind_noahmp)           :: PhotosynLeafTmp       ! temporary leaf photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp)           :: NitrogenFoliageFac    ! foliage nitrogen adjustment factor (0 to 1)
    real(kind=kind_noahmp)           :: CarboxylRateMax       ! maximum rate of carbonylation [umol co2/m2/s]
    real(kind=kind_noahmp)           :: MPE                   ! prevents overflow for division by zero
    real(kind=kind_noahmp)           :: RLB                   ! boundary layer resistance [s m2 / umol]
    real(kind=kind_noahmp)           :: TC                    ! foliage temperature [deg C]
    real(kind=kind_noahmp)           :: CS                    ! co2 concentration at leaf surface [Pa]
    real(kind=kind_noahmp)           :: KC                    ! co2 Michaelis-Menten constant [Pa]
    real(kind=kind_noahmp)           :: KO                    ! o2 Michaelis-Menten constant [Pa]
    real(kind=kind_noahmp)           :: B,C                   ! coefficients of quadratic equation for GS
    real(kind=kind_noahmp)           :: PPF                   ! absorb photosynthetic photon flux [umol photons/m2/s]
    real(kind=kind_noahmp)           :: WC                    ! Rubisco limited photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp)           :: WJ                    ! light limited photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp)           :: WE                    ! export limited photosynthesis [umol co2/m2/s]
    real(kind=kind_noahmp)           :: CP                    ! co2 compensation point [Pa]
    real(kind=kind_noahmp)           :: CI                    ! internal co2 [Pa]
    real(kind=kind_noahmp)           :: AWC                   ! intermediate calculation for wc
    real(kind=kind_noahmp)           :: J                     ! electron transport [umol co2/m2/s]
    real(kind=kind_noahmp)           :: CF                    ! [s m2/umol] -> [s/m]
    real(kind=kind_noahmp)           :: GBW                   ! leaf boundary layer conductance [umol/m2/s]
    real(kind=kind_noahmp)           :: GS                    ! leaf stomatal conductance [umol/m2/s]
    real(kind=kind_noahmp)           :: VPD                   ! leaf-to-canopy-air vapor pressure deficit [kPa]
    real(kind=kind_noahmp)           :: T1                    ! 1.6*An/cs [umol/m2/s]
    real(kind=kind_noahmp), parameter :: VPDMin = 100.0       ! minimum vapor pressure deficit [Pa]
    real(kind=kind_noahmp)           :: T                     ! temporary var
! local statement functions
    real(kind=kind_noahmp)           :: F1                    ! generic temperature response (statement function)
    real(kind=kind_noahmp)           :: F2                    ! generic temperature inhibition (statement function)
    real(kind=kind_noahmp)           :: AB                    ! used in statement functions
    real(kind=kind_noahmp)           :: BC                    ! used in statement functions
    F1(AB, BC) = AB**( (BC - 25.0) / 10.0 )
    F2(AB)     = 1.0 + exp( (-2.2e05 + 710.0 * (AB + 273.16)) / (8.314 * (AB + 273.16)) )

! --------------------------------------------------------------------
    associate(                                                                        &
              PressureAirRefHeight    => noahmp%forcing%PressureAirRefHeight         ,& ! in,  air pressure [Pa] at reference height
              TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight      ,& ! in,  air temperature [K] at reference height
              SoilTranspFacAcc        => noahmp%water%state%SoilTranspFacAcc         ,& ! in,  accumulated soil water transpiration factor (0 to 1)
              IndexGrowSeason         => noahmp%biochem%state%IndexGrowSeason        ,& ! in,  growing season index (0=off, 1=on)
              NitrogenConcFoliage     => noahmp%biochem%state%NitrogenConcFoliage    ,& ! in,  foliage nitrogen concentration [%]
              NitrogenConcFoliageMax  => noahmp%biochem%param%NitrogenConcFoliageMax ,& ! in,  foliage nitrogen concentration when f(n)=1 [%]
              QuantumEfficiency25C    => noahmp%biochem%param%QuantumEfficiency25C   ,& ! in,  quantum efficiency at 25c [umol co2 / umol photon]
              CarboxylRateMax25C      => noahmp%biochem%param%CarboxylRateMax25C     ,& ! in,  maximum rate of carboxylation at 25c [umol co2/m**2/s]
              CarboxylRateMaxQ10      => noahmp%biochem%param%CarboxylRateMaxQ10     ,& ! in,  change in maximum rate of carboxylation for each 10C temp change
              PhotosynPathC3          => noahmp%biochem%param%PhotosynPathC3         ,& ! in,  C3 photosynthetic pathway indicator: 0. = c4, 1. = c3
              SlopeConductMedlyn      => noahmp%biochem%param%SlopeConductMedlyn     ,& ! in,  Medlyn slope parameter g1 [kPa^0.5]
              ConductanceLeafMinMedlyn=> noahmp%biochem%param%ConductanceLeafMinMedlyn,& ! in, Medlyn minimum leaf conductance g0 [umol/m2/s]
              Co2MmConst25C           => noahmp%energy%param%Co2MmConst25C           ,& ! in,  co2 michaelis-menten constant at 25c [Pa]
              O2MmConst25C            => noahmp%energy%param%O2MmConst25C            ,& ! in,  o2 michaelis-menten constant at 25c [Pa]
              Co2MmConstQ10           => noahmp%energy%param%Co2MmConstQ10           ,& ! in,  q10 for Co2MmConst25C
              O2MmConstQ10            => noahmp%energy%param%O2MmConstQ10            ,& ! in,  q10 for ko25
              TemperatureCanopy       => noahmp%energy%state%TemperatureCanopy       ,& ! in,  vegetation temperature [K]
              VapPresSatCanopy        => noahmp%energy%state%VapPresSatCanopy        ,& ! in,  canopy saturation vapor pressure at TV [Pa]
              PressureVaporCanAir     => noahmp%energy%state%PressureVaporCanAir     ,& ! in,  canopy air vapor pressure [Pa]
              PressureAtmosO2         => noahmp%energy%state%PressureAtmosO2         ,& ! in,  atmospheric o2 pressure [Pa]
              PressureAtmosCO2        => noahmp%energy%state%PressureAtmosCO2        ,& ! in,  atmospheric co2 pressure [Pa]
              ResistanceLeafBoundary  => noahmp%energy%state%ResistanceLeafBoundary  ,& ! in,  leaf boundary layer resistance [s/m]
              VegFrac                 => noahmp%energy%state%VegFrac                 ,& ! in,  greeness vegetation fraction
              RadPhotoActAbsSunlit    => noahmp%energy%flux%RadPhotoActAbsSunlit     ,& ! in,  average absorbed par for sunlit leaves [W/m2]
              RadPhotoActAbsShade     => noahmp%energy%flux%RadPhotoActAbsShade      ,& ! in,  average absorbed par for shaded leaves [W/m2]
              ResistanceStomataSunlit => noahmp%energy%state%ResistanceStomataSunlit ,& ! out, sunlit leaf stomatal resistance [s/m]
              ResistanceStomataShade  => noahmp%energy%state%ResistanceStomataShade  ,& ! out, shaded leaf stomatal resistance [s/m]
              PhotosynLeafSunlit      => noahmp%biochem%flux%PhotosynLeafSunlit      ,& ! out, sunlit leaf photosynthesis [umol co2/m2/s]
              PhotosynLeafShade       => noahmp%biochem%flux%PhotosynLeafShade        & ! out, shaded leaf photosynthesis [umol co2/m2/s]
             )
! ----------------------------------------------------------------------

    ! initialization
    MPE = 1.0e-6

    ! initialize ResistanceStomata=maximum value and photosynthesis=0 because will only do calculations
    ! for RadPhotoActAbs  > 0, in which case ResistanceStomata <= maximum value and photosynthesis >= 0
    CF = PressureAirRefHeight / (8.314 * TemperatureAirRefHeight) * 1.0e06  ! unit conversion factor
    ResistanceStomataTmp = 1.0 / ConductanceLeafMinMedlyn * CF
    PhotosynLeafTmp      = 0.0
    if ( IndexShade == 0 ) RadPhotoActAbsTmp = RadPhotoActAbsSunlit / max(VegFrac,1.0e-6)  ! Sunlit case
    if ( IndexShade == 1 ) RadPhotoActAbsTmp = RadPhotoActAbsShade  / max(VegFrac,1.0e-6)  ! Shaded case

    ! only compute when there is radiation absorption
    if ( RadPhotoActAbsTmp > 0.0 ) then

       NitrogenFoliageFac = min(NitrogenConcFoliage/max(MPE, NitrogenConcFoliageMax), 1.0)
       TC                 = TemperatureCanopy - ConstFreezePoint
       PPF                = 4.6 * RadPhotoActAbsTmp
       J                  = PPF * QuantumEfficiency25C
       KC                 = Co2MmConst25C * F1(Co2MmConstQ10, TC)
       KO                 = O2MmConst25C * F1(O2MmConstQ10, TC)
       AWC                = KC * ( 1.0 + PressureAtmosO2 / KO )
       CP                 = 0.5 * KC / KO * PressureAtmosO2 * 0.21
       CarboxylRateMax    = CarboxylRateMax25C / F2(TC) * NitrogenFoliageFac * &
                            SoilTranspFacAcc * F1(CarboxylRateMaxQ10, TC)
       ! first guess ci
       CI  = 0.7 * PressureAtmosCO2 * PhotosynPathC3 + 0.4 * PressureAtmosCO2 * (1.0 - PhotosynPathC3)
       ! ResistanceLeafBoundary: s/m -> s m**2 / umol
       RLB = ResistanceLeafBoundary / CF
       GBW = 1.0 / RLB
       ! leaf-to-canopy-air vapor pressure deficit [Pa -> kPa]
       VPD = max(VapPresSatCanopy - min(PressureVaporCanAir,VapPresSatCanopy), VPDMin) * 0.001

       ! ci iteration
       do IndIter = 1, NumIter
          WJ = max(CI-CP, 0.0) * J / (CI + 2.0*CP) * PhotosynPathC3 + J * (1.0 - PhotosynPathC3)
          WC = max(CI-CP, 0.0) * CarboxylRateMax / (CI + AWC) * PhotosynPathC3 + &
               CarboxylRateMax * (1.0 - PhotosynPathC3)
          WE = 0.5 * CarboxylRateMax * PhotosynPathC3 + &
               4000.0 * CarboxylRateMax * CI / PressureAirRefHeight * (1.0 - PhotosynPathC3)
          PhotosynLeafTmp = min(WJ, WC, WE) * IndexGrowSeason
          CS = max(PressureAtmosCO2-1.37*RLB*PressureAirRefHeight*PhotosynLeafTmp, MPE)
          ! Medlyn closure with leaf surface VPD implicitly coupled through GBW:
          ! GS**2 + B*GS + C = 0, take the larger root
          if ( PhotosynLeafTmp > 0.0 ) then
             T1 = 1.6 * PhotosynLeafTmp * PressureAirRefHeight / CS
             B  = -( 2.0 * (ConductanceLeafMinMedlyn + T1) + (SlopeConductMedlyn * T1)**2 / (GBW * VPD) )
             C  = ConductanceLeafMinMedlyn * ConductanceLeafMinMedlyn + &
                  ( 2.0 * ConductanceLeafMinMedlyn + T1 * (1.0 - SlopeConductMedlyn * SlopeConductMedlyn / VPD) ) * T1
             GS = 0.5 * ( -B + sqrt(max(B*B - 4.0*C, 0.0)) )
             GS = max(GS, ConductanceLeafMinMedlyn)
          else
             GS = ConductanceLeafMinMedlyn
          endif
          ResistanceStomataTmp = 1.0 / GS
          CI   = max(CS-PhotosynLeafTmp*PressureAirRefHeight*1.65*ResistanceStomataTmp, 0.0)
       enddo

       ! ResistanceStomata:  s m**2 / umol -> s/m
       ResistanceStomataTmp = ResistanceStomataTmp * CF

    endif ! RadPhotoActAbsTmp > 0.0

    ! assign updated values
    ! Sunlit case
    if ( IndexShade == 0 ) then
       ResistanceStomataSunlit = ResistanceStomataTmp
       PhotosynLeafSunlit      = PhotosynLeafTmp
    endif
    ! Shaded case
    if ( IndexShade == 1 ) then
       ResistanceStomataShade  = ResistanceStomataTmp
       PhotosynLeafShade       = PhotosynLeafTmp
    endif

    end associate

  end subroutine ResistanceCanopyStomataMedlyn

end module ResistanceCanopyStomataMedlynMod
