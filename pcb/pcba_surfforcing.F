C $Header: /u/gcmpack/MITgcm/pkg/pcb/pcb01_surfforcing.F, Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

#ifdef USE_EXFIWR
#include "SEAICE_OPTIONS.h"
#include "EXF_OPTIONS.h"
#endif


CBOP
C     !ROUTINE: PCBa_SURFFORCING
C     !INTERFACE:
      SUBROUTINE PCBa_SURFFORCING( Ptr, SurfPCBa,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myIter, myTime, myThid )

      !DESCRIPTION:
      !==========================================================
      ! SUBROUTINE PCBa_SURFFORCING
      !   Calculate the PCBa air-sea fluxes 
      !
      ! MODIFICATION HISTORY
      !   XX Jul 2014 - hma - PFOS code inherited from Xianming Zhang and
      !                       adopted for PCBs.
      !   07 Aug 2014 - hma - add inputs from atmospheric deposition
      !   19 May 2015 - hma - Update USE_EXFIWR calls based on Hg
      !                       so you can use online ice, radiation, 
      !                       etc. when running with ECCOv4. 
      !   01 Jun 2015 - hma - Correct sign on pcbaevaso and use of
      !                       wind and ice variables.
      !   02 Jul 2015 - hma - Update air-sea exchange equations.
      !   23 Sept 2016 - ccw- Updated dynamic drag coefficient, 
      !                       unit errors in kw and Scw
      !   07 Oct 2016 - ccw - Salinity adjustment of Henry's law
      !                       constant
      !   07 Dec 2016 - ccw - Added SURFSALT, SURFTEMP, PCBaFUGA
      !                       diagnostics
      !==========================================================

      !USES:
      IMPLICIT NONE

      !----------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"

#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PCB_SIZE.h" 
#include "PCB.h"
#include "GCHEM.h"    

#ifdef USE_EXFIWR
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#include "EXF_FIELDS.h"
#endif


      !----------------------------------------------------------
      ! INPUT/OUTPUT PARAMETERS:
      !----------------------------------------------------------
      _RL  PTR (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr, 
     &             nSx, nSy, PTRACERS_num)       ! ocean PCBa concentration
      _RL  SurfPCBa (1-OLx:sNx+OLx,1-OLy:sNy+OLy)! ?
      _RL  FluxPCBa(1-OLx:sNx+OLx,1-OLy:sNy+OLy) ! air-sea PCBa fluxes
      _RL  pcbarivero                            ! input of PCBa from rivers
      _RL  pcbaevaso                             ! PCBa evasion
      _RL  pcbadepoo                             ! PCBa atmospheric deposition (mol/m2/sec)
      _RL  pcbaconco                             ! atmospheric concentration (mol/m3)
      _RL  myTime                                ! current time in simulation
      _RL  iceo                                  ! ice fraction
      _RL  windo                                 ! wind speed (m/s)          
      _RL  surfsalt                              ! surface water salinity (g/kg)
      _RL  surftemp                              ! surface water temperature (celsius)
      _RL  pcbafug                               ! fugacity between air and water

      
      INTEGER iMin  , iMax                       ! computation domain, 1st index bounds 
      INTEGER jMin  , jMax                       ! computation domain, 2nd  index bounds
      INTEGER bi    , bj                         ! current tile indices
      INTEGER myIter                             ! current iteration number
      INTEGER myThid                             ! my Thread Id number
      INTEGER i1,j1,month    
CEOP

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB

      !----------------------------------------------------------
      ! LOCAL VARIABLES:
      !----------------------------------------------------------

      ! i, j, k - Loop counters
      INTEGER i,j,k

      ! Solubility relation coefficients
      _RL etaA, rhoA
      _RL etaW
      _RL SchmidtAPCBa
      _RL SchmidtWPCBa
      _RL DiffusivityAir(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL DiffusivitySurfWater(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C hma, 02 jul 2015
C      _RL ka(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C      _RL kw(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C      _RL Kaexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C      _RL KawPCBa(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C      _RL KawPCBeff(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ka                                          ! single-phase air-side transfer velocity
      _RL kw                                          ! single-phase water-side transfer velocity
      _RL Kwexch                                      ! total water-side transfer velocity
      _RL F_evas, F_down, F_up                        ! Flux of PCB gas across the air-sea interface
C hma, 02 jul 2015 (these don't appear to be used anywhere)
C      _RL ttemp, ttemp2
      _RL temp1, temp2, temp3, temp4, temp5, temp6, temp7
      _RL doco                                        ! dissolved organic carbon concentration
      _RL C2K                                         ! conversion factor, Celcius to Kelvin 
      _RL R                                           ! universal gas constant 
      _RL TK0                                         ! temperature at standard conditions 
      _RL KH                                          ! dimensionlessnliquid-over-gas Henry's law constant 
      _RL KS                                          ! Setchenow constant (ccw)
C hma, 02 jul 2015
C      _RL KHcc                                        ! dimensionless liquid-over-gas Henry's law constant 
C      _RL KawD                                        ! dimensionless gas-over-liquid Henry's law constant 
      _RL eta_s                                       ! dynamic viscocity of water 
      _RL Ms                                          ! relative molecular mass of water
      _RL PHI                                         ! association factor of water
#ifdef ALLOW_DIAGNOSTICS      
      _RL pcbaconc_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)  ! atmospheric concentration (mol/m3)
      _RL pcbaevas_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)  ! evasion                   (mol/m2/sec)
      _RL pcbariver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy) ! river inputs              (mol/m2/sec)
      _RL pcbadep_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)   ! atmospheric deposition    (mol/m2/sec)
      _RL temp_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)      ! surface sewater temperature (C)
      _RL salt_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)      ! surface salinity (g/kg)
      _RL pcbafug_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)   ! fugacity
      _RL ASF_up_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)    ! upward flux (mol/m2/sec)
      _RL ASF_dn_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)    ! downward flux (mol/m2/sec)
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */

      !----------------------------------------------------------
      ! Set up constants
      !----------------------------------------------------------

      ! Conversion factor, Celcius to Kelvin
      C2K = 273.15d0     
 
      ! Universal gas constant (m^3 Pa mol^-1 K^-1 or J K^-1 mol^-1)
      R   = 8.314d0
 
      ! Temperature at standard conditions (Kelvin)
      TK0 = 298.15d0       

      ! Setchenow constant
      KS = 0.35d0
      !----------------------------------------------------------                         
      ! Initialize diagnostics
      !----------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            ! atmospheric concentration
            pcbaconc_diag(i,j)     = 0.d0
            
            ! evasion 
            pcbaevas_diag(i,j)     = 0.d0

            ! inputs from rivers
            pcbariver_diag(i,j)    = 0.d0

            ! inputs from atmospheric deposition
            pcbadep_diag(i,j)      = 0.d0

            ! surface seawater temperature
            temp_diag(i,j)         = 0.d0

            ! surface seawater salinity
            salt_diag(i,j)         = 0.d0

            ! fugacity of pcba
            pcbafug_diag(i,j)      = 0.d0

            ! upward flux of pcba (air-sea exchange)
            ASF_up_diag(i,j)       = 0.d0

            ! downward flux of pcba (air-sea exchange)
            ASF_dn_diag(i,j)       = 0.d0

        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */   

#ifdef USE_EXFIWR
       ! Get current month
        month=mod(int(myTime/2629800.0), 12) + 1
#endif /* USE_EXFIWR */   

      
        !--------------------------------------------------------
        ! Air-sea exchange, inputs from rivers, inputs from deposition
        !--------------------------------------------------------
        DO j=jMin,jMax
         DO i=iMin,iMax

#ifdef USE_EXFIWR
            ! use online ECCOv4 ice and wind information            
            iceo  = AREA(i,j,bi,bj)    ! ice fraction
            windo = wspeed(i,j,bi,bj)  ! wind speed (m/s)
#else /* USE_EXFIWR */
            ! read ice and wind data from offline file
            iceo  = fice(i,j,bi,bj)
            windo = wind(i,j,bi,bj)
#endif /* USE_EXFIWR */

            !---------------------------------------
            ! Are you in an ocean box?
            !   Ocean box when HFacC = 1
            !   Land  box when HFacC = 0
            !---------------------------------------
            IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN ! Yanxu's option

              !-------------------------------------
              ! Local variables
              !-------------------------------------

              ! Atmospheric PCB concentration (mol/m3)
              pcbaconco = max(pcbaconc(i,j,bi,bj),0.d0)

              ! river : kg/s    -> mol/m2/s
              pcbarivero = PCBaRiver(i,j,bi,bj)*recip_rA(i,j,bi,bj)
     &           *1000.d0/Mb

              ! Inputs from atmospheric deposition 
C Units on deposition file are mol/m2/sec (hma, 29 jun 2015)
C              ! kg/m2/s -> mol/m2/s
C              pcbadepoo = PCBaDEP(i,j,bi,bj) *1000.d0/Mb
              pcbadepoo = PCBaDEP(i,j,bi,bj)
#ifdef PCBAREFRACTORY
              ! Assume part of the deposition is refractory
     &                 * (1.0 _d 0 - f_depo_refrac) 
#endif /* PCBAREFRACTORY */

              ! DOC concentration (mmol C/m3)
              doco     = max(doc(i,j,1,bi,bj),0.d0)

      !=========================================================
      ! CALCULATE AIR-SEA EXCHANGE
      !
      ! References:
      ! Johnson, M. T. A numerical scheme to calculate temperature and 
      !   salinity dependent air-water transfer velocities for any gas.
      !   Ocean Science, 6, 913-932, 2010.
      !
      ! Sander, R: Compilation of Henry's law constant for inorganic 
      !   and organic species of potential importance in environmental
      !   chemistry, 1999.
      !   http://www.mpch-mainz.mpg.de/~sander/res/henry.html.
      !
      ! Tsilingiris, P. Thermophysical and transport properties of 
      !   humid air at temperature ranges between 0 and 100 C. Energ.
      !   Convers. Manage, 49, 1098-1110, 2008.
      !
      ! Tucker, W. A. and Nelken, L. H. Diffusion coefficients in air
      !   and water, in: Handbook of Chemical Property Estimatation
      !   Methods, American Chemical Society, Washington D. C., 1990.
      !=========================================================
                  
      !---------------------------------------------------------
      ! CALCULATE SCHMIDT NUMBER FOR PCBa
      !
      !   From Johnson (2010), Ocean Science:
      !
      !   Schmidt number in air:
      !   Sca = v_a / D_a = eta_a / ( rho_a * D_a )                [Eq. 23]
      !   
      !   dynamic viscosity:
      !   etaA = sv_0 + sv_1*t + sv_2*t^2 + sv_3*t^3 + sv_4*t^4    [Eq. 24]
      !
      !   density:
      !   rhoA = sd_0 + sd_1*t + sd_2*t^2 + sd_3*t^3               [Eq. 25]
      !
      !   diffusion coefficient in air:
      !   Da = 0.001 * (T^1.75) * (Mr^0.5) / ( (P*Va^1/3) + Vb^1/3 )^2   [Eq. 26]
      !
      !   Mr = ( Ma + Mb ) / (Ma * Mb)                             [Eq. 27]
      !
      !   where t is temperature in degrees Celcius, the units on etaA are
      !   kg m^-3, and the units on rhoA are kg m^-3. Values of sv_* and 
      !   sd_* are in Table A1 of Johnson (2010), which references 
      !   Tsilingiris (2008). P is the atmospheric pressure in atm. Va and
      !   Vb are the molar volumes of air and of the gas of interest, 
      !   respectively. Ma and Mb are the molecular masses of air and the
      !   gas of interest, respectively. T is temperautre in Kelvin.
      !---------------------------------------------------------

            k=1
            !---------------------------------------------------
            ! Dynamic viscosity of air (kg m^-3)
            !
            ! Eq. 24 from Johnson (2010)
            !---------------------------------------------------
            etaA =  sv_0
     &          +   sv_1 * theta(i,j,k,bi,bj)
     &          + ( sv_2 * theta(i,j,k,bi,bj)**2 )
     &          + ( sv_3 * theta(i,j,k,bi,bj)**3 )
     &          + ( sv_4 * theta(i,j,k,bi,bj)**4 )

            !---------------------------------------------------
            ! Density of air (kg m^-3)
            !
            ! Eq. 25 from Johnson (2010) 
            !---------------------------------------------------
            rhoA = sd_0
     &          +   sd_1 * theta(i,j,k,bi,bj)
     &          + ( sd_2 * theta(i,j,k,bi,bj)**2 )
     &          + ( sd_3 * theta(i,j,k,bi,bj)**3 )     
            
            !---------------------------------------------------
            ! Diffusion coefficient in air (cm^2 sec^-1)
            !
            ! Eq.26-27 from Johnson (2010) 
            !---------------------------------------------------
            temp1 = (Ma + Mb) / (Ma * Mb)

            temp2 = Va**(1.0d0/3.0d0) +Vb**(1.0d0/3.0d0)

            DiffusivityAir (i,j) = 1.0d-3*(theta(i,j,k,bi,bj)+C2K)
     &       **1.75d0*temp1**0.5d0/(Patm*temp2**2.0d0) 
         
            !---------------------------------------------------
            ! Schmidt number in air (units?)
            !
            ! Eq. 23 from Johnson (2010)
            !---------------------------------------------------
            SchmidtAPCBa = etaA/rhoA/DiffusivityAir(i,j) 
     
            !---------------------------------------------------
            ! Diffusivity of the gas of interest in water (cm^2 sec^-1)
            !
            ! Dw = (7.4E-8*T*sqrt(Phi*Ms)) / (etas*Vb^0.6) [Eq. 33]
            !
            ! where T the temperature in Kelvin, Ms the molecular mass
            ! of the solvent, etas the dynamic viscosity of hte solvent
            ! in cP and Phi the association factor of the solvent and Vb
            ! the liquid molar volume.
            !
            ! Reference: Johnson et al. 2010, Wilke and Chang (1955)
            !---------------------------------------------------
            C2K = 273.15d0     ! convert Celcius to Kelvin

            ! dynamic viscosity of the solvent (water)
            !
            ! Xianming Zhang estimate this value from Table 2 in Johnson (2010)
            eta_s = 1.219d0

            ! Relative molecular mass of the solvent (water)
            Ms = 18.01d0

            ! Association factor of the solvent (water)
            PHI = 2.6d0

            DiffusivitySurfWater (i,j)
     &      = (7.4d-8 * ( theta(i,j,k,bi,bj)+C2K ) * SQRT( PHI * Ms ) )
     &      / ( eta_s * (Vb**0.6d0) ) 
        
            !---------------------------------------------------
            ! Dynamic viscosity of water (centipoise cP)
            !
            !etaW = (t+246) / (0.05594*t^2+5.2842t+137.37)  [Eq.31]
            ! where t the water temperature in Celsium
            !
            ! Reference: Johnson et al (2010)
            !---------------------------------------------------
            temp3 = 0.05594d0*theta(i,j,k,bi,bj)**2
     &            + 5.2842d0*theta(i,j,k,bi,bj) + 137.37d0

            etaW = (theta(i,j,k,bi,bj)+246d0)
     &           /(temp3)
           
            !---------------------------------------------------
            ! Schmidt number in water 
            !
            ! Scw = etaw/rhow/Dw                        [Eq. 29]
            ! 
            ! where etaw is the dynamic viscosity of water (cP = 10^3
            ! kg m^-3), rhow is the density of water 1.03d3 (kg m^-3),and Dw
            ! is the Diffusivity of water (cm2 s^-1)
            !
            ! Ref:  Johnson et al.  (2010) 
            !---------------------------------------------------
            SchmidtWPCBa = etaW/1.03d3/DiffusivitySurfWater (i,j)*10.0d0
           
            !---------------------------------------------------
            ! Single-phase air side transfer velocity (m/sec)
            !
            ! ka = 1E-3 + u*/
            !      (13.3*sqrt(SCa)+CD^(-0.5)-5+ln(SCa)/(2K)) [Eq.15]
            ! where u* id hyr gtviyion velocity, Sca is the Schmidt number of 
            ! the compound in air, CD is the drag coefficient and K is 
            !the Karman constant 0.4
            !
            ! 10^3*CD = 0.61+0.063*u10       [Eq.11, Johnson 2010]
            ! where CD is the drag coefficient and u10 is the windspeed
            ! at 10m
            !
            !Reference: Smith, S. D.: Wind stress and heat flux overthe
            !ocean in gale force winds, J. Phys. Oceanogr., 10,709–726,
            !1980
            !--------------------------------------------------- 
            temp4 =  windo
     &           *sqrt(6.1d-4+6.3d-5*windo) !  u*, Johsnon 2010
            
            temp5 = (0.61d0+0.063d0*windo) /1.0d3 ! CD, Johnson 2010
           
            temp6 = 13.3d0*sqrt(SchmidtAPCBa)
     &           +temp5**(-0.5d0)-5
     &           +log(SchmidtAPCBa)/2.d0/0.4d0
            
            ka = (1.0d-3+temp4/temp6)
           
            !---------------------------------------------------
            ! Single-phaser water side transfer velocity (m/sec)
            !
            ! Eq. 28 from Johnson (2010)
            !--------------------------------------------------- 
            temp7 = 0.222d0*windo*windo
     &            + 0.333d0*windo
            
            kw = temp7*(SchmidtWPCBa/600.d0)**(-0.5d0)
            kw = kw /360000d0 ! conversion from cm/hr to m/s
            
            !--------------------------------------------------- 
            ! Temperature dependence of Henry's law constant
            ! 
            ! Reference: Sander et al. (1999), Eqn 31
            !--------------------------------------------------- 
            KH   = KH0 * exp( (-dH_AW / R) 
     &          * (1/(theta(i,j,k,bi,bj)+C2K) - 1/TK0 )  )

            !---------------------------------------------------
            ! Salinity adjustment
            ! KH = KHO * 10^ (KS*S)                [Eq.19]
            ! where KS the Setchenow constat=0.35 and S the salinity
            ! in Molar, dividing by 58.44 turns mol/kg into Molar
            !
            ! Reference: Johnson et al. 2010             
            !---------------------------------------------------
            KH = KH * 10**(KS*salt(i,j,k,bi,bj)/58.44d0)
               
C(hma, 02 jul 2015)
C            ! Dimensionless liquid-over-gas Henry's law constant 
C            KHcc = KH * R * ( theta(i,j,k,bi,bj)+C2K )
C
C            ! Dimensionless gas-over-liquid Henry's law constant
C            KawD = 1.d0 / KHcc
C
C            !--------------------------------------------------- 
C            ! Total transfer velocity, gas-phase (m sec^-1)
C            !
C            ! From Johnson et al. (2010):
C            !   Kaexch = (1/ka + KawD/kw)^-1           [Eq. 4]                   
C            !--------------------------------------------------- 
C            Kaexch (i,j) = (1. _d 0 - iceo )
C     &          / ( 1.d0 / ka(i, j) + KawD / kw(i,j)  )  


            !--------------------------------------------------- 
            ! Total water-side transfer velocity (m sec^-1)
            !
            ! From Johnson et al. (2010):
            !   Kwexch = (1/kw + 1/(ka*KH))^-1           [Eq. 3]      
            !---------------------------------------------------
            Kwexch = (1. _d 0 - iceo )
     &          / ( 1.d0 / kw + 1.d0/(KH*ka) )  
           

            ! Turn back on ice rafting (ewl, 10/4/16)
            !--------------------------------------------------- 
            ! Air-sea exchange enchanced in seasonal sea ice
            ! region due to turbulance and mixing. 
            !
            ! Reference: Loose et al., 2014, Ocean Science
            ! Adopted from Yanxu Zhang's Hg package.
            !--------------------------------------------------- 
            IF ( (iceo .gt. 0.01 .and. iceo .lt. 0.99)) THEN
              Kwexch = 2.d0 * Kwexch
            ENDIF
            
            ! Increase another factor of 2 when fraction of open
            ! water is between 0.2 and 0.6. Non-linear effect in
            ! Figure 10 of Loose et al., 2014.
            IF ( iceo .gt. 0.4 .and. iceo .lt. 0.8 ) THEN
              Kwexch = 2.d0 * Kwexch
            ENDIF

            !--------------------------------------------------- 
            ! Flux of gas across the air-water interface 
            ! (mol m^-2 sec^-1)
            !
            ! From Johnson (2010):
            !   F_evas = Kwexch * (Cair / KH - Cwater)
            !
            ! Original reference is Liss & Slater (1974)
            !---------------------------------------------------
C pcbaconco is already in units of mol/m3 (hma, 29 jun 2015)
C            pcbaevaso =  Kaexch(i,j)
C     &        *(pcbaconco*1.0d-9/Mb !ng/m3 to mol/m3
C     &          - PTR(i,j,1,bi,bj,ipcba)*KawD )  

C Use the same expression for F as Yanxu uses for chemistry, which
C considers the total transfer velocity from the water's perspective
c (hma, 02 jul 2015)
C            pcbaevaso =  Kaexch(i,j)
C     &        *(pcbaconco - PTR(i,j,1,bi,bj,ipcba)*KawD ) 
C
             F_evas = Kwexch * (pcbaconco / KH  
     &                 - max(PTR(i,j,1,bi,bj,ipcba), 0d0) ) 

            ! Add up and down flux
            F_down        = Kwexch * (pcbaconco / KH)
            F_up          = Kwexch * max(PTR(i,j,1,bi,bj,ipcba), 0d0)  

            !---------------------------------------------------
            ! Total flux of PCBa (mol/m2/s)
            !---------------------------------------------------
            fluxPCBa(i,j) = F_evas 
     &                      + pcbarivero  
 
            !---------------------------------------------------
            ! Temperature (Celsius)
            !---------------------------------------------------
            surftemp           = theta(i,j,k,bi,bj)

            !---------------------------------------------------
            ! Salinity (g/kg)
            !---------------------------------------------------
            surfsalt           = salt(i,j,k,bi,bj)

            !---------------------------------------------------
            !  Fugacity
            !---------------------------------------------------
            pcbafug            = ( pcbaconco * 1.0d3 ) /
     &             ( (max(PTR(i,j,1,bi,bj,ipcba),1.0d-30)) * KH )
            ! atm concentration in pg/m3 and seawater concnetration in
            ! pg/l (see Zhang et al. 2015.
            ! Environ.Sci.Technol,49,13760-13767
            ! Use min value of 1e-30 to avoid div by zero (ewl, 1/25/17)

            !---------------------------------------------------            
            ! Calculate diagnostics
            !--------------------------------------------------- 
#ifdef ALLOW_DIAGNOSTICS
            ! Air concentration (mol m^-3)
            pcbaconc_diag(i,j)  = pcbaconco
 
            ! From rivers (mol m^-2 sec^-1)
            pcbariver_diag(i,j) = pcbarivero

            ! Ocean-to-air evasion (mol m^-2 sec^-1) 
            pcbaevas_diag(i,j)  = -F_evas  

            ! Air sea exchange fluxes, upward and downward (mol m-2 sec-1)
            ASF_up_diag(i,j)    = F_up 
            ASF_dn_diag(i,j)    = F_down           

            ! Inputs from atmospheric deposition (mol m^-2 sec^-1)
            pcbadep_diag(i,j)   = pcbadepoo

            ! Temperature (Celsius)
            temp_diag(i,j)      = surftemp

            ! Salinity (g/kg)
            salt_diag(i,j)      = surfsalt

            !  Fugacity
            pcbafug_diag(i,j)   = pcbafug
#endif       
           ELSE
             fluxPCBa(i,j) = 0. _d 0
           ENDIF

         ENDDO
        ENDDO

        !-------------------------------------------------------
        ! Update surface tendencies
        !-------------------------------------------------------
        DO j=jMin,jMax
          DO i=iMin,iMax
           SURFPCBa(i,j)=
     &         HFacC(i,j,1,bi,bj)*FluxPCBa(i,j)*recip_drF(1)
          ENDDO
         ENDDO

        !-------------------------------------------------------
        ! Write diagnostics
        !-------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
        ! Atmospheric concentration (mol/m3)
        WRITE(diagname,'(A8)') 'PCBaCONC'
        CALL DIAGNOSTICS_FILL( pcbaconc_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )

        ! Inputs from rivers (mol/m2/sec)
        WRITE(diagname,'(A8)') 'PCBaRIVR'
        CALL DIAGNOSTICS_FILL( pcbariver_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )

        ! Atmospheric deposistion (mol/m2/sec)
        WRITE(diagname,'(A8)') 'PCBaDEPO'
        CALL DIAGNOSTICS_FILL( pcbadep_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )

        ! Atmospheric deposistion (mol/m2/sec)
        WRITE(diagname,'(A8)') 'PCBaEVAS'
        CALL DIAGNOSTICS_FILL( pcbaevas_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )     

       ! Surface water Temperature (Celsius)
       WRITE(diagname,'(A8)') 'SURFTEMP'
       CALL DIAGNOSTICS_FILL( temp_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid ) 

       ! Surface water salinity (g/kg)
       WRITE(diagname,'(A8)') 'SURFSALT'
       CALL DIAGNOSTICS_FILL( salt_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )

       ! PCBa fugacity
       WRITE(diagname,'(A8)') 'PCBaFUGA'
       CALL DIAGNOSTICS_FILL( pcbafug_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )

       ! PCBa Upward flux
       WRITE(diagname, '(A8)') 'PCBaFLUP'
       CALL DIAGNOSTICS_FILL( ASF_up_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )

       ! PCBa Downward Flux
       WRITE(diagname, '(A8)') 'PCBaFLDN'
       CALL DIAGNOSTICS_FILL( ASF_dn_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )


#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_PCB */
#endif /* ALLOW_PTRACERS */

      RETURN 
      END




        
        

        
        
        
