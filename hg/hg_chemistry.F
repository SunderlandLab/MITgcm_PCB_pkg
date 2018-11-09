C $Header: /u/gcmpack/MITgcm/pkg/hg/hg_chemistry.F,v 1.10 2010/03/16 00:14:47 jmc Exp $
C $Name: checkpoint64g $

C To calculate the chemistry of Hg species in the ocean, yxzhang, 5/30/13
C

#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"
#ifdef USE_EXFIWR
#include "EXF_OPTIONS.h"
#include "SEAICE_OPTIONS.h"
#endif

CStartOfInterFace
      SUBROUTINE HG_CHEMISTRY( Ptr,
     &                         Gchem,    
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)

C     /==========================================================\
C     | SUBROUTINE HG_CHEMISTRY                                  |
C     | o Calculate the changes to Hg species because            |
C     |   photo- and bio-chemical reactions, and food web        |
C     |==========================================================|
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#include "HG_SIZE.h"
#include "HG.h"
#include "chem_hg_Rates.h"
#include "chem_hg_Parameters.h"
#ifdef USE_EXFIWR
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#include "EXF_FIELDS.h"
#include "FFIELDS.h"
#endif

C     == Routine arguments ==
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)    
#ifdef ALLOW_DIAGNOSTICS      
      _RL  hg2redn_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  hg2redb_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  hg2redp_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  hg0oxdb_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  hg0oxdp_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  hg0oxdd_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef METHY_MODE
      _RL   hg2metb_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  mmhgdmeb_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  mmhgdmep_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  mmhgmetb_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  dmhgdmep_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
      CHARACTER*8 diagname
#endif
      _RL,SAVE :: meansuncos(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      LOGICAL,SAVE :: FIRST
      DATA FIRST /.TRUE./
      INTEGER nl
      INTEGER i,j,k,l,m,bi, bj, imin, imax, jmin, jmax
      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG
C     == Local variables ==
      _RL radz, radz0, ecz1, ecz2
      _RL suncos
      _RL hg0_old, hg2_old, hg0_new, hg2_new
      _RL mmhg_old,dmhg_old,mmhg_new,dmhg_new
      _RL CONCN(NVAR)
      _RL concn_scale, small_num, small_num2
      _RL inv_concn_scale, inv_concn_scale_dt
      PARAMETER ( small_num = 1.0 _d -8 )
      PARAMETER ( small_num2 = 1.0 _d -24 )     
      PARAMETER ( concn_scale = 1.0 _d 12 )    ! transfer mol/m3 to fM
      PARAMETER ( inv_concn_scale = 1.0 _d -12 )

C     == Parameters
C     -- Radiation extinction
      _RL EC_w    ! m-1
      _RL EC_doc  ! m-1
      _RL EC_chla ! m-1
      data EC_w /0.0145 _d 0/
C     data EC_doc /0.654 _d 0/ ! from Soerensen et al. 2010, too large
      data EC_doc /0.026 _d 0/
      data EC_chla /3.1 _d 1/

C     -- Reaction base, scale factors
C     -- Inorganic species, from Soerensen 2010
      _RL k_hg2_hg0_rad  ! s-1
      _RL k_hg2_hg0_bio
      _RL k_hg0_hg2_rad
      _RL k_hg0_hg2_dak
      _RL k_hg0_hg2_bio

      _RL k_hg2_hg0_rad_base
      _RL k_hg2_hg0_bio_base
      _RL k_hg0_hg2_rad_base
      _RL k_hg0_hg2_dak_base
      _RL k_hg0_hg2_bio_base
      _RL k_hg2_red_fraction      
      data k_hg2_hg0_rad_base / 1.56d-6/
C     Sorensen 2010, based on production 
C     data k_hg2_hg0_bio_base / 4.51d-10/  s-1 (mgCm-2d-1 npp)-1
C     based on remin, through entire water column, s-1(mgCm-3d-1 fremin)-1
#ifdef ISHG2LABILE
      data k_hg2_hg0_bio_base / 2.48d-7/ !/ 6.2d-7 /
#else
      data k_hg2_hg0_bio_base / 6.2d-7 /
#endif
C     assume Hg0:Hg2 in deep ocean = 1:4      
      data k_hg0_hg2_bio_base / 9.9d-7 /
C     Soerensen uses 6.64d-6 for k_hg0_hg2_rad_base
      data k_hg0_hg2_rad_base / 4.7d-6 /
C     Sorensen 2010, constant
C     now scaled based on monthly mean rad
C     data k_hg0_hg2_dak_base / 1d-7/
      data k_hg0_hg2_dak_base / 1.3d-9/
#ifndef VARI_RED_POOL
#ifdef ISHG2LABILE
      data k_hg2_red_fraction / 1.0d0 / !/4.0 _d -1/
#else
      data k_hg2_red_fraction /4.0 _d -1/
#endif
#endif

C     -- Methylation and demethylation
      _RL k_hg2_mmhg_bio_base
      _RL k_mmhg_hg2_rad_base
      _RL k_mmhg_hg0_bio_base
      _RL k_mmhg_hg0_che      
      _RL k_mmhg_dmhg_bio_base
      _RL k_dmhg_mmhg_rad_base
      _RL k_rad2par
      
#ifdef METHY_MODE
      !anls data k_hg2_mmhg_bio_base / 4.0d-7 / *(1-pe)*CHL  mg/m3         !2.3d-9 to 5d-7
      !anls data k_mmhg_hg2_rad_base / 3.0d-6 / *PARvis      E/m2/d
      !anls data k_mmhg_hg0_bio_base / 4.6d-8 / *CHL         mg/m3         !4.0d-3 / 86400, (0.1-3) before vs039 =0.3d-3
      !anls data k_mmhg_dmhg_bio_base /1.0d-7 / *CHL         mg/m3
      !anls data k_dmhg_mmhg_rad_base /8.7d-9 / *PARvis      E/m2/d
      ! In MITgcm, averagely, chl (mg/m3) = 0.6 * fremin
      ! (mmolC/m3/d), 
      ! Soerensen estimates remin = (1-pe)*chl
      ! Globally, pe has a value of approximately 20%, so this
      ! means: remin = 0.8 * chl          
      !data k_rad2par /4.32d-1/  !W/m2 -> E/m2/d       

      data k_hg2_mmhg_bio_base  / 2.0d-7 /      
      data k_mmhg_hg2_rad_base  / 1.94d-5 /      
      !data k_mmhg_hg0_bio_base  / 2.76d-8 / ! now obselete
      data k_mmhg_hg0_che       / 9.2d-10 /                  
      data k_mmhg_dmhg_bio_base / 6.0d-8 /
      data k_dmhg_mmhg_rad_base / 3.76d-9 /
#else
      data k_hg2_mmhg_bio_base / 0.0 _d 0 /      
      data k_mmhg_hg2_rad_base /0.0 _d 0 /      
      !data k_mmhg_hg0_bio_base /0.0 _d 0/
      data k_mmhg_hg0_che     /0.0 _d 0 /                  
      data k_mmhg_dmhg_bio_base /0.0 _d 0/
      data k_dmhg_mmhg_rad_base / 0.0 _d 0/
#endif

C     ==local variables      
      _RL iceo,nppo,chlo,fremino,doco,o2o,dfeo
      _RL radswo,radswo0,meank,tco     
      _RL myTimeo
      
C     == Start calculation
C     ==time variables
      inv_concn_scale_dt = inv_concn_scale / deltaTclock

C     == calculate monthly mean suncos
#ifndef USE_EXFIWR
      IF ( mod( myTime, 2629800.0 ) .lt. 0.1 .or. FIRST) THEN
       FIRST = .FALSE.   
       nl = int(2629800 / deltaTclock)
       DO j=jMin,jMax
        DO i=iMin,iMax
         IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN
           meansuncos(i,j)=0. _d 0  
           DO l=1,nl
            myTimeo=myTime+dble(l-1)*deltaTclock  
            CALL hg_suncos(i,j,bi,bj,myTimeo,suncos)  
            meansuncos(i,j)=meansuncos(i,j)+suncos
           ENDDO
           meansuncos(i,j)=meansuncos(i,j)/dble(nl)   
         ENDIF                                   
        ENDDO
       ENDDO
      ENDIF
#endif
      
      DO j=jMin,jMax
       DO i=iMin,iMax
       IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN
#ifdef USE_EXFIWR
        iceo = AREA(i,j,bi,bj)
C       radswo0 = swdown(i,j,bi,bj)*(1d0-iceo)
C       Qsw   :: Net upward shortwave radiation in W/m^2
C              Qsw = - ( downward - ice and snow absorption - reflected)
C              > 0 for decrease in theta (ocean cooling)
C              Typical range: -350 < Qsw < 0  
C       If seaice package is enabled, ice fraction is already taken into
C       account.     
        radswo0 = -Qsw(i,j,bi,bj)
        radswo = radswo0
        tco = theta(i,j,1,bi,bj)
#else
        iceo = fice(i,j,bi,bj)
        radswo0 = max(radsw(i,j,bi,bj),0d0)*(1d0-iceo)  ! wm-2
C       impose a diurnal cycle
        CALL hg_suncos(i,j,bi,bj,myTime,suncos)
        radswo = radswo0 * suncos / (meansuncos(i,j) + 1. _d -32)
        ! tco = ???
#endif   
        dfeo = max(dfe(i,j,bi,bj),0d0) 
C     == Radiation attenuation calculation
        ecz1 = 0d0
        ecz2 = 0d0
C     == localize variables
                          
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
C     == localize variables
           nppo =    max(npp(i,j,k,bi,bj),0d0)
           fremino = max(Fremin(i,j,k,bi,bj),0d0)
           chlo =    max(chl(i,j,k,bi,bj),0d0)
           doco =    max(doc(i,j,k,bi,bj),0d0)
           o2o  =    max(o2(i,j,k,bi,bj),0d0)

           ecz2 = ecz1 + drF(k) * ( EC_w
     &           + EC_doc * doco * 1.2d-2   ! mmolC/m3 -> mg/L
     &           + EC_chla * chlo * 1.0d-3 ) ! mg chl/m3 -> mg/L
           meank = 1. _d 0 / (ecz2 - ecz1)
     &          * ( exp( -ecz1 ) - exp( -ecz2 ) )
           radz = radswo * meank
           radz0 = radswo0 * meank ! for dark oxidation
           ecz1 = ecz2                  

C     == Rate coefficent calculation
c     -- inorganic species
#ifdef VARI_RED_POOL
C          assume a linear relationship between the HgII reducible
C          pool and the DOC concentrations
C          based on Soerensen et al. (2013)
C          coast: DOC 110-150 uM, Hg0% = 9%
C          Bermuda: DOC 60-70 uM, Hg0% = 18%
C          reducible pool: (130 uM, 20%) (45 uM, 40%)
C          fraction(%)= -0.235*([DOC]-45)+40 
C          45 is the average surface ocean DOC, uM
#ifdef ISHG2LABILE
           k_hg2_red_fraction = 1.0d0
#else
           k_hg2_red_fraction = max((-0.235*(doco-45.0)+40.0)/100.0,0d0)
#endif
#endif

           RC3 = k_hg2_hg0_rad_base * radz * k_hg2_red_fraction
           RC4 = k_hg2_hg0_bio_base * fremino * 
     &                     k_hg2_red_fraction  ! mmolCm-3d-1
c          coast
           IF (doco .gt. 4.99 _d 2) THEN
                RC3 = RC3 * 0.5
                RC4 = RC4 * 0.5
           ENDIF     
 
           RC1 = k_hg0_hg2_rad_base * radz
           RC5 = k_hg0_hg2_bio_base * fremino  ! mmolCm-3d-1
           RC2 = k_hg0_hg2_dak_base * radz0

#ifdef METHY_MODE          
c     -- methylation/demethylation
           IF ( dfeo .gt. 5.0e-13 ) THEN
             RC6 = k_hg2_mmhg_bio_base * fremino * 
c    &             300.0 / ( o2o + small_num )   ! anerobic condition boosts it 
     &             dfeo / 2.0e-13             ! iron deposition enhances methylation
           ELSE
             RC6 = k_hg2_mmhg_bio_base * fremino
           ENDIF 
           RC8 = k_mmhg_hg2_rad_base * radz 
           !RC7 = k_mmhg_hg0_bio_base * fremino  ! this also happens at subsurface waters
           RC7 = k_mmhg_hg0_che *
     &            0.69*exp(-5457.3*(1.0/(273.15+tco)-1.0/(293.15)))
           ! because [MeHg]/[HgII] is proportional to fremin
           ! cold temperature helps reserve MMHg

c     -- transfer between mmhg and dmhg                      
           RC9 = k_mmhg_dmhg_bio_base * fremino
           RC10 = k_dmhg_mmhg_rad_base * radz
#else                      
           RC6 = 0.0 _d 0
           RC7 = 0.0 _d 0
           RC8 = 0.0 _d 0
           RC9 = 0.0 _d 0
           RC10 = 0.0 _d 0
#endif
           
C     == Local concentrations
           hg0_old = max( PTR(i,j,k,bi,bj,ihg0), 0d0)
           hg2_old = max( PTR(i,j,k,bi,bj,ihg2), 0d0)
#ifdef METHY_MODE
           mmhg_old = max( PTR(i,j,k,bi,bj,immhg), 0d0)
           dmhg_old = max( PTR(i,j,k,bi,bj,idmhg), 0d0)
#else
           mmhg_old = 0.0 _d 0
           dmhg_old = 0.0 _d 0
#endif

C     == Feed the KPP solver
C     -- Tracer concentrations
           CONCN(ind_Hg2_aq) = max( hg2_old  * concn_scale, 0. _d 0 )
           CONCN(ind_Hg0_aq) = max( hg0_old  * concn_scale, 0. _d 0 )
           CONCN(ind_MMHg) = max( mmhg_old * concn_scale, 0. _d 0 )
           CONCN(ind_DMHg) = max( dmhg_old * concn_scale, 0. _d 0 )
C     -- Reaction rates, Note: they behave like tracers
           CONCN(ind_RR10) = 0d0
           CONCN(ind_RR9 ) = 0d0
           CONCN(ind_RR8 ) = 0d0
           CONCN(ind_RR7 ) = 0d0
           CONCN(ind_RR6 ) = 0d0
           CONCN(ind_RR5 ) = 0d0
           CONCN(ind_RR4 ) = 0d0
           CONCN(ind_RR3 ) = 0d0
           CONCN(ind_RR2 ) = 0d0
           CONCN(ind_RR1 ) = 0d0

C      -- Safety check before calling KPP solver:
C         at least one tracer larger than 10^-8 fM
C         at least one reaction coefficient larger than 10^-32 s-1
           IF ( ( CONCN(ind_DMHg) .GT. small_num .OR. 
     &            CONCN(ind_MMHg) .GT. small_num .OR.
     &            CONCN(ind_Hg0_aq) .GT. small_num .OR. 
     &            CONCN(ind_Hg2_aq) .GT. small_num ) .AND.
     &          ( RC1 .GT. small_num2 .OR. 
     &            RC2 .GT. small_num2 .OR. 
     &            RC3 .GT. small_num2 .OR.
     &            RC4 .GT. small_num2 .OR.
     &            RC5 .GT. small_num2 .OR.
     &            RC6 .GT. small_num2 .OR.
     &            RC7 .GT. small_num2 .OR.
     &            RC8 .GT. small_num2 .OR.
     &            RC9 .GT. small_num2 .OR.
     &           RC10 .GT. small_num2 ) ) THEN
              CALL CHEM_HG_KPP( deltaTclock, CONCN )
           ENDIF

           IF ( isnan(CONCN(ind_DMHg)) .OR.
     &          isnan(CONCN(ind_MMHg)) .OR.
     &          isnan(CONCN(ind_Hg0_aq)) .OR.
     &          isnan(CONCN(ind_Hg2_aq)) ) THEN
              dmhg_new = dmhg_old
              mmhg_new = mmhg_old
              hg0_new  = hg0_old
              hg2_new  = hg2_old
           ELSE
              dmhg_new = CONCN(ind_DMHg) * inv_concn_scale
              mmhg_new = CONCN(ind_MMHg) * inv_concn_scale
              hg0_new  = CONCN(ind_Hg0_aq) * inv_concn_scale
              hg2_new  = CONCN(ind_Hg2_aq) * inv_concn_scale
           ENDIF

C Calculate diagnistics                                                  
#ifdef ALLOW_DIAGNOSTICS
            hg2redn_diag(i,j,k) = ( CONCN(ind_RR3) + CONCN(ind_RR4)   !net reduction, mol/m^3/s
     &                  - CONCN(ind_RR1) - CONCN(ind_RR2)
     &                  - CONCN(ind_RR5) )
     &                  * inv_concn_scale_dt
            hg2redb_diag(i,j,k) = CONCN(ind_RR4) * inv_concn_scale_dt
            hg2redp_diag(i,j,k) = CONCN(ind_RR3) * inv_concn_scale_dt
            hg0oxdb_diag(i,j,k) = CONCN(ind_RR5) * inv_concn_scale_dt
            hg0oxdp_diag(i,j,k) = CONCN(ind_RR1) * inv_concn_scale_dt
            hg0oxdd_diag(i,j,k) = CONCN(ind_RR2) * inv_concn_scale_dt
#ifdef METHY_MODE
             hg2metb_diag(i,j,k) = CONCN(ind_RR6) * inv_concn_scale_dt 
            mmhgdmeb_diag(i,j,k) = CONCN(ind_RR7) * inv_concn_scale_dt
            mmhgdmep_diag(i,j,k) = CONCN(ind_RR8) * inv_concn_scale_dt
            mmhgmetb_diag(i,j,k) = CONCN(ind_RR9) * inv_concn_scale_dt
            dmhgdmep_diag(i,j,k) = CONCN(ind_RR10) * inv_concn_scale_dt
#endif
#endif

C     == Accumulate tendency
           Gchem(i,j,k,bi,bj,ihg0)=Gchem(i,j,k,bi,bj,ihg0)+
     &                             (hg0_new-hg0_old)
     &                             /deltaTclock            
           Gchem(i,j,k,bi,bj,ihg2)=Gchem(i,j,k,bi,bj,ihg2)+
     &                             (hg2_new-hg2_old)
     &                             /deltaTclock                       
#ifdef METHY_MODE           
           Gchem(i,j,k,bi,bj,immhg)=Gchem(i,j,k,bi,bj,immhg)+
     &                              (mmhg_new-mmhg_old)
     &                             /deltaTclock     
           Gchem(i,j,k,bi,bj,idmhg)=Gchem(i,j,k,bi,bj,idmhg)+
     &                              (dmhg_new-dmhg_old)
     &                             /deltaTclock     
#endif

         ENDIF
        ENDDO
        
       endif
       ENDDO
      ENDDO
      
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HG2REDN '
        CALL DIAGNOSTICS_FILL( hg2redn_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'HG2RBIO '
        CALL DIAGNOSTICS_FILL( hg2redb_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'HG2RPHO '
        CALL DIAGNOSTICS_FILL( hg2redp_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'HG0OBIO '
        CALL DIAGNOSTICS_FILL( hg0oxdb_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'HG0OPHO '
        CALL DIAGNOSTICS_FILL( hg0oxdp_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'HG0ODAK '
        CALL DIAGNOSTICS_FILL( hg0oxdd_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#ifdef METHY_MODE
        WRITE(diagname,'(A8)') 'HG2METB '
        CALL DIAGNOSTICS_FILL( hg2metb_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'MMHGDMEB'
        CALL DIAGNOSTICS_FILL( mmhgdmeb_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'MMHGDMEP'
        CALL DIAGNOSTICS_FILL( mmhgdmep_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'MMHGMETB'
        CALL DIAGNOSTICS_FILL( mmhgmetb_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
        WRITE(diagname,'(A8)') 'DMHGDMEP'
        CALL DIAGNOSTICS_FILL( dmhgdmep_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif
#endif

#endif
#endif

       RETURN
       END
