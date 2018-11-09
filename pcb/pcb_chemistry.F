C $Header: /u/gcmpack/MITgcm/pkg/hg/hg_chemistry.F,v 1.10 2010/03/16 00:14:47 jmc Exp $
C $Name: checkpoint64g $

#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

#ifdef USE_EXFIWR
#include "EXF_OPTIONS.h"
#include "SEAICE_OPTIONS.h"
#endif /* USE_EXFIWR */

CStartOfInterFace
      SUBROUTINE PCB_CHEMISTRY( Ptr,
     &                         Gchem,    
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)

      !==========================================================
      ! SUBROUTINE PCB_CHEMISTRY                                  
      !   Calculate the changes to dissolved-phased PCBs due to  
      !   degradation.                                           
      !
      ! REVISION HISTORY
      !   14 Jan 2015 - H. Amos - Adopted copy of hg_chemistry.F
      !                           for PCB degradation.
      !   19 May 2015 - H. Amos - Update USE_EXFIWR calls based on Hg
      !                           so you can use online ice, radiation, 
      !                           etc. when running with ECCOv4. Add
      !                           call PCB_SUNCOS function.
      !   01 Oct 2015 - H. Amos - Add temperature dependence to
      !                           degradation.
      !==========================================================
      IMPLICIT NONE

      !--------------------------------------------------------
      ! Global variables
      !--------------------------------------------------------
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#include "PCB_SIZE.h"
#include "PCB.h"

#ifdef USE_EXFIWR
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#include "EXF_FIELDS.h"
#include "FFIELDS.h"
#endif /* USE_EXFIWR */

      !--------------------------------------------------------
      ! Routine arguments
      !--------------------------------------------------------
      _RL myTime
      _RL Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)    
      INTEGER myThid
      INTEGER myIter
#ifdef ALLOW_DIAGNOSTICS
      _RL pcbaDegrad_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr) ! PCBa lost via degradation      
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */
      _RL,SAVE :: meansuncos(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      LOGICAL,SAVE :: FIRST
      DATA FIRST /.TRUE./
      INTEGER nl
      INTEGER i,j,k,l,m,bi, bj, imin, imax, jmin, jmax
      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
      !--------------------------------------------------------
      ! Local variables
      !--------------------------------------------------------
      _RL radz    , radz0    , ecz1, ecz2
      _RL suncos
      _RL pcba_old, pcba_new
      _RL Kdeg_pcb, Kdeg_pcb_base
      _RL T_CORREC ! HM Amos, added 01 Oct 2015
      _RL etz1, etz2 ! cw, 03 Nov 2016

      !--------------------------------------------------------
      ! Parameters
      !--------------------------------------------------------

      ! Degradation
      _RL R                ! universal gas consntant [J/mol/K]
      _RL TK0               ! standard temperature    [K]
      _RL C2K              ! conversion factor, degrees C to K
      _RL EACT             ! activation energy       [J/mol]
      data R    /8.3145d0/
      data TK0   /298d0/
      data C2K  /273.15d0/
      data EACT /30000d0/  ! value selected by Friedman & Selin 
                           ! [2015] for consistency with Lamon
                           ! et al. [2009] and Gouin et al. [2013]
      

      ! Radiation extinction
      _RL EC_w    ! 1/m
      _RL EC_doc  ! 1/m
      _RL EC_chla ! 1/m
      data EC_w /0.0145 _d 0/
      data EC_doc /0.026 _d 0/
      data EC_chla /3.1 _d 1/

      !--------------------------------------------------------
      ! Local variables      
      !--------------------------------------------------------
      _RL iceo,nppo,chlo,fremino,doco,radswo,radswo0,meank     
      _RL myTimeo

      !--------------------------------------------------------
      ! Initialize diagnostics
      !--------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            pcbaDegrad_diag(i,j,k)=0d0
          ENDDO
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */

      !--------------------------------------------------------
      ! Calculate monthly mean suncos
      !--------------------------------------------------------
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
            CALL PCB_SUNCOS(i,j,bi,bj,myTimeo,suncos)  
            meansuncos(i,j)=meansuncos(i,j)+suncos
           ENDDO
           meansuncos(i,j)=meansuncos(i,j)/dble(nl)   
         ENDIF                                   
        ENDDO
       ENDDO
      ENDIF
#endif /* USE_EXFIWR */
      
      DO j=jMin,jMax
       DO i=iMin,iMax
       IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN

        !/////////////////////////////////////////////////////////////////////
        ! UPDATED 23 Aug 2016, CC Wagner
        ! to account for light
        ! attenuation in the ocean, degradations is now linearly
        ! dependent on light penetration, code was adapted from Yanxu
        ! Zhang's Hg-code and Soerensen et al. 2015
        !/////////////////////////////////////////////////////////////////////

        !------------------------------------------------------------------------------
        ! Compute baseline radiation at ocean surface
        !------------------------------------------------------------------------------
        ! radswo0 = -Qsw(i,j,bi,bj)  ! W/m2
        !         = - ( downward - ice and snow absorption - reflected)
        ! where Qsw is the net upward shortwave radiation in W/m2
        ! Typical range: -350 < Qsw < 0
        ! If sea ice package USE_EXFIWR is enabled absorption due to ice
        ! and
        ! snow
        ! absorption and reflection upward is already taken into account
#ifdef USE_EXFIWR
        iceo = AREA(i,j,bi,bj)
        radswo0 = -Qsw(i,j,bi,bj)
        radswo = radswo0
#else
        iceo = fice(i,j,bi,bj)
        radswo0 = max(radsw(i,j,bi,bj),0d0)*(1d0-iceo)  ! wm-2

        ! Impose a diurnal cycle
        CALL PCB_SUNCOS(i,j,bi,bj,myTime,suncos)
        radswo = radswo0 * suncos / (meansuncos(i,j) + 1. _d -32)
#endif /* USE_EXFIWR */

        !------------------------------------------------------------------------------
        ! Radiation attenuation calculation
        !------------------------------------------------------------------------------
        ! Ref: Soerensen et al. 2016, SI, Table S5

        !--------------------------------------
        ! Intiate variables eta*x1 and eta*x2
        !-------------------------------------
        etz1 = 0d0      !cwag, 1 Sep 2016
        etz2 = 0d0      !cwag, 1 Sep 2016

        !--------------------------------------------------------
        ! Localize variables
        !--------------------------------------------------------
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN

           !--------------------------------------------------------

           !--------------------------------------------------------
           ! Local concentrations 
           !--------------------------------------------------------
           chlo =    max(chl(i,j,k,bi,bj),0d0)   ! chlorophyll conc.
           doco =    max(doc(i,j,k,bi,bj),0d0)   ! DOC concentration
           fremino = max(Fremin(i,j,k,bi,bj), 0d0) ! remineralisation fraction
           pcba_old = max( PTR(i,j,k,bi,bj,ipcba), 0d0)

          !---------------------------------------------------------
          ! Compute radiation attenuation
          !--------------------------------------------------------
          ! adapted from Soerensen et al. 2015 a, Table S5 and Yanxu's
          ! Hg Code
          !       RAD    = 1/(x2-x1) *
          !       Ri/eta*(exp(-eta*x1)-exp(-eta*x2))
          !  <=>  RAD/Ri = 1/(eta*x2-eta*x1) *
          !  (exp(-eta*x1)-exp(-eta*x2))
          ! where RAD is the radiation at a specific depth, Ri (here
          ! radswo) is the base line radiation at the ocean surface,
          ! eta is the extinction coefficient for radiation, due to
          ! water, pigments and dissolved organic carbon
           etz2 = etz1 + drF(k) * ( EC_w
     &           + EC_doc * doco * 1.2d-2   ! mmol C/m3 -> mg C/L
     &           + EC_chla * chlo * 1.0d-3 ) ! mg Chl/m3 -> mg Chl/L
           meank = 1. _d 0 / (etz2 - etz1)
     &          * ( EXP( -etz1 ) - EXP( -etz2 ) ) ! equivalent to RAD/Ri
           radz = radswo * meank ! RAD at specific depth z [W/m2]
           ! in Soerensen et al. 2016 this corresponds to RAD*Frad,
           ! because
           ! radswo already accounts for sea ice cover and reflection

           etz1 = etz2  ! etz1 will be the box' surface for the next layer

!Updated to be consistent with Carey Friedman's GEOS-Chem simulation
!and Matt McCleod's BETR-GLOBAL models.
!01 October 2015
!
!           !--------------------------------------------------------
!           ! Temperature dependence of degradation
!           ! Kdeg_pcb has units of 1/sec
!           !--------------------------------------------------------
!           Kdeg_pcb = a_deg * EXP( b_deg * theta(i,j,k,bi,bj) )

           !--------------------------------------------------------
           ! Temperature dependence of degradation
           !   Kdeg = a* ( ln(2)/tHalf )
           !   where a is a temperature correction
           !--------------------------------------------------------
           ! Temperature dependence
           T_CORREC = EXP((-EACT / R) 
     &                   * (1d0/(theta(i,j,k,bi,bj)+C2K) - 1d0/TK0))

           ! Calculate degradation rate [1/sec]
           Kdeg_pcb_base = 0.693d0 / tHalf /3600d0    ! convert [1/hr] --> [1/sec]
           
           Kdeg_pcb = 100
     &              * ( 0.9d0 * (Kdeg_pcb_base/172d0)*radz        ! apply light attenuation
     &              +   0.1d0 * Kdeg_pcb_base * fremino)         ! apply remineralisation rate
           ! rad(z) in [W/m2] at specific depth

           Kdeg_pcb = T_CORREC * Kdeg_pcb  ! apply temperature dependence
           !--------------------------------------------------------
           ! Calculate new concentration of PCBa after 
           ! degradation (photolytic + biodegradation)
           !--------------------------------------------------------
           pcba_new = pcba_old - Kdeg_pcb * pcba_old * deltaTclock

           !--------------------------------------------------------
           ! Calculate diagnostics                                             
           !--------------------------------------------------------     
#ifdef ALLOW_DIAGNOSTICS
           ! bulk degradation (mol m^-3 sec^-1)
           pcbaDegrad_diag(i,j,k) = Kdeg_pcb * pcba_old 
#endif /* ALLOW_DIAGNOSTICS */

           !--------------------------------------------------------
           ! Accumulate tendency
           !--------------------------------------------------------
           Gchem(i,j,k,bi,bj,ipcba)=Gchem(i,j,k,bi,bj,ipcba)+
     &                             (pcba_new - pcba_old)
     &                             /deltaTclock            

         ENDIF
        ENDDO
        
       endif
       ENDDO
      ENDDO

        !--------------------------------------------------------      
        !  Write diagnostics
        !--------------------------------------------------------      
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'PCBaDEGR'
        CALL DIAGNOSTICS_FILL( pcbaDegrad_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_PCB */
#endif /* ALLOW_PTRACERS */

       RETURN
       END
