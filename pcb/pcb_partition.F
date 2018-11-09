C $Header: /MITgcm/pkg/pcb/pcb_partition.F 

C==============================================================
C OBJECTIVES
C   Calculate paritioning between dissolved and particulate 
C   PCB phases in the ocean. 
C
C REVISION HISTORY
C  From Xianming Zhang:
C    Based on the code of Hg species in the ocean by Yanxu Zhang.
C    To calculate the partition of PFOS species in the ocean
C    Local equlibrium between the anionic, neutral, and 
C    particulate species are assumed to be reached instantaneously
C
C  06 Aug 2015 - hma - PFOS code adopted for the simulation of 
C                      PCBs in the ocean. Comments added and 
C                      obsolete/commented-out code deleted to 
C                      clean things up.
C  06 Jul 2015 - hma - Move Koc and Kow to pcb_param.F
C  05 Oct 2015 - hma - Updated Kow-Koc relationship
C===============================================================

#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

CStartOfInterFace
      SUBROUTINE PCB_PARTITION( Ptr,
     &                          bi    , bj    ,
     &                          imin  , imax  , jmin  , jmax ,
     &                          myIter, myTime, myThid        )

       
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
#include "PCB_SIZE.h"
#include "PCB.h"

      !--------------------------------------------------------
      ! Routine arguments
      !--------------------------------------------------------
      INTEGER myIter
      _RL     myTime
      INTEGER myThid
      _RL  Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
#ifdef ALLOW_DIAGNOSTICS      
      _RL     pcbapart_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL     pcbppart_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)            
      CHARACTER*8 diagname
#endif /* ALLOW_DIAGNOSTICS */
      _RL     TOT_PCB
      _RL     FRAC_TotDissolved
      _RL     FRAC_Dissolved
      _RL     FRAC_POCBound
      _RL     poco, doco
      INTEGER i,j,k,bi, bj, imin, imax, jmin, jmax
      _RL     pcba_old, pcbp_old, pcba_new, pcbp_new
      _RL     logKow_T, logKoc

      !--------------------------------------------------------
      ! local Parameters
      !--------------------------------------------------------      
      _RL R, TK0, KS, C2K

      R    = 8.3145d0
      TK0  = 298d0
      C2K  = 273.15d0
      KS   = 0.35d0

C Kow moved to pcb_param.F (hma, 06 jul 2015)
C      _RL     Kow   ! octanol-water partitioning coefficient <-- put this in pcb_param.F?
C      _RL     pcba_deg   ! degradation of dissolved-phase PCB (hma, 14 Jan 2015)

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB

C     == Parameters
C     -- Partition coefficient
c      _RL Koc_pcb  ! L/kg
c      data  Koc_pcb /372./ 
c xz: declared as a common variable in PCB.h and give the value in pcb_param.F
c  xz: Koc for PFOS 372 L/kg, need to be changed for other PCBs
c #ifdef METHY_MODE 
c       _RL Kd_mmhg_poc  ! L/kg
c       data Kd_mmhg_poc /126491./
c #endif

      !--------------------------------------------------------
      ! Initialize diagnostics
      !--------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            pcbapart_diag(i,j,k) = 0d0  
            pcbppart_diag(i,j,k) = 0d0         
          ENDDO
        ENDDO
      ENDDO
#endif /* ALLOW_DIAGNOSTICS */

      !--------------------------------------------------------
      ! Start calculation
      !--------------------------------------------------------
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr

         
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN

         !---------------------------------------------------
         ! instantaneous equilibrium
         !---------------------------------------------------  
C *** unterminated #ifdef. yanxu doesn't have it PART_EQ, so I'm
C     just removiting it from the code (hma, 30 May 2015) 
C#ifdef PART_EQ             

           !------------------------------------------------
           ! Local variables
           !------------------------------------------------
           poco = max( poc(i,j,k,bi,bj), 0d0 )  ! POC concentration (units?)
           doco = max( doc(i,j,k,bi,bj), 0d0 )  ! DOC concentration (units?)

           !------------------------------------------------
           ! PCB concentrations before partitioning 
           !------------------------------------------------
           pcba_old = max(Ptr(i,j,k,bi,bj,ipcba),0d0)
           pcbp_old = max(Ptr(i,j,k,bi,bj,ipcbp),0d0)
           TOT_PCB = pcba_old + pcbp_old

           IF (TOT_PCB .gt. 0. _d 0) THEN

             !---------------------------------------------------------
             ! Partition coefficient Kow to Koc
             !---------------------------------------------------------
             !   logKow_T = logKow + (dh_OW/2.303 * R) * (1/T0 - 1/T)
             !   logKoc   = 0.88d0(+/-0.07) *logKow + 0.9 (+/-0.47)
             !   Koc      = 10^logKoc * 10^ (KS *[salt])
             ! where dh_OW is the enthalpy of the octanol water phase
             ! change, T0 =298K, T is the local seawater temperature,
             ! (+/-) are the confidence intervals, KS is the Setchenow
             ! constant and [salt] is the seawater salinity in mol/kg
             !
             ! Note: instead of Kow values by Hawker and connell (1988),
             !       as suggested by Sobek Kow values from Schenker et
             !       al. 2005 were chosen for internal model consistency
             !
             !Reference:
             !   Sobek et al. (2004), Particle-water partitioning of
             !   PCBs
             !   in the photic zone: a 25-month study in the open
             !   Baltic Sea, Environ Sci Technol
             !   Schenker et al. (2005).
             !---------------------------------------------------------
             ! temperature adjustment
             !--------------------------------------------------------
             !--------------------------------------------------------
             logKow_T  = logKow +(dh_OW/(2.303d0 *R)) * (1.0d0/TK0 -
     &                   1.0d0/(theta(i,j,k,bi,bj)+C2K))
             logKoc    = (0.88d0+0.07d0) * logKow_T + (0.90d0+0.47d0)

             !---------------------------------------------------------
             ! Salinity adjustment
             !---------------------------------------------------------
             Koc_pcb   = 10**logKoc
             Koc_pcb   = Koc_pcb * 10**(KS *
     &                   (salt(i,j,k,bi,bj)/58.44d0))
                         ! salinity conversion from g/kg to mol kg

             !---------------------------------------------------------
             ! Calculate fraction in the dissolved phase
             ! convert mmolC/m3 -> kgC/L
             !---------------------------------------------------------

              FRAC_Dissolved = 1. _d 0 
     &             /(1. _d 0 + Koc_pcb * poco
     &             * 1.2d-8)   ! mmol C/m3 -> kg C/L     

             !------------------------------------------------
             ! pH of seawater
             !
             ! Note: Currently not considering the real pH
             !       (xz)
             !------------------------------------------------
             pH(i,j,k,bi,bj) = 8. _d 0  
           
             !------------------------------------------------
             ! Particulate fraction (i.e., PCB bound to POC)
             !------------------------------------------------
             FRAC_POCBound = 1.d0 - FRAC_Dissolved

             !------------------------------------------------
             ! New PCB concentrations (mol/m3) in the dissolved
             ! and particle phases after consdering partitioning
             !------------------------------------------------ 
             pcba_new = TOT_PCB * FRAC_Dissolved
             pcbp_new = TOT_PCB * FRAC_POCBound

C             !------------------------------------------------        
C             ! Calculate degradation (photolytic + biodegradation)
C             ! for dissolved-phase PCBs (units of Ptr/sec)
C             ! 
C             ! Units of Kdeg_pcb are 1/hr, so factor of 3600d0
C             ! is to convert to 1/sec.
C             !------------------------------------------------
C             pcba_deg = Kdeg_pcb/3600d0 * pcba_old    

             !------------------------------------------------
             ! Update tracer array
             !------------------------------------------------
             Ptr(i,j,k,bi,bj,ipcba) = pcba_new  
C     &                                - pcba_deg * deltaTclock
             Ptr(i,j,k,bi,bj,ipcbp) = pcbp_new

#ifdef ALLOW_DIAGNOSTICS
             pcbapart_diag(i,j,k) = (pcba_old - pcba_new) / deltaTclock
             pcbppart_diag(i,j,k) = (pcbp_old - pcbp_new) / deltaTclock
#endif /* ALLOW_DIAGNOSTICS */            
           ENDIF


C      No partitioning  
#ifdef PART_NO
#endif /* PART_NO */

C      Slow dissolution or sorption
#ifdef PART_DI
#endif /* PART_DI */
         
         ENDIF !(HFacC(i,j,k,bi,bj).gt.0. _d 0) 
        ENDDO   !k
       ENDDO    !i
      ENDDO     !j

      !--------------------------------------------------------
      ! Write diagnostics
      !--------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS   
        WRITE(diagname,'(A8)') 'PCBaPART'
        CALL DIAGNOSTICS_FILL( pcbapart_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
     
        WRITE(diagname,'(A8)') 'PCBpPART'
        CALL DIAGNOSTICS_FILL( pcbppart_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )   
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_PCB */
#endif /* ALLOW_PTRACERS */

       RETURN
       END
