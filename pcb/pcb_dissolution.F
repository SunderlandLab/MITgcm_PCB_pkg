C $Header: /u/gcmpack/MITgcm/pkg/hg/hg_dissolution.F,v 1.10 2010/03/16 00:14:47 jmc Exp $
C $Name: checkpoint64g $


#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

CStartOfInterFace
      SUBROUTINE PCB_DISSOLUTION( Ptr,
     &                         Gchem,    
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)

      !===========================================================
      ! SUBROUTINE PCB_DISSOLUTION                                
      !   Calculate the dissolution of refractory PCB to PCBa
      !    
      ! REVISION HISTORY
      !   04 Feb 2015 - H. Amos - Copied from Yanxu Zhang's Hg
      !                           code and adapted for PCBs.                                               
      !===========================================================
      IMPLICIT NONE

      !-----------------------------------------------------------
      ! Global variables 
      !-----------------------------------------------------------
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

      !-----------------------------------------------------------
      ! Routine arguments
      !-----------------------------------------------------------
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)    
#ifdef ALLOW_DIAGNOSTICS      
#ifdef PCBAREFRACTORY       
      _RL  pcbardiss_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif 
#ifdef PCBPREFRACTORY       
      _RL  pcbprdiss_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif                  
      CHARACTER*8 diagname
#endif
      INTEGER i,j,k,l,bi, bj, imin, imax, jmin, jmax
      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB

      !-----------------------------------------------------------
      ! Gross tendencies
      !-----------------------------------------------------------
#ifdef PCBAREFRACTORY
      _RL pcbar_pcba
#endif 
#ifdef PCBPREFRACTORY
      _RL pcbpr_pcba
#endif      

      !-----------------------------------------------------------
      ! Parameters
      !-----------------------------------------------------------
#ifdef PCBAREFRACTORY       
      _RL k_pcbar_pcba
      data k_pcbar_pcba / 3.85 _d -8 /  !10%/ (30 day)
#endif
#ifdef PCBPREFRACTORY       
      _RL k_pcbpr_pcba
      data k_pcbpr_pcba / 3.85 _d -8 /  !10%/ (30 day)
      INTEGER nriver 
#endif            

      !-----------------------------------------------------------      
      ! Start calculation
      !-----------------------------------------------------------

      ! Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax  
#ifdef PCBAREFRACTORY     
            pcbardiss_diag(i,j,k)=0d0
#endif
#ifdef PCBPREFRACTORY     
            pcbprdiss_diag(i,j,k)=0d0
#endif
          ENDDO
        ENDDO
      ENDDO
#endif

#if defined(PCBAREFRACTORY) & defined(PCBPREFRACTORY) & defined(PCBPR_SLOWRELEASING)
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN

            !-----------------------------------------------------
            ! Tendency calcualtion                      
            !-----------------------------------------------------
#ifdef PCBAREFRACTORY
            pcbar_pcba = 0d0
            pcbar_pcba = pcbar_pcba + k_pcbar_pcba
     &           * max( Ptr(i,j,k,bi,bj,ipcbar), 0d0)
#endif
#ifdef PCBPREFRACTORY
            pcbpr_pcba = 0d0
            DO nriver=1,nrmax
              pcbpr_pcba = pcbpr_pcba + k_pcbpr_pcba 
     &           * max( Ptr(i,j,k,bi,bj,ipcbpr+nriver-1), 0d0)
            ENDDO
#endif

            !-----------------------------------------------------           
            ! Calculate diagnostics                                           
            !-----------------------------------------------------       
#ifdef ALLOW_DIAGNOSTICS
#ifdef PCBAREFRACTORY
            pcbardiss_diag(i,j,k) = pcbar_pcba
#endif 
#ifdef PCBPREFRACTORY
            pcbprdiss_diag(i,j,k) = pcbpr_pcba
#endif                            
#endif

           !------------------------------------------------------
           ! Accumulate tendency
           !------------------------------------------------------
#ifdef PCBAREFRACTORY
           Gchem(i,j,k,bi,bj,ipcba)=Gchem(i,j,k,bi,bj,ipcba)+pcbar_pcba
           Gchem(i,j,k,bi,bj,ipcbar)=
     &            Gchem(i,j,k,bi,bj,ipcbar) - k_pcbar_pcba *
     &            max(Ptr(i,j,k,bi,bj,ipcbar), 0d0)
#endif 
#ifdef PCBPREFRACTORY
           Gchem(i,j,k,bi,bj,ipcba)=Gchem(i,j,k,bi,bj,ipcba)+pcbpr_pcba
           DO nriver=1,nrmax
             Gchem(i,j,k,bi,bj,ipcbpr+nriver-1)=
     &            Gchem(i,j,k,bi,bj,ipcbpr+nriver-1) - k_pcbpr_pcba * 
     &            max(Ptr(i,j,k,bi,bj,ipcbpr+nriver-1), 0d0)
           ENDDO
#endif                        
         ENDIF
        ENDDO
       ENDDO
      ENDDO
#endif
      
        !---------------------------------------------------------
        ! write diagnostics
        !---------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
#ifdef PCBAREFRACTORY
        WRITE(diagname,'(A8)') 'PCBARDIS'
        CALL DIAGNOSTICS_FILL( pcbardiss_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif
#ifdef PCBPREFRACTORY
        WRITE(diagname,'(A8)') 'PCBPRDIS'
        CALL DIAGNOSTICS_FILL( pcbprdiss_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif       
#endif

#endif
#endif

       RETURN
       END
