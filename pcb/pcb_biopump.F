C $Header: /u/gcmpack/MITgcm/pkg/pcb/pcb_biopump.F 
C $Name: checkpoint64g $

#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

CStartOfInterFace
      SUBROUTINE PCB_BIOPUMP(  Ptr,Gchem,      
     &                        bi,bj,imin,imax,jmin,jmax,
     &                        myIter,myTime,myThid)

      !==========================================================
      ! SUBROUTINE PCB_BIOPUMP                                   
      !   Calculate the changes to PCBP through particulate 
      !   sinking                           
      !
      ! REVISION HISTORY
      !   03 Feb 2015 - H. Amos - Update behavior of PCB associated
      !                           with riverine particles based on
      !                           Hg code of Yanxu Zhang.    
      !==========================================================
      IMPLICIT NONE

      !---------------------------------------------------------
      ! Global variables
      !---------------------------------------------------------
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PCB_SIZE.h"
#include "PCB.h"

      !---------------------------------------------------------
      ! Routine arguments 
      !---------------------------------------------------------
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)    
#ifdef ALLOW_DIAGNOSTICS      
      _RL  pcbpsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      CHARACTER*8 diagname
#endif      
      _RL  Fout_PCBp(Nr)

#ifdef PCBPREFRACTORY
      INTEGER nriver
      _RL  whgr(nrmax)
C     data whgr /100.0 _d 0, 10.0 _d 0, 1.0 _d 0/
C     data whgr /3.0 _d 0, 3.0 _d 0, 3.0 _d 0/        
      data whgr /3.0 _d 0 /
C     data whgr /30.0 _d 0 /
#if nrmax > 3
      What is the sinking speed of riverine particles?
#endif             
      _RL  Fout_PCBpR(Nr,nrmax)
#ifdef ALLOW_DIAGNOSTICS      
      _RL  pcbprsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif      
#endif        
      _RL  poco, fpoco
      _RL  tiny
      _RL  SecPerDay
      data tiny /1d-32/
      data SecPerDay /8.64 d 4/  ! seconds per day
      INTEGER i,j,k,bi, bj, imin, imax, jmin, jmax      

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
      !---------------------------------------------------------
      ! Local variables 
      !---------------------------------------------------------

      !---------------------------------------------------------
      ! Initialize diagnostics
      !---------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            pcbpsink_diag(i,j,k)=0d0
c #ifdef HGPRIVER 
c            hgrsink_diag(i,j,k)=0d0
c #endif                        
c #ifdef METHY_MODE 
c            mmhgpsink_diag(i,j,k)=0d0
c #endif            
          ENDDO
        ENDDO
      ENDDO
#endif      
      
      DO j=jMin,jMax
       DO i=iMin,iMax

        !-------------------------------------------------------
        ! Zero flux
        !-------------------------------------------------------
        DO k=1,Nr
         Fout_PCBp(k) = 0. _d 0
#ifdef PCBPREFRACTORY
         DO nriver=1,nrmax 
          Fout_PCBpr(k,nriver) = 0. _d 0
         ENDDO
#endif                     
        ENDDO


        !-------------------------------------------------------
        ! Flux calculation
        !-------------------------------------------------------
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           ! localize variables
           poco  =  max(poc(i,j,k,bi,bj),0d0)
           fpoco = max(Fpoc(i,j,k,bi,bj), 0d0)

           IF (poco .gt. tiny .and. fpoco .gt. tiny) THEN   
             Fout_PCBp(k) = max(PTR(i,j,k,bi,bj,ipcbp),0d0)
     &          * fpoco / poco
     &          / SecPerDay
           ELSE
             Fout_PCBp(k) = 0. _d 0
           ENDIF 
     
           !----------------------------------------------------
           ! Calculate diagnostics
           !----------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
           pcbpsink_diag(i,j,k) = Fout_PCBp(k)  !mol/m^2/s   
#endif

#ifdef PCBPREFRACTORY
           DO nriver=1,nrmax
            Fout_PCBpr(k,nriver) = max(PTR(i,j,k,bi,bj,
     &            ipcbpr+nriver-1),0d0)
#ifdef HGPR_IRREDUCIBLE 
     &          * fpoco / poco
#elif defined(HGPR_SLOWRELEASING)
     &          * whgr(nriver)
#endif
     &          / SecPerDay

#ifdef WALSH_FILTER
#ifdef walsh_depth_0
C           do nothing
            IF ( 1 .eq. 0 ) THEN
#endif                
#ifdef walsh_depth_1
c           the particle sinking at river mouth is handled by the Walsh Typology filter, 4/18/2014
            IF ( pcbprriver(i,j,bi,bj,nriver) .gt. tiny ) THEN
#endif                
#ifdef walsh_depth_2
c           the particle sinking over region shallower than 55m depth is handled by the Walsh Typology filter, 4/19/2014
            IF ( HFacC(i,j,5,bi,bj) .lt. tiny ) THEN
#endif                
#ifdef walsh_depth_3
c           the particle sinking over region shallower than 185m depth is handled by the Walsh Typology filter, 4/19/2014
            IF ( HFacC(i,j,9,bi,bj) .lt. tiny ) THEN
#endif
#ifdef walsh_depth_4
c           the particle sinking over region shallower than 55m depth, or
C           at the river mouth is handled by the Walsh Typology filter, 4/19/2014
            IF ( HFacC(i,j,5,bi,bj) .lt. tiny .or. 
     &           pcbprriver(i,j,bi,bj,nriver) .gt. tiny) THEN
#endif                  
             Fout_PCBpr(k,nriver) = 0. _d 0
            ENDIF
#endif
           ENDDO
#ifdef ALLOW_DIAGNOSTICS
           DO nriver=1,nrmax
            pcbprsink_diag(i,j,k) = pcbprsink_diag(i,j,k)
     &           +Fout_PCBpr(k,nriver)  !mol/m^2/s
           ENDDO 
#endif
#endif                
                                                
         ENDIF
        ENDDO

        !--------------------------------------------------------
        ! Tendency calculation
        !--------------------------------------------------------
        IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN
         Gchem(i,j,1,bi,bj,ipcbp)=Gchem(i,j,1,bi,bj,ipcbp) 
     &                          - Fout_PCBp(1) * recip_drF(1)

#ifdef PCBPREFRACTORY
         DO nriver=1,nrmax
          Gchem(i,j,1,bi,bj,ipcbpr+nriver-1)=
     &     Gchem(i,j,1,bi,bj,ipcbpr+nriver-1) 
     &     - Fout_pcbpr(1,nriver) * recip_drF(1)
         ENDDO
#endif         

        ENDIF
  
        DO k=2,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
          Gchem(i,j,k,bi,bj,ipcbp)=Gchem(i,j,k,bi,bj,ipcbp) + 
     &             (Fout_PCBp(k-1) - Fout_PCBp(k)) * recip_drF(k)

#ifdef PCBPREFRACTORY
          DO nriver=1,nrmax     
           Gchem(i,j,k,bi,bj,ipcbpr+nriver-1)=
     &        Gchem(i,j,k,bi,bj,ipcbpr+nriver-1) + 
     &        (Fout_pcbpr(k-1,nriver) - Fout_pcbpr(k,nriver)) 
     &        * recip_drF(k)
          ENDDO
#endif          

         ENDIF
        ENDDO

       ENDDO
      ENDDO

      !----------------------------------------------------------
      ! Write diagnostics
      !----------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'PCBpSINK'
        CALL DIAGNOSTICS_FILL( PCBpsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )

#ifdef PCBPREFRACTORY
        WRITE(diagname,'(A8)') 'PCBPRSNK'
        CALL DIAGNOSTICS_FILL( pcbprsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif     

#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_PCB */
#endif /* ALLOW_PTRACERS */

       RETURN
       END
