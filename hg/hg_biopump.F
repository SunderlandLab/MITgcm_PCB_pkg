C $Header: /u/gcmpack/MITgcm/pkg/hg/hg_biopump.F,v 1.10 2010/03/16 00:14:47 jmc Exp $
C $Name: checkpoint64g $

C to calculate the sinking of Hg species along POM
C also include the methylation and demethylation
C

#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

CStartOfInterFace
      SUBROUTINE HG_BIOPUMP(  Ptr,Gchem,      
     &                        bi,bj,imin,imax,jmin,jmax,
     &                        myIter,myTime,myThid)

C     /==========================================================\
C     | SUBROUTINE HG_BIOPUMP                                   |
C     | o Calculate the changes to Hgp and MMHgp                |
C     |   through particulate sinking                           |
C     |==========================================================|
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "HG_SIZE.h"
#include "HG.h"

C     == Routine arguments ==
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num) 
#ifdef ALLOW_DIAGNOSTICS      
      _RL  hgpsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      CHARACTER*8 diagname
#endif      
      _RL  Fout_Hgp(Nr)
#ifdef HGPREFRACTORY
      INTEGER nriver
      _RL  whgr(nrmax)
C     data whgr /100.0 _d 0, 10.0 _d 0, 1.0 _d 0/
C     data whgr /3.0 _d 0, 3.0 _d 0, 3.0 _d 0/        
      data whgr /3.0 _d 0 /
C     data whgr /30.0 _d 0 /
#if nrmax > 3
      What is the sinking speed of riverine particles?
#endif             
      _RL  Fout_HgPR(Nr,nrmax)
#ifdef ALLOW_DIAGNOSTICS      
      _RL  hgprsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif      
#endif        
#ifdef METHY_MODE       
      _RL  Fout_MMHgp(Nr)
#ifdef ALLOW_DIAGNOSTICS      
      _RL  mmhgpsink_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif      
#endif      
      
      _RL  poco, fpoco
      _RL  tiny
      data tiny /1d-32/
      INTEGER i,j,k,bi, bj, imin, imax, jmin, jmax 

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax
            hgpsink_diag(i,j,k)=0d0
#ifdef HGPREFRACTORY 
            hgprsink_diag(i,j,k)=0d0
#endif                        
#ifdef METHY_MODE 
            mmhgpsink_diag(i,j,k)=0d0
#endif            
          ENDDO
        ENDDO
      ENDDO
#endif      
      
      DO j=jMin,jMax
       DO i=iMin,iMax

C       Zero flux
        DO k=1,Nr
         Fout_Hgp(k) = 0. _d 0
#ifdef HGPREFRACTORY
         DO nriver=1,nrmax 
          Fout_Hgpr(k,nriver) = 0. _d 0
         ENDDO
#endif                     
#ifdef METHY_MODE 
         Fout_MMHgp(k) = 0. _d 0
#endif            
        ENDDO

C       Flux calculation
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           ! localize variables
           poco =  max(poc(i,j,k,bi,bj),0d0)
           fpoco = max(Fpoc(i,j,k,bi,bj), 0d0)

           IF (poco .gt. tiny .and. fpoco .gt. tiny) THEN   
             Fout_Hgp(k) = max(PTR(i,j,k,bi,bj,ihgp),0d0)
     &          * fpoco / poco
     &          / 8.64 d 4
           ELSE
             Fout_Hgp(k) = 0. _d 0
           ENDIF 
     
C Calculate diagnistics
#ifdef ALLOW_DIAGNOSTICS
           hgpsink_diag(i,j,k) = Fout_Hgp(k)  !mol/m^2/s   
#endif

#ifdef HGPREFRACTORY
           DO nriver=1,nrmax
            Fout_Hgpr(k,nriver) = max(PTR(i,j,k,bi,bj,
     &            ihgpr+nriver-1),0d0)
#ifdef HGPR_IRREDUCIBLE 
     &          * fpoco / poco
#elif defined(HGPR_SLOWRELEASING)
     &          * whgr(nriver)
#endif
     &          / 8.64 d 4

#ifdef WALSH_FILTER
#ifdef walsh_depth_0
C           do nothing
            IF ( 1 .eq. 0 ) THEN
#endif                
#ifdef walsh_depth_1
c           the particle sinking at river mouth is handled by the Walsh Typology filter, 4/18/2014
            IF ( hgprriver(i,j,bi,bj,nriver) .gt. tiny ) THEN
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
     &           hgprriver(i,j,bi,bj,nriver) .gt. tiny) THEN
#endif                  
             Fout_Hgpr(k,nriver) = 0. _d 0
            ENDIF
#endif
           ENDDO
#ifdef ALLOW_DIAGNOSTICS
           DO nriver=1,nrmax
            hgprsink_diag(i,j,k) = hgprsink_diag(i,j,k)
     &           +Fout_Hgpr(k,nriver)  !mol/m^2/s
           ENDDO 
#endif
#endif                
                                                
#ifdef METHY_MODE
           IF (poco .gt. tiny .and. fpoco .gt. tiny) THEN
#ifdef FOODW_MODE
c    need to figure it out, 8/21/2013
           Fout_MMHgp(k) = max(PTR(i,j,k,bi,bj,immhgp),0d0)
     &          * fpoco / ( poco + tiny )
     &          / 8.64 d 4
#else
           Fout_MMHgp(k) = max(PTR(i,j,k,bi,bj,immhgp),0d0)
     &          * fpoco / ( poco + tiny )
     &          / 8.64 d 4
#endif
           ELSE
             Fout_MMHgp(k) = 0. _d 0
           ENDIF
#ifdef ALLOW_DIAGNOSTICS
           mmhgpsink_diag(i,j,k) = Fout_MMHgp(k)  !mol/m^2/s   
#endif
#endif
         ENDIF
        ENDDO

C----&---1---------2---------3---------4---------5---------6---------7--    
C       Tendency calculation
        IF (HFacC(i,j,1,bi,bj).gt.0. _d 0) THEN
         Gchem(i,j,1,bi,bj,ihgp)=Gchem(i,j,1,bi,bj,ihgp) 
     &                          - Fout_Hgp(1) * recip_drF(1)
#ifdef HGPREFRACTORY
         DO nriver=1,nrmax
          Gchem(i,j,1,bi,bj,ihgpr+nriver-1)=
     &     Gchem(i,j,1,bi,bj,ihgpr+nriver-1) 
     &     - Fout_Hgpr(1,nriver) * recip_drF(1)
         ENDDO
#endif         
#ifdef METHY_MODE         
         Gchem(i,j,1,bi,bj,immhgp)=Gchem(i,j,1,bi,bj,immhgp)
     &     - Fout_MMHgp(1) * recip_drF(1)
#endif         
        ENDIF
  
        DO k=2,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
          Gchem(i,j,k,bi,bj,ihgp)=Gchem(i,j,k,bi,bj,ihgp) + 
     &             (Fout_Hgp(k-1) - Fout_Hgp(k)) * recip_drF(k)
#ifdef HGPREFRACTORY
          DO nriver=1,nrmax     
           Gchem(i,j,k,bi,bj,ihgpr+nriver-1)=
     &        Gchem(i,j,k,bi,bj,ihgpr+nriver-1) + 
     &        (Fout_Hgpr(k-1,nriver) - Fout_Hgpr(k,nriver)) 
     &        * recip_drF(k)
          ENDDO
#endif          
#ifdef METHY_MODE     
          Gchem(i,j,k,bi,bj,immhgp)=Gchem(i,j,k,bi,bj,immhgp) + 
     &             (Fout_MMHgp(k-1) - Fout_MMHgp(k)) * recip_drF(k)
#endif     
         ENDIF
        ENDDO

       ENDDO
      ENDDO

c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'HGPSINK '
        CALL DIAGNOSTICS_FILL( hgpsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#ifdef HGPREFRACTORY
        WRITE(diagname,'(A8)') 'HGPRSINK'
        CALL DIAGNOSTICS_FILL( hgprsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif     
#ifdef METHY_MODE
        WRITE(diagname,'(A8)') 'MHGPSINK'
        CALL DIAGNOSTICS_FILL( mmhgpsink_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif       
#endif

#endif
#endif

       RETURN
       END
