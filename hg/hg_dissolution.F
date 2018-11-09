C $Header: /u/gcmpack/MITgcm/pkg/hg/hg_dissolution.F,v 1.10 2010/03/16 00:14:47 jmc Exp $
C $Name: checkpoint64g $

C To calculate the chemistry of Hg species in the ocean, yxzhang, 5/30/13
C

#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

CStartOfInterFace
      SUBROUTINE HG_DISSOLUTION( Ptr,
     &                         Gchem,    
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)

C     /==========================================================\
C     | SUBROUTINE HG_DISSOLUTION                                |
C     | o Calculate the dissolution of HgR to Hg2                |
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

C     == Routine arguments ==
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)    
#ifdef ALLOW_DIAGNOSTICS      
#ifdef HG2REFRACTORY       
      _RL  hg2rdiss_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif 
#ifdef HGPREFRACTORY       
      _RL  hgprdiss_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif                  
      CHARACTER*8 diagname
#endif
      INTEGER i,j,k,l,bi, bj, imin, imax, jmin, jmax
      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C     -- Gross tendencies
#ifdef HG2REFRACTORY
      _RL hg2r_hg2
#endif 
#ifdef HGPREFRACTORY
      _RL hgpr_hg2
#endif      

C     == Parameters
#ifdef HG2REFRACTORY       
      _RL k_hg2r_hg2
      data k_hg2r_hg2 / 3.85 _d -8 /  !10%/ (30 day)
C      data k_hg2r_hg2 / 0.0 _d 0 /
#endif
#ifdef HGPREFRACTORY       
      _RL k_hgpr_hg2
      data k_hgpr_hg2 / 3.85 _d -8 /  !10%/ (30 day)
      INTEGER nriver 
#endif            
      
C     == Start calculation
C Initialize diagnostics
#ifdef ALLOW_DIAGNOSTICS
      DO k=1,Nr
        DO j=jMin,jMax
          DO i=iMin,iMax  
#ifdef HG2REFRACTORY     
            hg2rdiss_diag(i,j,k)=0d0
#endif
#ifdef HGPREFRACTORY     
            hgprdiss_diag(i,j,k)=0d0
#endif
          ENDDO
        ENDDO
      ENDDO
#endif

#if defined(HG2REFRACTORY) & defined(HGPREFRACTORY) & defined(HGPR_SLOWRELEASING)
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
C     == Tendency calcualtion                      
#ifdef HG2REFRACTORY
            hg2r_hg2 = 0d0
            hg2r_hg2 = hg2r_hg2 + k_hg2r_hg2
     &           * max( Ptr(i,j,k,bi,bj,ihg2r), 0d0)
#endif
#ifdef HGPREFRACTORY
            hgpr_hg2 = 0d0
            DO nriver=1,nrmax
              hgpr_hg2 = hgpr_hg2 + k_hgpr_hg2 
     &           * max( Ptr(i,j,k,bi,bj,ihgpr+nriver-1), 0d0)
            ENDDO
#endif
           
C Calculate diagnistics                                                  
#ifdef ALLOW_DIAGNOSTICS
#ifdef HG2REFRACTORY
            hg2rdiss_diag(i,j,k) = hg2r_hg2
#endif 
#ifdef HGPREFRACTORY
            hgprdiss_diag(i,j,k) = hgpr_hg2
#endif                            
#endif

C     == Accumulate tendency
#ifdef HG2REFRACTORY
           Gchem(i,j,k,bi,bj,ihg2)=Gchem(i,j,k,bi,bj,ihg2)+hg2r_hg2
           Gchem(i,j,k,bi,bj,ihg2r)=
     &            Gchem(i,j,k,bi,bj,ihg2r) - k_hg2r_hg2 *
     &            max(Ptr(i,j,k,bi,bj,ihg2r), 0d0)
#endif 
#ifdef HGPREFRACTORY
           Gchem(i,j,k,bi,bj,ihg2)=Gchem(i,j,k,bi,bj,ihg2)+hgpr_hg2
           DO nriver=1,nrmax
             Gchem(i,j,k,bi,bj,ihgpr+nriver-1)=
     &            Gchem(i,j,k,bi,bj,ihgpr+nriver-1) - k_hgpr_hg2 * 
     &            max(Ptr(i,j,k,bi,bj,ihgpr+nriver-1), 0d0)
           ENDDO
#endif                        
         ENDIF
        ENDDO
       ENDDO
      ENDDO
#endif
      
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
#ifdef HG2REFRACTORY
        WRITE(diagname,'(A8)') 'HG2RDISS'
        CALL DIAGNOSTICS_FILL( hg2rdiss_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif
#ifdef HGPREFRACTORY
        WRITE(diagname,'(A8)') 'HGPRDISS'
        CALL DIAGNOSTICS_FILL( hgprdiss_diag(1-Olx,1-Oly,1), diagname,
     &                         0,Nr,2,bi,bj,myThid )
#endif       
#endif

#endif
#endif

       RETURN
       END
