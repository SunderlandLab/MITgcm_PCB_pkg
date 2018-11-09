C $Header: /u/gcmpack/MITgcm/pkg/hg/hgr_surfforcing.F,v 1.3 2006/08/06 01:29:12 jmc Exp $
C $Name: checkpoint64g $

#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"
#ifdef PCBPREFRACTORY

CStartOfInterFace
      SUBROUTINE PCBpR_SURFFORCING( Ptr, surfpcbpr,
     I           bi,bj,iMin,iMax,jMin,jMax,
     I           myIter, myTime, myThid )
      !===========================================================
      ! SUBROUTINE PCBpR_SURFFORCING
      !   External inputs of refractory PCBp from rivers. 
      !    
      ! REVISION HISTORY
      !   04 Feb 2015 - H. Amos - Copied from Yanxu Zhang's Hg
      !                           code and adapted for PCBp.                               
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
#include "PCB_SIZE.h"
#include "PCB.h"

      !-----------------------------------------------------------
      ! Routine arguments
      !-----------------------------------------------------------
      INTEGER myIter, myThid
      _RL myTime
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  surfpcbpr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nrmax)
      
#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB
      _RL  Fluxpcbpr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nrmax)
      _RL  pcbprrivero(nrmax)
#ifdef ALLOW_DIAGNOSTICS
      _RL  pcbprriver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagname

      !-----------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------
      INTEGER I,J,nriver  ! I, J - loop counters

      !-----------------------------------------------------------
      ! Initialize diagnostics
      !-----------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            pcbprriver_diag(i,j)=0d0
        ENDDO
      ENDDO
#endif
      
      DO j=jMin,jMax
        DO i=iMin,iMax
          IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN
            !-----------------------------------------------------
            ! Localize variables
            !-----------------------------------------------------

            ! Inputs from rivers
            ! Convert kg/s -> mol/m2/s
            DO nriver=1,nrmax
             pcbprrivero(nriver) = pcbprriver(i,j,bi,bj,nriver)
     &           *recip_rA(i,j,bi,bj)*1000.d0/Mb
                           
             Fluxpcbpr(i,j,nriver)=pcbprrivero(nriver)
#ifdef ALLOW_DIAGNOSTICS             
             pcbprriver_diag(i,j)=pcbprriver_diag(i,j)+
     &           Fluxpcbpr(i,j,nriver)
#endif             
            ENDDO
          ENDIF
        ENDDO
      ENDDO

        !---------------------------------------------------------
        ! update surface tendencies
        !---------------------------------------------------------
        DO j=jMin,jMax
          DO i=iMin,iMax
           DO nriver=1,nrmax
            SURFpcbpr(i,j,nriver)=
     &         HFacC(i,j,1,bi,bj)*Fluxpcbpr(i,j,nriver)*recip_drF(1)
           ENDDO
          ENDDO
         ENDDO
         
        !---------------------------------------------------------
        ! write diagnostics
        !---------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
        WRITE(diagname,'(A8)') 'PCBPRRIV'
        CALL DIAGNOSTICS_FILL( pcbprriver_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )         
#endif
         
#endif
#endif
#endif
        RETURN
        END
#endif

