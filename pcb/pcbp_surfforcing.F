C $Header: /u/gcmpack/MITgcm/pkg/hg/hg2_surfforcing.F,v 1.3 2006/08/06 01:29:12 jmc Exp $
C $Name: checkpoint64g $

#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

CStartOfInterFace
      SUBROUTINE pcbp_SURFFORCING( PTR, surfpcbp,
     &           bi,bj,iMin,iMax,jMin,jMax,
     &           myIter, myTime, myThid )

      !==========================================================
      ! SUBROUTINE PCBp_SURFFORCING                               
      !   External inputs of particulate PCB from atmospheric 
      !   deposition.
      !
      ! REVISION HISTORY
      !   ?? Jul 2014 - H. Amos - PFOS code inherited from Xianming 
      !                           Zhang and adopted for PCBs.
      !   07 Aug 2014 - H. Amos - add inputs from atmospheric 
      !                           deposition 
      !   03 Feb 2015 - H. Amos - Update river inputs to mimic 
      !                           Yanxu Zhang's Hg river code.                             
      !==========================================================
      IMPLICIT NONE

      !----------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PCB_SIZE.h"
#include "PCB.h"

      !----------------------------------------------------------
      ! Routine arguments 
      !----------------------------------------------------------
      INTEGER myIter, myThid
      _RL myTime
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  surfpcbp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  Fluxpcbp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  pcbprivero  ! riverine discharge
      _RL  pcbpdepoo   ! atmospheric deposition
#ifdef ALLOW_DIAGNOSTICS
       ! Diagnostics      
      _RL  pcbpriver_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy) ! riverine discharge
      _RL  pcbpdep_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)   ! atmospheric deposition
      CHARACTER*8 diagname
#endif      
      INTEGER iMin,iMax,jMin,jMax, bi, bj

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB

      !----------------------------------------------------------
      ! Local variables
      !----------------------------------------------------------

      ! I, J - Loop counters
      INTEGER I,J

      !----------------------------------------------------------
      ! Initialize diagnostics
      !----------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS
      DO j=jMin,jMax
        DO i=iMin,iMax
            ! inputs from rivers
            pcbpriver_diag(i,j) = 0d0

            ! inputs from atmospheric deposition
            pcbpdep_diag(i,j)   = 0d0
        ENDDO
      ENDDO
#endif
      
      DO j=jMin,jMax
        DO i=iMin,iMax

          !------------------------------------------------------
          ! Are you in an ocean box?
          !   Ocean box when HFacC = 1
          !   Land  box when HFacC = 0
          !------------------------------------------------------
          IF (HFacC(i,j,1,bi,bj).GT.0.5) THEN

             !---------------------------------------------------
             ! Local variables
             !---------------------------------------------------

             ! Inputs from rivers
             ! kg/s -> mol/m2/s
             pcbprivero = PCBpriver(i,j,bi,bj)
     &             *recip_rA(i,j,bi,bj)*1000.d0/Mb
                           
             ! Inputs from  atmospheric deposition 
C Units of input file are already in mol/m2/sec
C             ! kg/m2/s -> mol/m2/s
C             pcbpdepoo = PCBpDEP(i,j,bi,bj)*1000.d0/Mb
             pcbpdepoo = PCBpDEP(i,j,bi,bj)

             ! Total particulate PCB inputs (mol/m2/s)
             Fluxpcbp(i,j) = pcbpdepoo + pcbprivero

#ifdef ALLOW_DIAGNOSTICS
             ! Inputs from rivers (mol/m2/s)             
             pcbpriver_diag(i,j) = pcbprivero
 
             ! Inputs from  atmospheric deposition (mol/m2/s)
             pcbpdep_diag(i,j)   = pcbpdepoo
#endif             
          ENDIF
        ENDDO
      ENDDO

        !--------------------------------------------------------
        ! Update surface tendencies
        !--------------------------------------------------------
        DO j=jMin,jMax
          DO i=iMin,iMax
           SURFPCBP(i,j)=
     &         HFacC(i,j,1,bi,bj)*Fluxpcbp(i,j)*recip_drF(1)
          ENDDO
         ENDDO
 
        !--------------------------------------------------------        
        ! Write diagnostics
        !--------------------------------------------------------
#ifdef ALLOW_DIAGNOSTICS

        ! Inputs from rivers (mol/m2/s)
        WRITE(diagname,'(A8)') 'PCBpRIVR'
        CALL DIAGNOSTICS_FILL( pcbpriver_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid )        

        ! Atmospheric deposition (mol/m2/s)
        WRITE(diagname,'(A8)') 'PCBpDEPO'
        CALL DIAGNOSTICS_FILL( pcbpdep_diag(1-Olx,1-Oly), diagname,
     &                         0,1,2,bi,bj,myThid ) 
#endif
         
#endif
#endif

        RETURN
        END
