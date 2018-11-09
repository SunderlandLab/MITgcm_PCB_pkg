C $Header: /u/gcmpack/MITgcm/pkg/hg/hg_foodweb.F,v 1.10 2010/03/16 00:14:47 jmc Exp $
C $Name: checkpoint64g $

C To calculate the chemistry of Hg species in the ocean, yxzhang, 5/30/13

#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

CStartOfInterFace
      SUBROUTINE HG_FOODWEB(   Ptr,
     &                         Gchem,
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)

C     /==========================================================\
C     | SUBROUTINE HG_FOODWEB                                    |
C     | o Calculate the changes to Hg species because            |
C     |   of foodweb dynamics                                    |
C     |   based on Asif Quesi 2013 paper                         |
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
#include "HG_SIZE.h"
#include "HG.h"

C     == Routine arguments ==
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL  Gchem(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      INTEGER bi, bj, imin, imax, jmin, jmax
#ifdef FOODW_MODE
C     == Local variables ==            
      INTEGER i,j,k
      INTEGER np,nz
      _RL total_mmhg
      _RL vcfo, vcfppo(npmax)  !volume concentration factor, m3/m3 
      _RL fphyo(npmax) !fraction of MMhg on phytoplankton
      _RL total_fphyo, phytopo(npmax), zoopo(nzmax), doco
      _RL gro, cppo, cwo    !temporary variable
      _RL phygrowo(npmax),zoograzeo(npmax,nzmax),zoomorto(nzmax)
      _RL grazeo,deatho, phymmhgo(npmax), zoommhgo(nzmax)
      _RL active(npmax)     !active species
      _RL tiny
      PARAMETER(tiny=1d-32)
      
C     == Phytoplankton variables
      _RL svr(npmax)  !surface area to volume ratios, um-1
      _RL fcarbon     !carbon content of wet phytoplankton
      _RL ae(nzmax)   !assimilation efficiency     
      DATA svr /7.5 _d -2  ,3 _d -1, 3 _d 0/   ! um-1, 3/r, three species with radius of 40 um, 10 um and 1 um, Queshi 2013
      DATA fcarbon /1 _d -1/
      DATA ae / 6 _d -1, 6 _d -1 / ! from Asif codes

#ifdef ALLOW_DIAGNOSTICS      
      CHARACTER*8   diagName
      _RL  mmhgbc_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,npmax)  ! MMHG bioaccumulation from sea water
      _RL  mmhggz_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &                  Nr,npmax,nzmax)  ! MMHG graze
      _RL  mmhgmt_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nzmax)   ! release of MMHg from zooplankton after death
      _RL  mmhgfp_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)         ! release of MMHg via floopy feeding
#endif

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C     == Partition MMHg onto phytoplankton
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           ! Localize variables
           doco = max(doc(i,j,k,bi,bj),0. _d 0)  !!mmol C/m3 -> umol C/L (uM)
           vcfo = 1.26 _d 4 * dexp(-5 _d -3 * doco)
           DO np=1,npmax
             vcfppo(np)=vcfo*svr(np)/0.94
           ENDDO
           cwo = max(Ptr(i,j,k,bi,bj,immhg),0 _d 0) ! MMHg seawater concentration, mol/m3
c
           DO np=1,npmax
               phymmhgo(np) = max(Ptr(i,j,k,bi,bj,iPhy+np-1),0. _d 0) !mol/m3
               phytopo(np) = max(phytoP(i,j,k,bi,bj,np),0. _d 0)    !mmol C/m3
               phygrowo(np) = max(phygrow(i,j,k,bi,bj,np),0. _d 0)  !mmol C/m3/s
               IF (phytopo(np) .gt. tiny .and. 
     &             phygrowo(np) .gt. tiny) THEN
                   active(np) = 1.0  ! active means update
               ELSE
                   active(np) = 0.0  ! non-active means keeps c not change
               ENDIF
           ENDDO
             
           ! total mmhg concentration
           total_mmhg = cwo
           DO np=1,npmax
             IF ( active(np) .gt. 0.5 ) THEN
               total_mmhg = total_mmhg + phymmhgo(np)
             ENDIF
           ENDDO
           
           ! calculate partitioning
           total_fphyo=1. _d 0
           DO np=1,npmax
             IF ( active(np) .gt. 0.5 ) THEN
               gro = phygrowo(np) / phytopo(np) * 8.6400 _d 4  !specific growth rate, s-1 -> d-1
               cppo = 1.0 _d 6 * vcfppo(np) / gro !* cwo       !ng g ww-1, note: cwo is canceled in the next equation
               fphyo(np) = cppo * 1.0 _d -9 / 2.006 _d 2      ! mol g ww-1
     &              * phytopo(np) * 1d-3 * 1.2 _d 1 / fcarbon  ! g ww /m3            
!     &              / cwo                                   ! mol/m3, note: cwo is canceled w/ the above equation
               total_fphyo = total_fphyo + fphyo(np)  
             ENDIF               
           ENDDO

#ifdef ALLOW_DIAGNOSTICS
           DO np = 1,npmax
            IF ( active(np) .gt. 0.5 ) THEN
             mmhgbc_diag(i,j,k,np) =  
     &           ( total_mmhg/(total_fphyo+tiny) 
     &           - Ptr(i,j,k,bi,bj,immhg) )
     &           / deltaTclock
            ELSE
             mmhgbc_diag(i,j,k,np) = 0. _d 0
            ENDIF
           ENDDO
#endif
           Ptr(i,j,k,bi,bj,immhg)=total_mmhg/(total_fphyo+tiny)
           DO np=1,npmax
            IF ( active(np) .gt. 0.5 ) THEN
               Ptr(i,j,k,bi,bj,iPhy+np-1)=total_mmhg*fphyo(np)
     &                  / (total_fphyo+tiny)
            ENDIF        
           ENDDO
C
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      
C     == Biomagnification
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           ! Localize variables
           DO np=1,npmax
               phytopo(np) = max(phytoP(i,j,k,bi,bj,np),0. _d 0)   !mmol C/m3
               phymmhgo(np) = max(Ptr(i,j,k,bi,bj,iPhy+np-1),0. _d 0) !mol/m3
           ENDDO
           DO nz=1,nzmax
               zoommhgo(nz) = max(Ptr(i,j,k,bi,bj,iZoo+nz-1),0. _d 0) !mol/m3
               zoopo(nz) = max(zooP(i,j,k,bi,bj,nz),0. _d 0)       !mmolC/m3
               zoomorto(nz) = max(zoomort(i,j,k,bi,bj,nz),0. _d 0) !mmolC/m3/s 
               DO np=1,npmax
                 zoograzeo(np,nz)=max(zoograze(i,j,k,bi,bj,np,nz), !mmolC/m3/s
     &                                0. _d 0)
               ENDDO 
           ENDDO

           ! Graze
           mmhgfp_diag(i,j,k) = 0. _d 0
           DO np=1,npmax
           DO nz=1,nzmax
             grazeo = zoograzeo(np,nz) / phytopo(np)                 ! /s
             grazeo = min(grazeo,1.3888 _d -4)                       ! cap it with 1/2 time step -1
             grazeo = grazeo * phymmhgo(np)                          ! mol/m3/s
             Gchem(i,j,k,bi,bj,iPhy+np-1) = 
     &                Gchem(i,j,k,bi,bj,iPhy+np-1) -
     &                grazeo
             Gchem(i,j,k,bi,bj,iZoo+nz-1) = 
     &                Gchem(i,j,k,bi,bj,iZoo+nz-1) +
     &                grazeo * ae(nz)
             Gchem(i,j,k,bi,bj,immhg) = 
     &                Gchem(i,j,k,bi,bj,immhg)
     &                + grazeo * (1.0 _d 0 - ae(nz))
#ifdef ALLOW_DIAGNOSTICS
             mmhggz_diag(i,j,k,np,nz) = grazeo
             mmhgfp_diag(i,j,k) = mmhgfp_diag(i,j,k) + grazeo * 
     &                (1.0 _d 0 - ae(nz))  
#endif
           ENDDO
           ENDDO

           ! release via death
           DO nz=1,nzmax
             deatho = zoomorto(nz) / zoopo(nz)
             deatho = min(deatho,2.7776 _d -4)
             deatho = deatho * zoommhgo(nz)                          ! mol/m3/s
             Gchem(i,j,k,bi,bj,iZoo+nz-1)=Gchem(i,j,k,bi,bj,iZoo+nz-1)
     &                - deatho
             Gchem(i,j,k,bi,bj,immhg)=Gchem(i,j,k,bi,bj,immhg)
     &                + deatho
#ifdef ALLOW_DIAGNOSTICS
             mmhgmt_diag(i,j,k,nz) = deatho
#endif
           ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO

c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        DO np=1,npmax
         WRITE(diagName,'(A6,I2.2)') 'MMHGBC',np
         CALL DIAGNOSTICS_FILL( mmhgbc_diag(1-Olx,1-Oly,1,np), diagname,
     &                          0,Nr,2,bi,bj,myThid )
         DO nz=1,nzmax
          WRITE(diagName,'(A4,I2.2,I2.2)') 'MHGZ',np,nz
          CALL DIAGNOSTICS_FILL( mmhggz_diag(1-Olx,1-Oly,1,np,nz),
     &                          diagname,
     &                          0,Nr,2,bi,bj,myThid )
         ENDDO
        ENDDO
        DO nz=1,nzmax
         WRITE(diagName,'(A6,I2.2)') 'MHMORT',nz
         CALL DIAGNOSTICS_FILL( mmhgmt_diag(1-Olx,1-Oly,1,nz), diagname,
     &                          0,Nr,2,bi,bj,myThid )
        ENDDO
        WRITE(diagName,'(A8)') 'MHHGFLOP'
        CALL DIAGNOSTICS_FILL( mmhgfp_diag(1-Olx,1-Oly,1), diagname,
     &                          0,Nr,2,bi,bj,myThid ) 
#endif

#endif
#endif
#endif
       RETURN
       END
