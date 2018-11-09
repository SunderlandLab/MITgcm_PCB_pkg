C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_calc_tendency.F,v 1.5 2013/06/10 02:52:57 jmc Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: GCHEM_CALC_TENDENCY

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_CALC_TENDENCY(
     I     myTime, myIter, myThid )

C !DESCRIPTION:
C In the case of GCHEM_SEPARATE_FORCING not defined,
C this subroutine computes the tracer tendencies due to a
C bio-geogchemistry or ecosystem model and stores them on an array
C gchemTendency, that will be incorporated into regular timestepping in
C in ptracers_intergrate.F
C The current example uses the CFC package, but his is meant to
C be replaced by anything that the user provides.
C Add Hg support, yxzhang, 5/16/13
C Add PCB package (H. M. Amos, 11 May 2015)

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM.h"
#include "GCHEM_FIELDS.h"

C !INPUT PARAMETERS: ===================================================
C  myTime               :: current time
C  myIter               :: current timestep
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_GCHEM
# ifndef GCHEM_SEPARATE_FORCING
C !LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  iTracer              :: ptracer number
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER i,j
      INTEGER iTracer
      INTEGER bi,bj,iMin,iMax,jMin,jMax,k

C     gchemTendency is re-initialized here
      DO iTracer = 1, PTRACERS_numInUse
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            gchemTendency(i,j,k,bi,bj,iTracer) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C
C     Here is the place for code to compute bio-geochemical
C     tendency terms (sometimes referred to as source-minus-sink
C     terms). The tendencies are stored on gchemTendency, as show
C     in the CFC-example.
C
C     loop over tiles
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
C     define horizontal loop ranges
        iMin = 1
        iMax = sNx
        jMin = 1
        jMax = sNy
c       DO iTracer = 1, PTRACER_numInUse
c        DO k = 1, Nr
c         DO j = 1-OLy, sNy+OLy
c          DO i = 1-OLx, sNx+OLx
c           gchemTendency(i,j,k,bi,bj,iTracer) = your specific model
c          ENDDO
c         ENDDO
c        ENDDO
c       ENDDO

C------------------------
C chemical forcing      |
C------------------------
#ifdef ALLOW_CFC
        IF ( useCFC ) THEN
         CALL CFC11_FORCING(
     I                       pTracer     (1-OLx,1-OLy,1,bi,bj,1),
     U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,1),
     I                       bi, bj, iMin, iMax, jMin, jMax,
     I                       myTime, myIter, myThid )
         CALL CFC12_FORCING(
     I                       pTracer     (1-OLx,1-OLy,1,bi,bj,2),
     U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,2),
     I                       bi, bj, iMin, iMax, jMin, jMax,
     I                       myTime, myIter, myThid )
        ENDIF
#endif /* ALLOW_CFC */

c yxzhang, 5/16/13
#ifdef ALLOW_HG
ccccccccccccccccccccccccc
c chemical forcing      c
ccccccccccccccccccccccccc
ccccccccccccccccccccccccccc HG cccccccccccccccccccccccccccccccc
c
c Partitioning, this should be put as the first step
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing Hg partition '
        _END_MASTER( myThid )
      ENDIF
        CALL HG_PARTITION(  pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

c Foodweb
#ifdef FOODW_MODE
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing Hg foodweb'
        _END_MASTER( myThid )
      ENDIF
        CALL HG_FOODWEB(    pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )
#endif

c Atmospheric forcings
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing Hg forcing '
        _END_MASTER( myThid )
      ENDIF
        CALL HG_FORCING(    pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

c Chemistry, redox + methylation + demethylation
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing Hg chemistry '
        _END_MASTER( myThid )
      ENDIF
        CALL HG_CHEMISTRY(  pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

c Chemistry, riverine particle dissolution
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing Hg dissolution '
        _END_MASTER( myThid )
      ENDIF
        CALL HG_DISSOLUTION(  pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

c Biological pump, particulate matter sinking
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing Hg biopump '
        _END_MASTER( myThid )
      ENDIF
        CALL HG_BIOPUMP(    pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

cccccccccccccccccccccccccc END HG cccccccccccccccccccccccccccccccccc
#endif /* ALLOW_HG */
c yxzhang

c hma, 11 May 2015
#ifdef ALLOW_PCB
ccccccccccccccccccccccccc
c chemical forcing      c
ccccccccccccccccccccccccc
ccccccccccccccccccccccccccc PCB cccccccccccccccccccccccccccccccc
c
c Partitioning, this should be put as the first step
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing PCB partition '
        _END_MASTER( myThid )
      ENDIF

       CALL PCB_PARTITION(  pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

c PCBs don't have a foodweb interaction at this point
c (hma, 11 May 2015)
c
c Foodweb
c#ifdef FOODW_MODE
c      IF ( debugLevel .GE. debLevZero ) THEN
c        _BEGIN_MASTER( myThid )
c        WRITE(standardmessageunit,'(A)')
c     &  'doing Hg foodweb'
c        _END_MASTER( myThid )
c      ENDIF
c        CALL HG_FOODWEB(    pTracer     (1-OLx,1-OLy,1,1,1,1),
c     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
c     &                      bi, bj, iMin, iMax, jMin, jMax,
c     &                      myIter, myTime, myThid )
c#endif

c Atmospheric forcings
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing PCB forcing '
        _END_MASTER( myThid )
      ENDIF

        CALL PCB_FORCING(   pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

c Chemistry, redox 
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing PCB chemistry '
        _END_MASTER( myThid )
      ENDIF

        CALL PCB_CHEMISTRY( pTracer     (1-OLx,1-OLy,1,1,1,1) ,
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

c Chemistry, riverine particle dissolution
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing PCB dissolution '
        _END_MASTER( myThid )
      ENDIF

        CALL PCB_DISSOLUTION(pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

c Biological pump, particulate matter sinking
      IF ( debugLevel .GE. debLevZero ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardmessageunit,'(A)')
     &  'doing PCB biopump '
        _END_MASTER( myThid )
      ENDIF

        CALL PCB_BIOPUMP(   pTracer     (1-OLx,1-OLy,1,1,1,1),
     &                      gchemTendency(1-OLx,1-OLy,1,1,1,1),
     &                      bi, bj, iMin, iMax, jMin, jMax,
     &                      myIter, myTime, myThid )

cccccccccccccccccccccccccc END PCB cccccccccccccccccccccccccccccccccc
#endif /* ALLOW_PCB */
c hma, end


C     end of tile-loops
       ENDDO
      ENDDO
# endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

      RETURN
      END
