C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/darwin_write_pickup.F,v 1.3 2013/05/15 14:42:54 jahn Exp $
C $Name:  $

#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_WRITE_PICKUP

C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_WRITE_PICKUP( permPickup,
     I                                suff, myTime, myIter, myThid )

C !DESCRIPTION:
C     Writes PAR_day array (needed for a restart) to a pickup file

C !USES: ===============================================================
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_IO.h"
#ifdef ALLOW_MONOD
#include "MONOD_FIELDS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C     permPickup :: write a permanent pickup
C     suff       :: suffix for pickup file (eg. ckptA or 0000000010)
C     myTime     :: Current time in simulation
C     myIter     :: Current iteration number in simulation
C     myThid     :: My Thread Id number
      LOGICAL permPickup
      CHARACTER*(*) suff
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_DARWIN
C add more ALLOWs here for other fields in pickup
#if defined(ALLOW_PAR_DAY) || (defined(ALLOW_MONOD) && defined(GEIDER) && !defined(DYNAMIC_CHL))

C     !LOCAL VARIABLES:
C     == Local variables ==
      CHARACTER*(MAX_LEN_FNAM) fn
      LOGICAL glf
      _RL     timList(1)
      INTEGER prec, irec, ifld, np
      INTEGER listDim, nWrFlds
      PARAMETER( listDim = 2+npmax )
      CHARACTER*(8) wrFldList(listDim)
      CHARACTER*(MAX_LEN_MBUF) msgBuf

c     IF ( DARWIN_pickup_write_mdsio ) THEN
        prec = precFloat64
        WRITE(fn,'(A,A)') 'pickup_darwin.',suff

        ifld = 0
        irec = 0

#ifdef ALLOW_PAR_DAY
C     record number < 0 : a hack not to write meta files now:
        ifld = ifld + 1
        irec = irec - 1
        CALL WRITE_REC_3D_RL( fn, prec, nR,
     &                        PARday(1-OLx,1-OLy,1,1,1,1),
     &                        irec, myIter, myThid )
        IF (ifld.LE.listDim) wrFldList(ifld) = 'PARday1 '
        ifld = ifld + 1
        irec = irec - 1
        CALL WRITE_REC_3D_RL( fn, prec, nR,
     &                        PARday(1-OLx,1-OLy,1,1,1,2),
     &                        irec, myIter, myThid )
        IF (ifld.LE.listDim) wrFldList(ifld) = 'PARday2 '
#endif

#if defined(ALLOW_MONOD) && defined(GEIDER) && !defined(DYNAMIC_CHL)
        DO np = 1, npmax
        ifld = ifld + 1
        irec = irec - 1
        print*,'prec',np,prec
        CALL WRITE_REC_3D_RL( fn, prec, nR,
     &                        Chl_phy(1-OLx,1-OLy,1,1,1,np),
     &                        irec, myIter, myThid )
        IF (ifld.LE.listDim) THEN
         wrFldList(ifld)(1:6) = 'ChlPhy'
         wrFldList(ifld)(7:8) = PTRACERS_ioLabel(np)
        ENDIF
        ENDDO
#endif

C--------------------------
        nWrFlds = ifld
        IF ( nWrFlds.GT.listDim ) THEN
          WRITE(msgBuf,'(2A,I5,A)') 'DARWIN_WRITE_PICKUP: ',
     &     'trying to write ',nWrFlds,' fields'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A,I5,A)') 'DARWIN_WRITE_PICKUP: ',
     &     'field-list dimension (listDim=',listDim,') too small'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DARWIN_WRITE_PICKUP (list-size Pb)'
        ENDIF

#ifdef ALLOW_MDSIO
C     uses this specific S/R to write (with more informations) only meta
C     files
        glf  = globalFiles
        timList(1) = myTime
        irec = ABS(irec)
        print*,'prec',irec,prec
        CALL MDS_WR_METAFILES( fn, prec, glf, .FALSE.,
     &                         0, 0, nR, ' ',
     &                         nWrFlds, wrFldList,
     &                         1, timList, oneRL,
     &                         irec, myIter, myThid )
#endif /* ALLOW_MDSIO */
C--------------------------

c     ENDIF /* DARWIN_pickup_write_mdsio */

#endif /*  ALLOW_PAR_DAY  */
#endif /*  ALLOW_DARWIN  */

      RETURN
      END
