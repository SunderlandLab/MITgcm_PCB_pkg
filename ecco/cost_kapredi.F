C $Header: /u/gcmpack/MITgcm/pkg/ecco/cost_kapredi.F,v 1.7 2014/10/18 18:15:44 gforget Exp $
C $Name:  $

#include "ECCO_OPTIONS.h"


      subroutine cost_kapredi(
     I                           myiter,
     I                           mytime,
     I                           mythid
     &                         )

C     o==========================================================o
C     | subroutine cost_kapredi                                    |
C     | o REDI coefficient adjustment penalization                 |
C     o==========================================================o

      implicit none

c     == global variables ==

#ifdef ALLOW_KAPREDI_COST_CONTRIBUTION
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_REDIREDI
# include "REDIREDI.h"
#endif

#include "ecco_cost.h"
#include "CTRL_SIZE.h"
#include "ctrl.h"
#include "ctrl_dummy.h"
#include "optim.h"
#endif

c     == routine arguments ==

      integer myiter
      _RL     mytime
      integer mythid

c     == local variables ==

#ifdef ALLOW_KAPREDI_COST_CONTRIBUTION

      integer bi,bj
      integer i,j,k
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer nrec
      integer irec
      integer ilfld

      _RL fctile
      _RL fcthread
      _RL tmpx

      logical doglobalread
      logical ladinit

      character*(80) fnamefld

      character*(MAX_LEN_MBUF) msgbuf

c     == external functions ==

      integer  ilnblnk
      external ilnblnk

c     == end of interface ==

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
      jmin = 1
      jmax = sny
      imin = 1
      imax = snx

c--   Read state record from global file.
      doglobalread = .false.
      ladinit      = .false.

      irec = 1

      if (optimcycle .ge. 0) then
        ilfld = ilnblnk( xx_kapredi_file )
        write(fnamefld(1:80),'(2a,i10.10)')
     &       xx_kapredi_file(1:ilfld),'.',optimcycle
      endif

      fcthread = 0. _d 0

      call active_read_xyz( fnamefld, tmpfld3d, irec, doglobalread,
     &                       ladinit, optimcycle, mythid
     &        , xx_kapredi_dummy )

c--     Loop over this thread tiles.
        do bj = jtlo,jthi
          do bi = itlo,ithi

c--         Determine the weights to be used.

            fctile = 0. _d 0
            do k = 1,nr
            do j = jmin,jmax
              do i = imin,imax
                if (_hFacC(i,j,k,bi,bj) .ne. 0.) then
c                  tmpx = (tmpfld3d(i,j,k,bi,bj)-GM_isopycK)
                  tmpx = tmpfld3d(i,j,k,bi,bj)
                IF ( .NOT.ctrlSmoothCorrel3D ) THEN
                  fctile = fctile
     &                 + wkaprediFld(i,j,k,bi,bj)*cosphi(i,j,bi,bj)
     &                 *tmpx*tmpx
                ELSE !IF ( .NOT.ctrlSmoothCorrel3D ) THEN
                  fctile = fctile + tmpx*tmpx
                ENDIF !IF ( .NOT.ctrlSmoothCorrel3D ) THEN
                endif
              enddo
            enddo
            enddo

            objf_kapredi(bi,bj) = objf_kapredi(bi,bj) + fctile
            fcthread          = fcthread + fctile

#ifdef ECCO_VERBOSE
c--         Print cost function for each tile in each thread.
            write(msgbuf,'(a)') ' '
            call print_message( msgbuf, standardmessageunit,
     &                          SQUEEZE_RIGHT , mythid)
            write(msgbuf,'(a,i8.8,1x,i3.3,1x,i3.3)')
     &        ' cost_kapredi: irec,bi,bj          =  ',irec,bi,bj
            call print_message( msgbuf, standardmessageunit,
     &                          SQUEEZE_RIGHT , mythid)
            write(msgbuf,'(a,d22.15)')
     &        '               cost function (dT(0)) = ',
     &        fctile
            call print_message( msgbuf, standardmessageunit,
     &                      SQUEEZE_RIGHT , mythid)
#endif
          enddo
        enddo

#ifdef ECCO_VERBOSE
c--     Print cost function for all tiles.
        _GLOBAL_SUM_RL( fcthread , myThid )
        write(msgbuf,'(a)') ' '
        call print_message( msgbuf, standardmessageunit,
     &                      SQUEEZE_RIGHT , mythid)
        write(msgbuf,'(a,i8.8)')
     &    ' cost_kapredi:                       irec =  ',irec
        call print_message( msgbuf, standardmessageunit,
     &                      SQUEEZE_RIGHT , mythid)
        write(msgbuf,'(a,d22.15)')
     &    '                 global cost function value = ',
     &    fcthread
        call print_message( msgbuf, standardmessageunit,
     &                      SQUEEZE_RIGHT , mythid)
        write(msgbuf,'(a)') ' '
        call print_message( msgbuf, standardmessageunit,
     &                      SQUEEZE_RIGHT , mythid)
#endif

#endif

      return
      end


