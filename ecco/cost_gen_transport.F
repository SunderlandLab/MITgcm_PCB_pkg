C $Header: /u/gcmpack/MITgcm/pkg/ecco/cost_gen_transport.F,v 1.8 2014/10/18 18:15:44 gforget Exp $
C $Name:  $

#include "ECCO_OPTIONS.h"


      subroutine cost_gen_transport(
     I                     myiter,
     I                     mytime,
     I                     mythid
     &                   )

c     ==================================================================
c     SUBROUTINE cost_gen_transport
c     ==================================================================
c
c     o Evaluate cost function contribution for transports
c
c     ==================================================================
c     SUBROUTINE cost_ssh
c     ==================================================================

      implicit none

c     == global variables ==

#ifdef ALLOW_TRANSPORT_COST_CONTRIBUTION
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"

#include "ecco_cost.h"
#include "ctrl.h"
#include "optim.h"
#include "DYNVARS.h"
#ifdef ALLOW_PROFILES
#include "profiles.h"
#endif
#endif

c     == routine arguments ==

      integer myiter
      _RL     mytime
      integer mythid

#ifdef ALLOW_TRANSPORT_COST_CONTRIBUTION
c     == local variables ==

      integer bi,bj
      integer itlo,ithi
      integer jtlo,jthi
      integer irec
      integer ilps
      logical doglobalread
      logical ladinit
      _RL     ftmp
      character*(max_len_mbuf) msgbuf

c     == external functions ==

      integer  ilnblnk
      external ilnblnk

c     == end of interface ==

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

c--   Loop over records for the second time.
      do irec = 1, ndaysrec

        ftmp = 0. _d 0
        do bj = jtlo,jthi
          do bi = itlo,ithi
             ftmp = ftmp + transpbar(irec,bi,bj)
          enddo
        enddo
        _GLOBAL_SUM_RL( ftmp , myThid )

cph        write(msgbuf,'(A,I,2(X,P2E22.14))')
cph     &       'ph-cost-transport FS: day, 1-model, 2-obs ',
cph     &       irec, ftmp, transpobs(irec)
cph        call print_message( msgbuf, standardmessageunit,
cph     &       SQUEEZE_RIGHT , mythid)

        if ( ftmp.NE.0. .AND. wtransp(irec).NE.0. ) then
           objf_transp = objf_transp +
     &          wtransp(irec)*( ftmp - transpobs(irec) )**2
           num_transp = num_transp + 1. _d 0
        endif

      enddo

#endif /* ifdef ALLOW_TRANSPORT_COST_CONTRIBUTION */

      end
