C $Header: /u/gcmpack/MITgcm/pkg/ecco/cost_drift.F,v 1.12 2014/10/18 18:15:44 gforget Exp $
C $Name:  $

#include "ECCO_OPTIONS.h"

      subroutine cost_drift( myiter, mytime, mythid )

c     ==================================================================
c     SUBROUTINE cost_drift
c     ==================================================================
c
c     o Evaluate cost function contribution of the t and S difference
c       between the first and the last year.
c
c     started: from the "old" code
c
c              Elisabeth Remy eremy@ucsd.edu   july 31 2001
c
c     ==================================================================
c     SUBROUTINE cost_drift
c     ==================================================================

      implicit none

c     == global variables ==

#ifdef ALLOW_DRIFT_COST_CONTRIBUTION
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

#include "cal.h"
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

#ifdef ALLOW_DRIFT_COST_CONTRIBUTION
c     == external functions ==
      integer  ilnblnk
      external ilnblnk

c     == local variables ==
      _RS        one_rs
      parameter( one_rs = 1. )

      integer bi,bj
      integer i,j,k
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer irec
      integer iltheta
      integer ilsalt
      integer nf, nl, nfmin
      integer minrec

      _RL fctilet
      _RL fctiles
      _RL fcthread_tdrift
      _RL fcthread_tdrifs
      _RL errtannurescal
      _RL errsannurescal

      character*(80) fnametheta
      character*(80) fnamesalt

      logical doglobalread
      logical ladinit

      character*(MAX_LEN_MBUF) msgbuf

      _RL diagnosfld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)

c     == end of interface ==

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
      jmin = 1
      jmax = sny
      imin = 1
      imax = snx

c--   Read tiled data.
      doglobalread = .false.
      ladinit      = .false.

c--   Rescale error read from data.err to convert to error/year
c--   (e.g. 0.08 degC year based on Levitus decadal drift 0.8 degC)
      errtannurescal = 0.16
      errsannurescal = 0.16

      if (optimcycle .ge. 0) then
        iltheta = ilnblnk( tbarfile )
        write(fnametheta(1:80),'(2a,i10.10)')
     &    tbarfile(1:iltheta),'.',optimcycle
        ilsalt = ilnblnk( sbarfile )
        write(fnamesalt(1:80),'(2a,i10.10)')
     &    sbarfile(1:ilsalt),'.',optimcycle
      endif

      fcthread_tdrift = 0. _d 0
      fcthread_tdrifs = 0. _d 0

      do bj = jtlo,jthi
        do bi = itlo,ithi
          do k = 1,nr
            do j = jmin,jmax
              do i = imin,imax
                Tfmean(i,j,k,bi,bj) = 0.0
                Sfmean(i,j,k,bi,bj) = 0.0
                Tlmean(i,j,k,bi,bj) = 0.0
                Slmean(i,j,k,bi,bj) = 0.0
              enddo
            enddo
          enddo
        enddo
      enddo

      nf = 0
      nl = 0
c--   Number of full years
      nfmin = MAX(INT(FLOAT(nmonsrec)/12.),1)
c--   Prevent code from crashing if integrated for less than a year
      minrec = MIN(nmonsrec,12)

c--   Loop over records.
      do irec = 1,minrec

c--     Read time averages and the monthly mean data.
        call active_read_xyz( fnametheta, tbar, irec,
     &                        doglobalread, ladinit,
     &                        optimcycle, mythid,
     &                        xx_tbar_mean_dummy )

        call active_read_xyz( fnamesalt, sbar, irec,
     &                        doglobalread, ladinit,
     &                        optimcycle, mythid,
     &                        xx_sbar_mean_dummy )

        nf = nf + 1
        do bj = jtlo,jthi
          do bi = itlo,ithi
            do k = 1,nr
              do j = jmin,jmax
                do i = imin,imax
                  Tfmean(i,j,k,bi,bj) = Tfmean(i,j,k,bi,bj) +
     &                  tbar(i,j,k,bi,bj)
                  Sfmean(i,j,k,bi,bj) = Sfmean(i,j,k,bi,bj) +
     &                  sbar(i,j,k,bi,bj)
                enddo
              enddo
            enddo
          enddo
        enddo

      enddo

      do irec = nmonsrec-minrec+1, nmonsrec

c--     Read time averages and the monthly mean data.
        call active_read_xyz( fnametheta, tbar, irec,
     &                        doglobalread, ladinit,
     &                        optimcycle, mythid,
     &                        xx_tbar_mean_dummy )

        call active_read_xyz( fnamesalt, sbar, irec,
     &                        doglobalread, ladinit,
     &                        optimcycle, mythid,
     &                        xx_sbar_mean_dummy )

        nl = nl + 1

        do bj = jtlo,jthi
          do bi = itlo,ithi
            do k = 1,nr
              do j = jmin,jmax
                do i = imin,imax
                  Tlmean(i,j,k,bi,bj) = Tlmean(i,j,k,bi,bj) +
     &                  tbar(i,j,k,bi,bj)
                  Slmean(i,j,k,bi,bj) = Slmean(i,j,k,bi,bj) +
     &                  sbar(i,j,k,bi,bj)
                enddo
              enddo
            enddo
          enddo
        enddo

      enddo

      do bj = jtlo,jthi
        do bi = itlo,ithi

c--       Loop over the model layers
          fctiles = 0. _d 0
          do k = 1,nr
            do j = jmin,jmax
             do i = imin,imax
              if (_hFacC(i,j,k,bi,bj) .ne. 0.) then
                     fctiles = fctiles +
     &                  (wsaltLev(i,j,k,bi,bj)*cosphi(i,j,bi,bj)/
     &          ((nfmin*errsannurescal)**2)*
     &          (Slmean(i,j,k,bi,bj)/nl - Sfmean(i,j,k,bi,bj)/nf)*
     &          (Slmean(i,j,k,bi,bj)/nl - Sfmean(i,j,k,bi,bj)/nf))
c
                     diagnosfld3d(i,j,k,bi,bj) =
     &                  (wsaltLev(i,j,k,bi,bj)*cosphi(i,j,bi,bj)/
     &          ((nfmin*errsannurescal)**2)*
     &          (Slmean(i,j,k,bi,bj)/nl - Sfmean(i,j,k,bi,bj)/nf)*
     &          (Slmean(i,j,k,bi,bj)/nl - Sfmean(i,j,k,bi,bj)/nf))
c
                     if ( wsaltLev(i,j,k,bi,bj)*cosphi(i,j,bi,bj)
     &                 .ne. 0. )
     &                 num_sdrift(bi,bj) = num_sdrift(bi,bj) + 1. _d 0
              else
                     diagnosfld3d(i,j,k,bi,bj) = 0.
              endif
             enddo
            enddo
          enddo
          fcthread_tdrifs    = fcthread_tdrifs    + fctiles
          objf_sdrift(bi,bj) = objf_sdrift(bi,bj) + fctiles

        enddo
      enddo
      CALL WRITE_REC_XYZ_RL( 'DiagnosCost_DriftSalt',
     &               diagnosfld3d, 1, optimcycle, mythid )

      do bj = jtlo,jthi
        do bi = itlo,ithi

          fctilet = 0. _d 0
          do k = 1,nr
            do j = jmin,jmax
             do i = imin,imax
              if (_hFacC(i,j,k,bi,bj) .ne. 0.) then
                     fctilet = fctilet +
     &                  (wthetaLev(i,j,k,bi,bj)*cosphi(i,j,bi,bj)/
     &          ((nfmin*errtannurescal)**2)*
     &          (Tlmean(i,j,k,bi,bj)/nl - Tfmean(i,j,k,bi,bj)/nf)*
     &          (Tlmean(i,j,k,bi,bj)/nl - Tfmean(i,j,k,bi,bj)/nf))
c
                     diagnosfld3d(i,j,k,bi,bj) =
     &                  (wthetaLev(i,j,k,bi,bj)*cosphi(i,j,bi,bj)/
     &          ((nfmin*errtannurescal)**2)*
     &          (Tlmean(i,j,k,bi,bj)/nl - Tfmean(i,j,k,bi,bj)/nf)*
     &          (Tlmean(i,j,k,bi,bj)/nl - Tfmean(i,j,k,bi,bj)/nf))
c
                     if ( wthetaLev(i,j,k,bi,bj)*cosphi(i,j,bi,bj)
     &                 .ne. 0. )
     &                 num_tdrift(bi,bj) = num_tdrift(bi,bj) + 1. _d 0
              else
                     diagnosfld3d(i,j,k,bi,bj) = 0.
              endif
             enddo
            enddo
          enddo
          fcthread_tdrift    = fcthread_tdrift    + fctilet
          objf_tdrift(bi,bj) = objf_tdrift(bi,bj) + fctilet

        enddo
      enddo
      CALL WRITE_REC_XYZ_RL( 'DiagnosCost_DriftTheta',
     &               diagnosfld3d, 1, optimcycle, mythid )

#endif /* ALLOW_DRIFT_COST_CONTRIBUTION */

      RETURN
      END
