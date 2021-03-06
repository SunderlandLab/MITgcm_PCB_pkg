C $Header: /u/gcmpack/MITgcm/pkg/atm2d/put_ocnvars.F,v 1.6 2010/04/30 20:25:59 jscott Exp $
C $Name:  $

#include "ctrparam.h"
#include "ATM2D_OPTIONS.h"

C     !INTERFACE:
      SUBROUTINE PUT_OCNVARS(  myTime, myIter, myThid )
C     *==========================================================*
C     | Put fluxes for the ocean where needed: MPI mode, into    |
c     | common block for passing; single-proc, into FFIELDS.h    |
C     *==========================================================*
        IMPLICIT NONE

C     === Global Atmosphere Variables ===
#include "ATMSIZE.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "THSICE_VARS.h"
#include "ATM2D_VARS.h"

#ifdef ATM2D_MPI_ON
#  include "OCNSIZE.h"
#  include "OCNVARS.h"
#  include "OCNIDS.h"
#else
#  include "FFIELDS.h"
#  include "SURFACE.h"
#  ifdef ALLOW_DIC
#    include "DIC_VARS.h"
#  endif
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime - current simulation time (ocean model time)
C     myIter - iteration number (ocean model)
C     myThid - Thread no. that called this routine.
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     LOCAL VARIABLES:
      INTEGER i,j

      DO j=1,sNy
        DO i=1,sNx

#ifdef ATM2D_MPI_ON

C OCNVARS common set from ATM2D common
          atmSLPr_ocn(i,j)= pass_slp(i,j)
          HeatFlux_ocn(i,j)= pass_qnet(i,j)
          qShortWave_ocn(i,j)= pass_solarnet(i,j)
          TauX_ocn(i,j)= pass_fu(i,j)
          TauY_ocn(i,j)= pass_fv(i,j)
C         convert FW to MITGCM units kg/m2/s for passing
          FWFlux_ocn(i,j)= (pass_precip(i,j) + pass_evap(i,j)
     &                     + pass_runoff(i,j)) * rhoConstFresh
          SaltFlx_ocn(i,j) = sFluxFromIce(i,j)
#  ifdef ATMOSPHERIC_LOADING
          sIceMass_ocn(i,j) = pass_sIceLoad(i,j)
#  endif /* ATMOSPHERIC_LOADING */
          IF ( ocnCpl_exchange_DIC ) THEN
            aCO2_ocn(i,j) = pass_pCO2(i,j)
            wSpeed_ocn(i,j) = pass_wspeed(i,j)
            sIceFrac_ocn(i,j) = iceMask(i,j,1,1)
          ENDIF

#else
C FFIELDS.h common block set from ATM2D common

          pLoad(i,j,1,1)= pass_slp(i,j)
          Qnet(i,j,1,1)= pass_qnet(i,j)
          Qsw(i,j,1,1)= pass_solarnet(i,j)
          fu(i,j,1,1)= pass_fu(i,j)
          fv(i,j,1,1)= pass_fv(i,j)
          EmPmR(i,j,1,1)= (pass_precip(i,j) + pass_evap(i,j)
     &                     + pass_runoff(i,j)) * rhoConstFresh
          saltFlux(i,j,1,1)= sFluxFromIce(i,j)

          CALL EXCH_UV_XY_RS(fu,fv,.TRUE.,myThid)
          _EXCH_XY_RS(EmPmR, myThid )   !for NL Free surf?
          _EXCH_XY_RS(saltFlux, myThid )
          _EXCH_XY_RS(Qnet, myThid )
          _EXCH_XY_RS(Qsw, myThid )

#  ifdef ATMOSPHERIC_LOADING
          sIceLoad(i,j,1,1) = pass_sIceLoad(i,j)
          _EXCH_XY_RS(pLoad, myThid )
          _EXCH_XY_RS( sIceLoad, myThid )
C          phi0surf(i,j,1,1) = pass_slp(i,j)*recip_rhoConst
C     &                 + gravity*seaIceMass(i,j,1,1)*recip_rhoConst
#  endif /* ATMOSPHERIC_LOADING */

#  ifdef ALLOW_DIC
C          icemask should be imported in DIC code
           AtmospCO2(i,j,1,1)=  pass_pCO2(i,j)
           wind(i,j,1,1)= pass_wspeed(i,j)
           _EXCH_XY_RL(AtmospCO2 ,  myThid )
           _EXCH_XY_RL(wind ,  myThid )
#  endif /* ALLOW_DIC */

#endif /* ATM2D_MPI_ON */

        ENDDO
      ENDDO

      RETURN
      END

