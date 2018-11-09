C $Header: /u/gcmpack/MITgcm/pkg/pcb/pcb_diagnostics_init.F,v 1.0 
C $Name:  $

#include "CPP_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"
#include "PCB_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_PCB

CBOP
C==================================================================
C ROUTINE HG_DIAGNOSTICS_INIT
C INTERFACE:
      SUBROUTINE PCB_DIAGNOSTICS_INIT( myThid )
C
C DESCRIPTION:
C   Routine to initialize PCB diagnostics
C
C MODIFICATION HISTORY 
C   15 Jul 2014 - hma - Comment out PCBi tracer. Rename PCBn to
C                       PCBa.    
C   03 Feb 2015 - hma - Update rivers based on Yanxu Zhang's Hg 
C                       code.
C   07 Dec 2016 - ccw - Added SURFSALT, SURFTEMP, PCBaFUGA
C==================================================================

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "PTRACERS_SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GCHEM.h"
#include "PCB_SIZE.h"
#include "PCB.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid -  Number of this instance of HG_DIAGNOSTICS_INIT
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf      - Informational/error meesage buffer
c     CHARACTER*(MAX_LEN_MBUF) msgBuf

      INTEGER       diagNum
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

      INTEGER      np

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Add diagnostics to the (long) list

C--   Add diagnostic PCBiDEPO, PCBaCONC, PCBaEVAS, PCBiREDN, and PCBPSINK

      diagName = 'PCBaSPRY'
      WRITE(diagTitle, '(A)') 'PCBa seaspray '
      diagUnits = 'mol/m^2/s          '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c Remove PCBi tracer. (hma, 15 Jul 2014)
c      diagName = 'PCBiSPRY'
c      WRITE(diagTitle, '(A)') 'PCBi seaspray '
c      diagUnits = 'mol/m^2/s          '
c      diagCode  = 'SM P    M1      '
c      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaCONC'
      WRITE(diagTitle, '(A)') 'Atmospheric PCBa concentration'
C      diagUnits = 'ng/m^3          ' ! updated units 29 jun 2015 (hma)
      diagUnits = 'mol/m^3          '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaEVAS'
      WRITE(diagTitle, '(A)') 'PCBa evasion flux from ocean'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaFLUP'
      WRITE(diagTitle, '(A)') 'PCBa flux from ocean'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
          
      diagName = 'PCBaFLDN'
      WRITE(diagTitle, '(A)') 'PCBa flux from atmosphere'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c Remove PCBi tracer. (hma, 15 Jul 2014)
c      diagName = 'PCBiDEPO'
c      WRITE(diagTitle, '(A)') 'Atmospheric PCBi deposition flux'
c      diagUnits = 'mol/m^2/s       '
c      diagCode  = 'SM P    M1      '
c      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c add diagnostics SURFTEMP, SURFSALT, PCBaFUGA (ccw 7 dec 2016)

      diagName = 'SURFTEMP'
      WRITE(diagTitle, '(A)') 'surface water temperature'
      diagUnits = 'Celsius        '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'SURFSALT'
      WRITE(diagTitle, '(A)') 'surface water salinity'
      diagUnits = 'g/kg            '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaFUGA'
      WRITE(diagTitle, '(A)') 'PCBa fugacity ratio fatm over focean'
      diagUnits = 'dimensionless   '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaDEPO'
      WRITE(diagTitle, '(A)') 'Atmospheric PCBa deposition flux'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBpDEPO'
      WRITE(diagTitle, '(A)') 'Atmospheric PCBp deposition flux'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBaRIVR'
      WRITE(diagTitle, '(A)') 'PCBa riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c Remove PCBi tracer. (hma, 15 Jul 2014)
c      diagName = 'PCBiRIVR'
c      WRITE(diagTitle, '(A)') 'PCBi riverine runoff'
c      diagUnits = 'mol/m^2/s       '
c      diagCode  = 'SM P    M1      '
c      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBpRIVR'
      WRITE(diagTitle, '(A)') 'PCBp riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 

      ! Mimic Hg for PCB river inputs
      ! (H. Amos, 03 Feb 2015) 
#ifdef PCBAREFRACTORY
      diagName = 'PCBaRRIV'
      WRITE(diagTitle, '(A)') 'Refractory riverine PCBa riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )          

      diagName = 'PCBaRDIS'
      WRITE(diagTitle, '(A)') 'Refractory river PCBa dissolution flux'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif  
     
#ifdef PCBAPREFRACTORY
      diagName = 'PCBPRSNK'
      WRITE(diagTitle, '(A)') 'Refrac PCBp sinking flux at level bottom'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBPRRIV'
      WRITE(diagTitle, '(A)') 'Refrac. PCBp riverine runoff'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )          
          
      diagName = 'PCBPRDIS'
      WRITE(diagTitle, '(A)') 'Refractory river PCBp dissolution flux'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif           

      diagName = 'PCBaPART'
      WRITE(diagTitle, '(A)') 'PCBa partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c Remove PCBi tracer. (hma, 15 Jul 2014)
c      diagName = 'PCBiPART'
c      WRITE(diagTitle, '(A)') 'PCBi partitioned onto particles'
c      diagUnits = 'mol/m^3/s       '
c      diagCode  = 'SM P    MR      '
c      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName = 'PCBpPART'
      WRITE(diagTitle, '(A)') 'PCBp partitioned onto particles'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName = 'PCBpSINK'
      WRITE(diagTitle, '(A)') 'PCBp sinking flux at level bottom'
      diagUnits = 'mol/m^2/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )         

      diagName = 'PCBaDEGR'
      WRITE(diagTitle, '(A)') 'PCBa degradation flux'
      diagUnits = 'mol/m^3/s       '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )         


#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END     
C============================================================================
#endif
#endif
