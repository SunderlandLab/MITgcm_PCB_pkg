C $Header: /u/gcmpack/MITgcm/pkg/pcb/pcb_readparms.F,Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "PCB_OPTIONS.h"

CBOP
C     !ROUTINE: PCB_READPARMS
C     !INTERFACE:
      SUBROUTINE PCB_READPARMS( myThid )

C     !DESCRIPTION: \bv

      !==========================================================
      ! SUBROUTINE PCB_READPARMS
      !   Initialize and read pcb package parameters
      !
      ! REVISION HISTORY
      !   03 Feb 2015 - H. Amos - Update river inputs of PCBs to
      !                           mimic Yanxu Zhang's Hg code.  
      !==========================================================
C     \ev

      !----------------------------------------------------------
      ! USES:
      !----------------------------------------------------------
      IMPLICIT NONE

      !----------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PCB_SIZE.h"
#include "PCB.h"

C Commented out to mirror Yanxu's updated Hg code in ECCOv4
C (hma, 20 May 2015)
C #ifdef ALLOW_PTRACERS
C # include "PTRACERS_SIZE.h"
C # include "PTRACERS_PARAMS.h"
C #endif

      !----------------------------------------------------------
      ! INPUT/OUTPUT PARAMETERS:
      !----------------------------------------------------------

      ! Routine arguments 
      !   myThid    :: My Thread Id. number
      INTEGER myThid
CEOP

#ifdef ALLOW_PCB

      ! Local variables
      !   msgBuf    :: Informational/error message buffer
      !   iUnit     :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

      NAMELIST /PCB_FORCING/
     &          PCB_windFile, PCB_iceFile,
     &          PCB_forcingPeriod, PCB_forcingCycle,
c xz:the following lines are modified from hg code  
c PCBa-neutral PCB; PCBi-ionized PCB; PCBp-particle bounded PCB      
     &          PCB_chemPeriod, PCB_chemCycle,
     &          PCB_emisPeriod, PCB_emisCycle,
c Remove PCBi tracer. (hma, 15 Jul 2014)
c     &          PCBa_concFile, PCBi_depoFile,
c     &          PCBi_riverFile, PCBp_riverFile,
     &          PCBp_depoFile, PCBa_depoFile,
     &          PCBa_concFile,
     &          PCBp_riverFile,
     &          PCBa_riverFile,
#ifdef PCBAREFRACTORY
     &          PCBAR_riverFile,
#endif
#ifdef PCBPREFRACTORY
     &          PCBPR_riverFile,
#endif     
     &          chl_file, npp_file,
     &          doc_file, poc_file,
C Fdoc doesn't appear to be used anywhere and there isn't an input file
C (hma, 01 july 2015)
C     &          Fdoc_file, Fpoc_file, Fremin_file
     &          Fpoc_file, Fremin_file
#ifdef FOODW_MODE     
     &          , NO3_file,phytoP_file, zooP_file
#endif     
     
      _BEGIN_MASTER(myThid)

       ! Initialize file names
       PCB_windFile   = ' '
       PCB_iceFile    = ' '
       radsw_File     = ' '
       PCBa_concFile  = ' '
       PCBp_riverFile = ' '      
       PCBa_riverFile = ' '      
#ifdef PCBAREFRACTORY
       PCBaR_riverFile= ' '
#endif
#ifdef PCBPREFRACTORY
       PCBpR_riverFile= ' '
#endif       
       chl_file       = ' '
       npp_file       = ' '
       doc_file       = ' '
       poc_file       = ' '
C Fdoc doesn't appear to be used anywhere and there isn't an input file
C (hma, 01 july 2015)
C       Fdoc_file      = ' '
       Fpoc_file      = ' '
       Fremin_file    = ' '
#ifdef FOODW_MODE       
       NO3_file       = ' '
       phytoP_file    = ' '
       zooP_file      = ' '
#endif       
      
       ! Default periodic forcing to same as for physics
       PCB_forcingPeriod = externForcingPeriod
       PCB_forcingCycle  = externForcingCycle
       PCB_emisPeriod = externForcingPeriod
       PCB_emisCycle  = externForcingCycle       

       ! Added based on hg code default periodic forcing 
       ! to same as for physics (Xianming)
       PCB_chemPeriod = externForcingPeriod
       PCB_chemCycle  = externForcingCycle       

      WRITE(msgBuf,'(A)') ' PCB_READPARMS: opening data.pcb'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      CALL OPEN_COPY_DATA_FILE( 'data.pcb', 'PCB_READPARMS',
     O                          iUnit, myThid )

      !----------------------------------------------------------
      ! Read parameters from open data file
      !----------------------------------------------------------

      ! Forcing filenames and parameters
      READ(UNIT=iUnit,NML=PCB_FORCING)

      WRITE(msgBuf,'(A)')
     &   ' PCB_READPARMS: finished reading data.pcb'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      ! Close the open data file
      CLOSE(iUnit)

c==============the following used in cfc but not in hg
c      IF ( atmPCB_timeOffset .EQ. UNSET_RL ) THEN
C-    Assume that pcbTime=0 corresponds to the beginning of the 1rst record
C     time-period. This is consistent with 1rst record value = time-averaged
C     atmos-PCB over time period: pcbTime= 0 to pcbTime= 1 x atmPCB_recSepTime
C-    To reproduce previous time interpolation (formerly in pcb02_forcing.F),
C     set default timeOffset according to PTRACERS_Iter0 :
c        atmPCB_timeOffset = atmPCB_recSepTime
c     &                     - deltaTClock*PTRACERS_Iter0
c      ENDIF

c      PCB_monFreq = 0. _d 0
c      IF ( debugLevel.GE.debLevB ) PCB_monFreq = PTRACERS_monitorFreq
c==============


      _END_MASTER(myThid)

       ! Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_PCB */

      RETURN
      END
