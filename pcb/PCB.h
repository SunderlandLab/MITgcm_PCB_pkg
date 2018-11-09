C $Header: /u/gcmpack/MITgcm/pkg/pcb/PCB.h,Exp $
C $Name:  $

      !==========================================================
      ! PCB.h
      !      
      ! REVISION HISTORY
      !   13 Jan 2015 - H. Amos - tHalf and Kdeg_pcb added. These
      !                           are variabiles to treat degradation
      !                           of PCBs in seawater. 
      !   03 Feb 2015 - H. Amos - Update river inputs based on 
      !                           Yanxu Zhang's Hg river code. 
      !   02 Jul 2015 - hma - Update air-sea exchange equations.
      !   06 Jul 2015 - hma - Added Kow
      !   05 Oct 2015 - hma - Kow moved pcb_param.F
      !   06 Oct 2016 - ccw - rename  dH to dH_AW, add dH_OW
      !==========================================================

       COMMON /PCB_FIELDS/
     &              PCBaConc,
     &              PCBaDEP, PCBpDEP,  
     &              PCBaRiver, PCBpRiver,
#ifdef PCBAREFRACTORY
     &              PCBaRRiver,
#endif
#ifdef PCBPREFRACTORY
     &              PCBpRRiver,
#endif      
     &              pisVel, fice, radsw, wind, ph,
     &              chl, npp, doc, poc,
C Fdoc is obsolete (hma, 06 2015)
C     &              Fdoc, Fpoc, Fremin
     &              Fpoc, Fremin

c Disable sea spray (H. Amos, 02 Feb 2015)
c#ifdef ALLOW_SEASPRAY
c     &           , seaspraywtrflux
c#endif     
      
#ifdef FOODW_MODE     
     &              , NO3, phytoP, zooP
#endif
      _RL  PCBaConc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)   ! atmospheric concentration
      _RL  PCBaDEP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)    ! atmospheric deposition
      _RL  PCBpDEP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)    ! atmospheric deposition
      _RL  PCBaRiver(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! river inputs, dissolved PCB
      _RL  PCBpRiver(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! river inputs, particulate PCB
#ifdef PCBaREFRACTORY
      _RL  PCBaRRiver(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef PCBPREFRACTORY
      _RL  PCBpRRiver(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nrmax)
#endif                  
c      _RL  PCBaAtmDep(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) ! atmospheric deposition, dissolved PCB
c      _RL  PCBpAtmDep(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) ! atmospheric deposition, particulate PCB 
      _RL  pisvel(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)     ! piston velocity 
      _RL  fice  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)     ! sea ice fraction
      _RL  wind (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)      ! wind speed
      _RL  ph (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)     
      _RL  radsw (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)     ! solar radiation
      _RL  chl   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  ! chlorophyll
      _RL  npp   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  ! net primary productivity
      _RL  doc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  ! dissolved organic carbon  
      _RL  poc   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  ! particulate organic carbon
C Fdoc isn't being used anywhere, so comment out (hma, 02 jul 2015)
C      _RL  Fdoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  ! Sinking flux of DOC
      _RL  Fpoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  ! Sinking flux of POC
      _RL  Fremin(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  ! Flux of remineralized POC and DOC

c Disable sea spray for PCBs (H. Amos, 02 Feb 2015)
c#ifdef ALLOW_SEASPRAY
c      _RL  seaspraywtrflux (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c#endif        
      
      
#ifdef FOODW_MODE      
      _RL  NO3   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)      
      _RL  phytoP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RL  zooP  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
#endif  

C     time-record currently loaded (in temp arrays *[1])
      COMMON /PCB_LOAD_I/
     &  PCB_ldRec_forcing, PCB_ldRec_chem, PCB_ldRec_emis
      INTEGER PCB_ldRec_forcing(nSx,nSy), PCB_ldRec_chem(nSx,nSy)
      INTEGER     PCB_ldRec_emis(nSx,nSy)

       COMMON /PCB_CHEM/
      !---------------------------------------------------  
      ! Physicochemical properties
      !---------------------------------------------------  
     &       pKa  , Koc_pcb , Kiw, logKow, dH_AW, dH_OW,
c     &       tHalf, Kdeg_pcb,
     &       tHalf,  
     &       a_deg, b_deg         ! hma, 22 Feb 2015  
        _RL  pKa  , Koc_pcb , Kiw, logKow
        _RL  tHalf ! degradation half life
        _RL  a_deg, b_deg         ! hma, 22 Feb 2015
        _RL  dH_AW, dH_OW         ! ccw, 07 Oct 2016 
      !---------------------------------------------------  
      ! Schmidt number coefficients and other constants
      ! to calculate air-sea exchange
      !---------------------------------------------------
      COMMON /PCB_PARAMETER/
     &     sv_0, sv_1, sv_2, sv_3, sv_4,
     &     sd_0, sd_1, sd_2, sd_3,
     &     Va  , Vb  , Ma  , Mb  , Patm,
     &     KH0                
      _RL  sv_0, sv_1, sv_2, sv_3, sv_4
      _RL  sd_0, sd_1, sd_2, sd_3
      _RL  Va  , Vb  , Ma  , Mb  , Patm
      _RL  KH0

      COMMON /PCB_LOAD/
     &    wind0, wind1, ice0, ice1, radsw0, radsw1,
c Remove PCBi tracer. (hma, 15 Jul 2014)
c     &    pcbaconc0, pcbaconc1, pcbidepo0, pcbidepo1,
c     &    pcbiriver0, pcbiriver1, pcbpriver0, pcbpriver1,
     &    pcbaconc0, pcbaconc1,    
     &    pcbadepo0, pcbadepo1,
     &    pcbpdepo0, pcbpdepo1,
     &    pcbpriver0, pcbpriver1,
     &    pcbariver0, pcbariver1,
#ifdef PCBAREFRACTORY
     &    pcbarriver0, pcbarriver1,
#endif
#ifdef PCBPREFRACTORY
     &    pcbprriver0, pcbprriver1,
#endif        
     &    chl0, chl1, npp0, npp1, doc0, doc1,
C Fdoc is obsolete (hma, 06 jul 2015)
C     &    poc0, poc1, Fdoc0, Fdoc1, Fpoc0, Fpoc1,
     &    poc0, poc1, Fpoc0, Fpoc1,
     &    Fremin0, Fremin1
#ifdef FOODW_MODE     
     &    , NO30, NO31,phytoP0, phytoP1, zooP0, zoop1
#endif
      _RS wind0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS wind1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS ice0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS ice1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS radsw0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS radsw1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS pcbaconc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric concentration
      _RS pcbaconc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric concentration
      _RS pcbadepo0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric deposition
      _RS pcbadepo1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric deposition
      _RS pcbpdepo0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric deposition
      _RS pcbpdepo1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric deposition
c Remove PCBi tracer (hma, 15 Jul 2014)
c      _RS pcbidepo0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric deposition
c      _RS pcbidepo0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric deposition
c      _RS pcbidepo1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  ! atmospheric deposition
c      _RS pcbiriver0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) ! river inputs
c      _RS pcbiriver1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) ! river inputs
      _RS pcbpriver0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) ! river inputs
      _RS pcbpriver1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) ! river inputs
      _RS pcbariver0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) ! river inputs
      _RS pcbariver1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) ! river inputs
#ifdef PCBAREFRACTORY
      _RS pcbarriver0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS pcbarriver1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef PCBPREFRACTORY
      _RS pcbprriver0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nrmax)
      _RS pcbprriver1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nrmax) 
#endif           
      
      _RS chl0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS chl1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS npp0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS npp1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS doc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS doc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS poc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS poc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C Fdoc is obsolete (hma, 07 july 2015)
C      _RS Fdoc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RS Fdoc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fpoc0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fpoc1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fremin0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS Fremin1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
#ifdef FOODW_MODE      
      _RS NO30    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS NO31    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)      
      _RS phytoP0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RS phytoP1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,npmax)
      _RS zooP0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
      _RS zooP1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nzmax)
#endif      


      COMMON /PCB_FILENAMES/
C  PCB_windFile      :: file name of wind speeds
C  PCB_iceFile       :: file name of seaice fraction
C  PCB_forcingPeriod :: periodic forcing parameter specific for PCB (seconds)
C  PCB_forcingCycle  :: periodic forcing parameter specific for PCB (seconds)
C  PCBa_concFile     :: file name of PCBa concentration in atmosphere
C  PCBi_depoFile     :: file name of deposition flux from atmosphere - REMOVED (hma, 15 Jul 2014)

C The river files are just dummy files for now. They're actually PFC
C inputs just renamed to PCBs. PCB river inputs are negligible and will
C likely not be considered in the final PCB simulation. (hma, 14 Jul 2014)
C  PCBi_riverFile    :: file name of PCBi riverine runoff - REMOVED (hma, 15 Jul 2014)
C  PCBp_riverFile    :: file name of PCBp riverine runoff
C  PCBa_riverFile    :: file name of PCBa riverine runoff
C  PCBpR_riverFile   :: file name of PCBpR riverine runoff
C  PCBaR_riverFile   :: file name of PCBaR riverine runoff
C  PCBa_depoFile     :: file name of PCBa atmospheric deposition
C  PCBp_depoFile     :: file name of PCBp atmospheric deposition

C  radsw_file       :: file name of short-wave radiation
C  chl_file         :: file name of chl concentration
C  npp_file         :: file name of net primary production
C  doc_file         :: file name of dissolved organic carbon concentration
C  poc_file         :: file name of particualte organic carbon concentration
C  Fdoc_file        :: file name of sinking dissolved organic carbon
C  Fpoc_file        :: file name of sinking particualte organic carbon
C  Fremin_file      :: file name of particualte organic carbon remineralization rate
C  NO3_file         :: file name of NO3 concentration
C  phytoP_file      :: file name of phytoplankton concentration
C  zooP_file        :: file name of zooplankton concentration
     &        PCB_windFile, PCB_iceFile, radsw_File,
     &        PCB_forcingPeriod, PCB_forcingCycle,
     &        PCB_chemPeriod, PCB_chemCycle,
     &        PCB_emisPeriod, PCB_emisCycle,
c Remove PCBi tracer. (hma, 15 Jul 2014)
c     &        PCBa_concFile , PCBi_depoFile,
c     &        PCBi_riverFile, PCBp_riverFile,   
     &        PCBa_concFile  , 
     &        PCBp_riverFile , 
     &        PCBa_riverFile , 
#ifdef PCBAREFRACTORY
     &        PCBAR_riverFile,
#endif
#ifdef PCBPREFRACTORY
     &        PCBPR_riverFile,
#endif       
     &        PCBa_depoFile  , 
     &        PCBp_depoFile  ,
     &        chl_file       , npp_file ,
     &        doc_file       , poc_file ,
C Fdoc is obsolete (hma, 07 july 2015)
C     &        Fdoc_file      , Fpoc_file, Fremin_file 
     &        Fpoc_file, Fremin_file 
#ifdef FOODW_MODE     
     &        , NO3_file,phytoP_file, zooP_file
#endif
      CHARACTER*(MAX_LEN_FNAM) PCB_windFile
      CHARACTER*(MAX_LEN_FNAM) PCB_iceFile
      CHARACTER*(MAX_LEN_FNAM) PCBa_concFile
      CHARACTER*(MAX_LEN_FNAM) PCBp_riverFile
      CHARACTER*(MAX_LEN_FNAM) PCBa_riverFile
#ifdef PCBAREFRACTORY
      CHARACTER*(MAX_LEN_FNAM) PCBAR_riverFile
#endif
#ifdef PCBPREFRACTORY
      CHARACTER*(MAX_LEN_FNAM) PCBR_riverFile 
#endif              
      CHARACTER*(MAX_LEN_FNAM) PCBa_depoFile
      CHARACTER*(MAX_LEN_FNAM) PCBp_depoFile            
      CHARACTER*(MAX_LEN_FNAM) radsw_File
      CHARACTER*(MAX_LEN_FNAM) chl_File
      CHARACTER*(MAX_LEN_FNAM) npp_File
      CHARACTER*(MAX_LEN_FNAM) doc_File
      CHARACTER*(MAX_LEN_FNAM) poc_File
C Fdoc is obsolete (hma, 06 july 2015)
C      CHARACTER*(MAX_LEN_FNAM) Fdoc_File
      CHARACTER*(MAX_LEN_FNAM) Fpoc_File
      CHARACTER*(MAX_LEN_FNAM) Fremin_File
#ifdef FOODW_MODE      
      CHARACTER*(MAX_LEN_FNAM) NO3_File
      CHARACTER*(MAX_LEN_FNAM) phytoP_File
      CHARACTER*(MAX_LEN_FNAM) zooP_File
#endif
      _RL     PCB_forcingPeriod   ! data related to air-sea exchange (i.e., atmospheric
      _RL     PCB_forcingCycle    !   concentration and deposition)
      _RL     PCB_chemPeriod      ! data related to PCB chemistry and DARWIN
      _RL     PCB_chemCycle       !   plankton biogeochemistry
      _RL     PCB_emisPeriod      ! currently doesn't control anything, legacy
      _RL     PCB_emisCycle       !   code from Xianming's PFCs
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|



