C $Header: HG_SIZE.h,v 1.0 2013/07/11 13:06:31 yanxu Exp $
C $Name:  $

c HG_SIZE.h
c description: Size specification for hg model
c
c
c npmax = no of "functional groups" of phytoplankton
c nzmax = no of "functional groups" of zooplankton
c
c
#ifdef FOODW_MODE      
         INTEGER npmax
         INTEGER nzmax
         PARAMETER(npmax=3,nzmax=2)
#endif

      INTEGER iHg0,iHg2,iHgP
#ifdef METHY_MODE
      INTEGER iMMHg,iMMHgP
      INTEGER iDMHg
#endif
#ifdef FOODW_MODE      
C iZoo          :: index of first zooplankton
C iPhy          :: index of first phytoplankton
C remember to bring the fields in data.ptracers in the right order !
      INTEGER iPhy
      INTEGER iZoo
#endif

#ifdef HG2REFRACTORY
      INTEGER iHg2R
#endif

#ifdef HGPREFRACTORY
      INTEGER nrmax
      PARAMETER(nrmax=1)
      INTEGER iHgPR
#endif

      PARAMETER (iHg0  = 1)
      PARAMETER (iHg2  = 2)
      PARAMETER (iHgP  = 3)

#ifdef HG2REFRACTORY
      PARAMETER (iHg2R = 4)
#endif

#ifdef HGPREFRACTORY
#ifdef HG2REFRACTORY
      PARAMETER (iHgPR  = 5)
#else
      PARAMETER (iHgPR  = 4)
#endif
#endif
      
#ifdef METHY_MODE
#ifdef HGPREFRACTORY
      PARAMETER (iMMHg  = iHgPR + nrmax)
#else
#ifdef HG2REFRACTORY
      PARAMETER (iMMHg  = iHg2R + 1)
#else
      PARAMETER (iMMHg  = 4)
#endif
#endif           

      PARAMETER (iMMHgP =iMMHg+1)
      PARAMETER (iDMHg  =iMMHgP+1) 
#ifdef FOODW_MODE                              
      PARAMETER (iPhy  =iDMHg+1)
      PARAMETER (iZoo  =iPhy  +npmax)
#endif
#endif
