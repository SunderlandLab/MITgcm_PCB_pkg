C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C                                                                  
C Initialization File                                              
C                                                                  
C Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  
C       (http://www.cs.vt.edu/~asandu/Software/KPP)                
C KPP is distributed under GPL, the general public licence         
C       (http://www.gnu.org/copyleft/gpl.html)                     
C (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           
C (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            
C     With important contributions from:                           
C        M. Damian, Villanova University, USA                      
C        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
C                                                                  
C File                 : chem_hg_Initialize.f                      
C Time                 : Thu Mar 27 11:36:36 2014                  
C Working directory    : /n/home05/yxzhang/kpp/kpp-2.2.3/mehg      
C Equation file        : chem_hg.kpp                               
C Output root filename : chem_hg                                   
C                                                                  
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C                                                                  
C Initialize - function to initialize concentrations               
C   Arguments :                                                    
C                                                                  
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE Initialize ( )

      IMPLICIT NONE
      INCLUDE 'chem_hg_Parameters.h'


      INCLUDE 'chem_hg_Global.h'
      INCLUDE 'chem_hg_Rates.h'
      INTEGER i
      REAL*8 x

      CFACTOR = 1.000000d+00

      x = (0.)*CFACTOR
      DO i = 1, NVAR
        VAR(i) = x
      END DO

      x = (0.)*CFACTOR
      DO i = 1, NFIX
        FIX(i) = x
      END DO

C constant rate coefficients                                       
C END constant rate coefficients                                   

C concentrations
C end concentrations

C INLINED initializations                                          
        Autonomous = 1 
        STEPSTART=STEPMIN
C End INLINED initializations                                      

      RETURN
      END

C End of Initialize function                                       
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


