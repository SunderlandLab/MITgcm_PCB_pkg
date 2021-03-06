C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C                                                                  
C Sparse Jacobian Data Structures File                             
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
C File                 : chem_hg_JacobianSP.f                      
C Time                 : Thu Mar 27 11:36:37 2014                  
C Working directory    : /n/home05/yxzhang/kpp/kpp-2.2.3/mehg      
C Equation file        : chem_hg.kpp                               
C Output root filename : chem_hg                                   
C                                                                  
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




C Sparse Jacobian Data                                             

      BLOCK DATA JACOBIAN_SPARSE_DATA

      INCLUDE 'chem_hg_Sparse.h'
      INTEGER i

      DATA LU_IROW /
     *  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,
     *  7,  7,  8,  8,  9,  9, 10, 10, 11, 11, 12, 12,
     * 12, 13, 13, 13, 14, 14, 14 / 


      DATA LU_ICOL /
     *  1, 13,  2, 13,  3, 14,  4, 14,  5, 13,  6, 14,
     *  7, 12,  8, 12,  9, 12, 10, 11, 11, 12, 11, 12,
     * 14, 12, 13, 14, 12, 13, 14 / 


      DATA LU_CROW /
     *  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23,
     * 26, 29, 32 / 


      DATA LU_DIAG /
     *  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 24,
     * 27, 31, 32 / 


      END


