      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NX=3200)
C
C.....Declare arrays.
      DIMENSION XA(NX),AR(NX,NX),WCH(NX),ZR(NX,NX) 
      DIMENSION FV1(NX),FV2(NX)
      DATA PI/3.141592653589793D0/
      DATA NPRIN/0/
      DATA NWRIT/20/
      INTEGER :: L
      LOGICAL :: MORLJ  


      ZMA   = 242271.817297
      ZMB   = 242271.817297
      L     = 1000
      RMAX  = 100.0D0
      RMIN  = 5.0D0
      C6    = 6870.1
      C8    = 0.0D0
      C10   = 0.0D0
      C12   = 709538221.031
      XMIN  = 6.5
      XMAX  = 25.0
      DX    = (XMAX-XMIN)/DBLE(NX)
      ZMU=ZMA*ZMB/(ZMA+ZMB)
      DHART = C6**2/(4.0*C12)
      R0    = (C12/C6)**(1.0/6.0)*2.0**(1.0/6.0)
      BETA  = 6.0/R0

      DO I=1,NX
       X = XMIN + (I-1)*DX
       CALL VSUB_M(X,V1,R0,BETA,DHART,L,ZMU)
       CALL VSUB_LJ(X,V2,R0,C6,C8,C10,C12,L,ZMU)
       CALL VSUB_MLR3(X,V3,L,ZMU)
       WRITE(*,*)X,V2,V3
      ENDDO
      END
C
      SUBROUTINE VSUB_MLR3(X,V,L,MU)
C     **************************************************************
C.....SUBROUTINE TO COMPUTE POTENTIAL ENERGY AT SEPERATION X (IN   
C.....ATOMIC UNITS). THE VALUE OF THE POTENTIAL IS RETURNED IN V,  
C.....AGAIN IN ATOMIC UNITS.THE POTENTIAL ENERGY IS ASSUMED TO     
C.....TO BE A MLR3 POTNETIAL AS DEFINDED IN JCP 132 094105 2010    
C.....FOR GROUND STATE CEASIUM DIMERS                              
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION PHI(19)

      BRtoAng = 0.52917721080000002
      R = X*BrtoAng
!     scipy.constants.physical_constants["Bohr radius"][0]*1.0E10
      CM_1toHar = 4.5563352527603877e-06
!     De = De*100.0/(scipy.constants.Rydberg*2.0)

      PHI(1)  = -2.48715027
      PHI(2)  = -0.7250527
      PHI(3)  = -1.988236
      PHI(4)  = -0.671755
      PHI(5)  = -1.202501
      PHI(6)  = -0.020151
      PHI(7)  = -0.5414
      PHI(8)  =  0.2298
      PHI(9)  =  1.3964
      PHI(10) =  0.687
      PHI(11) = -8.655
      PHI(12) =  1.73
      PHI(13) =  32.2
      PHI(14) = -2.66
      PHI(15) = -61.0
      PHI(16) =  6.1
      PHI(17) =  65.6
      PHI(18) = -2.0
      PHI(19) = -28.0
 
      Re      = 4.6479723
      De      = 3649.847
      C6      = 3.315E7
      C8      = 1.29962E9
      C12     = 5.136E10
      MCS     = 132.905451933
      RMIN    = 3.1
      RMAX    = 20.1
      DR      = 0.001
      RREF    = 5.47
      P       = 5.0
      M       = 6.0
      Q       = 4.0
      A       = 1.79

      ULRR    = C6/R**6 + C8/R**8 + C12/R**12
      ULRRe   = C6/Re**6 + C8/Re**8 + C12/Re**12
      YM      = (R**M-RREF**M)/(R**M+RREF**M)
      YQ      = (R**Q-RREF**Q)/(R**Q+RREF**Q)


      PHIMLR3 = 0.0
      DO I=1,19
       PHIMLR3 = PHIMLR3 + PHI(I)*YQ!**I
      ENDDO
      PHIMLR3 = PHIMLR3*(1.0-YM)
      PHIINF  = 0.5*((2.0*De*Re**6)/(C6))
      PHIMLR3 = PHIMLR3 + YM*PHIINF
      YPA     = (X**P-Re**P)/(X**P+A*Re**P)
      V       = De*(1.0-(ULRR/ULRRe)*EXP(-PHIMLR3*YPA))**2 - De
      V = V*4.5563352527603877e-06
      V = V + DBLE(L)*(DBLE(L)+1.0)/(2.0*MU*X**2)
      RETURN 
      END
  
      SUBROUTINE VSUB_M(X,V,R0,BETA,DHART,L,MU)
C     **************************************************************
C.....SUBROUTINE TO COMPUTE POTENTIAL ENERGY AT SEPERATION X (IN   
C.....ATOMIC UNITS). THE VALUE OF THE POTENTIAL IS RETURNED IN V,  
C.....AGAIN IN ATOMIC UNITS.THE POTENTIAL ENERGY IS ASSUMED TO     
C..... BE A MORSE CURVE:V(X) = DHART * (EXP(-BETA*(X-R0)) - 1)**2  
C..... THE FOLLOWING DATA CONSTANTS REPRESENT :                    
C.....    DHART : THE DISSOCIATION ENERGY FOR THE GIVEN MOLECULE   
C.....    BETA  : THE EXPONENTIAL PARAMETER                        
C.....    R0    : THE BOND LENGTH AT EQUILIBRIUM FOR GIVEN MOLECULE
C     **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION :: MU
      INTEGER :: L
C      DATA DHART/0.2250073497D0/
C      DATA BETA/1.1741D0/
C      DATA R0/1.7329D0/
C.....Set morse curve parameters.
      ARG=X-R0
      ARG=-BETA*ARG
      EF=DEXP(ARG)
      EF=1.0D0-EF
      EF=EF*EF
      V=DHART*EF-DHART
      V = V + DBLE(L)*(DBLE(L)+1.0)/(2.0*MU*X**2)
      RETURN
      END
CC
      SUBROUTINE VSUB_LJ(X,V,R0,C6,C8,C10,C12,L,MU)
C     **************************************************************
C.....SUBROUTINE TO COMPUTE POTENTIAL ENERGY AT SEPERATION X (IN   
C.....ATOMIC UNITS). THE VALUE OF THE POTENTIAL IS RETURNED IN V,  
C.....AGAIN IN ATOMIC UNITS.THE POTENTIAL ENERGY IS ASSUMED TO     
C..... BE A MORSE CURVE:V(X) = DHART * (EXP(-BETA*(X-R0)) - 1)**2  
C..... THE FOLLOWING DATA CONSTANTS REPRESENT :                    
C.....    DHART : THE DISSOCIATION ENERGY FOR THE GIVEN MOLECULE   
C.....    BETA  : THE EXPONENTIAL PARAMETER                        
C.....    R0    : THE BOND LENGTH AT EQUILIBRIUM FOR GIVEN MOLECULE
C     **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION :: MU
      INTEGER :: L
C      DATA DHART/0.2250073497D0/
C      DATA BETA/1.1741D0/
C      DATA R0/1.7329D0/
C.....Set morse curve parameters.
      V =C12/X**12-C10/X**6-C8/X**6-C6/X**6
      V = V + DBLE(L)*(DBLE(L)+1.0)/(2.0*MU*X**2)
      RETURN
      END
C
C      
