C      FGHEVEN
C      F. Gogtas, G.G. Balint-Kurti and C.C. Marston,
C      QCPE Program No. 647, (1993).
C
C     ***************************************************************
C     This program solves one dimensional Schrodinger equation for
C     bound state eigenvalues and eigenfunctions corresponding to a
C     potential V(x).
C
C     The Fourier Grid Hamiltonian method is used in the program.
C     This method is fully described in:
C     C.C. Marston and G.G.Balint-Kurti, J.Chem. Phys.,91,3571(1989).
C     and
C     G.G. Balint-Kurti, R.N. Dixon and C.C. Marston, Internat. Rev. 
C     Phys. Chem.,11, 317 (1992).
C
C     The program uses an even number of grid points. The Hamiltonian
C     matrix elements are calculated using analytic formula described
C     in the second of the above references.
C
C     The analytical Hamiltonian expression given in the reference
C     contains a small error.  The formula should read: 
C
C                H(i,j) = {(h**2)/(4*m*(L**2)} *
C                          [(N-1)(N-2)/6 + (N/2)] + V(Xi),   if i=j
C
C                H(i,j) = {[(-1)**(i-j)] / m } *
C                          { h/[2*L*sin(pi*(i-j)/N)]}**2 ,     if i#j
C
C     The eigenvalues of the Hamiltonian matrix which lie below the
C     asymptotic (large x) value of V(x) correspond to the bound state
C     energies. The corresponding eigenvectors are the representation
C     of the bound state wavefunctions on a regular one dimensional
C     grid.
C
C     The user must supply a self contained subroutine VSUB(X,V)
C     which is called with a value of X and returns the value of the
C     potential in V.
C
C     As supplied the program uses a Morse potential with parameters
C     chosen 
C     to represent the HF molecule. 
C
C     The Program should work on any computer
C     Programming language used                    : Fortran77
C     Memory required to execute with typical data : 1 Mbyte
C
C     To run the program on a SUN workstation (or any Unix based system)
C     First type:
C        f77 -o fgheven fgheven.f
C     to create an executable program (fgheven).
C     Then type:
C        fgheven > fgheven.out &   
C     to run the program.  
C     
C     The dimensions of the arrays generally depend on the number of
C     grid
C     points used.  This number NX is set in a PARAMETER statement.
C     The following arrays represent  :
C
C        XA  : X-coordinates i.e. position on a grid.
C
C        AR  : The Hamiltonian Matrix
C         
C        WCH : The eigenvalues for respective energy levels.
C
C        ZR  : The eigenvectors (X,Y); where
C                                      X : Wavefunction.
C                                      Y : Energy level.       
C    The folowing data constans represent :
C        ZMA, ZMB : Mass of atoms.
C        R0       : Equilibrium bond separation
C        NPRIN    : 0 or 1 for eigenvalue or eigenvector respectively
C        NWRIT    : Number of eigenvales/eigenvectors to be printed
C        RMIN     : Starting point of grid
C        RMAX     : End point of grid
C        ZL       : Grid length
C        DX       : Grid spacings
C
C       All quantities in au (unless otherwise stated)
C     ****************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NX=128)
C
C.....Declare arrays.
      DIMENSION XA(NX),AR(NX,NX),WCH(NX),ZR(NX,NX) 
      DIMENSION FV1(NX),FV2(NX)
C
C.....Variable data input.
      DATA ZMA/1836.9822D0/
      DATA ZMB/34629.61319D0/
      DATA R0/1.7329D0/
      DATA PI/3.141592653589793D0/
      DATA NPRIN/1/
      DATA NWRIT/10/
C
C....Test that NX is even
      ITEST=MOD(NX,2)
      IF(ITEST.NE.0) THEN
        WRITE(6,*)' **** NX MUST BE EVEN-FATAL ERROR  ****'
        STOP
      END IF
C
C.....Set up grid
      WRITE(6,*)'Grid paremeters:'
      WRITE(6,*)'   Number of grid poins    = ',NX
      RMIN=R0*0.2D0
      RMAX=R0*3.0D0
      ZL=(RMAX-RMIN)
      WRITE(6,*)'   Grid length             = ',ZL
      DX=ZL/DFLOAT(NX)
      WRITE(6,*)'   Grid spacings           = ',DX
C.....Print mass of atoms A and B
      WRITE(6,*)'Mass of atoms:'
      WRITE(6,*)'   Atomic mass of atom A      = ',ZMA
      WRITE(6,*)'   Atomic mass of atom B      = ',ZMB
C.....Compute reduced mass.
      ZMU=ZMA*ZMB/(ZMA+ZMB)
      WRITE(6,*)'   Reduced mass of AB      = ',ZMU
C
C.....Compute contants:
      PSQ=PI*PI
      CONST1=PSQ/(ZMU*(ZL**2))
      NFACT1=(NX-1)*(NX-2)
      NFACT2=(NX-2)/2
      CONST2=CONST1*(DFLOAT(NFACT1)/6.0D0 +1.0D0+DFLOAT(NFACT2))
      DARG=PI/DFLOAT(NX)
C
C.....Now compute Hamiltonian matrix:
      X=RMIN
      DO 003 I=1,NX
         XA(I)=X
       DO 004 J=1,I
          IJD=(I-J)
        IF(IJD.EQ.0)THEN
          AR(I,J)=CONST2
        ELSE
          RATIO=1.0D0/DSIN(DARG*DFLOAT(IJD))
          AR(I,J)=((-1)**IJD)*CONST1*(RATIO**2)
        END IF
 004  CONTINUE
C.....Find the potential value at x
      CALL VSUB(X,VV)
C.....Add the potential value when the kronecker delta function
C.....equals to one, i.e. when I and J are equal
        AR(I,I)=AR(I,I)+VV
        X=X+DX
 003  CONTINUE
C.....Now fill out Hamiltonian matrix.
      DO 006 I=1,NX
       DO 007 J=1,I
        AR(J,I)=AR(I,J)
007    CONTINUE
006   CONTINUE
C.....Now call eigenvalue solver.
      call flush(6)
      CALL RS(NX,NX,AR,WCH,NPRIN,ZR,FV1,FV2,IERR)
C.....Print out eigenvalues and eigenvectors.
      WRITE(6,12)
12    FORMAT(/)
      WRITE(6,*)' The first',NWRIT,' energy levels for AB molecule '
      DO 009 I=1,NWRIT  
       WRITE(6,*)' ENERGY LEVEL NO ',I,' EIGENVALUE= ',WCH(I)
 009  CONTINUE
      WRITE(6,12)
      WRITE(6,*)' THE CORRESPONDING EIGENFUNCTIONS ARE:'
      DO 010 I=1,NWRIT  
       WRITE(6,*)' ENERGY LEVEL NO ',I,' EIGENVALUE= ',WCH(I)
       IF (NPRIN.EQ.1) THEN
       DO 011 J=1,NX
         WRITE(6,*)' R = ',XA(J),' WAVEFUNCTION=',ZR(J,I)
 011   CONTINUE
       ENDIF
 010  CONTINUE
      STOP
      END
C
C
      SUBROUTINE VSUB(X,V)
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
      DATA DHART/0.2250073497D0/
      DATA BETA/1.1741D0/
      DATA R0/1.7329D0/
C.....Set morse curve parameters.
      ARG=X-R0
      ARG=-BETA*ARG
      EF=DEXP(ARG)
      EF=1.0D0-EF
      EF=EF*EF
      V=DHART*EF
      RETURN
      END
C
C      
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C     ****************************************************************
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                     
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)         
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC MATRIX.
C
C     ON INPUT : 
C
C        N,M  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM 
C        DIMENSION STATEMENT,
C        N  IS THE ORDER OF THE MATRIX  A,
C        A  CONTAINS THE REAL SYMMETRIC MATRIX,
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT :
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C    *******************************************************************
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     ********** FIND EIGENVALUES ONLY **********
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
C     ********** LAST CARD OF RS **********
      END
C
C
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NM,N),D(N),E(N),E2(N)
C     ***************************************************************
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT :
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT :
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C 
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ****************************************************************
C
      DO 100 I = 1, N
  100 D(I) = A(I,I)
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- ********** 
      DO 300 II = 1, N 
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) ********** 
         DO 120 K = 1, L 
  120    SCALE = SCALE + DABS(A(I,K)) 
C
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290
C
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE 
            H = H + A(I,K) * A(I,K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = A(I,L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 270
         F = 0.0D0
C
        DO 240 J = 1, L 
            G = 0.0D0
C     ********** FORM ELEMENT OF A*U ********** 
            DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
C     ********** FORM ELEMENT OF P **********
  220       E(J) = G / H
            F = F + E(J) * A(I,J)
  240    CONTINUE
C
         H = F / (H + H)
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L 
            F = A(I,J)
            G = E(J) - H * F
            E(J) = G
C
            DO 260 K = 1, J
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE 
C
  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)
C
  290    H = D(I)
         D(I) = A(I,I)
         A(I,I) = H
  300 CONTINUE
C 
      RETURN
C     ********** LAST CARD OF TRED1 **********
      END
C
C
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NM,N),D(N),E(N),Z(NM,N) 
C     **************************************************************** 
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A 
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT :
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM 
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE 
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT :
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX, 
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C 
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX 
C          PRODUCED IN THE REDUCTION,
C
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
C 
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     *************************************************************
      DO 100 I = 1, N 
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE 
C
      IF (N .EQ. 1) GO TO 320
C     ********** FOR I=N STEP -1 UNTIL 2 DO -- ********** 
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(Z(I,K))
C
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
C
  140    DO 150 K = 1, L 
            Z(I,K) = Z(I,K) / SCALE 
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
C 
         F = Z(I,L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0D0
C 
         DO 240 J = 1, L 
            Z(J,I) = Z(I,J) / H
            G = 0.0D0
C     ********** FORM ELEMENT OF A*U ********** 
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     ********** FORM ELEMENT OF P **********
  220       E(J) = G / H 
            F = F + E(J) * Z(I,J)
  240    CONTINUE
C  
         HH = F / (H + H) 
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L 
            F = Z(I,J) 
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
C 
  290    D(I) = H 
  300 CONTINUE
C 
  320 D(1) = 0.0D0
      E(1) = 0.0D0
C     ********** ACCUMULATION OF TRANSFORMATION MATRICES **********
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0D0) GO TO 380
C
         DO 360 J = 1, L
            G = 0.0D0
C 
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0D0
         IF (L .LT. 1) GO TO 500
C
         DO 400 J = 1, L
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  400    CONTINUE
C                
  500 CONTINUE
C
      RETURN 
C     ********** LAST CARD OF TRED2 **********
      END
C
C
      SUBROUTINE TQLRAT(N,D,E2,IERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION D(N),E2(N) 
      REAL*8 MACHEP
C     *****************************************************************
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C 
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC 
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT :
C
C        N IS THE ORDER OF THE MATRIX,
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX'
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT :   
C 
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES, 
C
C        IERR IS SET TO 
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     **************************************************************
C 
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C 
C
      MACHEP = 2.D0**(-26)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C 
      F = 0.0D0
      B = 0.0D0
      C = 0.0D0
      E2(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
         H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))
         IF (B .GT. H) GO TO 105
         B = H                   
         C = B * B
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
  105    DO 110 M = L, N 
            IF (E2(M) .LE. C) GO TO 120
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT 
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
         WRITE(6,*)' **** FATAL ERROR IN TQLRAT **** ' 
         WRITE(6,*)' **** FALLEN THROUGH BOTTOM OF LOOP 110 *** '
         STOP
C 
  120    IF (M .EQ. L) GO TO 210 
  130    IF (J .EQ. 30) GO TO 1000 
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1 
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = DSQRT(P*P+1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C 
         F = F + H
C     ********** RATIONAL QL TRANSFORMATION ********** 
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L 
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- ********** 
         DO 200 II = 1, MML 
            I = M - II
            P = G * H 
            R = P + E2(I)
            E2(I+1) = S * R 
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0D0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE.DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
C     ********** ORDER EIGENVALUES ********** 
         IF (L .EQ. 1) GO TO 250
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- *********
         DO 230 II = 2, L 
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L 
 1001 RETURN 
C     ********** LAST CARD OF TQLRAT **********
      END 
C 
C  
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION D(N),E(N),Z(NM,N)
      REAL*8 MACHEP
C     ************************************************************** 
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C 
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS 
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C 
C     ON INPUT :
C 
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL 
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C 
C        N IS THE ORDER OF THE MATRIX,
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C  
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS 
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C 
C      ON OUTPUT :
C 
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN  
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT 
C          UNORDERED FOR INDICES 1,2,...,IERR-1,   
C 
C        E HAS BEEN DESTROYED,   
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE, 
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ***************************************************************
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.

      MACHEP = 2.D0**(-26)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C 
      DO 100 I = 2, N 
  100 E(I-1) = E(I) 
C
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0 
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))
         IF (B .LT. H) B = H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
         DO 110 M = L, N
            IF (DABS(E(M)) .LE. B) GO TO 120
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT 
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE 
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L)) 
         R = DSQRT(P*P+1.0D0)
         D(L) = E(L) / (P + DSIGN(R,P))
         H = G - D(L)
         DO 140 I = L1, N
  140    D(I) = D(I) - H 
C
         F = F + H 
C     ********** QL TRANSFORMATION ********** 
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         MML = M - L 
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML 
            I = M - II
            G = C * E(I)
            H = C * P 
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150 
            C = E(I) / P
            R = DSQRT(C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R 
            C = 1.0D0 / R
            GO TO 160 
  150       C = P / E(I)
            R = DSQRT(C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     ********** FORM VECTOR **********
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
  200    CONTINUE
         E(L) = S * P
         D(L) = C * P 
         IF (DABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 300 II = 2, N 
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P 
  280    CONTINUE
  300 CONTINUE
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ********** 
 1000 IERR = L 
 1001 RETURN
C     ********** LAST CARD OF TQL2 **********
      END
