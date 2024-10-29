!*****************************************************************************
!***  Subroutine for the calculation of axis transformation                ***
!***  --> Jacobi-Method (cyclic)                                         ***
!***      vgl. z.B. I.N. Bronstein, K.A. Semendjajew, "Taschenbuch der     ***
!***                Mathematik", 25. Auflage                               ***
!***  --> Input: inertia tensor (symmetric 3x3 matrix)                     ***
!***      Output: matrix of orthonormal eigenvectors                       ***
!*****************************************************************************

SUBROUTINE HATsym3x3(A)
  INTEGER, PARAMETER :: MaxIterationen=25
  INTEGER i,Iterationen,p,q,Zyklus
  REAL(8) C,COSphi,SINphi,Skalierung,TANphi
  REAL(8), DIMENSION(1:3,1:3) ::  A,T,X

  Skalierung = MAXVAL( MAXVAL(A,DIM=1) )
  IF (Skalierung.EQ.0.0) THEN ;  Skalierung = 1.0 ;  END IF
  A = A/Skalierung

  X = 0.0
  DO i=1,3,1
    X(i,i) = 1.0
  END DO

  Iterationen = 0
  DO WHILE (      ( A(1,2)**2+A(1,3)**2+A(2,3)**2 .GT. EPSILON(A) )&
            .AND. (Iterationen .LE. MaxIterationen) )
    DO Zyklus=1,3,1
      SELECT CASE(Zyklus)
        CASE(1) ;  p = 1 ;  q = 2
        CASE(2) ;  p = 1 ;  q = 3
        CASE(3) ;  p = 2 ;  q = 3
      END SELECT

      IF (A(p,q).NE.0.0) THEN
        C = (A(q,q)-A(p,p)) / (2.0*A(p,q))
        IF ( C.GE.SQRT(HUGE(C)) ) THEN
          TANphi = 1.0 / (2.0*C)
        ELSE
          TANphi = 1.0 / (C + SIGN(SQRT(1+C**2),C))
        END IF
        COSphi =  1.0 / SQRT(1+TANphi**2)
        SINphi = COSphi*TANphi

        T = 0.0
        DO i=1,3,1
          T(i,i) = 1.0
        END DO
        T(p,p) = COSphi
        T(q,q) = COSphi
        T(p,q) = SINphi
        T(q,p) = -SINphi

        A = MATMUL( MATMUL(TRANSPOSE(T),A) , T )
        ! Symmetrie erhalten --> oberes Dreiek spiegeln
        A(2,1) = A(1,2)
        A(3,1) = A(1,3)
        A(3,2) = A(2,3)

        X = MATMUL(X,T)
      END IF
    END DO

    Iterationen = Iterationen+1
  END DO

  A = X
END SUBROUTINE HATsym3x3
