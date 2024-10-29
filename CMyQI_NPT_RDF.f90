! gfortran -O3 -fdefault-integer-8 CMyQI_NPT_RDF.f90

  PROGRAM ZellFIQA
  IMPLICIT NONE
  CHARACTER(26), PARAMETER :: Version = '   (Version April 2021)'

  ! variables from input file
  INTEGER AnzKomp, BilderMax, EquiNVTSchritte, EquiNPTSchritte, FilmAb, MaxAnzLadungen, MaxAnzDipole,&
          MaxAnzQuadrupole, MaxAnzSites, N, Output, OutputVis, SimSchritte
  INTEGER, ALLOCATABLE, DIMENSION(:) :: AnzLadungen, AnzDipole, AnzQuadrupole, AnzSites,&
                                        NKomp
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Farbe
  LOGICAL Film, NPT
  REAL(8) dt, epsilonRF, rc, rho, T, P, VMasse
  REAL(8), DIMENSION(1:3) :: SimBoxRatio
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: absLadung, absMy, absQ, etaLB, epsilonSite,&
                                          IBodyDummy, mLadung, mSite, sigmaSite, xiLB, MyGes
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: eMyBody, eQBody, rLadungBody, rDipolBody,&
                                            rQuadrupolBody, rSiteBody

  ! variables for programm
  CHARACTER(128) Dateiname, Programmname, Text, xMal
  INTEGER AnzZellNachbarn, BildNr, dZN, i, IARGC, IntSchritt, iZ, j, jZ,&
          MaxIterationen, NmaxZelle, Richtung, Si, Sj, x, xi, xj,&
          y, yi, yj, z, zi, zj, ZT, k, WWCount
  INTEGER, DIMENSION(1:3) :: kfzZellen, Zellen
  INTEGER, ALLOCATABLE, DIMENSION(:) :: FHG, Komp
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ZellNachbar, NSchale
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: NZelle
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: Liste
  REAL(8) AtomareMasseDim, AvogadroDim, betaRot, betaRotMax, betaRotMin,&
          betaTrans, betaTransMax, betaTransMin, BoltzmannDim,&
	  Cos2Ti, Cos2Tj, CosGij, CosTi, CosTj, CosTi3, CosTj3,&
	  delta, DKorr, DKorr_n1, dr2, Drittel,&
          epsRFInvrc3, EChargeDim, EpsilonVakuum,&
	  Invdr1, Invdr2, drMol2, dtInv2, dtInv4,&
          epsilonRefDim, mRefDim, E2, My2, MyFaktor, MyQFaktor, UMyRF,&
          MySelbstTerm, PartialRijInvdr1, PartialTiInvdr1, PartialTjInvdr1,&
          PartialGij, pi, QFaktor, qKorr, rc2, LJTerm6, LJTerm12, RFFaktor,&
          sigmaRefDim, sumDruck, sumIw2, sumMv2, sumUpot, sumUpotKorrMy,&
          sumUpotMy, sumUpotE, sumUpotMyQ, sumUpotQ, sumVirialE, sumVirialKorrMy, sumVirialMy,&
          sumVirialMyQ, sumVirialQ, tau, tau1, tau2, TermKlammer, TermU, TermU2,&
          UpotKorrLJ, UPotE, UpotKorrMy, UpotMy, UpotMyQ, UpotQ, UpotSite,&
          VirialKorrLJ, VirialE, VirialKorrMy, ViriralE, VirialMy, VirialMyQ, VirialQ,&
          VirialSite, vKorr, wKorr, Zufall, V0, V1, V2, PAktuell, sumrho, MaxSiteCOM
  REAL(8), DIMENSION(1:3) :: dD, D_n1, dr, drE, drMol, dv, eiXej, eXrij,&
                             Fij, InvLSimBox, InvLZelle, LSimBox, LSimBoxAlt, LZelle,&
                             Mij, Mji, wBody, ZweiInvLSimBox, MRF
  REAL(8), DIMENSION(0:3) :: qKonst, qAlt
  REAL(8), ALLOCATABLE, DIMENSION(:) :: dtInv2mKomp, mKomp, VolumenSchale
  REAL(8), DIMENSION(1:3,1:3) :: A,AT
  REAL(8), DIMENSION(1:4,1:3) :: Q4x3
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: D, Disp, F, IBody, InvIBody, M, q, r, v, F_v
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: drDipol, drQuadrupol, drSite, drLadung,&
                                            eMy, MyRFSpace, eQ, FDipol, FQuadrupol, FSite, FLadung, FSite_v
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: epsilon24, sigma2, Invsigma

  ! variables for spherical lrc
  INTEGER :: NShells, Ki, Kj, MeanIndex, TotalStep, NSMean
  REAL(8) :: drShells, drShells05, VierDPi, factorF, UCorrShells, FCorrTemp, FCorrShells, SumEkin, factorU, &
             UCorrTemp, UCorrSum, deltaShells, realk, UpotTot, UPotMin, UPotMax, sigma6, sumv, &
             rcmax, rhol, rhov, Dmin, Dmax, R0, D0, r10, r90, rdash, rlow, &
             rlowInv, rlowInv2, rdashInv, rdashInv2, factorP, PNCorrTemp, PNCorrShells, PTCorrShells, &
             PTCorrTemp, VirialTemp, FijTemp, vsquare, dr_nor, drMol_nor, rlow2, rdash2
  REAL(8), DIMENSION(1:3) :: SystemCOM, SystemCenter, TestVir, ksi_i, ksi_j, rss_tan, rcom_tan
  REAL(8), ALLOCATABLE, DIMENSION(:) :: VShells, RShells, ksi, rhoShellsAve, RShells2, ksi2, &
  PShells_N, PShells_T, PNShells_Mean, PTShells_Mean, UShells_Mean, FShells_Mean
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: rhoShells, rhoShellsTemp, rhoShellsT
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: rhoShellsMean
  INTEGER, ALLOCATABLE, DIMENSION(:) :: PartShells, lowerS, upperS, interS
  CHARACTER(128) :: InputFile


!*****************************************************************************
!***  Prelude: read data; "zero-th" step                                   ***
!*****************************************************************************
  INCLUDE 'Preparations.f90'
  INCLUDE 'ReadConfig.f90'
!*******   end of  prelude   ************************************************

!*****************************************************************************
!***  One time step with dt (translation/rotation)                         ***
!***  1) Leapfrog Integrator (Part I)                                      ***
!***     with FIQA (Fincham's Implicit Quaternion Algorithm                ***
!***  2) Assign molecules to cells                                         ***
!***  3) Calculation of forces and momenta for present configuration       ***
!***  4) Leapfrog Integrator (Part II)                                     ***
!***     with FIQA (Fincham's Implicit Quaternion Algorithm                ***
!*****************************************************************************

  ! Initialisations
  betaRotMax = 1.0
  betaRotMin = 1.0
  betaTransMax = 1.0
  betaTransMin = 1.0
  Disp = 0.0
  MaxIterationen = 0
  sumDruck = 0.0
  sumrho   = 0.0
  sumUpot = 0.0
  sumUpotE= 0.0
  sumUpotMy = 0.0
  sumUpotKorrMy = 0.0
  sumUpotMyQ = 0.0
  sumUpotQ = 0.0
  sumVirialE = 0.0
  sumVirialMy = 0.0
  sumVirialKorrMy = 0.0
  sumVirialMyQ = 0.0
  sumVirialQ = 0.0
  TotalStep = -1


  DO IntSchritt=(-EquiNVTSchritte-EquiNPTSchritte+1),SimSchritte,1   ! time step = t

    TotalStep = TotalStep + 1

    ! Align system center and spherical phase
    SystemCOM = 0.0
    DO i=1, N, 1
      SystemCOM = SystemCOM + r(:,i)
    END DO
    SystemCOM = SystemCOM/N
    DO i=1, N, 1
      r(:,i) = r(:,i) - SystemCOM + SystemCenter
      r(:,i) = r(:,i) - LSimBox*FLOOR(0.5*ZweiInvLSimBox*r(:,i))
    END DO

    ! 1) Leapfrog Integrator (Part I)
    !    with FIQA (Fincham's Implicit Quaternion Algorithm)
    ! 1a) Translation
    vKorr = 2.0 - 1.0/betaTrans
    DO i=1,N,1
      v(:,i) = vKorr*v(:,i) + dtInv2mKomp(Komp(i))*F(:,i)
      dr = dt*v(:,i)
      r(:,i) = r(:,i)+dr
      Disp(:,i) = Disp(:,i)+dr
      r(:,i) = r(:,i) - LSimBox*FLOOR((r(:,i)+delta)*InvLSimBox)   ! PBC
    END DO

    ! 1b) Rotation
    DKorr = 2.0 - 1.0/betaRot
    DKorr_n1 = 3.0 - 2.0/betaRot
    DO i=1,N,1
      A(1,1) = q(0,i)*q(0,i) + q(1,i)*q(1,i) - q(2,i)*q(2,i) - q(3,i)*q(3,i)
      A(1,2) = 2.0 * ( q(1,i)*q(2,i) + q(0,i)*q(3,i) )
      A(1,3) = 2.0 * ( q(1,i)*q(3,i) - q(0,i)*q(2,i) )
      A(2,1) = 2.0 * ( q(1,i)*q(2,i) - q(0,i)*q(3,i) )
      A(2,2) = q(0,i)*q(0,i) - q(1,i)*q(1,i) + q(2,i)*q(2,i) - q(3,i)*q(3,i)
      A(2,3) = 2.0 * ( q(2,i)*q(3,i) + q(0,i)*q(1,i) )
      A(3,1) = 2.0 * ( q(1,i)*q(3,i) + q(0,i)*q(2,i) )
      A(3,2) = 2.0 * ( q(2,i)*q(3,i) - q(0,i)*q(1,i) )
      A(3,3) = q(0,i)*q(0,i) - q(1,i)*q(1,i) - q(2,i)*q(2,i) + q(3,i)*q(3,i)
      wBody =  InvIBody(:,Komp(i)) * MATMUL(A,D(:,i))
      Q4x3(:,1) = (/-q(1,i),+q(0,i),+q(3,i),-q(2,i)/)
      Q4x3(:,2) = (/-q(2,i),-q(3,i),+q(0,i),+q(1,i)/)
      Q4x3(:,3) = (/-q(3,i),+q(2,i),-q(1,i),+q(0,i)/)
      qKonst = q(:,i) + dtInv4*MATMUL(Q4x3,wBody)

      D_n1 = DKorr_n1*D(:,i) + dt*M(:,i)
      qAlt = q(:,i)
      q(:,i) = q(:,i) + dtInv2*MATMUL(Q4x3,wBody)   ! initial guess

      x = 0
      DO WHILE ( MAXVAL(ABS(q(:,i)-qAlt)) .GT. 100.0*EPSILON(q) )
        qAlt = q(:,i)
        A(1,1) = q(0,i)*q(0,i) + q(1,i)*q(1,i) - q(2,i)*q(2,i) - q(3,i)*q(3,i)
        A(1,2) = 2.0 * ( q(1,i)*q(2,i) + q(0,i)*q(3,i) )
        A(1,3) = 2.0 * ( q(1,i)*q(3,i) - q(0,i)*q(2,i) )
        A(2,1) = 2.0 * ( q(1,i)*q(2,i) - q(0,i)*q(3,i) )
        A(2,2) = q(0,i)*q(0,i) - q(1,i)*q(1,i) + q(2,i)*q(2,i) - q(3,i)*q(3,i)
        A(2,3) = 2.0 * ( q(2,i)*q(3,i) + q(0,i)*q(1,i) )
        A(3,1) = 2.0 * ( q(1,i)*q(3,i) + q(0,i)*q(2,i) )
        A(3,2) = 2.0 * ( q(2,i)*q(3,i) - q(0,i)*q(1,i) )
        A(3,3) = q(0,i)*q(0,i) - q(1,i)*q(1,i) - q(2,i)*q(2,i) + q(3,i)*q(3,i)
        wBody =  InvIBody(:,Komp(i)) * MATMUL(A,D_n1)
        Q4x3(:,1) = (/-q(1,i),+q(0,i),+q(3,i),-q(2,i)/)
        Q4x3(:,2) = (/-q(2,i),-q(3,i),+q(0,i),+q(1,i)/)
        Q4x3(:,3) = (/-q(3,i),+q(2,i),-q(1,i),+q(0,i)/)
        q(:,i) = qKonst + dtInv4*MATMUL(Q4x3,wBody)
        x = x+1
      END DO
      qKorr = 1.0 / SQRT(DOT_PRODUCT(q(:,i),q(:,i)))
      q(:,i) = qKorr*q(:,i)
      MaxIterationen = MAX( MaxIterationen , x )

      D(:,i) = DKorr*D(:,i) + dtInv2*M(:,i)
    END DO

    ! 1c) Volume
    IF ((NPT).AND.(IntSchritt.GE.(-1*EquiNPTSchritte))) THEN
      LSimBoxAlt=LSimBox
      V1 = V1 + V2
      V0 = V0 + V1
      rho = REAL(N)/V0
      LSimBox(1) = ( (N*SimBoxRatio(1)**2.0)&
                    /(rho*SimBoxRatio(2)*SimBoxRatio(3)) ) ** (1.0/3.0)
      LSimBox(2) = ( (N*SimBoxRatio(2)**2.0)&
                    /(rho*SimBoxRatio(1)*SimBoxRatio(3)) ) ** (1.0/3.0)
      LSimBox(3) = ( (N*SimBoxRatio(3)**2.0)&
                    /(rho*SimBoxRatio(1)*SimBoxRatio(2)) ) ** (1.0/3.0)
      DO i=1,N,1
        r(:,i) = r(:,i) * LSimBox / LSimBoxAlt
      END DO
      LZelle = LSimBox/Zellen
      IF (MINVAL(LZelle).LT.rc) THEN
        OPEN(5,file=TRIM(Dateiname)//'_Erg.dat',POSITION='APPEND',STATUS='OLD')
          WRITE(5,*) 'Error: Side-length of sub-cells lower than rc.'
	  WRITE(5,*) 'Suggestion: Reduce initial density or choose lower rc.'
        CLOSE(5)
	STOP
      END IF
      InvLZelle = 1.0 / LZelle
      InvLSimBox = 1.0 / LSimBox
      ZweiInvLSimBox = 2.0 * InvLSimBox
    END IF

    ! 2) Assign molecules to cells
    NZelle = 0
    DO i=1,N,1
      xi = INT(r(1,i)*InvLZelle(1))
      yi = INT(r(2,i)*InvLZelle(2))
      zi = INT(r(3,i)*InvLZelle(3))
      NZelle(xi,yi,zi) = NZelle(xi,yi,zi) + 1
      IF (NZelle(xi,yi,zi).GT.NmaxZelle) THEN
        DEALLOCATE(Liste)
        NmaxZelle = NmaxZelle + 1
        ALLOCATE( Liste(1:NmaxZelle,0:Zellen(1)-1,&
                        0:Zellen(2)-1,0:Zellen(3)-1) )
        NZelle = 0
        DO j=1,i,1
          xj = INT(r(1,j)*InvLZelle(1))
          yj = INT(r(2,j)*InvLZelle(2))
          zj = INT(r(3,j)*InvLZelle(3))
          NZelle(xj,yj,zj) = NZelle(xj,yj,zj) + 1
          Liste(NZelle(xj,yj,zj),xj,yj,zj) = j
        END DO
      ELSE
        Liste(NZelle(xi,yi,zi),xi,yi,zi) = i
      END IF
    END DO

    ! 3) Calculation of forces and momenta for time step (t+dt)
    ! Preparation: transform Site,My,Q-vectors and eMy,eQ
    !               ins raumfeste System
    DO i=1,N,1
      AT(1,1) = q(0,i)*q(0,i) + q(1,i)*q(1,i) - q(2,i)*q(2,i) - q(3,i)*q(3,i)
      AT(1,2) = 2.0 * ( q(1,i)*q(2,i) - q(0,i)*q(3,i) )
      AT(1,3) = 2.0 * ( q(1,i)*q(3,i) + q(0,i)*q(2,i) )
      AT(2,1) = 2.0 * ( q(1,i)*q(2,i) + q(0,i)*q(3,i) )
      AT(2,2) = q(0,i)*q(0,i) - q(1,i)*q(1,i) + q(2,i)*q(2,i) - q(3,i)*q(3,i)
      AT(2,3) = 2.0 * ( q(2,i)*q(3,i) - q(0,i)*q(1,i) )
      AT(3,1) = 2.0 * ( q(1,i)*q(3,i) - q(0,i)*q(2,i) )
      AT(3,2) = 2.0 * ( q(2,i)*q(3,i) + q(0,i)*q(1,i) )
      AT(3,3) = q(0,i)*q(0,i) - q(1,i)*q(1,i) - q(2,i)*q(2,i) + q(3,i)*q(3,i)
      DO Si=1,AnzSites(Komp(i)),1
        drSite(:,Si,i) = MATMUL( AT,rSiteBody(:,Si,Komp(i)) )
      END DO
      DO Si=1,AnzLadungen(Komp(i)),1
        drLadung(:,Si,i) = MATMUL( AT,rLadungBody(:,Si,Komp(i)) )
      END DO
      DO Si=1,AnzDipole(Komp(i)),1
        drDipol(:,Si,i) = MATMUL( AT,rDipolBody(:,Si,Komp(i)) )
        eMy(:,Si,i) = MATMUL( AT,eMyBody(:,Si,Komp(i)) )
      END DO
      DO Si=1,AnzQuadrupole(Komp(i)),1
        drQuadrupol(:,Si,i) = MATMUL( AT,rQuadrupolBody(:,Si,Komp(i)) )
        eQ(:,Si,i) = MATMUL( AT,eQBody(:,Si,Komp(i)) )
      END DO
      MyRFSpace(:,Komp(i),i) = MATMUL( AT,MyGes(:,Komp(i)) )
    END DO

    ! Calculation of pairwise interactions (--> Site,My,Q-Kräfte)
    F = 0.0
    M = 0.0

    FSite = 0.0
    FLadung = 0.0
    FDipol = 0.0
    FQuadrupol = 0.0
    UMyRF = 0.0
    UpotSite = 0.0
    UpotE = 0.0
    UpotMy = 0.0
    UpotMyQ = 0.0
    UpotQ = 0.0
    VirialSite = 0.0
    VirialE = 0.0
    VirialMy = 0.0
    VirialMyQ = 0.0
    VirialQ = 0.0
    WWCount=0

    INCLUDE 'ShellInit.f90'

    DO zi=0,Zellen(3)-1,1
      DO yi=0,Zellen(2)-1,1
        DO xi=0,Zellen(1)-1,1
          IF (    (MIN(xi,yi,zi).LT.ZT)&
              .OR.(MAXVAL( (/xi,yi,zi/)/(Zellen-ZT) ).EQ.1) ) THEN

    ! block for outer part --> with MIC
    DO iZ=1,NZelle(xi,yi,zi),1
      i = Liste(iZ,xi,yi,zi)

      DO jZ=(iZ+1),NZelle(xi,yi,zi),1
        j = Liste(jZ,xi,yi,zi)
        drMol = r(:,i) - r(:,j)
        drMol2 = DOT_PRODUCT(drMol,drMol)
        IF (drMol2.LT.rc2) THEN
          INCLUDE 'Force_ij.f90'
        END IF
      END DO   ! jZ (same cell)

      DO dZN=1,AnzZellNachbarn,1
        x = xi+ZellNachbar(1,dZN)
        y = yi+ZellNachbar(2,dZN)
        z = zi+ZellNachbar(3,dZN)
        xj = x - Zellen(1)*FLOOR((REAL(x,KIND=8))/REAL(Zellen(1),KIND=8))
        yj = y - Zellen(2)*FLOOR((REAL(y,KIND=8))/REAL(Zellen(2),KIND=8))
        zj = z - Zellen(3)*FLOOR((REAL(z,KIND=8))/REAL(Zellen(3),KIND=8))
        DO jZ=1,NZelle(xj,yj,zj),1
          j = Liste(jZ,xj,yj,zj)
          drMol = r(:,i) - r(:,j)
          drMol = drMol - LSimBox*INT(ZweiInvLSimBox*drMol)   ! MIC
          drMol2 = DOT_PRODUCT(drMol,drMol)
          IF (drMol2.LT.rc2) THEN
            INCLUDE 'Force_ij.f90'
          END IF
        END DO   ! jZ (neighbour cell)
      END DO   ! dZN := increment for neighbour cell
    END DO   ! iZ

          ELSE   ! distinction between outer part (above), inner part (below)

    ! block for inner part --> ganz ohne MIC
    DO iZ=1,NZelle(xi,yi,zi),1
      i = Liste(iZ,xi,yi,zi)

      DO jZ=(iZ+1),NZelle(xi,yi,zi),1
        j = Liste(jZ,xi,yi,zi)
        drMol = r(:,i) - r(:,j)
        drMol2 = DOT_PRODUCT(drMol,drMol)
        IF (drMol2.LT.rc2) THEN
          INCLUDE 'Force_ij.f90'
        END IF
      END DO   ! jZ (same cell)

      DO dZN=1,AnzZellNachbarn,1
        xj = xi+ZellNachbar(1,dZN)
        yj = yi+ZellNachbar(2,dZN)
        zj = zi+ZellNachbar(3,dZN)
        DO jZ=1,NZelle(xj,yj,zj),1
          j = Liste(jZ,xj,yj,zj)
          drMol = r(:,i) - r(:,j)
          drMol2 = DOT_PRODUCT(drMol,drMol)
          IF (drMol2.LT.rc2) THEN
            INCLUDE 'Force_ij.f90'
          END IF
        END DO   ! jZ (neighbour cell)
      END DO   ! dZN := increment for neighbour cell
    END DO   ! iZ

          END IF   ! inner/outer part
        END DO   ! xi
      END DO   ! yi
    END DO   ! zi
    UpotSite = UpotSite/6.0   ! since in interaction-block epsilon24 instead of espilon4
    UMyRF=UMyRF*epsRFInvrc3
    UpotKorrMy = (UMyRF + MySelbstTerm) / N
    VirialKorrMy = (rho/N) * UMyRF

    ! Transform site-forces to resulting forces/momenta
    DO i=1,N,1
      DO Si=1,AnzSites(Komp(i)),1
        F(:,i) = F(:,i) + FSite(:,Si,i)
        M(1,i) = M(1,i) + drSite(2,Si,i)*FSite(3,Si,i)&
                        - drSite(3,Si,i)*FSite(2,Si,i)
        M(2,i) = M(2,i) + drSite(3,Si,i)*FSite(1,Si,i)&
                        - drSite(1,Si,i)*FSite(3,Si,i)
        M(3,i) = M(3,i) + drSite(1,Si,i)*FSite(2,Si,i)&
                        - drSite(2,Si,i)*FSite(1,Si,i)
      END DO
      DO Si=1,AnzLadungen(Komp(i)),1
        F(:,i) = F(:,i) + FLadung(:,Si,i)
        M(1,i) = M(1,i) + drLadung(2,Si,i)*FLadung(3,Si,i)&
                        - drLadung(3,Si,i)*FLadung(2,Si,i)
        M(2,i) = M(2,i) + drLadung(3,Si,i)*FLadung(1,Si,i)&
                        - drLadung(1,Si,i)*FLadung(3,Si,i)
        M(3,i) = M(3,i) + drLadung(1,Si,i)*FLadung(2,Si,i)&
                        - drLadung(2,Si,i)*FLadung(1,Si,i)
      END DO
      DO Si=1,AnzDipole(Komp(i)),1
        F(:,i) = F(:,i) + FDipol(:,Si,i)
        M(1,i) = M(1,i) + drDipol(2,Si,i)*FDipol(3,Si,i)&
                        - drDipol(3,Si,i)*FDipol(2,Si,i)
        M(2,i) = M(2,i) + drDipol(3,Si,i)*FDipol(1,Si,i)&
                        - drDipol(1,Si,i)*FDipol(3,Si,i)
        M(3,i) = M(3,i) + drDipol(1,Si,i)*FDipol(2,Si,i)&
                        - drDipol(2,Si,i)*FDipol(1,Si,i)
      END DO
      DO Si=1,AnzQuadrupole(Komp(i)),1
        F(:,i) = F(:,i) + FQuadrupol(:,Si,i)
        M(1,i) = M(1,i) + drQuadrupol(2,Si,i)*FQuadrupol(3,Si,i)&
                        - drQuadrupol(3,Si,i)*FQuadrupol(2,Si,i)
        M(2,i) = M(2,i) + drQuadrupol(3,Si,i)*FQuadrupol(1,Si,i)&
                        - drQuadrupol(1,Si,i)*FQuadrupol(3,Si,i)
        M(3,i) = M(3,i) + drQuadrupol(1,Si,i)*FQuadrupol(2,Si,i)&
                        - drQuadrupol(2,Si,i)*FQuadrupol(1,Si,i)
      END DO
    END DO

    ! 4) Leapfrog Integrator (Part II)
    !    with FEQA (Fincham's Explicit Quaternion Algorithm
    ! 4a) Translation
    sumMv2 = 0.0
    DO i=1,N,1
      v(:,i) = v(:,i) + dtInv2mKomp(Komp(i))*F(:,i)
      sumMv2 = sumMv2 + mKomp(Komp(i))*DOT_PRODUCT(v(:,i),v(:,i))
    END DO
    betaTrans = SQRT(3.0*N*T/sumMv2)
    v = betaTrans*v

    ! 4b) Rotation
    sumIw2 = 0.0
    DO i=1,N,1
      D(:,i) = D(:,i) + dtInv2*M(:,i)

      A(1,1) = q(0,i)*q(0,i) + q(1,i)*q(1,i) - q(2,i)*q(2,i) - q(3,i)*q(3,i)
      A(1,2) = 2.0 * ( q(1,i)*q(2,i) + q(0,i)*q(3,i) )
      A(1,3) = 2.0 * ( q(1,i)*q(3,i) - q(0,i)*q(2,i) )
      A(2,1) = 2.0 * ( q(1,i)*q(2,i) - q(0,i)*q(3,i) )
      A(2,2) = q(0,i)*q(0,i) - q(1,i)*q(1,i) + q(2,i)*q(2,i) - q(3,i)*q(3,i)
      A(2,3) = 2.0 * ( q(2,i)*q(3,i) + q(0,i)*q(1,i) )
      A(3,1) = 2.0 * ( q(1,i)*q(3,i) + q(0,i)*q(2,i) )
      A(3,2) = 2.0 * ( q(2,i)*q(3,i) - q(0,i)*q(1,i) )
      A(3,3) = q(0,i)*q(0,i) - q(1,i)*q(1,i) - q(2,i)*q(2,i) + q(3,i)*q(3,i)
      wBody =  InvIBody(:,Komp(i)) * MATMUL(A,D(:,i))

      sumIw2 = sumIw2 + SUM( IBody(:,Komp(i))*wBody*wBody )
    END DO
    IF (sumIw2.EQ.0.0) THEN   ! avoid zero devision for 1CLJ
      betaRot = 1.0
    ELSE
      betaRot = SQRT( SUM(FHG*NKomp)*T/sumIw2 )
    END IF
    D = betaRot*D

    ! 4c) Volume
    PAktuell = T*rho&
              + rho*(VirialSite+VirialE+VirialMy+VirialMyQ+VirialQ)/(3.0*N)&
              + VirialKorrLJ*rho*rho + VirialKorrMy
    IF ((NPT).AND.(IntSchritt.GE.(-1*EquiNPTSchritte))) THEN
      V2 = (PAktuell - P) * 0.5 * dt * dt / VMasse
      V1 = V1 + V2
    END IF

    ! check of thermostat: output in "..._Sim.dat"
    betaRotMax = MAX(betaRot,betaRotMax)
    betaRotMin = MIN(betaRot,betaRotMin)
    betaTransMax = MAX(betaTrans,betaTransMax)
    betaTransMin = MIN(betaTrans,betaTransMin)

    sumDruck =  sumDruck + PAktuell
    sumUpot  =  sumUpot&
              + (UpotSite+UpotE+UpotMy+UpotMyQ+UpotQ)/N&
              + UCorrSum/(2*N) + UpotKorrMy + UpotKorrLJ
    sumrho   = sumrho + rho
    IF (IntSchritt == 0) THEN
      rhoShellsAve = 0
    END IF
    IF (IntSchritt > 0) THEN
      ! Temperature
      sumv = 0.0
      DO i=1, N, 1
        vsquare = v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i)
        sumv = sumv + vsquare
      END DO
      SumEkin = SumEkin + 0.5*sumv
      UpotTot  = UpotTot&
                + (UpotSite+UpotE+UpotMy+UpotMyQ+UpotQ)/N&
                + UCorrSum/(2*N) + UpotKorrLJ + UpotKorrMy
      IF (IntSchritt .EQ. SimSchritte/2) THEN
        UPotMin = UpotTot/IntSchritt
        UPotMax = UPotTot/IntSchritt
      END IF
      IF (MOD(IntSchritt,Output).EQ.0) THEN
        ! error by method after Lotfi et al. Ref. 44
        IF (IntSchritt .GE. SimSchritte/2) THEN
          IF (UpotTot/IntSchritt < UPotMin) THEN
            UPotMin = UpotTot/IntSchritt
          ELSE IF (UpotTot/IntSchritt > UPotMax) THEN
            UPotMax = UpotTot/IntSchritt
          END IF
        END IF
        ! running average of potential energy
        OPEN(12,FILE=TRIM(Dateiname)//'_RAV.csv',POSITION='APPEND',STATUS='OLD')
          WRITE(12,'(I10,";",E20.10,";",E20.10)') IntSchritt, UpotTot/IntSchritt,&
         MAX( (UPotTot/IntSchritt-UPotMin),(UPotMax-UPotTot/IntSchritt) )
        CLOSE(12)
        ! average density and current corrections per shell
        OPEN(13,FILE=TRIM(Dateiname)//'_Shells.csv',STATUS='REPLACE')
          WRITE(13,'(A,A)') '            R;        rho;          U_corr;        '&
                            'F_corr;        PiN_corr;        PiT_corr'
          DO i=1,NShells,1
            WRITE(13,'(F12.6,";",F12.6,";",F12.6,";",F12.6,";",F12.6,";",F12.6)') RShells(i),&
                      rhoShellsAve(i)/(IntSchritt), &
                      UShells_Mean(i), &
                      FShells_Mean(i), &
                      PNShells_Mean(i), &
                      PTShells_Mean(i)
          END DO
        CLOSE(13)
        ! averages of normal and tangential part of virial and pressure per shell
        OPEN(49,FILE=TRIM(Dateiname)//'_Pres.csv',STATUS='REPLACE')
          WRITE(49,'(A,A)') '            R;           Pi_N;           Pi_T;           &
            p_N;          p_T'
          DO i=1,NShells,1
            WRITE(49,'(F12.8,";",F12.8,";",F14.8,";",F14.8,";",F14.8)') &
                    RShells(i), &
                    -PShells_N(i)/(VShells(i)*(TotalStep+1)), &
                    -PShells_T(i)/(VShells(i)*(TotalStep+1)), &
                    T*rhoShellsAve(i)/(IntSchritt) + PShells_N(i)/(VShells(i)*(TotalStep+1)), &
                    T*rhoShellsAve(i)/(IntSchritt) + PShells_T(i)/(VShells(i)*(TotalStep+1))
           END DO
        CLOSE(49)
        ! positions, orientations, velocities, momenta of particles; can be used for restart
        OPEN(10,FILE=TRIM(Dateiname)//'_POS.dat',STATUS='REPLACE')
          DO i=1,N,1
            WRITE(10,'(13F12.8)') r(:,i), q(:,i), v(:,i), D(:,i)
          END DO
        CLOSE(10)
      END IF
    END IF

    ! Output of results: Average over "Output"-steps
    IF (MOD(IntSchritt,Output).EQ.0) THEN
      OPEN(5,file=TRIM(Dateiname)//'_Erg.dat',POSITION='APPEND',STATUS='OLD')
      WRITE(5,'(I15,4E20.10)') IntSchritt,SUM(SUM(Disp*Disp,Dim=1))/N,&
        sumUpot/Output,sumDruck/Output,sumrho/Output
      CLOSE(5)
      sumDruck = 0.0
      sumUpot  = 0.0
      sumrho   = 0.0
    END IF

    ! Output of visualisaition if desired
    IF (       Film&
        .AND. ( BildNr .LT. BilderMax )&
        .AND. ( IntSchritt .GE. FilmAb )&
        .AND. ( MOD((IntSchritt-FilmAb),OutputVis) .EQ. 0 ) ) THEN
      OPEN(4,file=TRIM(Dateiname)//'.vim',POSITION='APPEND',STATUS='OLD')
      WRITE(4,'(A,F10.3)') '# ',MAXVAL(LSimBox)
      DO i=1,N,1
        WRITE(4,'(A,I2,3I4,4I5)') &
          '! ',Komp(i),INT(r(:,i)/MAXVAL(LSimBox) * 999),INT(q(:,i)*999)
      END DO
      WRITE(4,*)
      BildNr = BildNr+1
      CLOSE(4)
    END IF

    ! Output after end of equilibration
    IF (IntSchritt.EQ.0) THEN
      OPEN(7,file=TRIM(Dateiname)//'_Sim.dat',POSITION='APPEND',STATUS='OLD')
        WRITE(7,'(A,I20,A)')&
          'max. Iterations   = ',MaxIterationen,'  (Equilibration)'
        WRITE(7,'(A,F10.8,A,F10.8,A)')&
          'beta-Faktoren     =     ',MIN(betaRotMin,betaTransMin),&
          ' --- ',MAX(betaRotMax,betaTransMax),'  (Equilibration)'
      CLOSE(7)
      betaRotMax = 1.0
      betaRotMin = 1.0
      betaTransMax = 1.0
      betaTransMin = 1.0
      Disp = 0.0
      MaxIterationen = 0
    END IF

  END DO   ! IntStep ;  time step = t+dt
!*******   end of integration step   ****************************************


  IF (Film) THEN
    OPEN(4,file=TRIM(Dateiname)//'.vim',POSITION='APPEND',STATUS='OLD')
    WRITE(4,'(A)') '##'
    CLOSE(4)
  END IF

  OPEN(7,file=TRIM(Dateiname)//'_Sim.dat',POSITION='APPEND',STATUS='OLD')
    WRITE(7,'(A,F10.8,A,F10.8,A)')&
      'beta-Factors      =     ',MIN(betaRotMin,betaTransMin),&
      ' --- ',MAX(betaRotMax,betaTransMax),'  (Simulation)'
    WRITE(7,'(A,I20,A)') 'max. Iterations   = ',MaxIterationen,'  (Simulation)'
    WRITE(7,'(A,I20,A)') 'NmaxCell          = ',NmaxZelle,'  (End)'
  CLOSE(7)


  CONTAINS
!*****************************************************************************
!***  Functions and Subroutines                                            ***
!*****************************************************************************
    INCLUDE 'HATsym3x3.f90'
    INCLUDE 'TIFunctions.f90'
    INCLUDE 'ShellFunctions.f90'
!*******   Functions and Subroutines   ***************************************

END PROGRAM ZellFIQA
