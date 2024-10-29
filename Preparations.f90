!*****************************************************************************
!***  1) Physical constants                                                ***
!***  2) Read input data from input file and write to revision file        ***
!***  3) Coordinate transformation                                         ***
!***  4) Calculation of necessary quantities and constants                 ***
!***  5) Preparations for long range corrections                           ***
!***  6) Preparation of output files                                       ***
!***  7) Allocate undefined arrays                                         ***
!*****************************************************************************

! 1) Physical constants
!    http://physics.nist.gov/constants
!    Peter J. Mohr and Barry N. Taylor: "CODATA Recommended Values of
!    the Fundamental Physical Constants", 2002
!*****************************************************************************
AtomareMasseDim = 1.660539E-27   ! [atomic_mass] = kg
AvogadroDim = 6.02214E+23        ! [Avogadro] = 1/mol
BoltzmannDim = 1.38065E-23       ! [Boltzmann] = J/K
EChargeDim =  1.60217733E-19     ! [elementary charge] = C
EpsilonVakuum = 8.8541878E-12    ! [electrical conductivity for vacuum] = C /Vm

Drittel = 1.0 / 3.0
pi = 2.0*ACOS(0.0)
VierDPi = 4 * Drittel * pi

! 2) Read input data from input file and write to revision file
!*****************************************************************************
i = IARGC()
IF (i.EQ.1) THEN
  CALL GETARG(0,Programmname)
  CALL GETARG(1,Dateiname)
  OPEN(1,FILE=TRIM(Dateiname)//'.dat',IOSTAT=i,STATUS='OLD')
  IF (i.EQ.0) THEN
    OPEN(2,FILE=TRIM(Dateiname)//'_Rev.dat',STATUS='REPLACE')
    READ(1,*,ERR=377)
    WRITE(2,'(5A)')&
      '! >',TRIM(ADJUSTL(Programmname)),' ',TRIM(ADJUSTL(Dateiname)),Version

    READ(1,*,ERR=377)
    WRITE(2,*)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A20,D15.0)',ERR=377) Text,epsilonRefDim
    WRITE(2,'(A,F15.8,A)') 'eps/kB_Ref        = ',epsilonRefDim,'  K'
    epsilonRefDim = epsilonRefDim*BoltzmannDim
    READ(1,'(A20,D15.0)',ERR=377) Text,mRefDim
    WRITE(2,'(A,F15.8,A)')&
      'mRef              = ',mRefDim,'  atomic mass unit'
    mRefDim = mRefDim*AtomareMasseDim
    READ(1,'(A20,D15.0)',ERR=377) Text,sigmaRefDim
    WRITE(2,'(A,F15.8,A)') 'sigmaRef          = ',sigmaRefDim,'  Angstroem'
    sigmaRefDim = sigmaRefDim*1.0E-10

    READ(1,*,ERR=377)
    WRITE(2,*)
    READ(1,*,ERR=377)
    WRITE(2,'(A)')&
      '!*****************   State point and Simulation   *****************'
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A20,I15)',ERR=377) Text,AnzKomp
    WRITE(2,'(A,I15)') 'NumComp           = ',AnzKomp
    ALLOCATE( NKomp(1:AnzKomp) )
    DO i=1,AnzKomp,1
      READ(1,'(A20,I15)',ERR=377) Text,NKomp(i)
      WRITE(2,'(A,I2,A,I15)') 'NComp(',i,')         = ',NKomp(i)
    END DO
    N = SUM(NKomp)
    READ(1,'(A20,D15.0)',ERR=377) Text,rho
    WRITE(2,'(A,F15.8,A)') 'rho               = ',rho,'  mol/l'
    !rho = rho*1.0E+03 * (sigmaRefDim**3.0)*AvogadroDim   !  reduced
    READ(1,'(A20,D15.0)',ERR=377) Text,T
    WRITE(2,'(A,F15.8,A)') 'T                 = ',T,'  K'
    ! T = T * BoltzmannDim/epsilonRefDim   !  reduced
    READ(1,'(A20,L15)',ERR=377) Text,NPT
    WRITE(2,'(A,L15,A)') 'NPT               = ',NPT,'  (T=true, F=false)'
    READ(1,'(A20,D15.0)',ERR=377) Text,P
    WRITE(2,'(A,F15.8,A)') 'P                 = ',P,'  MPa' 
    ! P = P * 1.0E+6 * (sigmaRefDim**3.0) / (epsilonRefDim)   !  reduced
    IF (.NOT.NPT) P = 0.0
    READ(1,'(A20,D15.0)',ERR=377) Text,VMasse
    WRITE(2,'(A,F15.8,A)') 'VMass             = ',VMasse,'  reduced' !  reduced
    V0=N/rho
    V1=0.0
    V2=0.0

    READ(1,*,ERR=377)
    WRITE(2,*)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A20,A15)',ERR=377) Text,InputFile
    InputFile = ADJUSTL(TRIM(InputFile))
    READ(1,'(A20,D15.0)',ERR=377) Text,dt
    WRITE(2,'(A,F15.8,A)') 'dt                = ',dt,'  fs'
    ! dt = dt*1.0E-15 * SQRT(epsilonRefDim/mRefDim)/sigmaRefDim   !  reduced
    READ(1,'(A20,D15.0)',ERR=377) Text,rc
    WRITE(2,'(A,F15.8,A)') 'rc                = ',rc,'  sigma'
    READ(1,'(A20,I15)',ERR=377) Text,EquiNVTSchritte
    WRITE(2,'(A,I15)') 'EquiNVTSteps      = ',EquiNVTSchritte
    READ(1,'(A20,I15)',ERR=377) Text,EquiNPTSchritte
    WRITE(2,'(A,I15,A)') 'EquiNPTSteps      = ',EquiNPTSchritte,'   =0, if NPT=F'
    IF (.NOT.NPT) EquiNPTSchritte=0
    READ(1,'(A20,I15)',ERR=377) Text,SimSchritte
    WRITE(2,'(A,I15)') 'SimSteps          = ',SimSchritte
    READ(1,'(A20,D15.0)',ERR=377) Text,SimBoxRatio(1)
    WRITE(2,'(A,F15.8)') 'BoxRatio(x)       = ',SimBoxRatio(1)
    READ(1,'(A20,D15.0)',ERR=377) Text,SimBoxRatio(2)
    WRITE(2,'(A,F15.8)') 'BoxRatio(y)       = ',SimBoxRatio(2)
    READ(1,'(A20,D15.0)',ERR=377) Text,SimBoxRatio(3)
    WRITE(2,'(A,F15.8)') 'BoxRatio(z)       = ',SimBoxRatio(3)
    READ(1,'(A20,I15)',ERR=377) Text,Output
    WRITE(2,'(A,I15)') 'Output            = ',Output
    READ(1,'(A20,I15)',ERR=377) Text,NShells
    WRITE(2,'(A,I15)') 'NShells           = ',NShells

    OPEN(5,file=TRIM(Dateiname)//'_Erg.dat',STATUS='REPLACE')
      WRITE(5,*) EquiNVTSchritte
      WRITE(5,*) EquiNPTSchritte
      WRITE(5,*) SimSchritte
      WRITE(5,*) Output
      WRITE(5,*) sigmaRefDim*1.0E+10
      WRITE(5,*) epsilonRefDim/BoltzmannDim
    CLOSE(5)

    READ(1,*,ERR=377)
    WRITE(2,*)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A20,L15)',ERR=377) Text,Film
    WRITE(2,'(A,L15,A)') 'Film              = ',Film,'  (T=true, F=false)'
    READ(1,'(A20,I15)',ERR=377) Text,BilderMax
    WRITE(2,'(A,I15)') 'ImagesMax         = ',BilderMax
    READ(1,'(A20,I15)',ERR=377) Text,FilmAb
    WRITE(2,'(A,I15)') 'FilmStart         = ',FilmAb
    READ(1,'(A20,I15)',ERR=377) Text,OutputVis
    WRITE(2,'(A,I15)') 'OutputVis         = ',OutputVis

    READ(1,*,ERR=377)
    WRITE(2,*)
    READ(1,*,ERR=377)
    WRITE(2,'(A)')&
      '!*****************   Description of molecules   *****************'
    READ(1,'(A20,I15)',ERR=377) Text,MaxAnzSites
    WRITE(2,'(A,I15)') 'MaxNumSites       = ',MaxAnzSites
    READ(1,'(A20,I15)',ERR=377) Text,MaxAnzLadungen
    WRITE(2,'(A,I15)') 'MaxNumCharges     = ',MaxAnzLadungen
    READ(1,'(A20,I15)',ERR=377) Text,MaxAnzDipole
    WRITE(2,'(A,I15)') 'MaxNumDipoles     = ',MaxAnzDipole
    READ(1,'(A20,I15)',ERR=377) Text,MaxAnzQuadrupole
    WRITE(2,'(A,I15)') 'MaxNumQuadrupoles = ',MaxAnzQuadrupole

    READ(1,*,ERR=377)
    WRITE(2,*)
    READ(1,'(A20,D15.0)',ERR=377) Text,epsilonRF
    WRITE(2,'(A,E15.8)') 'eps_ReactionField = ',epsilonRF

    READ(1,*,ERR=377)
    WRITE(2,*)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    ALLOCATE( etaLB(1:AnzKomp,1:AnzKomp) )
    ALLOCATE( xiLB(1:AnzKomp,1:AnzKomp) )
    etaLB = 1.0   ! to ensure etaLB(i,i)=1.0 
    xiLB = 1.0   ! to ensure xiLB(i,i)=1.0 
    DO j=1,AnzKomp-1,1
      DO i=j+1,AnzKomp,1
        READ(1,'(A20,2D15.0)',ERR=377) Text,etaLB(i,j),xiLB(i,j)
        WRITE(2,'(A,I2,A,I2,A,2F15.8)')&
          '(eta,xi) (',i,',',j,')  = ',etaLB(i,j),xiLB(i,j)
        etaLB(j,i) = etaLB(i,j)   ! symmetric
        xiLB(j,i) = xiLB(i,j)   ! symmetric
      END DO
      READ(1,*,ERR=377)
      WRITE(2,*)
    END DO

    ALLOCATE( AnzSites(1:AnzKomp) )
    ALLOCATE( AnzLadungen(1:AnzKomp) )
    ALLOCATE( AnzDipole(1:AnzKomp) )
    ALLOCATE( AnzQuadrupole(1:AnzKomp) )
    ALLOCATE( IBodyDummy(1:3,1:AnzKomp) )
    IBodyDummy = 0.0   
    ALLOCATE( rSiteBody(1:3,1:MaxAnzSites,1:AnzKomp) )
    rSiteBody = 0.0   
    ALLOCATE( rLadungBody(1:3,1:MaxAnzLadungen,1:AnzKomp) )
    rLadungBody = 0.0   
    ALLOCATE( rDipolBody(1:3,1:MaxAnzDipole,1:AnzKomp) )
    rDipolBody = 0.0   
    ALLOCATE( rQuadrupolBody(1:3,1:MaxAnzQuadrupole,1:AnzKomp) )
    rQuadrupolBody = 0.0   
    ALLOCATE( epsilonSite(1:MaxAnzSites,1:AnzKomp) )
    epsilonSite = 0.0   
    ALLOCATE( mSite(1:MaxAnzSites,1:AnzKomp) )
    mSite = 0.0    
    ALLOCATE( mLadung(1:MaxAnzLadungen,1:AnzKomp) )
    mLadung = 0.0    
    ALLOCATE( MyGes(1:3,1:AnzKomp) )
    MyGes = 0.0    
    ALLOCATE( sigmaSite(1:MaxAnzSites,1:AnzKomp) )
    sigmaSite = 0.0    
    ALLOCATE( absLadung(1:MaxAnzLadungen,1:AnzKomp) )
    absLadung = 0.0    
    ALLOCATE( eMyBody(1:3,1:MaxAnzDipole,1:AnzKomp) )
    eMyBody = 0.0    
    ALLOCATE( absMy(1:MaxAnzDipole,1:AnzKomp) )
    absMy = 0.0   
    ALLOCATE( eQBody(1:3,1:MaxAnzQuadrupole,1:AnzKomp) )
    eQBody = 0.0   
    ALLOCATE( absQ(1:MaxAnzQuadrupole,1:AnzKomp) )
    absQ = 0.0   

    IF (AnzKomp.EQ.1) THEN
      READ(1,*,ERR=377)
      WRITE(2,*)
    ELSE
      
    END IF
    DO i=1,AnzKomp,1
      READ(1,'(A)',ERR=377) Text
      WRITE(2,'(A)') TRIM(Text)
      READ(1,'(A)',ERR=377) Text
      WRITE(2,'(A)') TRIM(Text)
      READ(1,'(A20,I15)',ERR=377) Text,AnzSites(i)
      WRITE(2,'(A,I15)') 'NumSites          = ',AnzSites(i)
      READ(1,'(A20,I15)',ERR=377) Text,AnzLadungen(i)
      WRITE(2,'(A,I15)') 'NumCharges        = ',AnzLadungen(i)
      READ(1,'(A20,I15)',ERR=377) Text,AnzDipole(i)
      WRITE(2,'(A,I15)') 'NumDipoles        = ',AnzDipole(i)
      READ(1,'(A20,I15)',ERR=377) Text,AnzQuadrupole(i)
      WRITE(2,'(A,I15)') 'NumQuadrupoles    = ',AnzQuadrupole(i)
      IF (     (AnzSites(i).GT.MaxAnzSites)&
          .OR. (AnzLadungen(i).GT.MaxAnzLadungen)&
          .OR. (AnzDipole(i).GT.MaxAnzDipole)&
          .OR. (AnzQuadrupole(i).GT.MaxAnzQuadrupole) ) THEN
        PRINT *, '>',TRIM(ADJUSTL(Programmname)),' ',TRIM(ADJUSTL(Dateiname))
        PRINT '(A,3I2,A)',&
          ' Enter for "MaxNum-(Sites,Charges,Dipoles,Quadrupoles)" = ',&
          MAX(AnzSites(i),MaxAnzSites),&
          MAX(AnzLadungen(i),MaxAnzLadungen),&
          MAX(AnzDipole(i),MaxAnzDipole),&
          MAX(AnzQuadrupole(i),MaxAnzQuadrupole)
        STOP
      END IF

      READ(1,*,ERR=377)
      WRITE(2,*)
      DO Si=1,AnzSites(i),1
        READ(1,'(A20,3D15.0)',ERR=377) Text,rSiteBody(1:3,Si,i)
        WRITE(2,'(A,I1,A,3F15.8,A)')&
          'rSite(1:3,',Si,')      = ',rSiteBody(1:3,Si,i),'  Angstroem'
      END DO
      DO Si=1,AnzLadungen(i),1
        READ(1,'(A20,3D15.0)',ERR=377) Text,rLadungBody(1:3,Si,i)
        WRITE(2,'(A,I1,A,3F15.8,A)')&
          'rCharge(1:3,',Si,')    = ',rLadungBody(1:3,Si,i),'  Angstroem'
      END DO
      DO Si=1,AnzDipole(i),1
        READ(1,'(A20,3D15.0)',ERR=377) Text,rDipolBody(1:3,Si,i)
        WRITE(2,'(A,I1,A,3F15.8,A)')&
          'rDipole(1:3,',Si,')     = ',rDipolBody(1:3,Si,i),'  Angstroem'
      END DO
      DO Si=1,AnzQuadrupole(i),1
        READ(1,'(A20,3D15.0)',ERR=377) Text,rQuadrupolBody(1:3,Si,i)
        WRITE(2,'(A,I1,A,3F15.8,A)')&
          'rQuadrupole(1:3,',Si,') = ',rQuadrupolBody(1:3,Si,i),'  Angstroem'
      END DO

      READ(1,*,ERR=377)
      WRITE(2,*)
      DO Si=1,AnzSites(i),1
        READ(1,'(A20,3D15.0)',ERR=377)&
          Text,epsilonSite(Si,i),mSite(Si,i),sigmaSite(Si,i)
        WRITE(2,'(A,I1,A,3F15.8,A)') '(eps/kB,m,sig)(',Si,') = ',&
          epsilonSite(Si,i),mSite(Si,i),sigmaSite(Si,i),&
          '  (K,atomic mass unit,Angstroem)'
      END DO
      DO Si=1,AnzLadungen(i),1
        READ(1,'(A20,2D15.0)',ERR=377)&
          Text,absLadung(Si,i),mLadung(Si,i)
        WRITE(2,'(A,I1,A,2F15.8,A)') 'abs,mCharge(',Si,')    = ',&
          absLadung(Si,i),mLadung(Si,i),&
          '  [absLadung,mLadung] = wieviele Protonenladung, atomic mass unit'
      END DO
      DO Si=1,AnzDipole(i),1
        READ(1,'(A20,4D15.0)',ERR=377)&
          Text,eMyBody(1:3,Si,i),absMy(Si,i)
        WRITE(2,'(A,I1,A,4F15.8,A)') 'eMy(1:3);absMy (',Si,')= ',&
          eMyBody(1:3,Si,i),absMy(Si,i),'  (scalar;Debyes)'
        dr2 = DOT_PRODUCT(eMyBody(1:3,Si,i),eMyBody(1:3,Si,i))
        IF ( dr2.NE.0.0 ) THEN
          eMyBody(1:3,Si,i) = eMyBody(1:3,Si,i) / SQRT(dr2)   ! normalised
        END IF
      END DO
      DO Si=1,AnzQuadrupole(i),1
        READ(1,'(A20,4D15.0)',ERR=377)&
          Text,eQBody(1:3,Si,i),absQ(Si,i)
        WRITE(2,'(A,I1,A,4F15.8,A)') 'eQ(1:3);absQ (',Si,')  = ',&
          eQBody(1:3,Si,i),absQ(Si,i),'  (Angstroem;Buckinghams)'
        dr2 = DOT_PRODUCT(eQBody(1:3,Si,i),eQBody(1:3,Si,i))
        IF ( dr2.NE.0.0 ) THEN
          eQBody(1:3,Si,i) = eQBody(1:3,Si,i) / SQRT(dr2)   ! normalised
        END IF
      END DO

      READ(1,*,ERR=377)
      WRITE(2,*)
      READ(1,'(A20,3D15.0)',ERR=377) Text,IBodyDummy(1:3,i)
      WRITE(2,'(A,3F15.8,A)') 'I_Dummy(1:3)      = ',&
        IBodyDummy(1:3,i),'  atomic mass unit * Angstroem^2'
      IBodyDummy(1:3,i) = IBodyDummy(1:3,i) * (AtomareMasseDim*1.0E-20)

      READ(1,*,ERR=377)
      WRITE(2,*)
    END DO

    ! reduce all quantities
    IBodyDummy = IBodyDummy / (mRefDim*sigmaRefDim**2.0)
    rSiteBody = rSiteBody*1.0E-10 / sigmaRefDim
    rLadungBody = rLadungBody*1.0E-10 / sigmaRefDim
    rDipolBody = rDipolBody*1.0E-10 / sigmaRefDim
    rQuadrupolBody = rQuadrupolBody*1.0E-10 / sigmaRefDim
    epsilonSite = epsilonSite*BoltzmannDim / epsilonRefDim
    mSite = mSite*AtomareMasseDim / mRefDim
    mLadung = mLadung*AtomareMasseDim / mRefDim
    sigmaSite = sigmaSite*1.0E-10 / sigmaRefDim
    absLadung = absLadung * EChargeDim / (SQRT ( 4.*2.0*ACOS(0.0)*EpsilonVakuum * epsilonRefDim * sigmaRefDim) )
    absMy =  absMy / (1.0E+24 * SQRT(10.0*epsilonRefDim*(sigmaRefDim**3.0)) )
    absQ =  absQ / (1.0E+34 * SQRT(10.0*epsilonRefDim*(sigmaRefDim**5.0)) )

    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    READ(1,'(A)',ERR=377) Text
    WRITE(2,'(A)') TRIM(Text)
    ALLOCATE( Farbe(1:MAXVAL(AnzSites),1:AnzKomp) )
    Farbe = 1   
    DO i=1,AnzKomp,1
       WRITE(xMal,*) AnzSites(i)
       READ(1,'(A20,'//TRIM(ADJUSTL(xMal))//'I15)',ERR=377) Text,Farbe(1:AnzSites(i),i)
       WRITE(2,'(A,I1,A,I2,A,'//TRIM(ADJUSTL(xMal))//'I15)')&
          'Colour(1:',AnzSites(i),',',i,')     = ',Farbe(1:AnzSites(i),i)
    END DO

    CLOSE(1)   ! Input file
    CLOSE(2)   ! Revision file
    GOTO 5867
      377 CONTINUE
      PRINT *, '>',TRIM(ADJUSTL(Programmname)),' ',TRIM(ADJUSTL(Dateiname))
      PRINT *, 'An error occurred while reading the input file. "'&
               ,TRIM(Dateiname),'.dat".'
      PRINT *, 'Please check syntax.'
      STOP
    5867  CONTINUE
  ELSE
    PRINT *, '>',TRIM(ADJUSTL(Programmname)),' ',TRIM(ADJUSTL(Dateiname))
    PRINT *, 'Input file "',TRIM(Dateiname),'.dat" could not be opened.'
    PRINT *, 'Check existence/name.'
    STOP
  END IF
ELSE
  PRINT *, '>',TRIM(ADJUSTL(Programmname)),' ',TRIM(ADJUSTL(Dateiname))
  PRINT *, 'There has to be EXACTLY one input parameter.'
  STOP
END IF


! Variables for temperature control
sumMv2 = 0.0
betaTrans = 1.0
x = 0
sumIw2 = 0.0
betaRot = 1.0


! 3) Coordinate transformation 
!*****************************************************************************
ALLOCATE( FHG(1:AnzKomp) )
ALLOCATE( mKomp(1:AnzKomp) )
mKomp = 0.0   ! Initialisation
ALLOCATE( IBody(1:3,1:AnzKomp) )
ALLOCATE( InvIBody(1:3,1:AnzKomp) )

OPEN(3,FILE=TRIM(Dateiname)//'_HAS.dat',STATUS='REPLACE')
WRITE(3,'(A)') 'Transformation of site coordinates to system in center of mass'
WRITE(3,*)
WRITE(3,'(A)') 'All quantities have the following units:'
WRITE(3,'(A)') '   [r] = Angstroem'
WRITE(3,'(A)') '   [I] = atomic mass unit * Angstroem^2'

DO i=1,AnzKomp,1
  WRITE(3,*)
  WRITE(3,*)
  WRITE(3,*)
  WRITE(3,'(A,I2,A)') 'Molecule (',i,') before transformation'
  DO Si=1,AnzSites(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'rSite(1:3,',Si,')      = ',&
      rSiteBody(1:3,Si,i) * sigmaRefDim / 1.0E-10
  END DO
  DO Si=1,AnzLadungen(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'rCharge(1:3,',Si,')    = ',&
      rLadungBody(1:3,Si,i) * sigmaRefDim / 1.0E-10
  END DO
  DO Si=1,AnzDipole(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'rDipole(1:3,',Si,')     = ',&
      rDipolBody(1:3,Si,i) * sigmaRefDim / 1.0E-10
  END DO
  DO Si=1,AnzQuadrupole(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'rQuadrupole(1:3,',Si,') = ',&
      rQuadrupolBody(1:3,Si,i) * sigmaRefDim / 1.0E-10
  END DO
  DO Si=1,AnzDipole(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'eMy(1:3,',Si,')        = ',eMyBody(1:3,Si,i)
  END DO
  DO Si=1,AnzQuadrupole(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'eQ(1:3,',Si,')         = ',eQBody(1:3,Si,i)
  END DO

  ! translate coordinate system to center of mass

  ! Total mass of component i:
  mKomp(i) = SUM( mSite(1:AnzSites(i),i) ) + SUM( mLadung(1:AnzLadungen(i),i) )
  dr = 0.0

  DO Si=1,AnzSites(i),1
    dr(1) = dr(1) + rSiteBody(1,Si,i) * mSite(Si,i)
    dr(2) = dr(2) + rSiteBody(2,Si,i) * mSite(Si,i)
    dr(3) = dr(3) + rSiteBody(3,Si,i) * mSite(Si,i)
  END DO
  DO Si=1,AnzLadungen(i),1
    dr(1) = dr(1) + rLadungBody(1,Si,i) * mLadung(Si,i)
    dr(2) = dr(2) + rLadungBody(2,Si,i) * mLadung(Si,i)
    dr(3) = dr(3) + rLadungBody(3,Si,i) * mLadung(Si,i)
  END DO
  dr(1) = dr(1)/mKomp(i)
  dr(2) = dr(2)/mKomp(i)
  dr(3) = dr(3)/mKomp(i)

  DO Si=1,AnzSites(i),1
    rSiteBody(:,Si,i) = rSiteBody(:,Si,i) - dr(:)
  END DO
  DO Si=1,AnzLadungen(i),1
    rLadungBody(:,Si,i) = rLadungBody(:,Si,i) - dr(:)
  END DO
  DO Si=1,AnzDipole(i),1
    rDipolBody(:,Si,i) = rDipolBody(:,Si,i) - dr(:)
  END DO
  DO Si=1,AnzQuadrupole(i),1
    rQuadrupolBody(:,Si,i) = rQuadrupolBody(:,Si,i) - dr(:)
  END DO

  A=0.0
  DO Si=1,AnzSites(i),1
     A(1,1) = A(1,1) + ( rSiteBody(2,Si,i)**2.0&
                     +   rSiteBody(3,Si,i)**2.0&
                       ) * mSite(Si,i)
     A(1,2) = A(1,2) - (  rSiteBody(1,Si,i)&
                     *    rSiteBody(2,Si,i)&
                       ) * mSite(Si,i)
     A(1,3) = A(1,3) - (  rSiteBody(1,Si,i)&
                     *    rSiteBody(3,Si,i)&
                       ) * mSite(Si,i)
     A(2,2) = A(2,2) + ( rSiteBody(1,Si,i)**2.0&
                     +   rSiteBody(3,Si,i)**2.0&
                       ) * mSite(Si,i)
     A(2,3) = A(2,3) - ( rSiteBody(2,Si,i)&
                     *   rSiteBody(3,Si,i)&
                       ) * mSite(Si,i)
     A(3,3) = A(3,3) + ( rSiteBody(1,Si,i)**2.0&
                     +   rSiteBody(2,Si,i)**2.0&
                        ) * mSite(Si,i)
  END DO
  DO Si=1,AnzLadungen(i),1
     A(1,1) = A(1,1) + ( rLadungBody(2,Si,i)**2.0&
                     +   rLadungBody(3,Si,i)**2.0&
                       ) * mLadung(Si,i)
     A(1,2) = A(1,2) - (  rLadungBody(1,Si,i)&
                     *    rLadungBody(2,Si,i)&
                       ) * mLadung(Si,i)
     A(1,3) = A(1,3) - (  rLadungBody(1,Si,i)&
                     *    rLadungBody(3,Si,i)&
                        ) * mLadung(Si,i)
     A(2,2) = A(2,2) + ( rLadungBody(1,Si,i)**2.0&
                     +   rLadungBody(3,Si,i)**2.0&
                       ) * mLadung(Si,i)
     A(2,3) = A(2,3) - ( rLadungBody(2,Si,i)&
                     *   rLadungBody(3,Si,i)&
                       ) * mLadung(Si,i)
     A(3,3) = A(3,3) + ( rLadungBody(1,Si,i)**2.0&
                     +   rLadungBody(2,Si,i)**2.0&
                        ) * mLadung(Si,i)
  END DO
  A(2,1) = A(1,2)
  A(3,1) = A(1,3)
  A(3,2) = A(2,3)

  ! Input: inertia tensor (symmetric 3x3 matrix) 
  ! Output: matrix of orthonormal eigenvectors 
  CALL HATsym3x3(A)
  DO Si=1,AnzSites(i),1
    rSiteBody(1:3,Si,i) = MATMUL( rSiteBody(1:3,Si,i),A )
  END DO
  DO Si=1,AnzLadungen(i),1
    rLadungBody(1:3,Si,i) = MATMUL( rLadungBody(1:3,Si,i),A )
  END DO
  DO Si=1,AnzDipole(i),1
    rDipolBody(1:3,Si,i) = MATMUL( rDipolBody(1:3,Si,i),A )
    eMyBody(1:3,Si,i) = MATMUL( eMyBody(1:3,Si,i),A )
  END DO
  DO Si=1,AnzQuadrupole(i),1
    rQuadrupolBody(1:3,Si,i) = MATMUL( rQuadrupolBody(1:3,Si,i),A )
    eQBody(1:3,Si,i) = MATMUL( eQBody(1:3,Si,i),A )
  END DO

  A=0.0
  DO Si=1,AnzSites(i),1
     A(1,1) = A(1,1) + ( rSiteBody(2,Si,i)**2.0&
                     +   rSiteBody(3,Si,i)**2.0&
                       ) * mSite(Si,i)
     A(1,2) = A(1,2) - (  rSiteBody(1,Si,i)&
                     *    rSiteBody(2,Si,i)&
                       ) * mSite(Si,i)
     A(1,3) = A(1,3) - (  rSiteBody(1,Si,i)&
                     *    rSiteBody(3,Si,i)&
                       ) * mSite(Si,i)
     A(2,2) = A(2,2) + ( rSiteBody(1,Si,i)**2.0&
                     +   rSiteBody(3,Si,i)**2.0&
                       ) * mSite(Si,i)
     A(2,3) = A(2,3) - ( rSiteBody(2,Si,i)&
                     *   rSiteBody(3,Si,i)&
                       ) * mSite(Si,i)
     A(3,3) = A(3,3) + ( rSiteBody(1,Si,i)**2.0&
                     +   rSiteBody(2,Si,i)**2.0&
                        ) * mSite(Si,i)
  END DO
  DO Si=1,AnzLadungen(i),1
     A(1,1) = A(1,1) + ( rLadungBody(2,Si,i)**2.0&
                     +   rLadungBody(3,Si,i)**2.0&
                       ) * mLadung(Si,i)
     A(1,2) = A(1,2) - (  rLadungBody(1,Si,i)&
                     *    rLadungBody(2,Si,i)&
                       ) * mLadung(Si,i)
     A(1,3) = A(1,3) - (  rLadungBody(1,Si,i)&
                     *    rLadungBody(3,Si,i)&
                        ) * mLadung(Si,i)
     A(2,2) = A(2,2) + ( rLadungBody(1,Si,i)**2.0&
                     +   rLadungBody(3,Si,i)**2.0&
                       ) * mLadung(Si,i)
     A(2,3) = A(2,3) - ( rLadungBody(2,Si,i)&
                     *   rLadungBody(3,Si,i)&
                       ) * mLadung(Si,i)
     A(3,3) = A(3,3) + ( rLadungBody(1,Si,i)**2.0&
                     +   rLadungBody(2,Si,i)**2.0&
                        ) * mLadung(Si,i)
  END DO
  A(2,1) = A(1,2)
  A(3,1) = A(1,3)
  A(3,2) = A(2,3)

  FHG(i) = 3
  DO Richtung=1,3,1
    IF (IBodyDummy(Richtung,i).LE.0.0) THEN
      IF (     A(Richtung,Richtung)&
          .LE. MAX(A(1,1),A(2,2),A(3,3))*100.0*EPSILON(IBody) ) THEN
        IBody(Richtung,i) = 0.0
        InvIBody(Richtung,i) = 0.0
        FHG(i) = FHG(i)-1
      ELSE
        IBody(Richtung,i) = A(Richtung,Richtung)
        InvIBody(Richtung,i) = 1.0/IBody(Richtung,i)
      END IF
    ELSE
      IBody(Richtung,i) = IBodyDummy(Richtung,i)
      InvIBody(Richtung,i) = 1.0/IBody(Richtung,i)
    END IF
  END DO

  WRITE(3,*)
  WRITE(3,'(A,I2,A)')&
    'Molecule (',i,') in center of mass system'
  DO Si=1,AnzSites(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'rSite(1:3,',Si,')      = ',&
      rSiteBody(1:3,Si,i) * sigmaRefDim / 1.0E-10
  END DO
  DO Si=1,AnzLadungen(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'rCharge(1:3,',Si,')    = ',&
      rLadungBody(1:3,Si,i) * sigmaRefDim / 1.0E-10
  END DO
  DO Si=1,AnzDipole(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'rDipole(1:3,',Si,')     = ',&
      rDipolBody(1:3,Si,i) * sigmaRefDim / 1.0E-10
  END DO
  DO Si=1,AnzQuadrupole(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'rQuadrupole(1:3,',Si,') = ',&
      rQuadrupolBody(1:3,Si,i) * sigmaRefDim / 1.0E-10
  END DO
  DO Si=1,AnzDipole(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'eMy(1:3,',Si,')        = ',eMyBody(1:3,Si,i)
  END DO
  DO Si=1,AnzQuadrupole(i),1
    WRITE(3,'(A,I1,A,3E20.10)') 'eQ(1:3,',Si,')         = ',eQBody(1:3,Si,i)
  END DO
  WRITE(3,*)
  DO Si=1,3,1
    WRITE(3,'(A,I1,A,3E20.10)') 'IBody(',Si,',1:3)      = ',&
      A(Si,1:3) * mRefDim*sigmaRefDim**2.0 / (AtomareMasseDim*1.0E-20)
  END DO
  WRITE(3,*)
  IF (IBodyDummy(1,i).LE.0.0) THEN; Text=''; ELSE; Text='  Dummy'; END IF
  WRITE(3,'(A,E20.10,A7)') 'I_xx              = ',&
    IBody(1,i) * mRefDim*sigmaRefDim**2.0 / (AtomareMasseDim*1.0E-20),Text
  IF (IBodyDummy(2,i).LE.0.0) THEN; Text=''; ELSE; Text='  Dummy'; END IF
  WRITE(3,'(A,E20.10,A7)') 'I_yy              = ',&
    IBody(2,i) * mRefDim*sigmaRefDim**2.0 / (AtomareMasseDim*1.0E-20),Text
  IF (IBodyDummy(3,i).LE.0.0) THEN; Text=''; ELSE; Text='  Dummy'; END IF
  WRITE(3,'(A,E20.10,A7)') 'I_zz              = ',&
    IBody(3,i) * mRefDim*sigmaRefDim**2.0 / (AtomareMasseDim*1.0E-20),Text
END DO
CLOSE(3)   ! HAS


! 4) Calculation of necessary quantities and constants  
!*****************************************************************************
dtInv2 = dt/2.0
ALLOCATE( dtInv2mKomp(1:AnzKomp) )
DO i=1,AnzKomp,1
  dtInv2mKomp(i) = dt / (2.0*mKomp(i))
END DO
dtInv4 = dt/4.0
ALLOCATE( epsilon24(1:MAXVAL(AnzSites),1:AnzKomp,&
                    1:MAXVAL(AnzSites),1:AnzKomp) )
DO i=1,AnzKomp,1
  DO Si=1,AnzSites(i),1
    DO j=1,AnzKomp,1
      DO Sj=1,AnzSites(j),1
        epsilon24(Si,i,Sj,j) =  24.0*xiLB(i,j)&
                              * SQRT(epsilonSite(Si,i)*epsilonSite(Sj,j))
      END DO
    END DO
  END DO
END DO
epsRFInvrc3 = -2.0*(epsilonRF-1.0) / ((rc**3.0)*(2.0*epsilonRF+1.0))

LSimBox(1) = ( (N*SimBoxRatio(1)**2.0)&
              /(rho*SimBoxRatio(2)*SimBoxRatio(3)) ) ** (1.0/3.0)
LSimBox(2) = ( (N*SimBoxRatio(2)**2.0)&
              /(rho*SimBoxRatio(1)*SimBoxRatio(3)) ) ** (1.0/3.0)
LSimBox(3) = ( (N*SimBoxRatio(3)**2.0)&
              /(rho*SimBoxRatio(1)*SimBoxRatio(2)) ) ** (1.0/3.0)

SystemCenter = LSimBox * 0.5

rc2 = rc**2.0
ALLOCATE( sigma2(1:MAXVAL(AnzSites),1:AnzKomp,&
                 1:MAXVAL(AnzSites),1:AnzKomp) )
ALLOCATE( Invsigma(1:MAXVAL(AnzSites),1:AnzKomp,&
                 1:MAXVAL(AnzSites),1:AnzKomp) )
DO i=1,AnzKomp,1
  DO Si=1,AnzSites(i),1
    DO j=1,AnzKomp,1
      DO Sj=1,AnzSites(j),1
        sigma2(Si,i,Sj,j)&
          =  (etaLB(i,j)*(sigmaSite(Si,i)+sigmaSite(Sj,j))/2.0) ** 2.0
        Invsigma(Si,i,Sj,j) = 2.0/(etaLB(i,j)*(sigmaSite(Si,i)+sigmaSite(Sj,j)))
      END DO
    END DO
  END DO
END DO
ZT = MAX( 1 , NINT((rho*(rc**3.0)/5.0)**(1.0/3.0)) )   ! ca. 5 molecules/cell

AnzZellNachbarn = ((2.0*ZT+1.0)**3.0 - 1.0) / 2.0
delta = MAXVAL(LSimBox) * 100.0*EPSILON(LSimBox)
InvLSimBox = 1.0/LSimBox
Zellen = FLOOR( LSimBox*ZT/rc )
IF (MINVAL(Zellen).LE.2*ZT) THEN
  PRINT *, '>',TRIM(ADJUSTL(Programmname)),' ',TRIM(ADJUSTL(Dateiname))
  PRINT *, 'Die Simulationsbox ist zu klein fuer den Zellenalgorithmus!'
  PRINT *, 'Abhilfe: Teilchenzahl N rauf; Cut-Off rc runter oder'
  PRINT *, '         SimBoxRatio angleichen.'
  STOP
END IF
ALLOCATE( ZellNachbar(1:3,1:AnzZellNachbarn) )
i = 0
ZN: DO z=-ZT,ZT,1   ! ZN := Neighbour cell
      DO y=-ZT,ZT,1
        DO x=-ZT,ZT,1
          IF ( (x.EQ.0).AND.(y.EQ.0).AND.(z.EQ.0) ) THEN
            EXIT ZN   
          END IF
          i = i+1
          ZellNachbar(:,i) = (/x,y,z/)
        END DO
      END DO
    END DO ZN

InvLZelle = Zellen/LSimBox
LZelle = LSimBox/Zellen
ZweiInvLSimBox = 2.0/LSimBox

NmaxZelle = CEILING(PRODUCT(LZelle)*rho)   ! first approximation



! 5) Long range corrections ;  MySelbstTerm
!*****************************************************************************
OPEN(12,file=TRIM(Dateiname)//'_RAV.csv',STATUS='REPLACE')
     WRITE(12,'(A10,2A20)') 'Step','Upot', 'Error'
CLOSE(12)

OPEN(30,file=TRIM(Dateiname)//'_tanh.csv',STATUS='REPLACE')
     WRITE(30,'(A10,4A20)') 'Step;', 'rho_v;', 'rho_l;', 'D0;', 'R0'
CLOSE(30)

IntSchritt = 0

SumEkin = 0.0
UPotMin = 0.0
UPotMax = 0.0
SystemCOM = 0.0

! Delta_k = drShells
drShells = MIN(LSimBox(1), LSimBox(2), LSimBox(3)) / 2 / (NShells + 1)
deltaShells =  rc/drShells

drShells05 = drShells/2

! time steps for averaging density *1000
NSMean = 50

ALLOCATE( VShells(1:NShells) )
ALLOCATE( RShells(1:NShells) )
ALLOCATE( RShells2(1:NShells) )
ALLOCATE( rhoShells(1:NShells, 1:AnzKomp) )
ALLOCATE( rhoShellsT(1:NShells, 1:AnzKomp) )
ALLOCATE( rhoShellsTemp(1:NShells, 1:AnzKomp) )
ALLOCATE( rhoShellsMean(1:NShells, 1:AnzKomp, 1:NSMean) )
ALLOCATE( rhoShellsAve(1:NShells) )
ALLOCATE( PartShells(1:N) )
ALLOCATE( lowerS(1:N) )
ALLOCATE( upperS(1:N) )
ALLOCATE( interS(1:N) )
rhoShellsAve = 0.0

ALLOCATE( PNShells_Mean(1:NShells) )
ALLOCATE( PTShells_Mean(1:NShells) )
ALLOCATE( UShells_Mean(1:NShells) )
ALLOCATE( FShells_Mean(1:NShells) )

ALLOCATE( ksi(1:N) )
ALLOCATE( ksi2(1:N) )

ALLOCATE( PShells_N(1:NShells) )
ALLOCATE( PShells_T(1:NShells) )

! radii and volumes of shells
DO i=1, NShells-1, 1
  RShells(i) = i * drShells
  RShells2(i) = RShells(i)*RShells(i)
  VShells(i) = VierDPi *( (RShells(i)+drShells05)**3 - (RShells(i)-drShells05)**3 )
END DO
RShells(NShells) = NShells * drShells
RShells2(NShells) = RShells(NShells)*RShells(NShells)
VShells(NShells) = V0 - VierDPi * (RShells(NShells)-drShells05)**3

rhoShells = 0.0
rhoShellsMean = 0.0
PShells_N = 0.0
PShells_T = 0.0


rcmax = 8.0

! homogeneous terms
UpotKorrLJ = 0.0
VirialKorrLJ = 0.0
DO i=1,AnzKomp,1
  DO j=1,AnzKomp,1
    DO Si=1,AnzSites(i),1
      DO Sj=1,AnzSites(j),1
        tau1 = SQRT(DOT_PRODUCT( rSiteBody(:,Si,i),rSiteBody(:,Si,i) ))
        tau2 = SQRT(DOT_PRODUCT( rSiteBody(:,Sj,j),rSiteBody(:,Sj,j) ))
        IF (tau1+tau2.GE.rc) THEN
          PRINT *, '>',TRIM(ADJUSTL(Programmname)),&
                   ' ',TRIM(ADJUSTL(Dateiname))
          PRINT *, 'Homogeneous corrections could not be computed.'
          PRINT *, 'Note: rc to small!'
          STOP
        END IF
        IF ( (tau1.EQ.0.0).AND.(tau2.EQ.0.0) ) THEN
          UpotKorrLJ = UpotKorrLJ&
            +   REAL(NKomp(i),KIND=8)/REAL(N,KIND=8)&
              * REAL(NKomp(j),KIND=8)/REAL(N,KIND=8)&
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TICCu(-6,rc,sigma2(Si,i,Sj,j))&
                 - TICCu(-3,rc,sigma2(Si,i,Sj,j)) )
          VirialKorrLJ = VirialKorrLJ&
            +   REAL(NKomp(i),KIND=8)/REAL(N,KIND=8)&
              * REAL(NKomp(j),KIND=8)/REAL(N,KIND=8)&
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TICCp(-6,rc,sigma2(Si,i,Sj,j))&
                 - TICCp(-3,rc,sigma2(Si,i,Sj,j)) )
        END IF
        IF ( (tau1.EQ.0.0).NEQV.(tau2.EQ.0.0) ) THEN
          tau = MAX(tau1,tau2)   
          UpotKorrLJ = UpotKorrLJ&
            +   REAL(NKomp(i),KIND=8)/REAL(N,KIND=8)&
              * REAL(NKomp(j),KIND=8)/REAL(N,KIND=8)&
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TICSu(-6,rc,sigma2(Si,i,Sj,j),tau)&
                 - TICSu(-3,rc,sigma2(Si,i,Sj,j),tau) )
          VirialKorrLJ = VirialKorrLJ&
            +   REAL(NKomp(i),KIND=8)/REAL(N,KIND=8)&
              * REAL(NKomp(j),KIND=8)/REAL(N,KIND=8)&
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TICSp(-6,rc,sigma2(Si,i,Sj,j),tau)&
                 - TICSp(-3,rc,sigma2(Si,i,Sj,j),tau) )
        END IF
        IF ( (tau1.NE.0.0).AND.(tau2.NE.0.0) ) THEN
          UpotKorrLJ = UpotKorrLJ&
            +   REAL(NKomp(i),KIND=8)/REAL(N,KIND=8)&
              * REAL(NKomp(j),KIND=8)/REAL(N,KIND=8)&
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TISSu(-6,rc,sigma2(Si,i,Sj,j),tau1,tau2)&
                 - TISSu(-3,rc,sigma2(Si,i,Sj,j),tau1,tau2) )
          VirialKorrLJ = VirialKorrLJ&
            +   REAL(NKomp(i),KIND=8)/REAL(N,KIND=8)&
              * REAL(NKomp(j),KIND=8)/REAL(N,KIND=8)&
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TISSp(-6,rc,sigma2(Si,i,Sj,j),tau1,tau2)&
                 - TISSp(-3,rc,sigma2(Si,i,Sj,j),tau1,tau2) )
        END IF
      END DO
    END DO
  END DO
END DO
UpotKorrLJ = 2.0*pi*UpotKorrLJ
VirialKorrLJ = (2.0/3.0)*pi*VirialKorrLJ

UpotTot = 0.0


MySelbstTerm = 0.0
DO i=1,AnzKomp,1
  DO j=1,AnzLadungen(i),1
    MyGes(1,i) = MyGes(1,i) + absLadung(j,i)*rLadungBody(1,j,i)
    MyGes(2,i) = MyGes(2,i) + absLadung(j,i)*rLadungBody(2,j,i)
    MyGes(3,i) = MyGes(3,i) + absLadung(j,i)*rLadungBody(3,j,i)
  END DO
  DO Si=1,AnzDipole(i),1
    MyGes(1,i) = MyGes(1,i) + absMy(Si,i) * eMyBody(1,Si,i)
    MyGes(2,i) = MyGes(2,i) + absMy(Si,i) * eMyBody(2,Si,i)
    MyGes(3,i) = MyGes(3,i) + absMy(Si,i) * eMyBody(3,Si,i)
  END DO
  MySelbstTerm =  MySelbstTerm&
                  + REAL(NKomp(i),KIND=8)*DOT_PRODUCT(MyGes(:,i),MyGes(:,i))
END DO
MySelbstTerm = 0.5*epsRFInvrc3*MySelbstTerm



! 6) Preparation of output files for visualisation, results, ...
!*****************************************************************************
IF (Film) THEN
  OPEN(4,file=TRIM(Dateiname)//'.vim',STATUS='REPLACE')
  DO i=1,AnzKomp,1
    WRITE(Text,*) i
    DO Si=1,AnzSites(i)
      WRITE(4,'(A,4F8.3,I3)') '~ '//TRIM(ADJUSTL(TEXT))//' LJ',&
        rSiteBody(:,Si,i),sigmaSite(Si,i),Farbe(Si,i)
    END DO
  END DO
  WRITE(4,*)
  CLOSE(4)   
  BildNr = 0
END IF

OPEN(5,file=TRIM(Dateiname)//'_Erg.dat',POSITION='APPEND',STATUS='OLD')
     WRITE(5,'(A15,4A20)') 'Nr.','Displacement','Upot-Total','Pressure-Total',&
    'Density-Total'
CLOSE(5)

OPEN(7,file=TRIM(Dateiname)//'_Sim.dat',STATUS='REPLACE')
  WRITE(7,'(A,E20.10,A)') 'atomic mass       = ',AtomareMasseDim,'  kg'
  WRITE(7,'(A,E20.10,A)') 'Avogadro          = ',AvogadroDim,'  1/mol'
  WRITE(7,'(A,E20.10,A)') 'Boltzmann         = ',BoltzmannDim,'  J/K'
  WRITE(7,'(A,E20.10,A)') 'elem. charge      = ',EChargeDim,'  C'
  WRITE(7,'(A,E20.10,A)') 'Vacuum perm.      = ',EpsilonVakuum,'  C/Vm'
  WRITE(7,*)
  WRITE(7,'(A,E20.10,A)') 'epsilon_reference = ',epsilonRefDim,'  J'
  WRITE(7,'(A,E20.10,A)') 'masse_reference   = ',mRefDim,'  kg'
  WRITE(7,'(A,E20.10,A)') 'sigma_reference   = ',sigmaRefDim,'  m'
  WRITE(7,*)
  WRITE(7,'(A,I20)') 'N components      = ',AnzKomp
  WRITE(7,'(A,I20)') 'N molecules       = ',N
  DO i=1,AnzKomp,1
    WRITE(7,'(A,I2,A,E20.10)')&
      'Molar Part x(',i,')  = ',REAL(NKomp(i),KIND=8)/REAL(N,KIND=8)
  END DO
  WRITE(7,'(A,E20.10)') 'rho               = ',rho
  WRITE(7,'(A,E20.10)') 'T                 = ',T
  WRITE(7,'(A,E20.10)') 'P                 = ',P
  WRITE(7,'(A,E20.10)') 'VMass             = ',VMasse
  WRITE(7,*)
  WRITE(7,'(A,E20.10)') 'dt                = ',dt
  WRITE(7,'(A,E20.10)') 'rc                = ',rc
  WRITE(7,'(A,I20)') 'EquiNVTSteps      = ',EquiNVTSchritte
  WRITE(7,'(A,I20)') 'EquiNPTSteps      = ',EquiNPTSchritte
  WRITE(7,'(A,I20)') 'SimSteps          = ',SimSchritte
  WRITE(7,'(A,I20)') 'Output            = ',Output
  WRITE(7,*)
  WRITE(7,'(A,E20.10)') 'UpotKorrLJ        = ',UpotKorrLJ
  WRITE(7,'(A,E20.10)') 'VirialKorrLJ      = ',VirialKorrLJ
  WRITE(7,*)
  WRITE(7,'(A,3E20.10)') 'LSimBox(x,y,z)    = ',LSimBox
  WRITE(7,'(A,3E20.10)') 'LCells(x,y,z)     = ',LZelle
  WRITE(7,'(A,3I20)') 'Cells(x,y,z)      = ',Zellen
  WRITE(7,*)
  WRITE(7,'(A,I20,A)') 'NmaxCell          = ',NmaxZelle,'  (Start)'
CLOSE(7)



! 7) ALLOCATE undefined arrays
!*****************************************************************************
ALLOCATE( D(1:3,1:N) )
ALLOCATE( Disp(1:3,1:N) )
ALLOCATE( drLadung(1:3,1:MAXVAL(AnzLadungen),1:N) )
ALLOCATE( drDipol(1:3,1:MAXVAL(AnzDipole),1:N) )
ALLOCATE( drQuadrupol(1:3,1:MAXVAL(AnzQuadrupole),1:N) )
ALLOCATE( drSite(1:3,1:MAXVAL(AnzSites),1:N) )
ALLOCATE( eMy(1:3,1:MAXVAL(AnzDipole),1:N) )
ALLOCATE( MyRFSpace(1:3,1:AnzKomp,1:N) )
ALLOCATE( eQ(1:3,1:MAXVAL(AnzQuadrupole),1:N) )
ALLOCATE( F(1:3,1:N) )
ALLOCATE( FDipol(1:3,MAXVAL(AnzDipole),1:N) )
ALLOCATE( FQuadrupol(1:3,MAXVAL(AnzQuadrupole),1:N) )
ALLOCATE( FSite(1:3,MAXVAL(AnzSites),1:N) )
ALLOCATE( FLadung(1:3,MAXVAL(AnzLadungen),1:N) )
ALLOCATE( Komp(1:N) )
ALLOCATE( Liste(1:NmaxZelle,0:Zellen(1)-1,0:Zellen(2)-1,0:Zellen(3)-1) )
ALLOCATE( M(1:3,1:N) )
ALLOCATE( NZelle(0:Zellen(1)-1,0:Zellen(2)-1,0:Zellen(3)-1) )
ALLOCATE( r(1:3,1:N) )
ALLOCATE( v(1:3,1:N) )

ALLOCATE( FSite_v(1:3,MAXVAL(AnzSites),1:N) )
ALLOCATE( F_v(1:3,1:N) )


!*****************************************************************************
