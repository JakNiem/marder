!*****************************************************************************
!***  Shell density and corrections per particle                           ***
!*****************************************************************************

! Determine current density in every time step
rhoShellsTemp = 0.0

DO i=1, N, 1
  F(:,i) = ( SystemCenter-r(:,i) )
  ksi2(i) = (r(1,i)-SystemCenter(1))*(r(1,i)-SystemCenter(1)) &
         + (r(2,i)-SystemCenter(2))*(r(2,i)-SystemCenter(2)) &
         + (r(3,i)-SystemCenter(3))*(r(3,i)-SystemCenter(3))
  ksi(i) = SQRT(ksi2(i))
  realk = ksi(i)/drShells
  lowerS(i) = MIN( INT(FLOOR((realk-deltaShells))), NShells) 
  interS(i) = MAX( INT(CEILING(ABS(realk-deltaShells))), 1 )
  upperS(i) = MIN( INT(CEILING(realk+deltaShells)), (NShells+1))
  k = NINT( realk )
  IF (k .GE. NShells) THEN
    rhoShellsTemp(NShells,Komp(i)) = rhoShellsTemp(NShells,Komp(i)) + 1
    PartShells(i) = NShells
  ELSE IF (k == 0) THEN
    rhoShellsTemp(1, Komp(i)) = rhoShellsTemp(1,Komp(i)) + 1
    PartShells(i) = 1
  ELSE 
    rhoShellsTemp(k, Komp(i)) = rhoShellsTemp(k,Komp(i)) + 1
    PartShells(i) = k
  END IF
END DO

DO i=1, NShells, 1
  rhoShellsTemp(i,:) = rhoShellsTemp(i,:) / VShells(i)
END DO

! average density over total simulation run
DO i=1, AnzKomp, 1
  rhoShellsAve = rhoShellsAve + rhoShellsTemp(:,i)
END DO

! densities of last 50.000 time steps 
MeanIndex = MOD( (FLOOR( TotalStep/1000.0 )), NSMean) + 1
IF (MOD(TotalStep, 1000) == 0) THEN
  rhoShellsMean(:,:,MeanIndex) = 0.0
END IF
rhoShellsMean(:,:,MeanIndex) = rhoShellsMean(:,:,MeanIndex) + rhoShellsTemp

! average density over last 50.000 time steps for tanh fit
IF (MOD( IntSchritt, 1000) == 0) THEN
  rhoShells = 0.0
  DO i=1, NSMean, 1
    rhoShells = rhoShells + 0.001/NSMean*rhoShellsMean(:,:,i)
  END DO

  INCLUDE 'RhoProfile.f90'
  INCLUDE 'ShellCorr.f90'

END IF

! Read corrections for each particle
UCorrSum = 0.0
DO i=1, N, 1    ! Loop over Particles
  UCorrSum = UCorrSum + UShells_Mean(PartShells(i))
  F(:,i) = F(:,i)/ksi(i) * FShells_Mean(PartShells(i))
  PShells_N(PartShells(i)) = PShells_N(PartShells(i)) - 0.5*PNShells_Mean(PartShells(i)) - VirialKorrLJ
  PShells_T(PartShells(i)) = PShells_T(PartShells(i)) - 0.25*PTShells_Mean(PartShells(i)) - VirialKorrLJ
END DO
