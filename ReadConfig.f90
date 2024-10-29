!*****************************************************************************
!***  1) Read initial configuration                                        ***
!***  2) Initial list for cell algorithm                                   ***  
!***  3) Zero-th step: forces/momenta = 0                                  ***
!***  4) Initial density profile                                           ***
!*****************************************************************************

! 1) Read initial configuration
!*****************************************************************************
ALLOCATE( q(0:3,1:N) )
OPEN(11,FILE=InputFile,STATUS='OLD')
  DO i = 1, N, 1
    READ(11, '(13F12.8)') r(1,i), r(2,i), r(3,i), q(0,i), q(1,i), q(2,i), q(3,i),&
      v(1,i), v(2,i), v(3,i), D(1,i), D(2,i), D(3,i)
    Komp(i) = 1
  END DO
CLOSE(11)


! 2) Initial list for cell algorithm
!*****************************************************************************
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


! 3) Zero-th step: forces/momenta = 0
!*****************************************************************************
F = 0.0
M = 0.0


! 4) Initial density profile and corrections
!*****************************************************************************

rhoShellsTemp = 0.0

DO i=1, N, 1
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
    PartShells(i) = MIN( k, NShells )
  ELSE ! k=0 is ignored
    rhoShellsTemp(k, Komp(i)) = rhoShellsTemp(k,Komp(i)) + 1
    PartShells(i) = k
  END IF
END DO
DO i=1, NShells, 1
  rhoShellsTemp(i,:) = rhoShellsTemp(i,:) / VShells(i)
END DO

rhoShells = rhoShellsTemp

DO i=1, NSMean, 1
  rhoShellsMean(:,:,i) = 1000*rhoShellsTemp
END DO

INCLUDE 'RhoProfile.f90'
INCLUDE 'ShellCorr.f90'
