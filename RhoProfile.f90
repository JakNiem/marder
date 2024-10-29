!*****************************************************************************
!*** Parameters for density profile                                        ***
!*****************************************************************************

! liquid and vapour density
rhol = 0.0
DO i=5, 14, 1
  rhol = rhol + 0.1*rhoShells(i,1)
END DO

DO i=15, NSHells, 1
  IF ( ABS(rhol-rhoShells(i,1)) < 0.1*rhol ) THEN
    rhol = ((i-1)*rhol + rhoShells(i,1))/i
  END IF
END DO

rhov = 0.0
DO i=(NShells-10), (NSHells-1), 1
  rhov = rhov + 0.1*rhoShells(i,1)
END DO

DO i=1, (NSHElls-20), 1
  j = NShells-10-i
  IF ( ABS(rhov-rhoShells(j,1)) < 0.2*rhov ) THEN
    rhov = ((10+i-1)*rhov + rhoShells(j,1))/(10+i)
  END IF
END DO

! D0 with 1090 (Baidakov et al.)
r10 = rhov + 0.1*(rhol-rhov)
r90 = rhov + 0.9*(rhol-rhov) 
DO i=1, (NSHells-10), 1
  IF ( rhoShells(i,1) > r90 ) THEN
    Dmin = RShells(i)
  END IF
END DO

DO i=1, (NSHElls-10), 1
  j = NShells-i
  IF ( rhoShells(j,1) < r10 ) THEN
    Dmax = RShells(j)
  END IF
END DO

D0 = (Dmax-Dmin)

! R0 
R0 = Dmin + 0.5*D0

! tanh density profile
DO i=1, NSHells, 1
  RhoShellsT(i,1) = RhoP(RShells(i), rhov, rhol, D0, R0)
END DO

OPEN(30,FILE=TRIM(Dateiname)//'_tanh.csv',POSITION='APPEND',STATUS='OLD')
  WRITE(30,'(I15,";",F12.8,";",F12.8,";",F12.8,";",F12.8)') IntSchritt, rhov, rhol, D0, R0
CLOSE(30)

! shift for artificial system
DO i=1, NShells, 1
  IF ( rhoShellsT(i,1) .GE. 1.02*rhov) THEN
    rhoShellsT(i,1) = rhoShellsT(i,1) - rhov
  ELSE
    rhoShellsT(i,1) = 0.0
  END IF
END DO

! homogeneous corrections
UpotKorrLJ = 0.0
VirialKorrLJ = 0.0 
DO i=1,AnzKomp,1
  DO j=1,AnzKomp,1
    DO Si=1,AnzSites(i),1
      DO Sj=1,AnzSites(j),1
        tau1 = SQRT(DOT_PRODUCT( rSiteBody(:,Si,i),rSiteBody(:,Si,i) ))
        tau2 = SQRT(DOT_PRODUCT( rSiteBody(:,Sj,j),rSiteBody(:,Sj,j) ))
        IF ( (tau1.EQ.0.0).AND.(tau2.EQ.0.0) ) THEN
          UpotKorrLJ = UpotKorrLJ&
              + rhov & 
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TICCu(-6,rc,sigma2(Si,i,Sj,j))&
                 - TICCu(-3,rc,sigma2(Si,i,Sj,j)) )
          VirialKorrLJ = VirialKorrLJ&
              + rhov &
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TICCp(-6,rc,sigma2(Si,i,Sj,j))&
                 - TICCp(-3,rc,sigma2(Si,i,Sj,j)) )
        END IF
        IF ( (tau1.EQ.0.0).NEQV.(tau2.EQ.0.0) ) THEN
          tau = MAX(tau1,tau2)   
          UpotKorrLJ = UpotKorrLJ&
              + rhov & 
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TICSu(-6,rc,sigma2(Si,i,Sj,j),tau)&
                 - TICSu(-3,rc,sigma2(Si,i,Sj,j),tau) )
          VirialKorrLJ = VirialKorrLJ&
              + rhov &
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TICSp(-6,rc,sigma2(Si,i,Sj,j),tau)&
                 - TICSp(-3,rc,sigma2(Si,i,Sj,j),tau) )
        END IF
        IF ( (tau1.NE.0.0).AND.(tau2.NE.0.0) ) THEN
          UpotKorrLJ = UpotKorrLJ&
              + rhov & 
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TISSu(-6,rc,sigma2(Si,i,Sj,j),tau1,tau2)&
                 - TISSu(-3,rc,sigma2(Si,i,Sj,j),tau1,tau2) )
          VirialKorrLJ = VirialKorrLJ&
              + rhov &
              * (epsilon24(Si,i,Sj,j)/6.0)&
              * (  TISSp(-6,rc,sigma2(Si,i,Sj,j),tau1,tau2)&
                 - TISSp(-3,rc,sigma2(Si,i,Sj,j),tau1,tau2) )
        END IF
      END DO
    END DO
  END DO
END DO
UpotKorrLJ = 2.0*pi*UpotKorrLJ 
VirialKorrLJ = 2.0*pi* VirialKorrLJ/3.0


