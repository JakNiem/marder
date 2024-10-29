!*****************************************************************************
!***  Average correction terms per shell                                   ***
!*****************************************************************************

! Corrections per Shell
UShells_Mean = 0.0
FShells_Mean = 0.0
PNShells_Mean = 0.0
PTShells_Mean = 0.0

DO i=1, N, 1                    ! Loop over Particles
  UCorrShells = 0.0             ! Corrections per particle
  FCorrShells = 0.0
  PNCorrShells = 0.0
  PTCorrShells = 0.0
  Ki = Komp(i)                  ! Component of current Particle
  DO Kj=1,AnzKomp,1             ! Loop over Components
    DO Si=1,AnzSites(Ki),1      ! Loop over Sites of current Component
      DO Sj=1,AnzSites(Kj),1    ! Loop over Sites of all Components
        tau1 = SQRT(DOT_PRODUCT( rSiteBody(:,Si,Ki),rSiteBody(:,Si,Ki) ))
        tau2 = SQRT(DOT_PRODUCT( rSiteBody(:,Sj,Kj),rSiteBody(:,Sj,Kj) ))
        sigma6 = sigma2(Si,Ki,Sj,Kj)*sigma2(Si,Ki,Sj,Kj)*sigma2(Si,Ki,Sj,Kj)
        factorU = -Pi*drShells*epsilon24(Si,Ki,Sj,Kj)*sigma6/(6*ksi(i))
        factorF = -factorU/ksi(i)
        factorP = 0.5*factorF/ksi(i)
        IF ( (tau1.EQ.0.0).AND.(tau2.EQ.0.0) ) THEN ! Center-Center
          DO j=1, (lowerS(i)), 1 ! Loop over Shells with smaller radius
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = ksi(i) - RShells(j) ! lower bound
              rlowInv = 1/rlow
              rlowInv2 = rlowInv * rlowInv
              rdashInv = 1/(RShells(j) + ksi(i))
              UCorrTemp = sigma6*0.2*(rdashInv**10 - rlowInv**10 ) &
                        - 0.5*(rdashInv**4 - rlowInv**4 )
              IF ( rlow < rcmax) THEN
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) ) ! upper bound
                rdashInv = 1/rdash
                rdashInv2 = rdashInv * rdashInv
                FCorrTemp = sigma6 &
                          * ( (6/5*rdash*rdash + ksi2(i) - RShells2(j)) * rdashInv2**6 &
                            - (6/5*rlow*rlow + ksi2(i) - RShells2(j)) * rlowInv2**6 ) &
                            - (1.5*rdash*rdash + ksi2(i) - RShells2(j)) * rdashInv2**3 &
                            + (1.5*rlow*rlow + ksi2(i) - RShells2(j)) * rlowInv2**3 
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = 1.5*sigma6 * (rdashInv2**4 - rlowInv2**4) &
                       - 3 * (rdashInv2 - rlowInv2) & 
                       + 2*(ksi2(i)-RShells2(j)) * ( 6/5*sigma6 * (rdashInv2**5 - rlowInv2**5) &
                       - 1.5 * (rdashInv2**2 - rlowInv2**2) ) & 
                       + (ksi2(i)-RShells2(j))**2 * ( sigma6 * ( rdashInv2**6 - rlowInv2**6) & 
                       - (rdashInv2**3 - rlowInv2**3) )
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = 6/5*sigma6 * (rdashInv2**5 - rlowInv2**5) &
                           - 1.5 * (rdashInv2**2 - rlowInv2**2) 
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF
              UCorrShells = UCorrShells + UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
            END IF
          END DO
          DO j=upperS(i), NSHells, 1 ! Loop over Shells with larger radius
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = RShells(j) - ksi(i)
              rlowInv = 1/rlow
              rlowInv2 = rlowInv * rlowInv
              rdashInv = 1/(RShells(j) + ksi(i))
              UCorrTemp = sigma6*0.2*(rdashInv**10 - rlowInv**10 ) &
                        - 0.5*(rdashInv**4 - rlowInv**4 )
              IF ( rlow < rcmax) THEN
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) )
                rdashInv = 1/rdash
                rdashInv2 = rdashInv * rdashInv
                FCorrTemp = sigma6 &
                          * ( (6/5*rdash*rdash + ksi2(i) - RShells2(j)) * rdashInv2**6 &
                            - (6/5*rlow*rlow + ksi2(i) - RShells2(j)) * rlowInv2**6 ) &
                            - (1.5*rdash*rdash + ksi2(i) - RShells2(j)) * rdashInv2**3 &
                            + (1.5*rlow*rlow + ksi2(i) - RShells2(j)) * rlowInv2**3 
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = 1.5*sigma6 * (rdashInv2**4 - rlowInv2**4) &
                            - 3 * (rdashInv2 - rlowInv2) & 
                            + 2*(ksi2(i)-RShells2(j)) * ( 6/5*sigma6 * (rdashInv2**5 - rlowInv2**5) &
                            - 1.5 * (rdashInv2**2 - rlowInv2**2) ) & 
                            + (ksi2(i)-RShells2(j))**2 * ( sigma6 * ( rdashInv2**6 - rlowInv2**6) & 
                            - (rdashInv2**3 - rlowInv2**3) )
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = 6/5*sigma6 * (rdashInv2**5 - rlowInv2**5) &
                            - 1.5 * (rdashInv2**2 - rlowInv2**2)
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF
              UCorrShells = UCorrShells + UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
            END IF
          END DO
          DO j=(interS(i)), (upperS(i)-1), 1 ! Loop over partly contributing Shells
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = rc
              rlowInv = 1/rlow
              rlowInv2 = rlowInv * rlowInv
              rdashInv = 1/(RShells(j) + ksi(i))
              UCorrTemp = sigma6*0.2*(rdashInv**10 - rlowInv**10 ) &
                        - 0.5*(rdashInv**4 - rlowInv**4 )
              IF (rlow<rcmax) THEN
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) )
                rdashInv = 1/rdash
                rdashInv2 = rdashInv * rdashInv
                FCorrTemp = sigma6 &
                        * ( (6/5*rdash*rdash + ksi2(i) - RShells2(j)) * rdashInv2**6 &
                          - (6/5*rlow*rlow + ksi2(i) - RShells2(j)) * rlowInv2**6 ) &
                          - (1.5*rdash*rdash + ksi2(i) - RShells2(j)) * rdashInv2**3 &
                          + (1.5*rlow*rlow + ksi2(i) - RShells2(j)) * rlowInv2**3 
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = 1.5*sigma6 * (rdashInv2**4 - rlowInv2**4) &
                          - 3 * (rdashInv2 - rlowInv2) & 
                          + 2*(ksi2(i)-RShells2(j)) * ( 6/5*sigma6 * (rdashInv2**5 - rlowInv2**5) &
                          - 1.5 * (rdashInv2**2 - rlowInv2**2) ) & 
                          + (ksi2(i)-RShells2(j))**2 * ( sigma6 * ( rdashInv2**6 - rlowInv2**6) & 
                          - (rdashInv2**3 - rlowInv2**3) )
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = 6/5*sigma6 * (rdashInv2**5 - rlowInv2**5) &
                          - 1.5 * (rdashInv2**2 - rlowInv2**2) 
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF
              UCorrShells = UCorrShells + UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
            END IF
          END DO
        END IF
        IF ( (tau1.EQ.0.0).NEQV.(tau2.EQ.0.0) ) THEN ! Center-Site
          tau = MAX(tau1,tau2)   
          DO j=1, (lowerS(i)), 1 ! Loop over Shells with smaller radius
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = ksi(i) - RShells(j)
              UCorrTemp = - SICSu(-3,RShells(j) + ksi(i),tau) &
                        + SICSu(-3,rlow,tau)
              UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
              IF ( rlow < rcmax ) THEN
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) )
                rdash2 = rdash*rdash
                rlow2 = rlow*rlow
                FCorrTemp = (rdash2 + ksi2(i) - RShells2(j)) &
                        * ( - CS(-3,rdash,tau) ) &
                        + (rlow2 + ksi2(i) - RShells2(j) ) &
                        * ( CS(-3,rlow,tau)) &
                        + 2*( SICSu(-3,rdash,tau) &
                            - SICSu(-3,rlow,tau) )
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = (rdash2 + ksi2(i) -RShells2(j))**2 &
                        * ( - CS(-3,rdash,tau) ) &
                        + (rlow2 + ksi2(i) -RShells2(j))**2 &
                        * ( CS(-3,rlow,tau) ) &
                        + 4 * (ksi2(i) - RShells2(j) + rdash2) &
                        * ( SICSu(-3,rdash,tau) ) &
                        - 4 * (ksi2(i) -RShells2(j) + rlow2) &
                        * ( SICSu(-3,rlow,tau) ) &
                        - 4 * ( rdash2*CS(-2,rdash,tau)/6 + 1/((rdash2-tau*tau)*12) ) &
                        + 4 * ( rlow2*CS(-2,rlow,tau)/6 + 1/((rlow2-tau*tau)*12) )
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = - CS(-3,rdash,tau)*rdash2 &
                        + CS(-3,rlow,tau)*rlow2 &
                        + 2 * ( SICSu(-3,rdash,tau) &
                              - SICSu(-3,rlow,tau) )
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF
            END IF
          END DO
          DO j=upperS(i), NSHells, 1 ! Loop over Shells with larger radius
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = RShells(j) - ksi(i)
              UCorrTemp = - SICSu(-3,RShells(j) + ksi(i),tau) &
                        + SICSu(-3,rlow,tau)
              UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
              IF ( rlow < rcmax ) THEN
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) )
                rdash2 = rdash*rdash
                rlow2 = rlow*rlow
                FCorrTemp = (rdash2 + ksi2(i) - RShells2(j)) &
                        * ( - CS(-3,rdash,tau) ) &
                        + (rlow2 + ksi2(i) - RShells2(j) ) &
                        * ( CS(-3,rlow,tau)) &
                        + 2*( SICSu(-3,rdash,tau) &
                            - SICSu(-3,rlow,tau) )
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = (rdash2 + ksi2(i) -RShells2(j))**2 &
                        * ( - CS(-3,rdash,tau) ) &
                        + (rlow2 + ksi2(i) -RShells2(j))**2 &
                        * ( CS(-3,rlow,tau) ) &
                        + 4 * (ksi2(i) - RShells2(j) + rdash2) &
                        * ( SICSu(-3,rdash,tau) ) &
                        - 4 * (ksi2(i) -RShells2(j) + rlow2) &
                        * ( SICSu(-3,rlow,tau) ) &
                        - 4 * ( rdash2*CS(-2,rdash,tau)/6 + 1/((rdash2-tau*tau)*12) ) &
                        + 4 * ( rlow2*CS(-2,rlow,tau)/6 + 1/((rlow2-tau*tau)*12) )
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = - CS(-3,rdash,tau)*rdash2 &
                        + CS(-3,rlow,tau)*rlow2 &
                        + 2 * ( SICSu(-3,rdash,tau) &
                              - SICSu(-3,rlow,tau) )
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF 
            END IF
          END DO
          DO j=(interS(i)), (upperS(i)-1), 1 ! Loop over partly contributing Shells
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = rc
              UCorrTemp = - SICSu(-3,RShells(j) + ksi(i),tau) &
                        + SICSu(-3,rlow,tau)
              UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
              IF ( rlow < rcmax ) THEN
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) )
                rdash2 = rdash*rdash
                rlow2 = rlow*rlow
                FCorrTemp = (rdash2 + ksi2(i) - RShells2(j)) &
                        * ( - CS(-3,rdash,tau) ) &
                        + (rlow2 + ksi2(i) - RShells2(j) ) &
                        * ( CS(-3,rlow,tau)) &
                        + 2*( SICSu(-3,rdash,tau) &
                            - SICSu(-3,rlow,tau) )
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = (rdash2 + ksi2(i) -RShells2(j))**2 &
                        * ( - CS(-3,rdash,tau) ) &
                        + (rlow2 + ksi2(i) -RShells2(j))**2 &
                        * ( CS(-3,rlow,tau) ) &
                        + 4 * (ksi2(i) - RShells2(j) + rdash2) &
                        * ( SICSu(-3,rdash,tau) ) &
                        - 4 * (ksi2(i) -RShells2(j) + rlow2) &
                        * ( SICSu(-3,rlow,tau) ) &
                        - 4 * ( rdash2*CS(-2,rdash,tau)/6 + 1/((rdash2-tau*tau)*12) ) &
                        + 4 * ( rlow2*CS(-2,rlow,tau)/6 + 1/((rlow2-tau*tau)*12) )
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = - CS(-3,rdash,tau)*rdash2 &
                        + CS(-3,rlow,tau)*rlow2 &
                        + 2 * ( SICSu(-3,rdash,tau) &
                              - SICSu(-3,rlow,tau) )
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF
            END IF
          END DO
        END IF
        IF ( (tau1.NE.0.0).AND.(tau2.NE.0.0) ) THEN ! Site-Site
          DO j=1, (lowerS(i)), 1 ! Loop over Shells with smaller radius
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = ksi(i) - RShells(j)
              UCorrTemp = - SISSu(-3,RShells(j) + ksi(i),tau1,tau2) &
                        + SISSu(-3,rlow,tau1,tau2)
              UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
              IF ( rlow < rcmax ) THEN
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) )
                rdash2 = rdash*rdash
                rlow2 = rlow*rlow
                FCorrTemp = -(rdash2 + ksi2(i) - RShells2(j)) &
                        * SS(-3,rdash,tau1,tau2) &
                        + (rlow2 + ksi2(i) - RShells2(j) ) &
                        * ( SS(-3,rlow,tau1,tau2)) &
                        + 2*( SISSu(-3,rdash,tau1,tau2) &
                            - SISSu(-3,rlow,tau1,tau2) )
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = (rdash2 + ksi2(i) -RShells2(j))**2 &
                           * ( - SS(-3,rdash,tau1,tau2) ) &
                           + (rlow2 + ksi2(i) -RShells2(j))**2 &
                           * ( SS(-3,rlow,tau1,tau2) ) &
                           + 4 * (ksi2(i) - RShells2(j) + rdash2) &
                           * ( SISSu(-3,rdash,tau1,tau2) ) &
                           - 4 * (ksi2(i) -RShells2(j) + rlow2) &
                           * ( SISSu(-3,rlow,tau1,tau2) ) &
                           - 4 * ( rdash2*SS(-2,rdash,tau1,tau2)/6 - SSLN(rdash,tau1,tau2) ) &
                           + 4 * ( rlow2*SS(-2,rlow,tau1,tau2)/6 - SSLN(rlow,tau1,tau2) )          
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = - SS(-3,rdash,tau1,tau2)*rdash2 &
                           + SS(-3,rlow,tau1,tau2)*rlow2 &
                           + 2 * ( SISSu(-3,rdash,tau1,tau2) &
                                 - SISSu(-3,rlow,tau1,tau2) )
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF
            END IF
          END DO
          DO j=upperS(i), NSHells, 1 ! Loop over Shells with larger radius
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = RShells(j) - ksi(i) 
              UCorrTemp = - SISSu(-3,RShells(j) + ksi(i),tau1,tau2) &
                              + SISSu(-3,rlow,tau1,tau2)
              UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
              IF ( rlow < rcmax ) THEN
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) )
                rdash2 = rdash*rdash
                rlow2 = rlow*rlow
                FCorrTemp = -(rdash2 + ksi2(i) - RShells2(j)) &
                        * SS(-3,rdash,tau1,tau2) &
                        + (rlow2 + ksi2(i) - RShells2(j) ) &
                        * ( SS(-3,rlow,tau1,tau2) ) &
                        + 2*( SISSu(-3,rdash,tau1,tau2) &
                            - SISSu(-3,rlow,tau1,tau2) )
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = (rdash2 + ksi2(i) -RShells2(j))**2 &
                        * ( - SS(-3,rdash,tau1,tau2) ) &
                        + (rlow2 + ksi2(i) -RShells2(j))**2 &
                        * ( SS(-3,rlow,tau1,tau2) ) &
                        + 4 * (ksi2(i) - RShells2(j) + rdash2) &
                        * ( SISSu(-3,rdash,tau1,tau2)) &
                        - 4 * (ksi2(i) -RShells2(j) + rlow2) &
                        * ( SISSu(-3,rlow,tau1,tau2)) &
                        - 4 * ( rdash2*SS(-2,rdash,tau1,tau2)/6 - SSLN(rdash,tau1,tau2) ) &
                        + 4 * ( rlow2*SS(-2,rlow,tau1,tau2)/6 - SSLN(rlow,tau1,tau2) )        
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = - SS(-3,rdash,tau1,tau2)*rdash2 &
                        + SS(-3,rlow,tau1,tau2)*rlow2 &
                        + 2 * ( SISSu(-3,rdash,tau1,tau2) &
                              - SISSu(-3,rlow,tau1,tau2) )
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF
            END IF
          END DO
          DO j=(interS(i)), (upperS(i)-1), 1 ! Loop over partly contributing Shells
            IF (rhoShellsT(j,Komp(i)) .NE. 0.0) THEN
              rlow = rc
              UCorrTemp = - SISSu(-3,RShells(j) + ksi(i),tau1,tau2) &
                          + SISSu(-3,rlow,tau1,tau2)
              UCorrShells = UCorrShells - 2*UCorrTemp*factorU*rhoShellsT(j,Komp(i))*RShells(j)
              IF (rlow < rcmax ) THEN 
                rdash = MIN( rcmax,(RShells(j) + ksi(i)) )
                rdash2 = rdash*rdash
                rlow2 = rlow*rlow
                FCorrTemp = -(rdash2 + ksi2(i) - RShells2(j)) &
                          * SS(-3,rdash,tau1,tau2) &
                          + (rlow2 + ksi2(i) - RShells2(j) ) &
                          * ( SS(-3,rlow,tau1,tau2) ) &
                          + 2*( SISSu(-3,rdash,tau1,tau2) &
                              - SISSu(-3,rlow,tau1,tau2) )
                FCorrShells = FCorrShells + FCorrTemp*factorF*rhoShellsT(j,Komp(i))*RShells(j)
                PNCorrTemp = (rdash2 + ksi2(i) -RShells2(j))**2 &
                          * ( - SS(-3,rdash,tau1,tau2) ) &
                          + (rlow2 + ksi2(i) -RShells2(j))**2 &
                          * ( SS(-3,rlow,tau1,tau2) ) &
                          + 4 * (ksi2(i) - RShells2(j) + rdash2) &
                          * ( SISSu(-3,rdash,tau1,tau2)) &
                          - 4 * (ksi2(i) -RShells2(j) + rlow2) &
                          * ( SISSu(-3,rlow,tau1,tau2)) &
                          - 4 * ( rdash2*SS(-2,rdash,tau1,tau2)/6 - SSLN(rdash,tau1,tau2) ) &
                          + 4 * ( rlow2*SS(-2,rlow,tau1,tau2)/6 - SSLN(rlow,tau1,tau2) )
                PNCorrShells = PNCorrShells + PNCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
                PTCorrTemp = - SS(-3,rdash,tau1,tau2)*rdash2 &
                        + SS(-3,rlow,tau1,tau2)*rlow2 &
                        + 2 * ( SISSu(-3,rdash,tau1,tau2) &
                              - SISSu(-3,rlow,tau1,tau2) )
                PTCorrShells = PTCorrShells + 4*ksi2(i)*PTCorrTemp*factorP*rhoShellsT(j,Komp(i))*RShells(j)
              END IF
            END IF
          END DO
        END IF
      END DO
    END DO
  END DO
  PTCorrShells = PTCorrShells - PNCorrShells
  PNShells_Mean(PartShells(i)) = PNShells_Mean(PartShells(i)) + PNCorrShells
  PTShells_Mean(PartShells(i)) = PTShells_Mean(PartShells(i)) + PTCorrShells
  UShells_Mean(PartShells(i)) = UShells_Mean(PartShells(i)) + UCorrShells
  FShells_Mean(PartShells(i)) = FShells_Mean(PartShells(i)) + FCorrShells
END DO

! Correction divided by number of particles in shell
DO i=1,NShells,1
  IF (RhoShellsTemp(i,1) .NE. 0.0) THEN
    PNShells_Mean(i) = PNShells_Mean(i) /(RhoShellsTemp(i,1)*VShells(i))
    PTShells_Mean(i) = PTShells_Mean(i) /(RhoShellsTemp(i,1)*VShells(i))
    UShells_Mean(i) = UShells_Mean(i) /(RhoShellsTemp(i,1)*VShells(i))
    FShells_Mean(i) = FShells_Mean(i) /(rhoShellsTemp(i,1)*VShells(i))
  END IF 
END DO
