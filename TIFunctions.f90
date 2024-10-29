!*****************************************************************************
!***  Functions for integrations in homogeneous lrc                        ***
!***  Note: must be called with sigma squared                              ***
!***  1) Integrals for potential energy                                    ***
!***  2) Integrals for pressure                                            ***
!*****************************************************************************

! 1) Integrals for potential energy
!*****************************************************************************
REAL(8) FUNCTION TICCu(n,rc,sigma2)
  INTEGER n
  REAL(8) rc,sigma2
  TICCu = -  ( rc**(2*n+3) )&
            /( sigma2**(n)*(2*n+3) )
END FUNCTION TICCu


REAL(8) FUNCTION TICSu(n,rc,sigma2,tau)
  INTEGER n
  REAL(8) rc,sigma2,tau
  TICSu = -  ( (rc+tau)**(2*n+3) - (rc-tau)**(2*n+3) ) * rc&
            /( 4*sigma2**(n)*tau*(n+1)*(2*n+3) )&
          +  ( (rc+tau)**(2*n+4) - (rc-tau)**(2*n+4) )&
            /( 4*sigma2**(n)*tau*(n+1)*(2*n+3)*(2*n+4) )
END FUNCTION TICSu


REAL(8) FUNCTION TISSu(n,rc,sigma2,tau1,tau2)
  INTEGER n
  REAL(8) rc,sigma2,tau1,tau2,tauMinus,tauPlus
  tauPlus = tau1+tau2
  tauMinus = tau1-tau2
  TISSu = -  (   (rc+tauPlus)**(2*n+4) - (rc+tauMinus)**(2*n+4)&
               - (rc-tauMinus)**(2*n+4) + (rc-tauPlus)**(2*n+4) ) * rc&
            /( 8*sigma2**(n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4) )&
          +  (   (rc+tauPlus)**(2*n+5) - (rc+tauMinus)**(2*n+5)&
               - (rc-tauMinus)**(2*n+5) + (rc-tauPlus)**(2*n+5) )&
            /( 8*sigma2**(n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4)*(2*n+5) )
END FUNCTION TISSu



! 2) Integrals for pressure  
!*****************************************************************************
REAL(8) FUNCTION TICCp(n,rc,sigma2)
  INTEGER n
  REAL(8) rc,sigma2
  TICCp = 2*n * TICCu(n,rc,sigma2)
END FUNCTION TICCp


REAL(8) FUNCTION TICSp(n,rc,sigma2,tau)
  INTEGER n
  REAL(8) rc,sigma2,tau
  TICSp = -  ( (rc+tau)**(2*n+2) - (rc-tau)**(2*n+2) ) * rc**(2)&
            /( 4*sigma2**(n)*tau*(n+1) )&
          - 3*TICSu(n,rc,sigma2,tau)
END FUNCTION TICSp


REAL(8) FUNCTION TISSp(n,rc,sigma2,tau1,tau2)
  INTEGER n
  REAL(8) rc,sigma2,tau1,tau2,tauMinus,tauPlus
  tauPlus = tau1+tau2
  tauMinus = tau1-tau2
  TISSp = -  (   (rc+tauPlus)**(2*n+3) - (rc+tauMinus)**(2*n+3)&
               - (rc-tauMinus)**(2*n+3) + (rc-tauPlus)**(2*n+3) ) * rc**(2)&
            /( 8*sigma2**(n)*tau1*tau2*(n+1)*(2*n+3) )&
          - 3*TISSu(n,rc,sigma2,tau1,tau2)
END FUNCTION TISSp
