!*****************************************************************************
!***  Functions for spherical lrc                                          ***
!*****************************************************************************

! Density profile
REAL(8) FUNCTION RhoP(r, rhov, rhol, D0, R0)

  REAL(8) r, rhov, rhol, D0, R0

  RhoP = 0.5*(rhol+rhov) - 0.5*(rhol-rhov)*TANH(2*(r-R0)/D0)

END FUNCTION RhoP

! h_cs(2l,r)
REAL(8) FUNCTION SICSu(l,r,tau) 

  INTEGER l
  REAL(8) r, tau

  SICSu = ( (r+tau)**(2*l+3) - (r-tau)**(2*l+3) ) &
          / ( 4*tau*(l+1)*(2*l+3) )

END FUNCTION SICSu

! s_cs
REAL(8) FUNCTION CS(l,r,tau)

  INTEGER l
  REAL(8) r, tau

  CS = ( (r+tau)**(2*l+2) - (r-tau)**(2*l+2) ) / ( 4*r*tau*(l+1) )

END FUNCTION CS

! h_ss(2l,r) 
REAL(8) FUNCTION SISSu(l,r,tau1,tau2)

  INTEGER l
  REAL(8) r, tau1, tau2, tauMinus, tauPlus
  tauPlus = tau1+tau2
  tauMinus = tau1-tau2

  SISSu =  ( (r+tauPlus)**(2*l+4) - (r+tauMinus)**(2*l+4) &
             - (r-tauMinus)**(2*l+4) + (r-tauPlus)**(2*l+4) ) &
            /( 8*tau1*tau2*(l+1)*(2*l+3)*(2*l+4) )

END FUNCTION SISSu

! s_ss
REAL(8) FUNCTION SS(l,r,tau1,tau2)

  INTEGER l
  REAL(8) r, tau1, tau2, tauMinus, tauPlus
  tauPlus = tau1+tau2
  tauMinus = tau1-tau2

  SS = ( (r+tauPlus)**(2*l+3)-(r+tauMinus)**(2*l+3) &
        -(r-tauMinus)**(2*l+3) + (r-tauPlus)**(2*l+3) ) &
        /( 8*tau1*tau2*r*(l+1)*(2*l+3) )

END FUNCTION SS

! last term of eq. (55)
REAL(8) FUNCTION SSLN(r,tau1,tau2)

  REAL(8) r, tau1, tau2, tauMinus, tauPlus
  tauPlus = tau1+tau2
  tauMinus = tau1-tau2

  SSLN = LOG( (r*r-tauPlus*tauPlus)/(r*r-tauMinus*tauMinus) )/( 48*tau1*tau2 )

END FUNCTION SSLN
