!*****************************************************************************
!***  Force calculation IF (drMol2.LT.rc2) ...                             ***
!*****************************************************************************

WWCount = WWCount + 1

DO Si=1,AnzSites(Komp(i)),1
  DO Sj =1,AnzSites(Komp(j)),1
    ! Site(i) -- Site(j)
    !*************************************************************************

    ! r_com = drMol
    ! r_ss = dr
    dr = drMol + drSite(:,Si,i) - drSite(:,Sj,j) 
    Invdr2 = 1.0/DOT_PRODUCT(dr,dr)
    LJTerm6 = sigma2(Si,Komp(i),Sj,Komp(j)) * Invdr2
    LJTerm6 = LJTerm6*LJTerm6*LJTerm6
    LJTerm12 = LJTerm6*LJTerm6

    ! f_ij/r_ss = FijTemp
    FijTemp = epsilon24(Si,Komp(i),Sj,Komp(j)) &
            * (LJTerm12+LJTerm12 - LJTerm6) * Invdr2
    Fij =  FijTemp * dr
    FSite(:,Si,i) = FSite(:,Si,i) + Fij
    FSite(:,Sj,j) = FSite(:,Sj,j) - Fij
    UpotSite =  UpotSite +  epsilon24(Si,Komp(i),Sj,Komp(j))&
                           *(LJTerm12-LJTerm6)
    VirialTemp = DOT_PRODUCT(drMol,Fij) 
    VirialSite = VirialSite + VirialTemp

    ! Virial decomposition for spherical interfaces
    ! part i 
    ksi_i = (SystemCenter-r(:,i))/ksi(i)
    dr_nor = DOT_PRODUCT(dr,ksi_i)
    drMol_nor = DOT_PRODUCT(drMol,ksi_i)
    PShells_N(PartShells(i)) = PShells_N(PartShells(i)) + 0.5*(FijTemp*dr_nor*drMol_nor)
    rss_tan = dr - dr_nor*ksi_i
    rcom_tan = drMol - drMol_nor*ksi_i
    PShells_T(PartShells(i)) = PShells_T(PartShells(i)) + 0.25* FijTemp*DOT_PRODUCT(rss_tan,rcom_tan) 
    ! part j 
    ksi_j = (SystemCenter-r(:,j))/ksi(j)
    dr_nor = DOT_PRODUCT(dr,ksi_j)
    drMol_nor = DOT_PRODUCT(drMol,ksi_j)
    PShells_N(PartShells(j)) = PShells_N(PartShells(j)) + 0.5*(FijTemp*dr_nor*drMol_nor)
    rss_tan = dr - dr_nor*ksi_j
    rcom_tan = drMol - drMol_nor*ksi_j
    PShells_T(PartShells(j)) = PShells_T(PartShells(j)) + 0.25* FijTemp*DOT_PRODUCT(rss_tan,rcom_tan) 


  END DO   ! Sj
END DO   ! Si


DO Si=1,AnzLadungen(Komp(i)),1
  DO Sj =1,AnzLadungen(Komp(j)),1
    ! Charge(i) -- Charge(j)
    !*************************************************************************
    dr = drMol + drLadung(:,Si,i) - drLadung(:,Sj,j)
    Invdr2 = 1.0/DOT_PRODUCT(dr,dr)
    Invdr1 = SQRT(Invdr2)
    E2 = absLadung(Si,Komp(i))*absLadung(Sj,Komp(j))
    Fij = E2*Invdr2*Invdr1*dr	 
    FLadung(:,Si,i) = FLadung(:,Si,i) + Fij
    FLadung(:,Sj,j) = FLadung(:,Sj,j) - Fij
    UpotE =  UpotE + E2*Invdr1
    VirialE = VirialE + DOT_PRODUCT(drMol,Fij)
  END DO   ! Sj
  
  DO Sj =1,AnzDipole(Komp(j)),1
    ! Charge(i) -- Dipole(j)
    !*************************************************************************
    dr = drMol + drLadung(:,Si,i) - drDipol(:,Sj,j)
    Invdr1 = 1.0/SQRT(DOT_PRODUCT(dr,dr))
    drE = dr * Invdr1
    CosTj = DOT_PRODUCT(eMy(1:3,Sj,j),drE)
    CosTj3 = 3.0 * CosTj
    TermU = absLadung(Si,Komp(i)) * absMy(Sj,Komp(j)) * Invdr1 * Invdr1
    TermU2 = TermU * Invdr1
 
    Fij = TermU2 * ( CosTj3 * drE - eMy(1:3,Sj,j) )
    FLadung(:,Si,i) = FLadung(:,Si,i) + Fij
    FDipol(:,Sj,j) = FDipol(:,Sj,j) - Fij
    eXrij(1) = eMy(2,Sj,j)*drE(3) - eMy(3,Sj,j)*drE(2)
    eXrij(2) = eMy(3,Sj,j)*drE(1) - eMy(1,Sj,j)*drE(3)
    eXrij(3) = eMy(1,Sj,j)*drE(2) - eMy(2,Sj,j)*drE(1)
    M(:,j) = M(:,j) - TermU * eXrij
    UpotE =  UpotE + TermU * CosTj
    VirialE = VirialE + DOT_PRODUCT(drMol,Fij)
  END DO   ! Sj
  
  DO Sj =1,AnzQuadrupole(Komp(j)),1
    ! Charge(i) -- Quadrupole(j)
    dr = drMol + drLadung(:,Si,i) - drQuadrupol(:,Sj,j)
    Invdr2 = 1.0/DOT_PRODUCT(dr,dr)
    Invdr1 = SQRT(Invdr2)    
    drE = dr * Invdr1
    CosTj = DOT_PRODUCT(eQ(1:3,Sj,j),drE)
    Cos2Tj = CosTj * CosTj
    TermU = 1.5 * absLadung(Si,Komp(i)) * absQ(Sj,Komp(j)) * Invdr1 * Invdr2

    Fij = TermU * Invdr1 * ( (5.0*Cos2Tj-1.0) * drE - 2.0*CosTj*eQ(1:3,Sj,j) )
    FLadung(:,Si,i) = FLadung(:,Si,i) + Fij
    FQuadrupol(:,Sj,j) = FQuadrupol(:,Sj,j) - Fij
    eXrij(1) = eQ(2,Sj,j)*drE(3) - eQ(3,Sj,j)*drE(2)
    eXrij(2) = eQ(3,Sj,j)*drE(1) - eQ(1,Sj,j)*drE(3)
    eXrij(3) = eQ(1,Sj,j)*drE(2) - eQ(2,Sj,j)*drE(1)        
    M(:,j) = M(:,j) - TermU * 2.0 * CosTj * eXrij
    UpotE =  UpotE + TermU * ( Cos2Tj - Drittel )
    VirialE = VirialE + DOT_PRODUCT(drMol,Fij)  
  END DO   ! Sj  
END DO   ! Si

DO Si=1,AnzDipole(Komp(i)),1
  DO Sj =1,AnzLadungen(Komp(j)),1
    ! Dipole(i) -- Charge(j)
    !*************************************************************************
    dr = drMol + drDipol(:,Si,i) - drLadung(:,Sj,j)
    Invdr1 = 1.0/SQRT(DOT_PRODUCT(dr,dr))
    drE = dr * Invdr1
    CosTi = DOT_PRODUCT(eMy(1:3,Si,i),drE)
    CosTi3 = 3.0 * CosTi
    TermU = absMy(Si,Komp(i)) * absLadung(Sj,Komp(j)) * Invdr1 * Invdr1
    TermU2 = TermU * Invdr1
  
    Fij = TermU2 * ( eMy(1:3,Si,i) - CosTi3 * drE )
    FDipol(:,Si,i) = FDipol(:,Si,i) + Fij
    FLadung(:,Sj,j) = FLadung(:,Sj,j) - Fij
    eXrij(1) = eMy(2,Si,i)*drE(3) - eMy(3,Si,i)*drE(2)
    eXrij(2) = eMy(3,Si,i)*drE(1) - eMy(1,Si,i)*drE(3)
    eXrij(3) = eMy(1,Si,i)*drE(2) - eMy(2,Si,i)*drE(1)    
    M(:,i) = M(:,i) + TermU * eXrij
    UpotMy =  UpotMy - TermU * CosTi
    VirialMy = VirialMy + DOT_PRODUCT(drMol,Fij)
  END DO   ! Sj

  DO Sj =1,AnzDipole(Komp(j)),1
    ! Dipole(i) -- Dipole(j)
    !*************************************************************************
    dr = drMol + drDipol(:,Si,i) - drDipol(:,Sj,j)
    Invdr1 = 1.0/SQRT(DOT_PRODUCT(dr,dr))
    My2 = absMy(Si,Komp(i))*absMy(Sj,Komp(j))
    MyFaktor =  My2*Invdr1*Invdr1*Invdr1
    CosTi = DOT_PRODUCT(eMy(1:3,Si,i),dr)*Invdr1
    CosTj = DOT_PRODUCT(eMy(1:3,Sj,j),dr)*Invdr1
    CosGij = DOT_PRODUCT(eMy(1:3,Si,i),eMy(1:3,Sj,j))
    TermU = MyFaktor * (CosGij - 3.0*CosTi*CosTj)
    PartialRijInvdr1 = -3.0*Invdr1*TermU * Invdr1
    PartialTiInvdr1 = -MyFaktor*3.0*CosTj * Invdr1
    PartialTjInvdr1 = -MyFaktor*3.0*CosTi * Invdr1
    PartialGij = MyFaktor

    Fij =  ( -PartialRijInvdr1 + (  CosTi*PartialTiInvdr1&
                                  + CosTj*PartialTjInvdr1 )*Invdr1 ) * dr&
         - PartialTiInvdr1 * eMy(1:3,Si,i)&
         - PartialTjInvdr1 * eMy(1:3,Sj,j)
    FDipol(:,Si,i) = FDipol(:,Si,i) + Fij
    FDipol(:,Sj,j) = FDipol(:,Sj,j) - Fij

    eiXej(1) = eMy(2,Si,i)*eMy(3,Sj,j) - eMy(3,Si,i)*eMy(2,Sj,j)
    eiXej(2) = eMy(3,Si,i)*eMy(1,Sj,j) - eMy(1,Si,i)*eMy(3,Sj,j)
    eiXej(3) = eMy(1,Si,i)*eMy(2,Sj,j) - eMy(2,Si,i)*eMy(1,Sj,j)
    eXrij(1) = eMy(2,Si,i)*dr(3) - eMy(3,Si,i)*dr(2)
    eXrij(2) = eMy(3,Si,i)*dr(1) - eMy(1,Si,i)*dr(3)
    eXrij(3) = eMy(1,Si,i)*dr(2) - eMy(2,Si,i)*dr(1)
    Mij = -PartialTiInvdr1*eXrij - PartialGij*eiXej
    M(:,i) = M(:,i) + Mij
    eXrij(1) = eMy(2,Sj,j)*dr(3) - eMy(3,Sj,j)*dr(2)
    eXrij(2) = eMy(3,Sj,j)*dr(1) - eMy(1,Sj,j)*dr(3)
    eXrij(3) = eMy(1,Sj,j)*dr(2) - eMy(2,Sj,j)*dr(1)
    Mji = -PartialTjInvdr1*eXrij + PartialGij*eiXej
    M(:,j) = M(:,j) + Mji

    UpotMy = UpotMy + TermU
    VirialMy = VirialMy + DOT_PRODUCT(drMol,Fij)
  END DO   ! Sj

  DO Sj =1,AnzQuadrupole(Komp(j)),1
    ! Dipole(i) -- Quadrupole(j)
    !*************************************************************************
    dr = drMol + drDipol(:,Si,i) - drQuadrupol(:,Sj,j)
    Invdr1 = 1.0/SQRT(DOT_PRODUCT(dr,dr))
    MyQFaktor =  1.5*absMy(Si,Komp(i))*absQ(Sj,Komp(j))&
               * Invdr1*Invdr1*Invdr1*Invdr1
    CosTi = DOT_PRODUCT(eMy(1:3,Si,i),dr)*Invdr1
    CosTj = DOT_PRODUCT(eQ(1:3,Sj,j),dr)*Invdr1
    CosGij = DOT_PRODUCT(eMy(1:3,Si,i),eQ(1:3,Sj,j))
    TermU = MyQFaktor * ( -CosTi*(5.0*CosTj*CosTj - 1.0) + 2.0*CosGij*CosTj)
    PartialRijInvdr1 = -4.0*Invdr1*TermU * Invdr1
    PartialTiInvdr1 = MyQFaktor * ( -5.0*CosTj*CosTj + 1.0 ) * Invdr1
    PartialTjInvdr1 = MyQFaktor * ( -10.0*CosTi*CosTj + 2.0*CosGij ) * Invdr1
    PartialGij = MyQFaktor*2.0*CosTj

    Fij =  ( -PartialRijInvdr1 + (  CosTi*PartialTiInvdr1&
                                  + CosTj*PartialTjInvdr1 )*Invdr1 ) * dr&
         - PartialTiInvdr1 * eMy(1:3,Si,i)&
         - PartialTjInvdr1 * eQ(1:3,Sj,j)
    FDipol(:,Si,i) = FDipol(:,Si,i) + Fij
    FQuadrupol(:,Sj,j) = FQuadrupol(:,Sj,j) - Fij

    eiXej(1) = eMy(2,Si,i)*eQ(3,Sj,j) - eMy(3,Si,i)*eQ(2,Sj,j)
    eiXej(2) = eMy(3,Si,i)*eQ(1,Sj,j) - eMy(1,Si,i)*eQ(3,Sj,j)
    eiXej(3) = eMy(1,Si,i)*eQ(2,Sj,j) - eMy(2,Si,i)*eQ(1,Sj,j)
    eXrij(1) = eMy(2,Si,i)*dr(3) - eMy(3,Si,i)*dr(2)
    eXrij(2) = eMy(3,Si,i)*dr(1) - eMy(1,Si,i)*dr(3)
    eXrij(3) = eMy(1,Si,i)*dr(2) - eMy(2,Si,i)*dr(1)
    Mij = -PartialTiInvdr1*eXrij - PartialGij*eiXej
    M(:,i) = M(:,i) + Mij
    eXrij(1) = eQ(2,Sj,j)*dr(3) - eQ(3,Sj,j)*dr(2)
    eXrij(2) = eQ(3,Sj,j)*dr(1) - eQ(1,Sj,j)*dr(3)
    eXrij(3) = eQ(1,Sj,j)*dr(2) - eQ(2,Sj,j)*dr(1)
    Mji = -PartialTjInvdr1*eXrij + PartialGij*eiXej
    M(:,j) = M(:,j) + Mji

    UpotMyQ = UpotMyQ + TermU
    virialMyQ = virialMyQ + DOT_PRODUCT(drMol,Fij)
  END DO   ! Sj
END DO   ! Si


DO Si=1,AnzQuadrupole(Komp(i)),1
  DO Sj =1,AnzLadungen(Komp(j)),1
    ! Quadrupole(i) -- Charge(j)
    dr = drMol + drQuadrupol(:,Si,i) - drLadung(:,Sj,j)
    Invdr2 = 1.0/DOT_PRODUCT(dr,dr)
    Invdr1 = SQRT(Invdr2)    
    drE = -dr * Invdr1
    CosTi = DOT_PRODUCT(eQ(1:3,Si,i),drE)
    Cos2Ti = CosTi * CosTi
    TermU = 1.5 * absQ(Si,Komp(i)) * absLadung(Sj,Komp(j)) * Invdr1 * Invdr2
    
    Fij = TermU * Invdr1 * ( 2.0*CosTi*eQ(1:3,Si,i) - (5.0*Cos2Ti-1.0) * drE )
    FQuadrupol(:,Si,i) = FQuadrupol(:,Si,i) + Fij
    FLadung(:,Sj,j) = FLadung(:,Sj,j) - Fij
    eXrij(1) = eQ(2,Si,i)*drE(3) - eQ(3,Si,i)*drE(2)
    eXrij(2) = eQ(3,Si,i)*drE(1) - eQ(1,Si,i)*drE(3)
    eXrij(3) = eQ(1,Si,i)*drE(2) - eQ(2,Si,i)*drE(1)        
    M(:,i) = M(:,i) - TermU * 2.0 * CosTi * eXrij 
    UpotQ =  UpotQ + TermU * ( Cos2Ti - Drittel )
    VirialQ = VirialQ + DOT_PRODUCT(drMol,Fij)
  END DO   ! Sj    

  DO Sj =1,AnzDipole(Komp(j)),1
    ! Quadrupole(i) -- Dipole(j)
    !*************************************************************************
    dr = drMol + drQuadrupol(:,Si,i) - drDipol(:,Sj,j)
    Invdr1 = 1.0/SQRT(DOT_PRODUCT(dr,dr))
    MyQFaktor =  1.5*absQ(Si,Komp(i))*absMy(Sj,Komp(j))&
               * Invdr1*Invdr1*Invdr1*Invdr1
    CosTi = DOT_PRODUCT(eQ(1:3,Si,i),dr)*Invdr1
    CosTj = DOT_PRODUCT(eMy(1:3,Sj,j),dr)*Invdr1
    CosGij = DOT_PRODUCT(eQ(1:3,Si,i),eMy(1:3,Sj,j))
    TermU = MyQFaktor * ( CosTj*(5.0*CosTi*CosTi - 1.0) - 2.0*CosGij*CosTi)
    PartialRijInvdr1 = -4.0*Invdr1*TermU * Invdr1
    PartialTiInvdr1 = MyQFaktor * ( 10.0*CosTi*CosTj - 2.0*CosGij ) * Invdr1
    PartialTjInvdr1 = MyQFaktor * ( 5.0*CosTi*CosTi - 1.0 ) * Invdr1
    PartialGij = -MyQFaktor*2.0*CosTi

    Fij =  ( -PartialRijInvdr1 + (  CosTi*PartialTiInvdr1&
                                  + CosTj*PartialTjInvdr1 )*Invdr1 ) * dr&
         - PartialTiInvdr1 * eQ(1:3,Si,i)&
         - PartialTjInvdr1 * eMy(1:3,Sj,j)
    FQuadrupol(:,Si,i) = FQuadrupol(:,Si,i) + Fij
    FDipol(:,Sj,j) = FDipol(:,Sj,j) - Fij

    eiXej(1) = eQ(2,Si,i)*eMy(3,Sj,j) - eQ(3,Si,i)*eMy(2,Sj,j)
    eiXej(2) = eQ(3,Si,i)*eMy(1,Sj,j) - eQ(1,Si,i)*eMy(3,Sj,j)
    eiXej(3) = eQ(1,Si,i)*eMy(2,Sj,j) - eQ(2,Si,i)*eMy(1,Sj,j)
    eXrij(1) = eQ(2,Si,i)*dr(3) - eQ(3,Si,i)*dr(2)
    eXrij(2) = eQ(3,Si,i)*dr(1) - eQ(1,Si,i)*dr(3)
    eXrij(3) = eQ(1,Si,i)*dr(2) - eQ(2,Si,i)*dr(1)
    Mij = -PartialTiInvdr1*eXrij - PartialGij*eiXej
    M(:,i) = M(:,i) + Mij
    eXrij(1) = eMy(2,Sj,j)*dr(3) - eMy(3,Sj,j)*dr(2)
    eXrij(2) = eMy(3,Sj,j)*dr(1) - eMy(1,Sj,j)*dr(3)
    eXrij(3) = eMy(1,Sj,j)*dr(2) - eMy(2,Sj,j)*dr(1)
    Mji = -PartialTjInvdr1*eXrij + PartialGij*eiXej
    M(:,j) = M(:,j) + Mji

    UpotMyQ = UpotMyQ + TermU
    virialMyQ = virialMyQ + DOT_PRODUCT(drMol,Fij)
  END DO   ! Sj



  DO Sj =1,AnzQuadrupole(Komp(j)),1
    ! Quadrupole(i) -- Quadrupole(j)
    !*************************************************************************
    dr = drMol + drQuadrupol(:,Si,i) - drQuadrupol(:,Sj,j)
    Invdr1 = 1.0/SQRT(DOT_PRODUCT(dr,dr))
    QFaktor = 0.75*absQ(Si,Komp(i))*absQ(Sj,Komp(j))&
             *Invdr1*Invdr1*Invdr1*Invdr1*Invdr1
    CosTi = DOT_PRODUCT(eQ(1:3,Si,i),dr)*Invdr1
    Cos2Ti = CosTi*CosTi
    CosTj = DOT_PRODUCT(eQ(1:3,Sj,j),dr)*Invdr1
    Cos2Tj = CosTj*CosTj
    CosGij = DOT_PRODUCT(eQ(1:3,Si,i),eQ(1:3,Sj,j))
    TermKlammer = (CosGij - 5.0*CosTi*CosTj)
    TermU =  QFaktor * (  1.0 - 5.0*(Cos2Ti+Cos2Tj) - 15.0*Cos2Ti*Cos2Tj&
                        + 2.0*TermKlammer*TermKlammer )
    PartialRijInvdr1 = -5.0*Invdr1*TermU * Invdr1
    PartialTiInvdr1 = -QFaktor * (  10.0*CosTi + 30.0*CosTi*Cos2Tj&
                                  + 20.0*CosTj*TermKlammer ) * Invdr1
    PartialTjInvdr1 = -QFaktor * (  10.0*CosTj + 30.0*Cos2Ti*CosTj&
                                  + 20.0*CosTi*TermKlammer ) * Invdr1
    PartialGij = QFaktor*4.0*TermKlammer

    Fij =  ( -PartialRijInvdr1 + (  CosTi*PartialTiInvdr1&
                                  + CosTj*PartialTjInvdr1 )*Invdr1 ) * dr&
         - PartialTiInvdr1 * eQ(1:3,Si,i)&
         - PartialTjInvdr1 * eQ(1:3,Sj,j)
    FQuadrupol(:,Si,i) = FQuadrupol(:,Si,i) + Fij
    FQuadrupol(:,Sj,j) = FQuadrupol(:,Sj,j) - Fij

    eiXej(1) = eQ(2,Si,i)*eQ(3,Sj,j) - eQ(3,Si,i)*eQ(2,Sj,j)
    eiXej(2) = eQ(3,Si,i)*eQ(1,Sj,j) - eQ(1,Si,i)*eQ(3,Sj,j)
    eiXej(3) = eQ(1,Si,i)*eQ(2,Sj,j) - eQ(2,Si,i)*eQ(1,Sj,j)
    eXrij(1) = eQ(2,Si,i)*dr(3) - eQ(3,Si,i)*dr(2)
    eXrij(2) = eQ(3,Si,i)*dr(1) - eQ(1,Si,i)*dr(3)
    eXrij(3) = eQ(1,Si,i)*dr(2) - eQ(2,Si,i)*dr(1)
    Mij = -PartialTiInvdr1*eXrij - PartialGij*eiXej
    M(:,i) = M(:,i) + Mij
    eXrij(1) = eQ(2,Sj,j)*dr(3) - eQ(3,Sj,j)*dr(2)
    eXrij(2) = eQ(3,Sj,j)*dr(1) - eQ(1,Sj,j)*dr(3)
    eXrij(3) = eQ(1,Sj,j)*dr(2) - eQ(2,Sj,j)*dr(1)
    Mji = -PartialTjInvdr1*eXrij + PartialGij*eiXej
    M(:,j) = M(:,j) + Mji

    UpotQ = UpotQ + TermU
    virialQ = virialQ + DOT_PRODUCT(drMol,Fij)
  END DO   ! Sj
END DO   ! Si


! ********* Reaction field contributions

MRF(1) = epsRFInvrc3*( MyRFSpace(3,Komp(i),i)*MyRFSpace(2,Komp(j),j)-&
                       MyRFSpace(2,Komp(i),i)*MyRFSpace(3,Komp(j),j) )
MRF(2) = epsRFInvrc3*( MyRFSpace(1,Komp(i),i)*MyRFSpace(3,Komp(j),j)-&
                       MyRFSpace(3,Komp(i),i)*MyRFSpace(1,Komp(j),j) )
MRF(3) = epsRFInvrc3*( MyRFSpace(2,Komp(i),i)*MyRFSpace(1,Komp(j),j)-&
                       MyRFSpace(1,Komp(i),i)*MyRFSpace(2,Komp(j),j) )

M(:,i) = M(:,i) + MRF
M(:,j) = M(:,j) - MRF		   
		   
UMyRF = UMyRF + DOT_PRODUCT( MyRFSpace(:,Komp(i),i), MyRFSpace(:,Komp(j),j))

!*****************************************************************************
