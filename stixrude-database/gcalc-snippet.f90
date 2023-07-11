!********************************************************************************
!====================== eduester18===============================
!===== volume function Stixrude and Lithgow-Bertelloni 2011, pirated from PERPLEX
!
!=====1 H0R = Helmoltz energy (F0, J/mol)
!=====2 NAT = number of atoms per formula unit (n)
!=====3 V0R = reference volume (V0)
!=====4 VPA = isothermal bulk modulus (K0, kbar)
!=====5 VPB = pressure derivative of the isothermal bulk modulus (K')
!=====6 D1  = Debye Temperature (Θ0, K)
!=====7 D2  = Gruneisen thermal parameter (γ0)
!=====8 D3  = Mie-Gruneisen exponent (q0)
!=====9 D4  = Shear strain derivative of the tensorial Gruneisen parameter (ηS0)
!=====10 VAA = Configurational (and magnetic) entropy (J/mol/K)
!===== SHM0 = Adiabatic shear modulus (μS0, kbar)
!===== SHMP = Pressure derivative of the adiabatic shear modulus (μS0')
!-----
      IF (VOLVO.EQ.12) THEN
      
      IZAP=0
      IZAP1=0    
!---  NAT number of atoms stored because of phases with CODE = +...
      NN0=NAT
!---            
      A0=VAA
      K0=VPA
      K01=VPB
!---

      NR9 = 9D0*NN0*R ! R = gas constant
      C1 = 9D0*K0*V0R
      C2 = K01/2D0-2D0
      C3 = 3D0*C1*C2
      C5 = 3D0*D2
      C7 = C5*(-2D0+6D0*D2-3D0*D3)
      QM1 = D3 - 1D0
      DT = T - T0
      

!---  'Newton–Raphson method' approximation:
!---          initial guess for volume
!             uses a 4 order expansion on
!             the polylog terms before 
!             differentiation and a second
!             order series on the resulting
!             derivative.


      VOLUM = V0R*(3D0*((2D0 - D3)*DT*D2*NR9 - 3D0*P*V0R) + C1)/ &
              (C1 - QM1*DT*D2*NR9)

      IF (VOLUM.LT.0D0) VOLUM = V0R

!---          run 'Newton–Raphson'
      FL = 1D9  ! bar
      ITIC = 0 
      TOL = 1D-6*P ! bar
      
      DO WHILE (DABS(FL).GT.TOL)
      

         ITIC = ITIC + 1
         VQ = (VOLUM/V0R)**D3
         GVQ = D2*VQ
         ! Birch-Murnaghan finite strain
         V23 = (V0R/VOLUM)**(2D0/3D0)
         F = 0.5D0*V23 - 0.5D0
         DF = -V23/VOLUM/3D0
         D2F = 5D0/9D0*V23/VOLUM**2D0
         ! Debye T/T
         ROOT = DSQRT(1D0+C7*F**2+2D0*C5*F)
         IF (ROOT**2D0.LT.0D0) GOTO 91 

         THT = D1/T*ROOT
         THT0 =  THT*T/T0

         ETHT  = DEXP(-THT )
         IF (1D0-ETHT.LE.0) GOTO 92 
         ETHT0 = DEXP(-THT0)
         LTHT  = DLOG(1D0 - ETHT )
         LTHT0 = DLOG(1D0 - ETHT0)
         ! diff(theta/T,v)
         DD2 = (C7*F+C5)/ROOT
         DD1 = D1*DF*DD2
         DTHT = DD1/T
         DTHT0 = DD1/T0
         ! diff(theta/T,v,v)
         DD3 = D1*(DF**2*(C7-C5**2)/ROOT**3 + D2F*DD2)
         D2THT  = DD3/T
         D2THT0 = DD3/T0

         PLGG   = plg(tht )
         PLGG0  = plg(tht0)
         DG  = THT **2*LTHT *DTHT
         DG0 = THT0**2*LTHT0*DTHT0
         D2G  = ((2D0*LTHT  + THT *ETHT /(1D0-ETHT ))*DTHT**2 + &
               THT *LTHT *D2THT )*THT
         D2G0 = ((2D0*LTHT0 + THT0*ETHT0/(1D0-ETHT0))*DTHT0**2 + &
               THT0*LTHT0*D2THT0)*THT0

         DFC = (C3*F+C1)*F*DF
         D2FC = (2D0*C3*F+C1)*DF**2+(C3*F+C1)*F*D2F

         DFT  = NR9*T /THT **3*(DG  -3D0/THT *PLGG *DTHT )
         DFT0 = NR9*T0/THT0**3*(DG0 -3D0/THT0*PLGG0*DTHT0)

         D2FT =  NR9*T /THT **3*(3D0/THT *(DTHT *(4D0/THT *PLGG *DTHT & 
                             - 2D0*DG ) - PLGG *D2THT ) + D2G )
         D2FT0 = NR9*T0/THT0**3*(3D0/THT0*(DTHT0*(4D0/THT0*PLGG0*DTHT0 &
                             - 2D0*DG0) - PLGG0*D2THT0) + D2G0)

         FL  = -DFC - DFT + DFT0 - P   ! P(V) bar

         DFL = -D2FC - D2FT + D2FT0    ! d P(V)/dV
!--- 'Newton–Raphson'-method:           
         VOLUM = VOLUM - FL/DFL
!---          
!                                 the limit on FL should be 
!                                 machine dependent
         IF (VOLUM.LE.0D0.OR.ITIC.GT.40.OR.DABS(FL).GT.1D40) GOTO 93 
      END DO
!                                 if everything is ok, now get 
!                                 helmoltz energy:
      GOTO 10 
!                                 if we get here, failed to converge
91    IF (IZAP.LT.100) THEN
         write (*,*) 'ERROR=91'
         write (*,1000) T,P,NAME(IP), IP
      ELSE IF (IZAP.EQ.100) THEN 
         write (*,1010)
      END IF 

92    IF (IZAP.LT.100) THEN
         write (*,*) 'ERROR=92'
         write (*,1000) T,P,NAME(IP), IP
      ELSE IF (IZAP.EQ.100) THEN 
         write (*,1010)
      END IF 
      
93    IF (IZAP.LT.100) THEN
         write (*,*) 'ERROR=93'
         write (*,1000) T,P,NAME(IP), IP
      ELSE IF (IZAP.EQ.100) THEN 
         write (*,1010)
      END IF             

      IZAP = IZAP + 1
!      use series expansion on failure
!      v = v0*(3d0*((2d0 - q)*dt*gamma0*nr9 - 3d0*p*v0) + c1) &
!            (c1 - qm1*dt*gamma0*nr9)
!      if (v.lt.0d0) v = v0
!      destabilize the phase:  
      GR  = 1D30
      GOTO 99  
!      
10    VQ = (VOLUM/V0R)**D3
      ! Volume (Birch-Murnaghan finite strain)
      F = 0.5D0*(V0R/VOLUM)**(2D0/3D0) - 0.5D0
      ROOT = DSQRT(1D0+C7*F**2+2D0*C5*F)
      THT = D1/T*ROOT
      THT0 =  THT*T/T0

!---  Helmholtz energy
      FF = H0R + C1*F**2*(0.5D0 + C2*F) &
          + NR9*(T/THT**3*plg(THT) -T0/THT0**3*plg(THT0))
!---  Legendre transformation G=F+PV
      GR = FF + P*VOLUM - T*VAA