thermo = [-794335.40000000002, -3.0000000000000000, -1.3670000000000000, 3275843.0000000000, 4.0155300000000000, 1140.7719999999999, 1.3746600000000000, 2.8351700000000002, 4.9710799999999997, 0.0000000000000000, 224.48914020000001, 40302696.428999998, 7.7650000000000219E-3, 938851.31331355765, 8.2479599999999991, -18.620182098000015, -6.3457399999999993, -9.3100910490000075, 4.1239799999999995, 67346.742060000004]
v0     = -thermo[3]
nr9    = thermo[11]
c1     = thermo[12]
c2     = thermo[13]
c3     = thermo[14]
aii    = thermo[15]
aiikk  = thermo[16]
aiikk2 = thermo[18]
aii2   = thermo[19]
nr9t0  = thermo[20]

# t1     = thermo(6,id)/t
# t2     = t/tr
# nr9t   = nr9*t
# # c                                 initial guess for volume:
# # c                                 taylor(diff(FTH,v),v=v0,1)
# # c                                 JADC Feb 26, 2008.
# # c                                 the dfth0 could be loaded as a
# # c                                 constant.
# tht    = t1
# tht0   = tht*t2
# gamma0 = thermo(7,id)
# k00    = thermo(4,id)
# k0p    = thermo(5,id)

# dfth   = nr9t*gamma0/v0*(3d0*plg(tht)/tht**3 - dlog(1d0-exp(-tht)))
# dfth0  = nr9t0*gamma0/v0*(3d0*plg(tht0)/tht0**3 - dlog(1d0-exp(-tht0)))
# # c                                 taylor(diff(FC,v),v=v0,2)
# # c     v       = (k00-dfth+dfth0-p)/k00*v0
# # c                                 taylor(diff(FC,v),v=v0,3)
# root = k00*((2d0+2d0*k0p)*(p+dfth-dfth0) + k00)

# if (root > 0d0) 
#       v = ((2d0+k0p)-dsqrt(root)/k00)*v0/(1d0+k0p)
#       if (v.lt.v0/1d1.or.v.gt.v0*1d1) v = v0
# else
#       v = v0
# end 

# itic = 0
# ibad = 4
# bad = .true.

# do

#       itic = itic + 1
# # c                                 f, and derivatives
#       v23 = (v0/v)**r23
#       f = 0.5d0*v23 - 0.5d0
#       df = -v23/v/3d0
#       d2f = r59*v23/v**2
# # c                                 cold part derivatives
#       dfc = (c3*f+c1)*f*df
#       d2fc = (2d0*c3*f+c1)*df**2+(c3*f+c1)*f*d2f
# # c                                 debye T/T (tht)
#       z  = 1d0+(aii+aiikk2*f)*f

#       if (z.lt.0d0.or.v/v0.gt.1d2.or.v/v0.lt.1d-2) exit

#       root = dsqrt(z)

#       tht   = t1*root
#       tht0  =  tht*t/tr
# # c                                 tht derivatives
#       a2f   = aii2+aiikk2*f
#       da    = a2f/root
#       dtht  = t1*da*df
#       d2tht = t1*((aiikk2/root-a2f**2/z**1.5d0)*df**2
# *               + da*d2f)

#       dtht0 = dtht*t2
#       d2tht0 = d2tht*t2
# # c                                 polylog functions:
#       fpoly   = 3d0*plg(tht )/tht**3
#       fpoly0  = 3d0*plg(tht0)/tht0**3
# # c                                 thermal part derivatives:
#       etht  = dexp(-tht )

#       if (1d0-etht.lt.0d0) exit

#       letht = dlog(1d0-etht)

#       dfth = (letht-fpoly)*nr9t*dtht/tht
#       d2fth = ((4d0*dtht**2/tht-d2tht)*(fpoly-letht)
# *         + dtht**2*etht/(1d0-etht))*nr9t/tht

#       etht0 = dexp(-tht0)

#       if (1d0-etht0.lt.0d0) exit

#       letht0 = dlog(1d0-etht0)

#       dfth0 = (letht0-fpoly0)*nr9t0*dtht0/tht0
#       d2fth0 = ((4d0*dtht0**2/tht0-d2tht0)*(fpoly0-letht0)
# *          + dtht0**2*etht0/(1d0-etht0))*nr9t0/tht0

#       f1  = -dfc - dfth + dfth0 - p

#       df1 = -d2fc - d2fth + d2fth0

#       dv = f1/df1

#       if (v - dv.lt.0d0) dv = v/2d0

#       v = v - dv

#       if (itic.gt.iopt(21).or.dabs(f1).gt.1d40) then
# c                                 allow bad result
#       if (dabs(f1/p).lt.zero) ibad = 5

#       exit

#       else if (dabs(dv/(1d0+v)).lt.nopt(51)) then

#       bad = .false.
#       exit

#       end if

# end do

# if (bad) then
# # c                                 if we get here, failed to converge
#       if (izap.le.iopt(1)) then

#       msg = 'STXGJI/'//names(id)

#       call conwrn (ibad,msg)

#       izap = izap + 1

#       if (izap.eq.iopt(1)) call warn (49,r,93,'STXGJI')

#       end if

#       if (ibad.eq.4) then 
# # c                                 destabilize the phase.
#       gstxgi  = 1d2*p
#       badend(id) = .true.
#       return

#       end if

# end if
# # c                                 everything is ok, now get
# # c                                 helmoltz energy:
# f = 0.5d0*(v0/v)**r23 - 0.5d0
# z = 1d0+(aii+aiikk2*f)*f
# root = dsqrt(z)
# # c                                 final estimate for tht
# tht   = t1*root
# tht0  = tht*t2
# # c                                 helmholtz enery
# a = thermo(1,id) + c1*f**2*(0.5d0 + c2*f)
# *    + nr9*(t/tht**3*plg(tht ) -tr/tht0**3*plg(tht0))

# gstxgi = a + p*v - t*thermo(10,id)
# # c                                 z = (theta/theta0)^2
# gamma = (2d0*f+1d0)*(aii+aiikk*f)/6d0/z
# etas = - gamma - thermo(17,id)/z*(2d0*f + 1d0)**2
# # c                                 thermal energy/V, based on
# # c                                 previous v estimate
# if (gamma.ne.0d0) then
#       ethv = (dfth0-dfth)/gamma
# else
#       ethv = 0d0
# end if
# # c                                 adiabatic shear modulus
# smu = (1d0 + 2d0*f)**(2.5d0)*(emod(1,id) + f*(thermo(21,id) + thermo(22,id)*f)) - etas*ethv