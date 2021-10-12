      subroutine gausl3 (n,xa,xb,wt,ab)
c
c weights and abscissas for nth order gaussian quadrature on (xa,xb).
c input arguments
c
c n  -the order desired
c xa -the left endpoint of the interval of integration
c xb -the right endpoint of the interval of integration
c output arguments
c ab -the n calculated abscissas
c wt -the n calculated weights
c
      implicit double precision (a-h,o-z)
c
cfj
      include '../include/rank.h'
cfj
      real  ab(n) ,wt(n),xa,xb
c
c machine dependent constants---
c  tol - convergence criterion for double precision iteration
c  pi  - given to 15 significant digits
c  c1  -  1/8                     these are coefficients in mcmahon"s
c  c2  -  -31/(2*3*8**2)          expansions of the kth zero of the
c  c3  -  3779/(2*3*5*8**3)       bessel function j0(x) (cf. abramowitz,
c  c4  -  -6277237/(3*5*7*8**5)   handbook of mathematical functions).
c  u   -  (1-(2/pi)**2)/4
c
      data tol/1.d-14/,pi/3.14159265358979/,u/.148678816357662/
      data c1,c2,c3,c4/.125,-.080729166666667,.246028645833333,
     1                -1.82443876720609 /
c
c maximum number of iterations before giving up on convergence
c
      data maxit /5/
c
c arithmetic statement function for converting integer to double
c
      dbli(i) = dble(float(i))
c
      ddif = .5d0*(dble(xb)-dble(xa))
      dsum = .5d0*(dble(xb)+dble(xa))
      if (n .gt. 1) go to 101
      ab(1) = 0.
      wt(1) = 2.*ddif
      go to 107
  101 continue
      nnp1 = n*(n+1)
      cond = 1./sqrt((.5+float(n))**2+u)
      lim = n/2
c
      do 105 k=1,lim
	 b = (float(k)-.25)*pi
	 bisq = 1./(b*b)
c
c rootbf approximates the kth zero of the bessel function j0(x)
c
	 rootbf = b*(1.+bisq*(c1+bisq*(c2+bisq*(c3+bisq*c4))))
c
c      initial guess for kth root of legendre poly p-sub-n(x)
c
	 dzero = cos(rootbf*cond)
	 do 103 i=1,maxit
c
	    dpm2 = 1.d0
	    dpm1 = dzero
c
c       recursion relation for legendre polynomials
c
	    do 102 nn=2,n
		dp = (dbli(2*nn-1)*dzero*dpm1-dbli(nn-1)*dpm2)/dbli(nn)
		dpm2 = dpm1
		dpm1 = dp
  102       continue
	    dtmp = 1.d0/(1.d0-dzero*dzero)
	    dppr = dbli(n)*(dpm2-dzero*dp)*dtmp
	    dp2pri = (2.d0*dzero*dppr-dbli(nnp1)*dp)*dtmp
	    drat = dp/dppr
c
c       cubically-convergent iterative improvement of root
c
	    dzeri = dzero-drat*(1.d0+drat*dp2pri/(2.d0*dppr))
	    ddum= dabs(dzeri-dzero)
	 if (ddum .le. tol) go to 104
	    dzero = dzeri
  103    continue
cfj
	 if(myrank .eq. 0) print 504
cfj
  504    format(1x,' in gausl3, convergence failed')
  104    continue
	 ddifx = ddif*dzero
	 ab(k) = dsum-ddifx
	 wt(k) = 2.d0*(1.d0-dzero*dzero)/(dbli(n)*dpm2)**2*ddif
	 i = n-k+1
	 ab(i) = dsum+ddifx
	 wt(i) = wt(k)
  105 continue
c
      if (mod(n,2) .eq. 0) go to 107
      ab(lim+1) = dsum
      nm1 = n-1
      dprod = n
      do 106 k=1,nm1,2
	 dprod = dbli(nm1-k)*dprod/dbli(n-k)
  106 continue
      wt(lim+1) = 2.d0/dprod**2*ddif
  107 return
      end
