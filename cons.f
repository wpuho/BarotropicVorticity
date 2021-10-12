      subroutine cons
c
      include '../include/param.h'
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
      include '../include/const.h'
      include '../include/fftcom.h'
c
      rad=6.37e6
      radsq=rad*rad
      grav=9.81
      omega=7.292e-5
      cp=1004.5
      rgas=287.
      radsq = rad**2
      pi=4.0*atan(1.0)
      capa=1.0/3.5
      flat=24.0
c
c  build pointer arrays for locating zonal and total wavenumber
c  values in the one-dimensional spherical harmonic arrays.
c
      call sortml (jtrun,mlmax,msort,lsort,mlsort)
c
      call make_list
c
      do 150 m =1,mlistnum
       mf=mlist(m)
cibm--
       rm=mf-1
       if (mf.eq.1)  rm= 0.0
       cim(m) = rm
cibm--
      do 150 l=mf,jtrun
       rl = l
       rlm= rl-1.0
       if (mf.eq.1)  rm= 0.0
       eps4(l,m)= rl*rlm/radsq
       if ( l.eq.1 ) then
        wdfac(1,m) = 0. ; wcfac(1,m) = 0.
       else
        wdfac(l,m) = 1.0/(radsq*eps4(l,m))
        wcfac(l,m) = cim(m)*wdfac(l,m)
       endif
  150 continue
c
c  gaussian quadrature weights and latitudes
c
      one= 1.0
      onem= -one
      call gausl3(my,onem,one,weight,sinl)
c
      my2= my/2
c
      do 180 j = 1, my2
      sinl(my+1-j)  = -sinl(j)
      weight(my+1-j)= weight(j)
      onocos(j)     = 1.0/(1.0-sinl(j)*sinl(j))
      onocos(my+1-j)= onocos(j)
      cosl(j)       = 1.0/sqrt(onocos(j))
      cosl(my+1-j)  = cosl(j)
  180 continue
c
c  define coriolis parameter for each latitude
c
c      do 190 j=1,my
c      cor(j)= 2.0*omega*sinl(j)
c 190 continue
 
      if(fplane)then
        d2r=pi/180.
        do j=1,my
          cor(j)=2.0*omega*sin(flat*d2r)
        enddo
      else
      do 190 j=1,my
      cor(j)= 2.0*omega*sinl(j)
  190 continue     
      endif
      
      basu = 0.
      basv = 0.
c
c
c  initialize ifax and trigs for rfftmlt routine
c
      call fftfax (nx,ifax,trigs)
c
cibm---beg
cibm ... if 'ibm_fft=1' set in fftcom.h
cibm ... check if the ESSL fft routine support the seq lenth of nx
cibm ... if not, rfftmlt routine will be used.
cibm ... contral logical 'lessl_fft' is in the common/ibm_fft/
c
         lessl_fft=.false.
         if (ibm_fft.eq.1) then
           nn = nx
           if (iand(nn,1).eq.1) go to 111   ! not an even number
           if (mod(nn,9).eq.0) then
             if (mod(nn/9,3).eq.0) go to 111   ! radix of 3**i, i>2
           endif
           if (mod(nn,5).eq.0) then
             if (mod(nn/5,5).eq.0) go to 111   ! radix of 5**i, i>1
           endif
           if (mod(nn,7).eq.0) then
             if (mod(nn/7,7).eq.0) go to 111   ! radix of 7**i, i>1
           endif
           if (mod(nn,11).eq.0) then
             if (mod(nn/11,11).eq.0) go to 111 ! radix of 11**i, i>1
           endif
           if (mod(nn,13).eq.0) go to 111   ! radix of 13
           if (mod(nn,17).eq.0) go to 111
           if (mod(nn,19).eq.0) go to 111
           if (mod(nn,23).eq.0) go to 111
           if (mod(nn,29).eq.0) go to 111
           lessl_fft=.true.
  111      continue
         endif
c
c  define associated legendre polynomials and their derivatives
c
      call lgndr (my2,jtrun,jtmax,sinl,poly,dpoly)
c
      return
      end
