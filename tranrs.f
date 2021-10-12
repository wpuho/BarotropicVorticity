      subroutine tranrs (jtrun,jtmax,nx,my,my_max,lev,poly,w,cc
     &                  ,wss,num,nsize)
c
c  subroutine to transform a scalar grid point field to spectral
c  coefficients
c
c *** input ***
c
c  jtrun: zonal wavenumber truncation limit
c  jtmax: maximum amount of zonal waves located in each pe
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  my_max: maximum amount of n-s grids located in each pe
c  lev: number of vertical levels to transform
c  poly: legendre polynomials
c  w: gaussian quadrature weights
c  cc: 3-dim input grid pt. field to be transformed
c  num: number of variables grouped together
c
c *** output ***
c
c  wss: spectral coefficient fields
c
c  **********************************
c
      include '../include/index.h'
      include '../include/paramt.h'
      include '../include/fftcom.h'
c
      dimension poly(jtrun,my/2,jtmax),w(my)
      dimension cc(nx+3,lev,num,my_max)
      dimension wss(lev,2,num,jtrun,jtmax)
c
csun  include '../include/paramt.h' .. change im,jm to nx,my
      dimension  gwk1(nx+3,lev,num,my_max)
c
      dimension wcc_fk (lev,2,num,jtmax,my_max*nsize)
      dimension twcc_fk(lev,2,num,jtmax*nsize,my_max)
c
      dimension wccSUM(lev*2*num,my/2)
      dimension wccDIF(lev*2*num,my/2)
cfong--
      dimension polyw(jm/2,jtr,jtr)
      logical wfirst
      data wfirst/.true./
      save wfirst, polyw
cfong--
c
      mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2
      lev2=lev*2
c
cfong--
      if( wfirst )then
      do j = 1, myhalf
      do m = 1, mlistnum
       mf=mlist(m)
      do l = mf, jtrun
       polyw(j,l,m) = poly(l,j,m)*w(j)
      end do
      end do
      end do
      wfirst = .false.
      endif
cfong--
c
c
c  fft for each guassian latitude of 2-d field
c
cfong
      if( length_fft .eq. 0 )then
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+3,nx,lev*jlistnum*num,-1)
      else
      do jj = 1, jlistnum
       call rfftmlt(cc(1,1,1,jj),gwk1,trigs,ifax,1,nx+3,nx,lev*num,-1)
      end do
      end if
cfong
c
cibm--beg

      mchk=iand(jtrun,3)

      do j =1, jlistnum
        do m=1,mchk
          mm= 2*m-1
          mp= mm+1
          mlst=nlist(m)
          do ii=1,num
            do k=1,lev
              twcc_fk(k,1,ii,mlst,j)=cc(mm,k,ii,j)
              twcc_fk(k,2,ii,mlst,j)=cc(mp,k,ii,j)
            enddo
          enddo
        enddo

        do m=mchk+1,jtrun,4
          mm= 2*m-1
          mp= mm+1
          mlst=nlist(m)
          mm1= 2*(m+1)-1
          mp1= mm1+1
          mlst1=nlist(m+1)
          mm2= 2*(m+2)-1
          mp2= mm2+1
          mlst2=nlist(m+2)
          mm3= 2*(m+3)-1
          mp3= mm3+1
          mlst3=nlist(m+3)
          do ii=1,num
            do k=1,lev
              twcc_fk(k,1,ii,mlst,j)=cc(mm,k,ii,j)
              twcc_fk(k,2,ii,mlst,j)=cc(mp,k,ii,j)
              twcc_fk(k,1,ii,mlst1,j)=cc(mm1,k,ii,j)
              twcc_fk(k,2,ii,mlst1,j)=cc(mp1,k,ii,j)
              twcc_fk(k,1,ii,mlst2,j)=cc(mm2,k,ii,j)
              twcc_fk(k,2,ii,mlst2,j)=cc(mp2,k,ii,j)
              twcc_fk(k,1,ii,mlst3,j)=cc(mm3,k,ii,j)
              twcc_fk(k,2,ii,mlst3,j)=cc(mp3,k,ii,j)
            enddo
          enddo
        enddo
      enddo
cibm--end

      call mpe_transpose_rs(twcc_fk,wcc_fk,lev*2*num,jtmax,my_max,nsize)

cibm--beg

      do m=1,mlistnum
         mf=mlist(m)

      do j=1,myhalf
      j1=jlist2(j)
      j2=jlist2(my-j+1)
      do k=1,lev2*num
       wccSUM(k,j)=(wcc_fk(k,1,1,m,j1)+wcc_fk(k,1,1,m,j2))
       wccDIF(k,j)=(wcc_fk(k,1,1,m,j1)-wcc_fk(k,1,1,m,j2))
      enddo
      enddo

      do l=mf,jtrun
      do k=1,lev2*num
        wss(k,1,1,l,m) = 0.
      enddo
      enddo

cfong--beg
c     i1=(jtrun-mf+1)/4
c     i2=(jtrun-mf+1-i1*4)/2
c     i3= jtrun-mf+1-i1*4-i2*2
c
c     do i=1,i1
c     l=mf+(i-1)*4
c     do j=1,myhalf
c     do k=1,lev2*num
c       wss(k,1,1,l  ,m)=wss(k,1,1,l  ,m)+poly(l  ,j,m)*w(j)*wccSUM(k,j)
c       wss(k,1,1,l+1,m)=wss(k,1,1,l+1,m)+poly(l+1,j,m)*w(j)*wccDIF(k,j)
c       wss(k,1,1,l+2,m)=wss(k,1,1,l+2,m)+poly(l+2,j,m)*w(j)*wccSUM(k,j)
c       wss(k,1,1,l+3,m)=wss(k,1,1,l+3,m)+poly(l+3,j,m)*w(j)*wccDIF(k,j)
c     enddo
c     enddo
c     enddo
c
c     do i=1,i2
c     l=mf+i1*4+(i-1)*2
c     do j=1,myhalf
c     do k=1,lev2*num
c       wss(k,1,1,l  ,m)=wss(k,1,1,l  ,m)+poly(l  ,j,m)*w(j)*wccSUM(k,j)
c       wss(k,1,1,l+1,m)=wss(k,1,1,l+1,m)+poly(l+1,j,m)*w(j)*wccDIF(k,j)
c     enddo
c     enddo
c     enddo
c
c     do i=1,i3
c     l=jtrun
c     do j=1,myhalf
c     do k=1,lev2*num
c       wss(k,1,1,l  ,m)=wss(k,1,1,l,m)+poly(l,j,m)*w(j)*wccSUM(k,j)
c     enddo
c     enddo
c     enddo
c
cfong-- cache blocking method
c   nb: multiple of 4
c
        nb = 32
c
        lchk = iand(jtrun-mf+1, 1)
        lle  = jtrun-lchk
c
        do ll = mf, lle, nb
        do kk = 1, lev2*num, nb
        do jj = 1, myhalf, nb

          do l = ll, min(ll+nb-1, lle), 2
          do k = kk, min(kk+nb-1, lev2*num), 2
            sa00 = wss(k  ,1,1,l  ,m)
            sa01 = wss(k  ,1,1,l+1,m)
            sa10 = wss(k+1,1,1,l  ,m)
            sa11 = wss(k+1,1,1,l+1,m)
          do j = jj, min(jj+nb-1, myhalf)
            sa00 = sa00 + polyw(j,l  ,m)*wccSUM(k,j)
            sa01 = sa01 + polyw(j,l+1,m)*wccDIF(k,j)
            sa10 = sa10 + polyw(j,l  ,m)*wccSUM(k+1,j)
            sa11 = sa11 + polyw(j,l+1,m)*wccDIF(k+1,j)
          enddo
            wss(k  ,1,1,l  ,m) = sa00
            wss(k  ,1,1,l+1,m) = sa01
            wss(k+1,1,1,l  ,m) = sa10
            wss(k+1,1,1,l+1,m) = sa11
          enddo
          enddo
        enddo
        enddo
        enddo
c
c odd number
c
        if( lchk.eq.1 )then
            l=jtrun
        do kk = 1, lev2*num, nb
        do jj = 1, myhalf, nb
          do k = kk, min(kk+nb-1, lev2*num), 2
            sa00 = wss(k  ,1,1,l,m)
            sa10 = wss(k+1,1,1,l,m)
          do j = jj, min(jj+nb-1, myhalf)
            sa00 = sa00 + polyw(j,l  ,m)*wccSUM(k,j)
            sa10 = sa10 + polyw(j,l  ,m)*wccSUM(k+1,j)
          enddo
            wss(k  ,1,1,l,m) = sa00
            wss(k+1,1,1,l,m) = sa10 
          enddo
        enddo
        enddo
        endif

cfong--end

      enddo
cibm---end

      return
      end
