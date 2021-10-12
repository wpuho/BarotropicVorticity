      subroutine tranrs1(jtrun,jtmax,nx,my,my_max,poly,w,r,s,nsize)
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
c  poly: legendre polynomials
c  w: gaussian quadrature weights
c  r: 2-dim input grid pt. field to be transformed
c
c *** output ***
c
c  s: spectral coefficient fields
c
c  **********************************
c
      include '../include/index.h'
      include '../include/paramt.h'
      include '../include/fftcom.h'
c
      dimension poly(jtrun,my/2,jtmax),w(my)
      dimension r(nx,my),s(jtrun,jtmax,2)
c
csun  include '../include/paramt.h' .. change im,jm to nx,my
      dimension gwk1(nx+3,1,6,my_max)
c
      dimension wcc_fk(jtmax,my_max*nsize,2)
      dimension twcc_fk(jtmax*nsize,my_max,2)
      dimension wss(jtrun,2)
      dimension cc(nx+3,my_max)
c
      dimension wccSUM(my,2)
      dimension wccDIF(my,2)
cfong--
      dimension polyw(jm/2,jtr,jtr)
      logical wfirst
      data wfirst/.true./
      save wfirst, polyw
cfong--
c
      mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2
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
c  put grid point fields into two dimensional horizontal array
c
c*** r1 start ***
      do 23 jj =1, jlistnum
      j=jlist1(jj)
!ocl repeat(im)
      do 23 i=1,nx
      cc(i,jj)= r(i,j)
   23 continue
c
c  fft for each guassian latitude of 2-d field
c
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+3,nx,jlistnum,-1)

c
      do j =1, jlistnum
!ocl repeat(jtr)
      do m=1,jtrun
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         twcc_fk(mlst,j,1)=cc(mm,j)
         twcc_fk(mlst,j,2)=cc(mp,j)
      enddo
      enddo

       call mpe_transpose_rs1(twcc_fk,wcc_fk,jtmax,my_max,2,nsize)

c*** r1  end  ***
cibm---beg

      do m=1,mlistnum
         mf=mlist(m)

      do l=mf,jtrun
        wss(l,1) = 0.
        wss(l,2) = 0.
      enddo

      i1=(jtrun-mf+1)/4
      i2=(jtrun-mf+1-i1*4)/2
      i3= jtrun-mf+1-i1*4-i2*2

      do j=1,myhalf
      j1=jlist2(j)
      j2=jlist2(my-j+1)
       wccSUM(j,1)=wcc_fk(m,j1,1)+wcc_fk(m,j2,1)
       wccDIF(j,1)=wcc_fk(m,j1,1)-wcc_fk(m,j2,1)
       wccSUM(j,2)=wcc_fk(m,j1,2)+wcc_fk(m,j2,2)
       wccDIF(j,2)=wcc_fk(m,j1,2)-wcc_fk(m,j2,2)
      enddo

      do i=1,i1
      l=mf+(i-1)*4
      do j=1,myhalf
      jj=my-j+1
        wss(l  ,1)=wss(l  ,1)+polyw(j,l,m)*wccSUM(j,1)
        wss(l+1,1)=wss(l+1,1)+polyw(j,l+1,m)*wccDIF(j,1)
        wss(l+2,1)=wss(l+2,1)+polyw(j,l+2,m)*wccSUM(j,1)
        wss(l+3,1)=wss(l+3,1)+polyw(j,l+3,m)*wccDIF(j,1)
        wss(l  ,2)=wss(l  ,2)+polyw(j,l,m)*wccSUM(j,2)
        wss(l+1,2)=wss(l+1,2)+polyw(j,l+1,m)*wccDIF(j,2)
        wss(l+2,2)=wss(l+2,2)+polyw(j,l+2,m)*wccSUM(j,2)
        wss(l+3,2)=wss(l+3,2)+polyw(j,l+3,m)*wccDIF(j,2)
      enddo
      enddo

      do i=1,i2
      l=mf+i1*4+(i-1)*2
      do j=1,myhalf
      jj=my-j+1
        wss(l  ,1)=wss(l  ,1)+polyw(j,l,m)*wccSUM(j,1)
        wss(l+1,1)=wss(l+1,1)+polyw(j,l+1,m)*wccDIF(j,1)
        wss(l  ,2)=wss(l  ,2)+polyw(j,l,m)*wccSUM(j,2)
        wss(l+1,2)=wss(l+1,2)+polyw(j,l+1,m)*wccDIF(j,2)
      enddo
      enddo

      do i=1,i3
      l=jtrun
      do j=1,myhalf
      jj=my-j+1
        wss(l  ,1)=wss(l  ,1)+polyw(j,l,m)*wccSUM(j,1)
        wss(l  ,2)=wss(l  ,2)+polyw(j,l,m)*wccSUM(j,2)
      enddo
      enddo

      do l=mf,jtrun
        s(l,m,1)=wss(l,1)
        s(l,m,2)=wss(l,2)
      enddo

      enddo
cibm---end

      return
      end
