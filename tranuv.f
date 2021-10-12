      subroutine tranuv (jtrun,jtmax,nx,my,my_max,lev,onocos,wcfac
     *                  ,wdfac,poly,dpoly,vor,div,ut,vt,nsize)
c
c  subroutine to transform vorticity and divergence to velocity
c  components
c
c *** input ***
c
c  jtrun: zonal wavenumber truncation limit
c  jtmax: maximum amount of zonal waves located in each pe
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  my_max: maximum amount of n-s grids located in each pe
c  lev: number of veritical levels to transform
c  onocos: 1.0/(cos(lat)**2)
c  wcfac: constants defined in cons
c  wdfac: constants defined in cons
c  poly: legendre polynomials
c  dpoly: d(poly)/d(sin(lat))
c  vor: spectral vorticity
c  div: spectral divergence
c
c *** output ***
c
c  ut: e-w velocity component
c  vt: n-s velocity component
c
c  ****************************************
c
      include '../include/index.h'
      include '../include/paramt.h'
      include '../include/fftcom.h'
c
      dimension onocos(my),wcfac(jtrun,jtmax),wdfac(jtrun,jtmax)
     &         ,poly(jtrun,my/2,jtmax),dpoly(jtrun,my/2,jtmax)
     &         ,vor(lev,2,jtrun,jtmax),div(lev,2,jtrun,jtmax)
     &         ,ut(nx,lev,my_max),vt(nx,lev,my_max)
c
csun  include '../include/paramt.h' .. change im,jm,mlm to nx,my,mlmax
      dimension  gwk1(nx+3,lev,6,my_max)
      dimension wcc_fk (lev,2,2,jtmax,my_max*nsize)
      dimension twcc_fk(lev,2,2,jtmax*nsize,my_max)
      dimension cc(nx+3,lev,2,my_max)
      dimension tcc(lev,2,2,my)
      dimension ws3(lev,2,2,jtrun)
      dimension ws4(lev,2,2,jtrun)
cibm---add
      dimension tc2(lev,2,2,my)
      dimension wc(jtrun,my/2),wd(jtrun,my/2)
      dimension coslr(jm)
      save coslr
c
      logical lfirst
      data lfirst/.true./
      save lfirst
c---
c
      myhalf=my/2
      lev2=lev*2
      mlx= (jtrun/2)*((jtrun+1)/2)
c
cfong--
      if (lfirst) then
        do j=1,my
          coslr(j)=1./onocos(j)
        enddo
        lfirst=.false.
      endif
cfong--

      do m=1,mlistnum
         mf=mlist(m)
cibm---add
        do j=1,myhalf
        do l=mf,jtrun
          wc(l,j)=wcfac(l,m)* poly(l,j,m)
          wd(l,j)=wdfac(l,m)*dpoly(l,j,m)*coslr(j)
        enddo
        enddo
cibm---

        do l=mf,jtrun
        do k=1,lev
          ws3(k,1,1,l) = +div(k,2,l,m)
          ws3(k,2,1,l) = -div(k,1,l,m)
          ws3(k,1,2,l) = +vor(k,2,l,m)
          ws3(k,2,2,l) = -vor(k,1,l,m)
          ws4(k,1,1,l) = +vor(k,1,l,m)
          ws4(k,2,1,l) = +vor(k,2,l,m)
          ws4(k,1,2,l) = -div(k,1,l,m)
          ws4(k,2,2,l) = -div(k,2,l,m)
        enddo
        enddo

        do k=1,lev*2*2*myhalf
         tcc(k,1,1,1)=0.
        enddo


cfong-- cache blocking method
c     do j=1,myhalf
c     do l=mf,jtrun
c     do k=1,lev*2*2
c     tcc(k,1,1,j  ) = tcc(k,1,1,j  )
c    &    + ws3(k,1,1,l)*wc(l,j) + ws4(k,1,1,l)*wd(l,j)
c     enddo
c     enddo
c     enddo

c nb : multiple of 4
c
        nb=32
        do jj=1,myhalf,nb
          jje= min(jj+nb-1, myhalf)
          ichk = iand(jje-jj+1, 1)
          jjn= jje-ichk

        do kk=1,lev2*2,nb
        do ll=mf,jtrun,nb

        do j = jj,jjn,2
          do k = kk,min(kk+nb-1,lev2*2),4
            sa00 = tcc(k  ,1,1,j  )
            sa10 = tcc(k+1,1,1,j  )
            sa20 = tcc(k+2,1,1,j  )
            sa30 = tcc(k+3,1,1,j  )
            sa01 = tcc(k  ,1,1,j+1)
            sa11 = tcc(k+1,1,1,J+1)
            sa21 = tcc(k+2,1,1,j+1)
            sa31 = tcc(k+3,1,1,j+1)
          do l = ll,min(ll+nb-1,jtrun)
            sa00 =sa00+ws3(k,1,1,l)*wc(l,j)+ws4(k,1,1,l)*wd(l,j)
            sa10 =sa10+ws3(k+1,1,1,l)*wc(l,j)+ws4(k+1,1,1,l)*wd(l,j)
            sa20 =sa20+ws3(k+2,1,1,l)*wc(l,j)+ws4(k+2,1,1,l)*wd(l,j)
            sa30 =sa30+ws3(k+3,1,1,l)*wc(l,j)+ws4(k+3,1,1,l)*wd(l,j)
            sa01 =sa01+ws3(k,1,1,l)*wc(l,j+1)+ws4(k,1,1,l)*wd(l,j+1)
            sa11 =sa11+ws3(k+1,1,1,l)*wc(l,j+1)+ws4(k+1,1,1,l)*wd(l,j+1)
            sa21 =sa21+ws3(k+2,1,1,l)*wc(l,j+1)+ws4(k+2,1,1,l)*wd(l,j+1)
            sa31 =sa31+ws3(k+3,1,1,l)*wc(l,j+1)+ws4(k+3,1,1,l)*wd(l,j+1)
          enddo
            tcc(k  ,1,1,j  ) = sa00
            tcc(k+1,1,1,j  ) = sa10
            tcc(k+2,1,1,j  ) = sa20
            tcc(k+3,1,1,j  ) = sa30
            tcc(k  ,1,1,j+1) = sa01
            tcc(k+1,1,1,J+1) = sa11
            tcc(k+2,1,1,j+1) = sa21
            tcc(k+3,1,1,j+1) = sa31
          enddo
        enddo
c
c odd number
c
        if( ichk.eq.1 )then
          j = jje
          do k = kk,min(kk+nb-1,lev2*2),4
            sa00 = tcc(k  ,1,1,j  )
            sa10 = tcc(k+1,1,1,j  )
            sa20 = tcc(k+2,1,1,j  )
            sa30 = tcc(k+3,1,1,j  )
          do l = ll,min(ll+nb-1,jtrun)
            sa00 = sa00+ws3(k,1,1,l)*wc(l,j)+ws4(k,1,1,l)*wd(l,j)
            sa10 = sa10+ws3(k+1,1,1,l)*wc(l,j)+ws4(k+1,1,1,l)*wd(l,j)
            sa20 = sa20+ws3(k+2,1,1,l)*wc(l,j)+ws4(k+2,1,1,l)*wd(l,j)
            sa30 = sa30+ws3(k+3,1,1,l)*wc(l,j)+ws4(k+3,1,1,l)*wd(l,j)
          enddo
            tcc(k  ,1,1,j  ) = sa00
            tcc(k+1,1,1,j  ) = sa10
            tcc(k+2,1,1,j  ) = sa20
            tcc(k+3,1,1,j  ) = sa30
          enddo
        endif

        enddo
        enddo
        enddo
cfong--end
      
        do l=mf,jtrun,2
        do k=1,lev
          ws3(k,1,1,l) = +div(k,2,l,m)
          ws3(k,2,1,l) = -div(k,1,l,m)
          ws3(k,1,2,l) = +vor(k,2,l,m)
          ws3(k,2,2,l) = -vor(k,1,l,m)
          ws4(k,1,1,l) = -vor(k,1,l,m)
          ws4(k,2,1,l) = -vor(k,2,l,m)
          ws4(k,1,2,l) = +div(k,1,l,m)
          ws4(k,2,2,l) = +div(k,2,l,m)
        enddo
        enddo
        do l=mf+1,jtrun,2
        do k=1,lev
          ws3(k,1,1,l) = -div(k,2,l,m)
          ws3(k,2,1,l) = +div(k,1,l,m)
          ws3(k,1,2,l) = -vor(k,2,l,m)
          ws3(k,2,2,l) = +vor(k,1,l,m)
          ws4(k,1,1,l) = +vor(k,1,l,m)
          ws4(k,2,1,l) = +vor(k,2,l,m)
          ws4(k,1,2,l) = -div(k,1,l,m)
          ws4(k,2,2,l) = -div(k,2,l,m)
        enddo
        enddo
c
        do k=1,lev*2*2*myhalf
         tc2(k,1,1,1)=0.
        enddo

cfong-- cache blocking method
c     do j=1,myhalf
c     do l=mf,jtrun
c     do k=1,lev*2*2
c     tc2(k,1,1,j ) = tc2(k,1,1,j )
c    &    + ws3(k,1,1,l)*wc(l,j) + ws4(k,1,1,l)*wd(l,j)
c     enddo
c     enddo
c     enddo
c nb : multiple of 4
c
        nb=32
        do jj=1,myhalf,nb
          jje= min(jj+nb-1, myhalf)
          ichk = iand(jje-jj+1, 1)
          jjn= jje-ichk

        do kk=1,lev2*2,nb
        do ll=mf,jtrun,nb

        do j = jj,jjn,2
          do k = kk,min(kk+nb-1,lev2*2),4
            sa00 = tc2(k  ,1,1,j  )
            sa10 = tc2(k+1,1,1,j  )
            sa20 = tc2(k+2,1,1,j  )
            sa30 = tc2(k+3,1,1,j  )
            sa01 = tc2(k  ,1,1,j+1)
            sa11 = tc2(k+1,1,1,J+1)
            sa21 = tc2(k+2,1,1,j+1)
            sa31 = tc2(k+3,1,1,j+1)
          do l = ll,min(ll+nb-1,jtrun)
            sa00 =sa00+ws3(k,1,1,l)*wc(l,j)+ws4(k,1,1,l)*wd(l,j)
            sa10 =sa10+ws3(k+1,1,1,l)*wc(l,j)+ws4(k+1,1,1,l)*wd(l,j)
            sa20 =sa20+ws3(k+2,1,1,l)*wc(l,j)+ws4(k+2,1,1,l)*wd(l,j)
            sa30 =sa30+ws3(k+3,1,1,l)*wc(l,j)+ws4(k+3,1,1,l)*wd(l,j)
            sa01 =sa01+ws3(k,1,1,l)*wc(l,j+1)+ws4(k,1,1,l)*wd(l,j+1)
            sa11 =sa11+ws3(k+1,1,1,l)*wc(l,j+1)+ws4(k+1,1,1,l)*wd(l,j+1)
            sa21 =sa21+ws3(k+2,1,1,l)*wc(l,j+1)+ws4(k+2,1,1,l)*wd(l,j+1)
            sa31 =sa31+ws3(k+3,1,1,l)*wc(l,j+1)+ws4(k+3,1,1,l)*wd(l,j+1)
          enddo
            tc2(k  ,1,1,j  ) = sa00
            tc2(k+1,1,1,j  ) = sa10
            tc2(k+2,1,1,j  ) = sa20
            tc2(k+3,1,1,j  ) = sa30
            tc2(k  ,1,1,j+1) = sa01
            tc2(k+1,1,1,J+1) = sa11
            tc2(k+2,1,1,j+1) = sa21
            tc2(k+3,1,1,j+1) = sa31
          enddo
        enddo
c
c odd number
c
        if( ichk.eq.1 )then
          j = jje
          do k = kk,min(kk+nb-1,lev2*2),4
            sa00 = tc2(k  ,1,1,j  )
            sa10 = tc2(k+1,1,1,j  )
            sa20 = tc2(k+2,1,1,j  )
            sa30 = tc2(k+3,1,1,j  )
          do l = ll,min(ll+nb-1,jtrun)
            sa00 = sa00+ws3(k,1,1,l)*wc(l,j)+ws4(k,1,1,l)*wd(l,j)
            sa10 = sa10+ws3(k+1,1,1,l)*wc(l,j)+ws4(k+1,1,1,l)*wd(l,j)
            sa20 = sa20+ws3(k+2,1,1,l)*wc(l,j)+ws4(k+2,1,1,l)*wd(l,j)
            sa30 = sa30+ws3(k+3,1,1,l)*wc(l,j)+ws4(k+3,1,1,l)*wd(l,j)
          enddo
            tc2(k  ,1,1,j  ) = sa00
            tc2(k+1,1,1,j  ) = sa10
            tc2(k+2,1,1,j  ) = sa20
            tc2(k+3,1,1,j  ) = sa30
          enddo
        endif

        enddo
        enddo
        enddo
cfong--end

        do j=1,myhalf
          jj=jlist2(j)
          jx=my-j+1
          j2=jlist2(jx)
          do k=1,lev*2*2
            wcc_fk(k,1,1,m,jj)=tcc(k,1,1,j)
            wcc_fk(k,1,1,m,j2)=tc2(k,1,1,j)
          enddo
        enddo
c
      enddo   ! end of big m loop


      call mpe_transpose_sr(wcc_fk,twcc_fk,lev*2*2,jtmax,my_max,nsize)

      do jj=1,jlistnum

      do i=1,(nx+3)*lev*2
       cc(i,1,1,jj)= 0.
      enddo
c
cfong--
      mchk=iand(jtrun,3)

      do m=1,mchk
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
      do k=1,lev
         cc(mm,k,1,jj)=twcc_fk(k,1,1,mlst,jj)
         cc(mp,k,1,jj)=twcc_fk(k,2,1,mlst,jj)
         cc(mm,k,2,jj)=twcc_fk(k,1,2,mlst,jj)
         cc(mp,k,2,jj)=twcc_fk(k,2,2,mlst,jj)
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
      do k=1,lev
         cc(mm,k,1,jj)=twcc_fk(k,1,1,mlst,jj)
         cc(mp,k,1,jj)=twcc_fk(k,2,1,mlst,jj)
         cc(mm,k,2,jj)=twcc_fk(k,1,2,mlst,jj)
         cc(mp,k,2,jj)=twcc_fk(k,2,2,mlst,jj)
         cc(mm1,k,1,jj)=twcc_fk(k,1,1,mlst1,jj)
         cc(mp1,k,1,jj)=twcc_fk(k,2,1,mlst1,jj)
         cc(mm1,k,2,jj)=twcc_fk(k,1,2,mlst1,jj)
         cc(mp1,k,2,jj)=twcc_fk(k,2,2,mlst1,jj)
         cc(mm2,k,1,jj)=twcc_fk(k,1,1,mlst2,jj)
         cc(mp2,k,1,jj)=twcc_fk(k,2,1,mlst2,jj)
         cc(mm2,k,2,jj)=twcc_fk(k,1,2,mlst2,jj)
         cc(mp2,k,2,jj)=twcc_fk(k,2,2,mlst2,jj)
         cc(mm3,k,1,jj)=twcc_fk(k,1,1,mlst3,jj)
         cc(mp3,k,1,jj)=twcc_fk(k,2,1,mlst3,jj)
         cc(mm3,k,2,jj)=twcc_fk(k,1,2,mlst3,jj)
         cc(mp3,k,2,jj)=twcc_fk(k,2,2,mlst3,jj)
      enddo
      enddo

      enddo
c
cfong
      if( length_fft .eq. 0 )then
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+3,nx,lev*jlistnum*2,1)
      else
      do jj = 1, jlistnum
       call rfftmlt(cc(1,1,1,jj),gwk1,trigs,ifax,1,nx+3,nx,lev*2,1)
      end do
      end if
cfong
c
      do 22 jj=1,jlistnum
      do 22 k=1,lev
      do 22 i=1,nx
      ut(i,k,jj)= cc(i,k,1,jj)
      vt(i,k,jj)= cc(i,k,2,jj)
   22 continue

      return
      end
