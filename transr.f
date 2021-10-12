      subroutine transr (jtrun,jtmax,nx,my,my_max,lev,poly,wss
     &                  ,cc,num,nsize)
c
c  subroutine to transform a spectral coefficient field to
c  grid point form
c
c *** input ***
c
c  jtrun: zonal wavenumber resolution limit
c  jtmax: maximum amount of zonal waves located in each pe
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  my_max: maximum amount of n-s grids located in each pe
c  lev: number of levels to transform
c  poly: legendre polynomials
c  wss: spectral coefficient array to transform
c  num: number of variables grouped together
c
c *** output ***
c
c  cc: 3-d output grid point fields
c
c  **************************************
c
      include '../include/index.h'
      include '../include/paramt.h'
      include '../include/fftcom.h'
c
      dimension poly(jtrun,my/2,jtmax)
      dimension cc(nx+3,lev,num,my_max),wss(lev,2,num,jtrun,jtmax)
c
csun  include '../include/paramt.h' .. change im,jm to nx,my
      dimension gwk1(nx+3,lev,num,my_max)
c
      dimension wcc_fk (lev,2,num,jtmax,my_max*nsize)
      dimension twcc_fk(lev,2,num,jtmax*nsize,my_max)
      dimension tcc(lev,2,num,my/2),tc2(lev,2,num,my/2)
      dimension ws2(lev,2,num,jtrun)
c
      mlx= (jtrun/2)*((jtrun+1)/2)
      myhalf=my/2
      lev2=lev*2
c
c-- do start
c
      do m=1,mlistnum
        mf=mlist(m)
        lmax=jtrun-mf+1
        lchk=iand(lmax,1)

      do l=mf,jtrun-1,2
      do k=1,lev2*num
      ws2(k,1,1,l) = wss(k,1,1,l,m)
      ws2(k,1,1,l+1) = -wss(k,1,1,l+1,m)
      enddo
      enddo
      if (lchk.eq.1) then
        l=jtrun
        do k=1,lev2*num
          ws2(k,1,1,l) = -wss(k,1,1,l,m)
        enddo
      endif
c............................  replaced by two dgemm (BLAS3 routine)
c     do k=1,lev2*num*myhalf
c      tcc(k,1,1,1)=0.
c      tc2(k,1,1,1)=0.
c     enddo
c
c     do j=1,myhalf
c       do l=mf,jtrun
c       do k=1,lev2*num
c         tcc(k,1,1,j)=tcc(k,1,1,j)+poly(l,j,m)*wss(k,1,1,l,m)
c         tc2(k,1,1,j)=tc2(k,1,1,j)+poly(l,j,m)*ws2(k,1,1,l)
c       enddo
c       enddo
c     enddo
cfong--beg
c     call dgemm('N','N',lev2*num,myhalf,lmax,
c    &             1.0, wss(1,1,1,mf,m),lev2*num,
c    &                 poly(mf,1,m),jtrun,
c    &             0.0,tcc,lev2*num)
c     call dgemm('N','N',lev2*num,myhalf,lmax,
c    &             1.0, ws2(1,1,1,mf),  lev2*num,
c    &                 poly(mf,1,m),jtrun,
c    &             0.0, tc2,lev2*num)
cfong--end
cfong-- cache blocking method
c
        do k=1,lev2*num*myhalf
         tcc(k,1,1,1)=0.
         tc2(k,1,1,1)=0.
        enddo
c
c  nb : multiple of 2
c
        nb=32
        jchk= iand(myhalf, 1)
        jje = myhalf-jchk
c
        do jj=1,jje,nb
        do kk=1,lev2*num,nb
        do ll=mf,jtrun,nb

        do j = jj,min(jj+nb-1,jje), 2
          do k = kk,min(kk+nb-1,lev2*num),2
            sa00=tcc(k  ,1,1,j)
            sa01=tcc(k  ,1,1,j+1)
            sa10=tcc(k+1,1,1,j)
            sa11=tcc(k+1,1,1,j+1)
            sb00=tc2(k  ,1,1,j)
            sb01=tc2(k  ,1,1,j+1)
            sb10=tc2(k+1,1,1,j)
            sb11=tc2(k+1,1,1,j+1)
            do l = ll,min(ll+nb-1,jtrun)
              sa00=sa00+poly(l,j,  m)*wss(k,1,1,l,m)
              sa01=sa01+poly(l,j+1,m)*wss(k,1,1,l,m)
              sa10=sa10+poly(l,j,  m)*wss(k+1,1,1,l,m)
              sa11=sa11+poly(l,j+1,m)*wss(k+1,1,1,l,m)
              sb00=sb00+poly(l,j,  m)*ws2(k,1,1,l)
              sb01=sb01+poly(l,j+1,m)*ws2(k,1,1,l)
              sb10=sb10+poly(l,j,  m)*ws2(k+1,1,1,l)
              sb11=sb11+poly(l,j+1,m)*ws2(k+1,1,1,l)
            enddo
            tcc(k,1,1,j  )=sa00
            tcc(k,1,1,j+1)=sa01
            tcc(k+1,1,1,j  )=sa10
            tcc(k+1,1,1,j+1)=sa11
            tc2(k,1,1,j)=sb00
            tc2(k,1,1,j+1)=sb01
            tc2(k+1,1,1,j)=sb10
            tc2(k+1,1,1,j+1)=sb11
          enddo
        enddo

        enddo
        enddo
        enddo
c
c odd number
c
        if( jchk.eq.1 )then
          j = myhalf
        do kk=1,lev2*num,nb
        do ll=mf,jtrun,nb
          do k = kk,min(kk+nb-1,lev2*num),2
            sa00=tcc(k,1,1,j)
            sa10=tcc(k+1,1,1,j)
            sb00=tc2(k,1,1,j)
            sb10=tc2(k+1,1,1,j)
            do l = ll,min(ll+nb-1,jtrun)
              sa00=sa00+poly(l,j,  m)*wss(k,1,1,l,m)
              sa10=sa10+poly(l,j,  m)*wss(k+1,1,1,l,m)
              sb00=sb00+poly(l,j,  m)*ws2(k,1,1,l)
              sb10=sb10+poly(l,j,  m)*ws2(k+1,1,1,l)
            enddo
            tcc(k,1,1,j  )=sa00
            tcc(k+1,1,1,j  )=sa10
            tc2(k,1,1,j)=sb00
            tc2(k+1,1,1,j)=sb10
          enddo
        enddo
        enddo
        endif

cfong--end


      do j=1, myhalf
      jj=jlist2(j)
      jx=my-j+1
      j2=jlist2(jx)
      do k=1,lev2*num
       wcc_fk(k,1,1,m,jj)=tcc(k,1,1,j)
       wcc_fk(k,1,1,m,j2)=tc2(k,1,1,j)
      enddo
      enddo

      enddo
c--do end

      call mpe_transpose_sr(wcc_fk,twcc_fk,lev*2*num,jtmax,my_max,nsize)

      do jj =1,jlistnum

      do ii=1,num
      do k=1,lev
      do i=1,nx+3
       cc(i,k,ii,jj)=0.
      enddo
      enddo
      enddo

cfong--
      mchk=iand(jtrun,3)

      do m=1,mchk
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
      do ii=1,num
      do k=1,lev
         cc(mm,k,ii,jj)=twcc_fk(k,1,ii,mlst,jj)
         cc(mp,k,ii,jj)=twcc_fk(k,2,ii,mlst,jj)
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
         cc(mm,k,ii,jj)=twcc_fk(k,1,ii,mlst,jj)
         cc(mp,k,ii,jj)=twcc_fk(k,2,ii,mlst,jj)
         cc(mm1,k,ii,jj)=twcc_fk(k,1,ii,mlst1,jj)
         cc(mp1,k,ii,jj)=twcc_fk(k,2,ii,mlst1,jj)
         cc(mm2,k,ii,jj)=twcc_fk(k,1,ii,mlst2,jj)
         cc(mp2,k,ii,jj)=twcc_fk(k,2,ii,mlst2,jj)
         cc(mm3,k,ii,jj)=twcc_fk(k,1,ii,mlst3,jj)
         cc(mp3,k,ii,jj)=twcc_fk(k,2,ii,mlst3,jj)
      enddo
      enddo
      enddo

      enddo

c
cfong
      if( length_fft .eq. 0 )then
      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+3,nx,lev*jlistnum*num,1)
      else
      do jj = 1, jlistnum
       call rfftmlt(cc(1,1,1,jj),gwk1,trigs,ifax,1,nx+3,nx,lev*num,1)
      end do
      end if
cfong
c
   20 continue

      return
      end
