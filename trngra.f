      subroutine trngra (jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,s
     *                  ,dlpl,dtpl,nsize)
c
c  subroutine to transform spectral terrain pressure to grid point
c  fields of zonal and meridional derivatives of terrain pressure
c
c *** input ***
c
c  jtrun: zonal wavenumber resolution limit
c  jtmax: maximum amount of zonal waves located in each pe
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  my_max: maximum amount of n-s grids located in each pe
c  cim: zonal wavenumber array
c  poly: legendre polynomials
c  dpoly: d(poly)/d(sin(lat))
c  s: spectral coefficient array
c  nsiz: pe number
c
c **** output ****
c
c  dlpl: d(pt)/d(longitude)
c  dtpl: d(pt)/d(sin(lat))
c
c ****************************************************
c
      include '../include/index.h'
      include '../include/paramt.h'
      include '../include/fftcom.h'
c
      dimension poly(jtrun,my/2,jtmax),dpoly(jtrun,my/2,jtmax)
      dimension s(jtrun,jtmax,2),dlpl(nx,my),dtpl(nx,my)
      dimension cim(jtmax)
c
csun  include '../include/paramt.h'
      dimension cc(nx+3,2,my_max)
c
      dimension gwk1(nx+3,1,6,my_max)
c
      dimension twcc_fk(my_max,jtmax*nsize,2)
      dimension twdd_fk(my_max,jtmax*nsize,2)
cibm---deleted
cibm     dimension wss(jtrun,jtmax,2)

      dimension wcu_fk(my_max*nsize,jtmax,2)
     &         ,wcv_fk(my_max*nsize,jtmax,2)
      dimension wcu_t(my,2),wcv_t(my,2)

      dimension ws3(jtrun,2,2)
      dimension ws4(jtrun,2,2)

c
      myhalf=my/2

      do m=1,mlistnum
         mf=mlist(m)

      do 55 j=1,my
      wcu_fk(j,m,1)= 0.0
      wcu_fk(j,m,2)= 0.0
      wcv_fk(j,m,1)= 0.0
      wcv_fk(j,m,2)= 0.0
   55 continue
c
      do l=mf,jtrun
      ws3(l,1,1) = +s(l,m,2)
      ws3(l,2,1) = -s(l,m,1)
      ws4(l,1,2) = -s(l,m,1)
      ws4(l,2,2) = -s(l,m,2)
      enddo

      do j=1,my*2
      wcu_t(j,1)=0.
      wcv_t(j,1)=0.
      enddo

      do l=mf,jtrun
*vocl loop,repeat(jmhalf)
      do j=1,myhalf
      wcu_t(j,1) = wcu_t(j,1) + ws3(l,1,1)*(cim(m)* poly(l,j,m))
      wcu_t(j,2) = wcu_t(j,2) + ws3(l,2,1)*(cim(m)* poly(l,j,m))
      wcv_t(j,1) = wcv_t(j,1) + ws4(l,1,2)*dpoly(l,j,m)
      wcv_t(j,2) = wcv_t(j,2) + ws4(l,2,2)*dpoly(l,j,m)
      enddo
      enddo

      do l=mf,jtrun,2
      ws3(l,1,1) = +s(l,m,2)
      ws3(l,2,1) = -s(l,m,1)
      ws4(l,1,2) = +s(l,m,1)
      ws4(l,2,2) = +s(l,m,2)
      enddo

      do l=mf+1,jtrun,2
      ws3(l,1,1) = -s(l,m,2)
      ws3(l,2,1) = +s(l,m,1)
      ws4(l,1,2) = -s(l,m,1)
      ws4(l,2,2) = -s(l,m,2)
      enddo

      do l=mf,jtrun
*vocl loop,repeat(jmhalf)
      do jj=myhalf+1,my
      j=my-jj+1
      wcu_t(jj,1) = wcu_t(jj,1) + ws3(l,1,1)*(cim(m)* poly(l,j,m))
      wcu_t(jj,2) = wcu_t(jj,2) + ws3(l,2,1)*(cim(m)* poly(l,j,m))
      wcv_t(jj,1) = wcv_t(jj,1) + ws4(l,1,2)*dpoly(l,j,m)
      wcv_t(jj,2) = wcv_t(jj,2) + ws4(l,2,2)*dpoly(l,j,m)
      enddo
      enddo

      do j=1,my
      jj=jlist2(j)
      wcu_fk(jj,m,1)=wcu_t(j,1)
      wcu_fk(jj,m,2)=wcu_t(j,2)
      wcv_fk(jj,m,1)=wcv_t(j,1)
      wcv_fk(jj,m,2)=wcv_t(j,2)
      enddo

      enddo

       call mpe_transpose_rs1(wcu_fk,twcc_fk,my_max,jtmax,2,nsize)
       call mpe_transpose_rs1(wcv_fk,twdd_fk,my_max,jtmax,2,nsize)

      do jj=1,jlistnum
      do i=1,nx+3
      cc(i,1,jj)= 0.
      cc(i,2,jj)= 0.
      enddo
      enddo

      do jj=1,jlistnum
      do m=1,jtrun
         mm= 2*m-1
         mp= mm+1
         mlst=nlist(m)
         cc(mm,1,jj)=twcc_fk(jj,mlst,1)
         cc(mp,1,jj)=twcc_fk(jj,mlst,2)
         cc(mm,2,jj)=twdd_fk(jj,mlst,1)
         cc(mp,2,jj)=twdd_fk(jj,mlst,2)
      enddo
      enddo

      call rfftmlt(cc,gwk1,trigs,ifax,1,nx+3,nx,jlistnum*2,1)

      do 22 jj=1,jlistnum
      j=jlist1(jj)
!ocl novrec
      do 22 i=1,nx
      dlpl(i,j)= -cc(i,1,jj)
      dtpl(i,j)= -cc(i,2,jj)
   22 continue
c
      return
      end
