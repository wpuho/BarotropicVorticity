      subroutine getrdy
c
      include '../include/param.h'
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
c
      dimension cc(nx+3,lev,3+ncld,my_max),wss(lev,2,3+ncld,jtrun,jtmax)
      dimension wks1(lev,2,jtrun,jtmax),tmp(nx,lev,my_max)
c
      radx=rad/1000.
      pi=4.0*atan(1.0)
      d2r=pi/180.
      r2d=1./d2r
c
      do 15 j=1,my
      xlon(1,j)=0.
      do 10 i=2,nx
        xlon(i,j)=xlon(1,j)+float(i-1)*360./nx
 10   continue
        xlat(j)=asin(sinl(j))*r2d
 15   continue

c-----ideal vortex
c-----../include/grid.h rvor(nx,lev,my_max),rdiv(nx,lev,my_max)
c-----../include/spec.h vorspc(lev,2,jtrun,jtmax),divspc(lev,2,jtrun,jtmax)                                     
c-----                  strspc(lev,2,jtrun,jtmax)
c-----
      call bogusuv(nx,my,lev,my_max,rad,ut,vt,rvor,rdiv,cosl
     &            ,tylat,tylon,vmax,rvm,b,ix,jy,xlon,xlat)

c-----transform rvor into spectral form
      call joinrs(cc,rvor,rdiv,dummy,dummy,nx,my_max,lev,jlistnum,2,ncld)
      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc
     *           ,wss,2,nsize)
      call ujoinrs(wss,vorspc,divspc,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
 

c-----compute stream function
      do 50 m=1,mlistnum
         mf=mlist(m)
      do 50 n=mf,jtrun
        if ( n.ne.1 ) then
        do 51 j = 1, 2
        do 51 k = 1, lev
          strspc(k,j,n,m) = -vorspc(k,j,n,m)/eps4(n,m)
  51    continue
        endif
  50  continue

c-----transform strspc into grid form 

      call joinsr(wss,strspc,strspc,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
      call transr(jtrun,jtmax,nx,my,my_max,lev,poly,wss,cc,2,nsize)
      call ujoinsr(cc,str,str,dummy,dummy,nx,my_max,lev,
     &             jlistnum,2,ncld)

c--------------------------------
c-----../include/grid.h u(nx,my),v(nx,my),vor(nx,my),div(nx,my),st(nx,my)          
c

      do jj = 1, jlistnum
        j=jlist1(jj)
        xx=rad/cosl(j)
      do i = 1, nx
        vor(i,j)=rvor(i,1,jj)
        u(i,j)=ut(i,1,jj)*xx  +basu(j)
        v(i,j)=vt(i,1,jj)*xx  +basv(j)
        st(i,j)=str(i,1,jj)
      enddo
      enddo
   
      call mpe_unify(vor,nx,my,2,mpe_double)
      call mpe_unify(u,nx,my,2,mpe_double)
      call mpe_unify(v,nx,my,2,mpe_double)
      call mpe_unify(st,nx,my,2,mpe_double)

      if(myrank.eq.0)print *,'vor=',vor(ix(1),jy(1))
      if(myrank.eq.0)print *,'u=',u(ix(1),jy(1)+5)
      if(myrank.eq.0)print *,'v=',v(ix(1),jy(1)+5)
      if(myrank.eq.0)print *,'st=',st(ix(1),jy(1))
c
c-------------------------------
c-----vortex track
c


c
c-----output of initail field
c
      
      call outflds(0.,1)
c
c
      return
      end

