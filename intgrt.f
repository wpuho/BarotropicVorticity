      subroutine intgrt
c
c  time integration use leapfog scheme
c
      include '../include/param.h'
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
      include '../include/semigrid.h'
c
      dimension wks(lev,2,jtrun,jtmax),wk(nx,lev,my_max),plnow(jtrun,jtmax,2)
      dimension cc(nx+3,lev,3+ncld,my_max),wss(lev,2,3+ncld,jtrun,jtmax)
      dimension vorten(lev,2,jtrun,jtmax),divten(lev,2,jtrun,jtmax)
     &,         fvt(nx,lev,my_max) 
      dimension rfvn(nx,lev,my_max),rfbt(nx,lev,my_max)
      
c
      logical forward
      data forward/.true./
      integer timef
c
      pi=4.0*atan(1.0)
      d2r=pi/180.
      r2d=1./d2r
      timef = 1
c
      iend=taue*3600./dt+0.001
      iout=tauo*3600./dt+0.001
      itrk=tauo*3600./dt+0.001
      nrec=1
c
c-----set initial field to old 
c
      do j=1,my
      do i=1,nx
        vorold(i,j)=vor(i,j)
      enddo
      enddo
c
      if(ktopo.eq.0)then
          topo=0.
      elseif(ktopo.eq.1)then
          nxmy8=nx*my*4
          open(66,file='/nwpr/gfs/ncust1/VORTEX/TAIWAN1.T720',
     &         access='direct'
     &,   form='unformatted',recl=nxmy8,status='unknown')
          call  reads(66,1,nx,my,topo)
          close(66)
          do j=1,my
          do i=1,nx
          topo(i,j)=topo(i,j)/grav
          enddo
          enddo
      endif

c-----gradient of terrain : topo(nx,my) into dlzl(nx,my),dtzl(nx,my) 
      call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,topo,plnow,nsize)
      call trngra ( jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,plnow,
     *              dlzl,dtzl,nsize )
c
c-----------------------
c
      if(myrank.eq.0)print *,'iend=',iend
c
      do 100 iter=1,iend
c
      taux=iter*dt/3600.
      tau=taux
c
      if(myrank.eq.0)print *,'forcast tau=',taux
c
c spectrum model 
c      
      do jj = 1, jlistnum
        j=jlist1(jj)
        xx=cosl(j)/rad
      do i = 1, nx
        rvor(i,1,jj)=vorold(i,j)
        rdiv(i,1,jj)=0.
        ut(i,1,jj)=u(i,j)*xx
        vt(i,1,jj)=v(i,j)*xx
      enddo
      enddo
c------------------------------------------------------------------------
      call joinrs(cc,rvor,rdiv,dummy,dummy,nx,my_max,lev,jlistnum,2,ncld) 
      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc
     *           ,wss,2,nsize)
      call ujoinrs(wss,vorspc,divspc,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
c
c-----gradient of vorticity : vor(nx,my) into dlpl(nx,my),dtpl(nx,my)
      vorten=0.
      divspc=0.
c
      call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,vor,plnow,nsize)
      call trngra ( jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,plnow,
     *              dlpl,dtpl,nsize ) 


c
c=====vortex vorticity tendency=====
      fvn =0.
      fbt =0.
c
      do jj = 1, jlistnum
        j=jlist1(jj)
      do i=1,nx
c
c-----earth vorticity advection
c
      if(fplane)then
      fbt(i,j)=0.
      else
      fbt(i,j)=-2.*omega*(vt(i,1,jj))
      endif
c
c-----relative vorticity advection
c
      fvn(i,j)=-(ut(i,1,jj))*dlpl(i,j)/cosl(j)**2
     &         -(vt(i,1,jj))*dtpl(i,j)      
c
c-----topographic beta effect
c
      ftn(i,j)=-((rvor(i,1,jj)+cor(j))/(hmean-topo(i,j)))
     &         *((ut(i,1,jj))*dlzl(i,j)/cosl(j)**2
     &         +(vt(i,1,jj))*dtzl(i,j))
c
      fvt(i,1,jj)=fbt(i,j)+fvn(i,j)
      rfbt(i,1,jj)=fbt(i,j)
      rfvn(i,1,jj)=fvn(i,j)
      enddo
      enddo  

c-----transform total vorticity tendency into vorten
      call joinrs(cc,fvt,rdiv,rfbt,rfvn,nx,my_max,lev,jlistnum,2,ncld) 
      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc
     *           ,wss,2,nsize)
      call ujoinrs(wss,vorten,divten,fbt,fvn,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)

c===================================
      if(forward)then
        dta=dt
      else
c-----leapfog time integration
        dta=2.*dt
      endif

      do 140 m=1,mlistnum
        mf=mlist(m)
      do 140 n=mf,jtrun
        if ( n.ne.1 ) then
          do 139 j = 1, 2
          do 139 k = 1, lev
            vorten(k,j,n,m)= dta*vorten(k,j,n,m) + vorspc(k,j,n,m)
 139      continue
        endif
 140  continue

c-----transform vorten into vorf(nx,my)
      divten=0.
      call joinsr(wss,vorten,divten,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
      call transr(jtrun,jtmax,nx,my,my_max,lev,poly,wss,cc,2,nsize)
      call ujoinsr(cc,rvor,rdiv,dummy,dummy,nx,my_max,lev,jlistnum,2,ncld)  
      do jj=1,jlistnum
        j=jlist1(jj)
      do i=1,nx
        vorf(i,j)=rvor(i,1,jj)
      enddo
      enddo
      call mpe_unify(vorf,nx,my,2,mpe_double)


c-------------------------------------------------
c-----robert time filter to get vorold
c
      if(.not. forward)then
      do j=1,my
      do i=1,nx
      vorold(i,j)=vor(i,j)+tfilt*(vorold(i,j)-2.0*vor(i,j)+vorf(i,j))
      enddo
      enddo
      endif

      forward=.false.
c------------------------------------------------
c-----transform vorf(nx,my) into vorspc for hdiffu

      do jj = 1, jlistnum
         j=jlist1(jj)
      do i = 1, nx
      rvor(i,1,jj)=vorf(i,j)
      rdiv(i,1,jj)=0.
      enddo
      enddo
      call joinrs(cc,rvor,rdiv,dummy,dummy,nx,my_max,lev,jlistnum,2,ncld) 
      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc
     *           ,wss,2,nsize)
      call ujoinrs(wss,vorspc,divspc,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)

c-------------------------------------------------
c-----horizontal diffusion for vor at time(t=t+1)
      call hdiffu
c
c-------------------------------------------------
c
c-----transform vorspc(output from hdiffu) into xx(nx,my)
     
      call joinsr(wss,vorspc,divspc,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
      call transr(jtrun,jtmax,nx,my,my_max,lev,poly,wss,cc,2,nsize)
      call ujoinsr(cc,rvor,rdiv,dummy,dummy,nx,my_max,lev,jlistnum,2,ncld)

      call tranuv ( jtrun,jtmax,nx,my,my_max,lev,onocos,wcfac,wdfac
     1            , poly,dpoly,vorspc,divspc,ut,vt,nsize )
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
      call joinsr(wss,strspc,strspc,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
      call transr(jtrun,jtmax,nx,my,my_max,lev,poly,wss,cc,2,nsize)
      call ujoinsr(cc,str,str,dummy,dummy,nx,my_max,lev,jlistnum,2,ncld)


c-----output at grid(nx,my)    
      do jj = 1, jlistnum
        j=jlist1(jj)
        xx=rad/cosl(j)
      do i = 1, nx
        vor(i,j)=rvor(i,1,jj)
        u(i,j)=ut(i,1,jj)*xx +basu(j)
        v(i,j)=vt(i,1,jj)*xx +basv(j)
        st(i,j)=str(i,1,jj)
      enddo
      enddo
      call mpe_unify(vor,nx,my,2,mpe_double)
      call mpe_unify(u,nx,my,2,mpe_double)
      call mpe_unify(v,nx,my,2,mpe_double)
      call mpe_unify(st,nx,my,2,mpe_double)


c-------------------------------------------------
c-------------------------------------------------
c-----vortex track & variables field output
c      
      call tracking(nx,my,vor,ix,jy)
      if(myrank.eq.0)print *,'after tracking ix,jy=',ix(1),jy(1)
      if(myrank.eq.0)print *,'vor=',vor(ix(1),jy(1))
      if(myrank.eq.0)print *,'u,v=',u(ix(1),jy(1)+5),v(ix(1),jy(1)+5)
c

      if(mod(iter,itrk).eq.0)then
        ntau=nrec 
        call tracking(nx,my,vor,ntau,ix,jy,xlat,xlon)
      endif

      if(mod(iter,iout).eq.0)then
        nrec=nrec+1
        call outflds(taux,nrec)
      endif
      if(mod(iter,iout).eq.0)then
        call mpe_unify(fvn,nx,my,2,mpe_double)
        call mpe_unify(fbt,nx,my,2,mpe_double)
        call writes(86,timef,nx,my,fvn)
        call writes(87,timef,nx,my,fbt)
        timef = timef +1
      endif
c
 100  continue
c
      if(myrank.eq.0)print *,' finish leapfor time integration'
c
      return
      end
