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
      dimension rvorten(nx,lev,my_max) 
c
      logical forward
      data forward/.true./
c
      pi=4.0*atan(1.0)
      d2r=pi/180.
      r2d=1./d2r
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
     &         access='direct',form='unformatted',recl=nxmy8,
     &         status='unknown')
          call reads(66,1,nx,my,topo)
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

      taux=iter*dt/3600.
      tau=taux
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

      call joinrs(cc,rvor,rdiv,dummy,dummy,nx,my_max,lev,
     &            jlistnum,2,ncld)
      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc
     *           ,wss,2,nsize)
      call ujoinrs(wss,vorspc,divspc,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
c------------------------------------------------------------------------
c
c-----gradient of vorticity : vor(nx,my) into dlpl(nx,my),dtpl(nx,my)

      vorten=0.
      divspc=0.
      call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,vor,plnow
     &             ,nsize)
      call trngra ( jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,plnow,
     *              dlpl,dtpl,nsize )

c
c=====vortex vorticity tendency=====
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
               -(vt(i,1,jj))*dtpl(i,j)
c
c-----topographic beta effect
c
      ftn(i,j)= -(vor(i,j)+cor(j))*(1/(hmean-topo(i,j)))*
     &          ((ut(i,1,jj))*dlzl(i,j)/cosl(j)**2 +
     &          (vt(i,1,jj))*dtzl(i,j))
      rvorten(i,1,jj) = fvn(i,j) +fbt(i,j)+ftn(i,j) 
      enddo
      enddo  
c-----transform total vorticity tendency into vorten
    
      call joinrs(cc,rvorten,rdiv,dummy,dummy,nx,my_max,lev,
     &            jlistnum,2,ncld)
      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc
     *           ,wss,2,nsize)
      call ujoinrs(wss,vorten,divten,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
c===================================
      if(forward)then
        dta=dt
      else
        dta=2.*dt
      endif
c-----leapfog time integration

      do 140 m=1,mlistnum
        mf=mlist(m)
      do 140 n=mf,jtrun
        if ( n.ne.1 ) then
          do 139 j = 1, 2
          do 139 k = 1, lev
          vorten(k,j,n,m)= dta*vorten(k,j,n,m) + vorspcold(k,j,n,m)
 139  continue
        endif
 140  continue
c-----transform vorten into vorf(nx,my)
      divten =0.

      call joinsr(wss,vorten,divten,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
      call transr(jtrun,jtmax,nx,my,my_max,lev,poly,wss,cc,2,nsize)
      call ujoinsr(cc,rvor,rdiv,dummy,dummy,nx,my_max,lev,jlistnum
     &    ,2, ncld)

      do jj = 1, jlistnum
        j=jlist1(jj)
      do i = 1, nx
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
   
c     call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,vorf,
c    *             vorspc,nsize)
      do jj = 1, jlistnum
        j=jlist1(jj)
      do i = 1, nx
        rvorten(i,1,jj) = vorold(i,j)
      enddo
      enddo

      call joinrs(cc,rvorten,rvorten,dummy,dummy,nx,my_max,lev,
     &            jlistnum,2,ncld)
      call tranrs(jtrun,jtmax,nx,my,my_max,lev,poly,weight,cc
     *           ,wss,2,nsize)
      call ujoinrs(wss,vorspcold,vorspcold,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)

      
      vorspc = vorten

c-------------------------------------------------
c-----horizontal diffusion for vor at time(t=t+1)
      call hdiffu
c-----transform vorspc(output from hdiffu) into xx(nx,my)
      call tranuv (jtrun,jtmax,nx,my,my_max,lev,onocos,wcfac,wdfac,ploy
     *             ,dpoly,vorspc,divspc,ut,vt,nsize )
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
      call joinsr(wss,strspc,vorspc,dummy,dummy,jtrun,jtmax,lev
     *           ,mlistnum,2,ncld)
      call transr(jtrun,jtmax,nx,my,my_max,lev,poly,wss,cc,2,nsize)
      call ujoinsr(cc,str,rvor,dummy,dummy,nx,my_max,lev,jlistnum,2,ncld)

c-----output at grid(nx,my)    
      vor =0.
      do jj = 1, jlistnum
        j=jlist1(jj)
        xx=rad/cosl(j)
      do i = 1, nx
        vor(i,j)=rvor(i,1,jj)
        u(i,j)=ut(i,1,jj)*xx
        v(i,j)=vt(i,1,jj)*xx
        st(i,j)=str(i,1,jj)
      enddo
      enddo
   
      call mpe_unify(vor,nx,my,2,mpe_double)
      call mpe_unify(u,nx,my,2,mpe_double)
      call mpe_unify(v,nx,my,2,mpe_double)
      call mpe_unify(st,nx,my,2,mpe_double)
      
c-------------------------------------------------
c-----vortex track & variables field output
c      
c     call tracking(nx,my,vor,ix,jy)
c
      if (myrank.eq.0) then
      if (amod(taux,tauo).eq. 0.)then
         print *,'ing ix iy tau', ix(1),jy(1),taux
         call outflds(taux,nrec) 
      endif
      endif
     
c
 100  continue
c
      if(myrank.eq.0)print *,' finish leapfor time integration'
c
      return
      end
