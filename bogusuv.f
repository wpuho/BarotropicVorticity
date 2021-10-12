      subroutine bogusuv(nx,my,lev,my_max,rad,ut,vt,rvor,rdiv,cosl
     &,                  tylat,tylon,vmax,rvm,b,ix,jy,xlon,xlat)
c
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
c
      dimension ut(nx,lev,my_max),vt(nx,lev,my_max),rvor(nx,lev,my_max)
     *,  rdiv(nx,lev,my_max),cosl(my),ix(3),jy(3),xlon(nx+1,my),xlat(my)
     
      real xdim , ydim ,distance , vofr, xx ,xdis ,ydis
      real dislat(3),dislon(3)
      integer ty

c
      radx=rad/1000.
      pi=4.0*atan(1.0)
      d2r=pi/180.
      r2d=1./d2r
      

 200 continue
c
      ix(1)=tylon/(360./nx)+1             
      jy(1)=(tylat-xlat(1))/(180./my)+1
      ix(2)=ix(1)+25
      jy(2)=jy(1)+25
      ix(3)=ix(1)
      jy(3)=jy(1)
      if(myrank.eq.0)print *,'vmax=',vmax,' rvm=',rvm,' b=',b
      if(myrank.eq.0)print *,'tylon=',tylon,' tylat=',tylat
      if(myrank.eq.0)print *,' ix=',ix(1),' jy=',jy(1)

      ut=0.
      vt=0.
      rvor=0.
c      
      do 101 ty =1,2
         if(ty.eq.1) then
             vmax = 50.
                b  = 1.1673007300730074
              rvm = 100.
         elseif(ty.eq.2) then
              vmax = 20.
                 b = 0.31907790779077905
               rvm = 100.
         endif

      do 100 jj=1,jlistnum
        j=jlist1(jj)
        xx = cosl(j)/rad   
        ydis =  radx* (xlat(j)-((jy(ty)-1.)*(180./my)+xlat(1)))*d2r
      do 100 i=1,nx
        xdis =  radx*cosl(j)*(xlon(i,j)-((ix(ty)-1)*360./nx))* d2r 
        distance =( xdis**2 + ydis**2 )**(0.5)
        if (ty.eq.1) then
            if     ( distance .lt. 300.) then
              vofr = 35*(distance/rvm)*exp((1-(distance/rvm)**0.9)/0.9)
            elseif (distance .ge. 300.)  then
               vofr = vmax*(distance/rvm)*exp((1-(distance/rvm)**b)/b)
            endif
        elseif (ty.eq.2) then
            if     ( distance .lt. 300.) then
              vofr = 35*(distance/rvm)*exp((1-(distance/rvm)**0.9)/0.9)
            elseif (distance .ge. 300.)  then
               vofr = vmax*(distance/rvm)*exp((1-(distance/rvm)**b)/b)
            endif
         endif
        
        rdiv(i,1,jj) = 0.
        ut(i,1,jj) = vofr *(ydis/distance)*xx *(-1)  +ut(i,1,jj)
        vt(i,1,jj) = vofr *(xdis/distance)*xx        +vt(i,1,jj)
        rvor(i,1,jj)=2.*(vmax/(rvm*1000.))*(1.-0.5*(distance/rvm)**(b))
     &               * exp((1-(distance/rvm)**b) /b )+rvor(i,1,jj)
 100  continue
 101  continue

c
      return
      end
