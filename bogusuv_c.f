      subroutine bogusuv(nx,my,lev,my_max,rad,ut,vt,rvor,rdiv,cosl
     &,                  tylat,tylon,vmax,rvm,b,ix,jy,xlon,xlat)
c
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
c
      dimension ut(nx,lev,my_max),vt(nx,lev,my_max),rvor(nx,lev,my_max)
     *,   rdiv(nx,lev,my_max),cosl(my),ix(3),jy(3),xlon(nx+1,my),xlat(my)
c
      radx=rad/1000.
      pi=4.0*atan(1.0)
      d2r=pi/180.
      r2d=1./d2r
c
      ix(1)=tylon/(360./nx)+1
      jy(1)=(tylat-xlat(1))/(180./my)+1
      ix(2)=ix(1)
      jy(2)=jy(1)
      ix(3)=ix(1)
      jy(3)=jy(1)
      if(myrank.eq.0)print *,'vmax=',vmax,' rvm=',rvm,' b=',b
      if(myrank.eq.0)print *,'tylon=',tylon,' tylat=',tylat
      if(myrank.eq.0)print *,' ix=',ix(1),' jy=',jy(1)
c      
      do 100 jj=1,jlistnum
        j=jlist1(jj)
		xx = cosl(j)/rad   
        ry = radx*(xlat(j)-tylat)*d2r
      do 100 i=1,nx
	    rx = radx*cosl(j)*(xlon(i,j)-tylon)* d2r 
        distance =(rx**2 + ry**2 )**(0.5)
	      
		  vorticity = vmax*(distance/rvm)*exp((1-(distance/rvm)**b)/b)
	  
          ut(i,1,jj) = vorticity*(ydis/distance)*xx *(-1)
              
          vt(i,1,jj)=.....
              
          rvor(i,1,jj)=.....
 100  continue
c
      return
      end
