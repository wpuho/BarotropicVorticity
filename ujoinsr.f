      subroutine ujoinsr(cc,r1,r2,r3,r4,nx,my_max,lev,jlistnum,num,ncld)
c
      dimension cc(*), r1(*), r2(*), r3(*), r4(*) 
c
      if(num .eq. 2) call ujoin2sr(cc,r1,r2,nx,my_max,lev
     &                            ,jlistnum,ncld)
      if(num .eq. 3) call ujoin3sr(cc,r1,r2,r3,nx,my_max,lev
     &                            ,jlistnum,ncld)
      if(num .eq. 4) call ujoin4sr(cc,r1,r2,r3,r4,nx,my_max,lev
     &                            ,jlistnum,ncld)
c
      return
      end
c
      subroutine ujoin2sr(cc,r1,r2,nx,my_max,lev,jlistnum,ncld)
c
      dimension cc(nx+3,lev,1+ncld,my_max)
      dimension r1(nx,lev,my_max)
      dimension r2(nx,lev*ncld,my_max)
c
      do 10 jj =1, jlistnum
      do 10 k=1,lev
      do 10 i=1,nx
        r1(i,k,jj)=cc(i,k,1,jj)
   10 continue
c
      do 20 jj =1, jlistnum
      do 20 n=1,ncld
      nk=(n-1)*lev
      do 20 k=1,lev
      kk=nk+k
      do 20 i=1,nx
        r2(i,kk,jj)=cc(i,k,1+n,jj)
   20 continue
c
      return
      end
c
      subroutine ujoin3sr(cc,r1,r2,r3,nx,my_max,lev,jlistnum,ncld)
c
      dimension cc(nx+3,lev,2+ncld,my_max)
      dimension r1(nx,lev,my_max)
      dimension r2(nx,lev,my_max)
      dimension r3(nx,lev*ncld,my_max)
c
      do 10 jj =1, jlistnum
      do 10 k=1,lev
      do 10 i=1,nx
        r1(i,k,jj)=cc(i,k,1,jj)
        r2(i,k,jj)=cc(i,k,2,jj)
   10 continue
c
      do 20 jj =1, jlistnum
      do 20 n=1,ncld
      nk=(n-1)*lev
      do 20 k=1,lev
      kk=nk+k
      do 20 i=1,nx
        r3(i,kk,jj)=cc(i,k,2+n,jj)
   20 continue
c
      return
      end
c
      subroutine ujoin4sr(cc,r1,r2,r3,r4,nx,my_max,lev,jlistnum,ncld)
c
      dimension cc(nx+3,lev,3+ncld,my_max)
      dimension r1(nx,lev,my_max)
      dimension r2(nx,lev,my_max)
      dimension r3(nx,lev,my_max)
      dimension r4(nx,lev*ncld,my_max)
c
      do 10 jj =1, jlistnum
      do 10 k=1,lev
      do 10 i=1,nx
        r1(i,k,jj)=cc(i,k,1,jj)
        r2(i,k,jj)=cc(i,k,2,jj)
        r3(i,k,jj)=cc(i,k,3,jj)
   10 continue
c
      do 20 jj =1, jlistnum
      do 20 n=1,ncld
      nk=(n-1)*lev
      do 20 k=1,lev
      kk=nk+k
      do 20 i=1,nx
        r4(i,kk,jj)=cc(i,k,3+n,jj)
   20 continue
c
      return
      end
