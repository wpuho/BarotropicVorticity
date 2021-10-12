      subroutine joinrs(cc,r1,r2,r3,r4,nx,my_max,lev,jlistnum,num,ncld)
c
      dimension cc(*), r1(*), r2(*), r3(*), r4(*) 
c
      if(num .eq. 2) call join2rs(cc,r1,r2,nx,my_max,lev,jlistnum,ncld)
      if(num .eq. 3) call join3rs(cc,r1,r2,r3,nx,my_max,lev,jlistnum,ncld)
c
      return
      end
c
      subroutine join2rs(cc,r1,r2,nx,my_max,lev,jlistnum,ncld)
c
      dimension cc(nx+3,lev,1+ncld,my_max)
      dimension r1(nx,lev,my_max)
      dimension r2(nx,lev*ncld,my_max)
c
      do 10 jj =1, jlistnum
      do 10 k=1,lev
      do 10 i=1,nx
      cc(i,k,1,jj)= r1(i,k,jj)
   10 continue
c
      do 20 jj =1, jlistnum
      do 20 n=1,ncld
      nk=(n-1)*lev
      do 20 k=1,lev
      kk=nk+k
      do 20 i=1,nx
      cc(i,k,1+n,jj)= r2(i,kk,jj)
   20 continue
c
      return
      end
c
      subroutine join3rs(cc,r1,r2,r3,nx,my_max,lev,jlistnum,ncld)
c
      dimension cc(nx+3,lev,2+ncld,my_max)
      dimension r1(nx+3,lev,my_max)
      dimension r2(nx+3,lev,my_max)
      dimension r3(nx+3,lev*ncld,my_max)
c
      do 10 jj =1, jlistnum
      do 10 k=1,lev
      do 10 i=1,nx
      cc(i,k,1,jj)= r1(i,k,jj)
      cc(i,k,2,jj)= r2(i,k,jj)
   10 continue
c
      do 20 jj =1, jlistnum
      do 20 n=1,ncld
      nk=(n-1)*lev
      do 20 k=1,lev
      kk=nk+k
      do 20 i=1,nx
      cc(i,k,2+n,jj)= r3(i,kk,jj)
   20 continue
c
      return
      end
