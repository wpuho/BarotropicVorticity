      subroutine ujoinrs(wss,s1,s2,s3,s4,jtrun,jtmax,lev
     &                  ,mlistnum,num,ncld)
c
      dimension wss(*), s1(*), s2(*), s3(*), s4(*) 
c
      if(num .eq. 2) 
     &  call ujoin2rs(wss,s1,s2,jtrun,jtmax,lev,mlistnum,ncld)
      if(num .eq. 3) 
     &  call ujoin3rs(wss,s1,s2,s3,jtrun,jtmax,lev,mlistnum,ncld)
      if(num .eq. 4) 
     &  call ujoin4rs(wss,s1,s2,s3,s4,jtrun,jtmax,lev,mlistnum,ncld)
c
      return
      end
c
      subroutine ujoin2rs(wss,s1,s2,jtrun,jtmax,lev,mlistnum,ncld)
c
      dimension wss (lev,2,1+ncld,jtrun,jtmax)
c
      dimension s1(lev,2,jtrun,jtmax)
      dimension s2(lev*ncld,2,jtrun,jtmax)
c
      do 10 m=1,mlistnum
      do 10 l=1,jtrun
      do 10 k = 1, lev*2
        s1(k,1,l,m) = wss(k,1,1,l,m) 
  10  continue
c
      do 20 m=1,mlistnum
      do 20 l=1,jtrun
      do 20 n=1,ncld
        nk=(n-1)*lev
      do 20 k = 1, lev
        kk=nk+k
        s2(kk,1,l,m) = wss(k,1,1+n,l,m)
        s2(kk,2,l,m) = wss(k,2,1+n,l,m)
  20  continue

c
      return
      end
c
      subroutine ujoin3rs(wss,s1,s2,s3,jtrun,jtmax,lev,mlistnum,ncld)
c
      dimension wss (lev,2,2+ncld,jtrun,jtmax)
c
      dimension s1(lev,2,jtrun,jtmax)
      dimension s2(lev,2,jtrun,jtmax)
      dimension s3(lev*ncld,2,jtrun,jtmax)
c
      do 10 m=1,mlistnum
      do 10 l=1,jtrun
      do 10 k = 1, lev*2
        s1(k,1,l,m) = wss(k,1,1,l,m) 
        s2(k,1,l,m) = wss(k,1,2,l,m)
  10  continue
c
      do 20 m=1,mlistnum
      do 20 l=1,jtrun
      do 20 n=1,ncld
        nk=(n-1)*lev
      do 20 k = 1, lev
        kk=nk+k
        s3(kk,1,l,m) = wss(k,1,2+n,l,m)
        s3(kk,2,l,m) = wss(k,2,2+n,l,m)
  20  continue
c
      return
      end
c
      subroutine ujoin4rs(wss,s1,s2,s3,s4,jtrun,jtmax,lev,mlistnum,ncld)
c
      dimension wss (lev,2,3+ncld,jtrun,jtmax)
c
      dimension s1(lev,2,jtrun,jtmax)
      dimension s2(lev,2,jtrun,jtmax)
      dimension s3(lev,2,jtrun,jtmax)
      dimension s4(lev*ncld,2,jtrun,jtmax)
c
      do 10 m=1,mlistnum
      do 10 l=1,jtrun
      do 10 k = 1, lev*2
        s1(k,1,l,m) = wss(k,1,1,l,m) 
        s2(k,1,l,m) = wss(k,1,2,l,m)
        s3(k,1,l,m) = wss(k,1,3,l,m)
  10  continue
c
      do 20 m=1,mlistnum
      do 20 l=1,jtrun
      do 20 n=1,ncld
        nk=(n-1)*lev
      do 20 k = 1, lev
        kk=nk+k
        s4(kk,1,l,m) = wss(k,1,3+n,l,m)
        s4(kk,2,l,m) = wss(k,2,3+n,l,m)
  20  continue
c
      return
      end
