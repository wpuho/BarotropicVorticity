      subroutine hdiffu
c
      include '../include/param.h'
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
c
c  compute diffusion coefficients
c
      dta=dt
      hfilt1=hfilt
c      hfilt2=hfilt
c      hfilt3=hfilt
c
c  diffuse vorticity and divergence fields
c
      k=1
c
      do 30 m=1,mlistnum
         mf=mlist(m)
      do 30 n=mf,jtrun
        c1=1.+dta*hfilt1*eps4(n,m)**2
c        c2=1.+dta*hfilt2*eps4(n,m)**2
c        c3=1.+dta*hfilt3*eps4(n,m)**2
c        c1=1.+dta*hfilt1*eps4(n,m)**3
c        c2=1.+dta*hfilt2*eps4(n,m)**3
c        c3=1.+dta*hfilt3*eps4(n,m)**3
        vorspc(k,1,n,m)=vorspc(k,1,n,m)/c1
        vorspc(k,2,n,m)=vorspc(k,2,n,m)/c1
c        divspc(k,1,n,m)=divspc(k,1,n,m)/c2
c        divspc(k,2,n,m)=divspc(k,2,n,m)/c2
c        phispc(k,1,n,m)=phispc(k,1,n,m)/c3
c        phispc(k,2,n,m)=phispc(k,2,n,m)/c3
 30   continue
c
      return
      end
