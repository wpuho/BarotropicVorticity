      subroutine make_list
c
c  m        : start point of Lugendre number
c  mlist    : index of Fourier number on PE
c  mlistnum : index of Fourier number on PE
c
c

      include '../include/param.h'
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'

      do mm=1,jtmax
         mlist(mm)=0
      enddo

      mlistnum=0
      do m=1,jtmax*nsize
         nlist(m)=0 ; ilist(m)=0
      enddo

c
c for spectral space, use folded cyclic allocation type
c
 
      m=jtrun
      do mm=1,jtmax+1,2
         do ipe=nsize,1,-1
            if(m.gt.0)then
             if(ipe-1 .eq. myrank) then
               mlistnum=mlistnum+1
               mlist(mm)=m
             endif
               n=(ipe-1)*jtmax+mm
               nlist(m)=n 
               if(ipe-1 .eq. myrank) ilist(m)=mm
               m=m-1
            endif
         enddo
         do ipe=1,nsize
            if(m.gt.0)then
             if(ipe-1 .eq. myrank) then
               mlistnum=mlistnum+1
               mlist(mm+1)=m
             endif
               n=(ipe-1)*jtmax+mm+1
               nlist(m)=n
               if(ipe-1 .eq. myrank) ilist(m)=mm+1
               m=m-1
            endif
         enddo
      enddo


      jlistnum=0

      do jj=1,my_max
         jlist1(jj)=0
      enddo

      do j=1,my
         jlist2(j)=0
      enddo

c
c for grid space, use cyclic allocation type
c

      j=1
      do jj=1,my_max
         do ipe=1,nsize
            if(j.le.my)then
             if(ipe-1 .eq. myrank) then
               jlistnum=jlistnum+1
               jlist1(jj)=j
             endif
               j2=(ipe-1)*my_max+jj
               jlist2(j)=j2
               j=j+1
            endif
         enddo
      enddo
c

c  Usage:
c        do ipe=1,npe
c        do mm=1,mlistnum(ipe)
c           mf=mlist(mm,ipe); m=nlist(mf)
c           ...
c        enddo
c        enddo

      return
      end
