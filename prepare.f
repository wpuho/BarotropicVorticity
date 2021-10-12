      subroutine prepare
c
      include '../include/param.h'
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
      include '../include/const.h'
      include '../include/grid.h'
c
      namelist/modlst/hmean,taui,taue,tauo,dt,hfilt,tfilt,ktopo
     *,  vmax,rvm,b,tylat,tylon,fplane
      data hmean /10000./
      data taui/0.0/, taue/48.0/, tauo/1.0/, dt/240.0/
      data hfilt/1.0e16/, ktopo/0/, tfilt/0.15/
      data vmax/38./, rvm/60./, b/0.9/
      data tylat/23.5/, tylon/125./, fplane/.false./
c
      open(unit=1,file='namlsts720',form='formatted')
      read(1,modlst,end=40)
 40   if(myrank.eq.0)print modlst
c
      call cons
c
c-----XXX.OUT for grads
      nxmy8=nx*my*4
      open(81,file='/nwpr/gfs/ncust1/VORTEX/output/2020b/windu_v.OUT',
     &access='direct',form='unformatted',recl=nxmy8
     &,convert='big_endian',status='unknown')

      open(83,file='/nwpr/gfs/ncust1/VORTEX/output/2020b/vor.OUT',
     &access='direct',form='unformatted',recl=nxmy8
     &,convert='big_endian',status='unknown')

      open(84,file='/nwpr/gfs/ncust1/VORTEX/output/2020b/st.OUT',
     &access='direct',form='unformatted',recl=nxmy8
     &,convert='big_endian',status='unknown')

      open(85,file='/nwpr/gfs/ncust1/VORTEX/output/2020b/tracking.OUT',
     &access='direct',form='unformatted',recl=4*2
     &,convert='big_endian',status='unknown')

      open(86,file='/nwpr/gfs/ncust1/VORTEX/output/2020b/fvn.OUT',
     &access='direct',form='unformatted',recl=nxmy8
     &,convert='big_endian',status='unknown')

      open(87,file='/nwpr/gfs/ncust1/VORTEX/output/2020b/fbt.OUT',
     &access='direct',form='unformatted',recl=nxmy8
     &,convert='big_endian',status='unknown')

      open(88,file='/nwpr/gfs/ncust1/VORTEX/output/2020b/trackxy.OUT',
     &status='unknown')

      write(88,"(5(A7))") 'times' ,'ix(1)','jy(1)','ix(1)','jy(2)'
c-----vortex output track
c     open(88,file='/...../XXX.dat'
c    &, form='formatted',status='unknown')
c
      return
      end
