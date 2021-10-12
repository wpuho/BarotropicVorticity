      subroutine outflds(taux,nrec)
c
c  output grid point value to a file
c
      include '../include/param.h'
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
      include '../include/const.h'
      include '../include/grid.h'
      integer record
      integer uwindrecord,vwindrecord
c
c  working array

      itau=taux
c     nrec = itau / tauo +1
      record = itau / tauo +1
      vwindrecord = 2*record
      uwindrecord = vwindrecord -1
c

      call  writes(81,uwindrecord,nx,my,u)
      call  writes(81,vwindrecord,nx,my,v)
      call  writes(83,record,nx,my,vor)
      call  writes(84,record,nx,my,st)
      write(85,rec=uwindrecord) ix(1)
      write(85,rec=vwindrecord) jy(1)

      write(88,"(5(I7))") record,ix(1),jy(1),ix(2),jy(2)
      return
      end
