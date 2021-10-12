      subroutine writes(nunit,nrec,nx,my,ary8)
      real*8 ary8(nx,my)
      real*4 ary4(nx,my)
      ary4=ary8
      write(nunit,rec=nrec)ary4
      return
      end
