      subroutine reads(nunit,nrec,nx,my,ary8)
      real*8 ary8(nx,my)
      real*4 ary4(nx,my)
      read(nunit,rec=nrec)ary4
      ary8=ary4
      return
      end
