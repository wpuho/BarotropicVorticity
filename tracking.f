      subroutine tracking(nx,my,vor,ix,jy)
c
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
c
      dimension vor(nx,my)
      dimension ix(3),jy(3)
      data ixyrange/5/
c
      do ty = 1,2
      ib=ix(ty)-ixyrange
      ie=ix(ty)+ixyrange
      jb=jy(ty)-ixyrange
      je=jy(ty)+ixyrange
      vormax=-99999.
      do 10 j=jb,je
      do 10 i=ib,ie
        if(vor(i,j).gt.vormax)then
          ix(ty)=i
          jy(ty)=j
          vormax=vor(i,j)
        endif
 10   continue
      enddo
      return
      end
