      program gsw
c
      include '../include/param.h'
      include '../include/mpe.h'
      include '../include/rank.h'
      include '../include/index.h'
      include '../include/const.h'
      include '../include/grid.h'
      include '../include/spec.h'
c
      call mpe_init(nsize,myrank)
c
      call prepare
c
      call getrdy
c
      call intgrt
c
      call mpe_finalize
c
      stop
      end
