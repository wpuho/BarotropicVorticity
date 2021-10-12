      subroutine mpe_finalize
c
      include 'mpif.h'
c
      call MPI_BARRIER(MPI_COMM_WORLD, IERR )
      call MPI_FINALIZE(IERR)
c
      return
      end
