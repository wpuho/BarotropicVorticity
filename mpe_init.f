      subroutine mpe_init(nsize, myrank)
c
      include 'mpif.h'
c
      call MPI_INIT( IERR )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, NSIZE , IERR )
      call MPI_COMM_RANK( MPI_COMM_WORLD, MYRANK, IERR )
c
      return
      end
