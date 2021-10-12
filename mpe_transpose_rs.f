      subroutine mpe_transpose_rs(sbuf,rbuf,lev,n,m,nsize)
c
      include 'mpif.h'
      real sbuf (lev,n,nsize,m)
      real swork(lev,n,m,nsize)
      real rbuf (lev,n,m*nsize)
c
      len_tr=m*n
c
      do j=1,m
      do ii=1,nsize
      do i=1,n
      do k=1,lev
        swork(k,i,j,ii)=sbuf(k,i,ii,j)
      enddo
      enddo
      enddo
      enddo
c
ccc add on 2002/6/27
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)
c
      call MPI_ALLTOALL( SWORK, LEN_TR*LEV, MPI_DOUBLE_PRECISION,
     &                   RBUF , LEN_TR*LEV, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, IERR )
c
      return
      end
