cfk   subroutine mpe_transpose_sr(sbuf1,rbuf1,lev,n,m,nsize)
      subroutine mpe_transpose_sr(sbuf,rbuf,lev,n,m,nsize)
c
      include 'mpif.h'
      real sbuf (lev,n,m*nsize)
      real rwork(lev,n,m,nsize)
      real rbuf (lev,n,nsize,m)
c
c
      len_tr=m*n
c
ccc add on 2002/6/27
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)
c
      call MPI_ALLTOALL( SBUF,  LEN_TR*LEV, MPI_DOUBLE_PRECISION,
     &                   RWORK, LEN_TR*LEV, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, IERR )
c
      do j=1,m
      do ii=1,nsize
      do i=1,n
      do k=1,lev
        rbuf(k,i,ii,j)=rwork(k,i,j,ii)
      enddo
      enddo
      enddo
      enddo
c
      return
      end
c
