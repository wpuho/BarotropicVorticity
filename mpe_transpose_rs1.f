      subroutine mpe_transpose_rs1(sbuf,rbuf,n,m,lev,nsize)
c
      include 'mpif.h'
      real sbuf (n,nsize,m,lev), rbuf (n,m*nsize,lev)
      real swork(lev,n,m,nsize), rwork(lev,n,m*nsize)
c
      len_tr=n*m

      do j=1,m
      do ii=1,nsize
      do i=1,n
      do k=1,lev
        swork(k,i,j,ii)=sbuf(i,ii,j,k)
      enddo
      enddo
      enddo
      enddo
c
ccc add on 2002/6/27
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)
c
      call MPI_ALLTOALL( SWORK, LEN_TR*LEV, MPI_DOUBLE_PRECISION,
     &                   RWORK, LEN_TR*LEV, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, IERR )
c
      do k=1,lev
      do j=1,m*nsize
      do i=1,n
        rbuf(i,j,k)=rwork(k,i,j)
      enddo
      enddo
      enddo
c
      return
      end
