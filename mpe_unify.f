      subroutine mpe_unify(a,n,m,idcmp,type)
c
      include 'mpif.h'
      include'../include/mpe.h'
      include'../include/rank.h'
      include'../include/index.h'
      integer n,m,idcmp,type
      dimension a(*)
c
      call MPI_COMM_SIZE( MPI_COMM_WORLD, NSIZE , IERR )
      if(idcmp .eq. 1) then
        my_max=n/nsize+1
        call mpe_unify1_r(a,n,m,my_max,jlistnum,jlist1,jlist2,nsize)
      else if(idcmp .eq. 2) then
        my_max=m/nsize+1
        if(type .eq. mpe_integer) then
          call mpe_unify2_i(a,n,m,my_max,jlistnum,jlist1,jlist2,nsize)
        else if(type .eq. mpe_double) then
          call mpe_unify2_r(a,n,m,my_max,jlistnum,jlist1,jlist2,nsize)
        else if(type .eq. mpe_logical) then
          call mpe_unify2_l(a,n,m,my_max,jlistnum,jlist1,jlist2,nsize)
        else
          write(6,*) 'mpe_unify: Argument(type) Error  RANK=',myrank
        endif
      else if(idcmp .eq. 3) then
        jtmax=m/nsize+1
        call mpe_unify3_r(a,n,m,jtmax,mlistnum,mlist,nlist,nsize)
      else
        write(6,*) 'mpe_unify: Argument(idcmp) Error  RANK=',myrank
      endif
c
      return
      end

      subroutine mpe_unify1_r(a,m,n,my_max,jlistnum,jlist1,jlist2,nsize)
c
      include 'mpif.h'
      dimension jlist1(my_max), jlist2(m)
c
      real a(m,n)
      real b1(n,my_max)
      real b2(n,my_max*nsize)
c
      do jj=1,jlistnum
        j1=jlist1(jj)
      do i=1,n
        b1(i,jj)=a(j1,i)
      enddo
      enddo
c
ccc add on 2002/6/27
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)
c
      call MPI_ALLGATHER( B1,N*MY_MAX,   MPI_DOUBLE_PRECISION,
     &                    B2,N*MY_MAX,   MPI_DOUBLE_PRECISION,
     &                     MPI_COMM_WORLD, IERR )
c
      do j=1,m
        jf=jlist2(j)
      do i=1,n
        a(j,i)=b2(i,jf)
      enddo
      enddo
c
      return
      end
c
      subroutine mpe_unify2_i(a,n,m,my_max,jlistnum,jlist1,jlist2,nsize)
c
      include 'mpif.h'
      dimension jlist1(my_max), jlist2(m)
c
      integer a(n,m)
      integer b1(n,my_max)
      integer b2(n,my_max*nsize)

      do jj=1,jlistnum
        j1=jlist1(jj)
      do i=1,n
        b1(i,jj)=a(i,j1)
      enddo
      enddo
c
ccc add on 2002/6/27
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)
c
      call MPI_ALLGATHER( B1,N*MY_MAX,   MPI_INTEGER,
     &                    B2,N*MY_MAX,   MPI_INTEGER,
     &                     MPI_COMM_WORLD, IERR )
c
      do j=1,m
        jf=jlist2(j)
      do i=1,n
        a(i,j)=b2(i,jf)
      enddo
      enddo

      return
      end

      subroutine mpe_unify2_r(a,n,m,my_max,jlistnum,jlist1,jlist2,nsize)
c
      include 'mpif.h'
      dimension jlist1(my_max), jlist2(m)
c
      real a(n,m)
      real b1(n,my_max)
      real b2(n,my_max*nsize)

      do jj=1,jlistnum
        j1=jlist1(jj)
      do i=1,n
        b1(i,jj)=a(i,j1)
      enddo
      enddo
c
ccc add on 2002/6/27
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)
c
      call MPI_ALLGATHER( B1,N*MY_MAX,   MPI_DOUBLE_PRECISION,
     &                    B2,N*MY_MAX,   MPI_DOUBLE_PRECISION,
     &                     MPI_COMM_WORLD, IERR )
c
      do j=1,m
        jf=jlist2(j)
      do i=1,n
        a(i,j)=b2(i,jf)
      enddo
      enddo

      return
      end

      subroutine mpe_unify2_l(a,n,m,my_max,jlistnum,jlist1,jlist2,nsize)
c
      include 'mpif.h'
      dimension jlist1(my_max), jlist2(m)
c
      dimension a(n,m)
      dimension b1(n,my_max)
      dimension b2(n,my_max*nsize)
      logical a, b1, b2

      do jj=1,jlistnum
        j1=jlist1(jj)
      do i=1,n
        b1(i,jj)=a(i,j1)
      enddo
      enddo
c
ccc add on 2002/6/27
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)
c
      call MPI_ALLGATHER( B1,N*MY_MAX,   MPI_LOGICAL,
     &                    B2,N*MY_MAX,   MPI_LOGICAL,
     &                     MPI_COMM_WORLD, IERR )
c
      do j=1,m
        jf=jlist2(j)
      do i=1,n
        a(i,j)=b2(i,jf)
      enddo
      enddo

      return
      end

      subroutine mpe_unify3_r(a,n,m,jtmax,mlistnum,mlist,nlist,nsize)
c
      include 'mpif.h'
      dimension mlist(jtmax), nlist(m)
c
      real a(n,m)
      real b1(n,jtmax)
      real b2(n,jtmax*nsize)

      do jj=1,mlistnum
        j1=mlist(jj)
      do i=1,n
        b1(i,jj)=a(i,j1)
      enddo
      enddo
c
ccc add on 2002/6/27
      call MPI_BARRIER(MPI_COMM_WORLD, IERR)
c
      call MPI_ALLGATHER( B1,N*JTMAX,   MPI_DOUBLE_PRECISION,
     &                    B2,N*JTMAX,   MPI_DOUBLE_PRECISION,
     &                     MPI_COMM_WORLD, IERR )
c
      do j=1,m
        jf=nlist(j)
      do i=1,n
        a(i,j)=b2(i,jf)
      enddo
      enddo

      return
      end
