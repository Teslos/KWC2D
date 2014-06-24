
      subroutine slapsol(ae,aw,an,as,ap,b,u,ni,nj,maxits,method)
c      implicit none
      include 'icslap.h'
      double precision ae(nelt),aw(nelt),an(nelt),as(nelt),ap(nelt),
     $ u(nelt)
      double precision a(nelt),af(nelt),b(nelt),s(nelt),
     $ x(nelt),r(nelt),p(nelt),tol,err
      integer ia(nelt),ja(nelt),li(nelt),itol,iunit,
     $ maxits
      integer i,j,k,ni,nj,ierr,ij,lda,n,nim,njm,nix,njy,its
      character(80) method
c
c     this is a test main program for the incomplete factorization conjugate
c     gradient algorithm which uses an incomplete lu decomposition.
c
      err = 0.0d0
      ierr = 0
      lda = 10 
      nix = ni-2
      njy = nj-2
      nim = ni-1
      njm = nj-1
      n = nix*njy
      iunit = 0
      itol = 2
      iunit = 6
      nel = 0
c
c     working array for conversion of the indices
c
      do 5 i=1,ni
        li(i) = (i-1)*nj
    5 continue
  
c
c     initialize the structure for the matrix.
c
      write(6,*) 'Initialize matrix'
      call inits(a,n,ia,ja)
c
c     setup the matrix and the right hand side.
c
      k = 0
c      write(6,*) 'Setup the matrix'
c      write(*,*) 'Size of the matrix',nix,njy
      do  j = 2,njm
      do  i = 2,nim
         ij = li(i) + j
c         write(*,*) 'i=',i,'j=',j
         k = k + 1 
c         write(*,*) 'k=',k
         if (k-1 .ge. 1) then
            call put(aw(ij),a,n,nel,ia,ja,k,k-1)
         end if
         if (k+1 .le. n) then
            call put(ae(ij),a,n,nel,ia,ja,k,k+1)
         end if   
         call put(ap(ij),a,n,nel,ia,ja,k,k)
         if( k-nix .ge. 1 ) then 
            call put(as(ij),a,n,nel,ia,ja,k,k-nix)
         end if   
         if( k+nix .le. n ) then 
            call put(an(ij),a,n,nel,ia,ja,k,k+nix)
         end if
c         b(i) = 6.0d0
      end do
      end do      
c      write(*,*) 'Length of the row=',ni
c      call put(1.0d0,a,n,nel,ia,ja,1,2)
c      call put(-1.0d0,a,n,nel,ia,ja,1,1+nix)
c      call put(2.0d0,a,n,nel,ia,ja,1,n)
c      call put(2.0d0,a,n,nel,ia,ja,n,1)
c      b(1) = 5.0d0
c      b(n) = 5.0d0
c      call put(4.0d0,a,n,nelt,ia,ja,1,1)
c
c     estimate the solution
c
      x(1) = 1.0d0
      do 20 i = 2,n
         x(i) = u(i)
   20 continue
c  
c     solve the system.
c
      maxits = 100
      tol = 0.00001d0
c      write(6,*) 'CALL iccglu'
      write(6,1050)(i,ia(i),ja(i),a(i),i=1,nel)
      CALL DS2Y( N, NEL, IA, JA, A, ISYM )
c      WRITE(6,1040) (K,IA(K),JA(K),A(K),K=1,NEL)
      call dcpplt(n,nel,ia,ja,a,isym,6)
c      call iccglu(a,n,ia,ja,af,b,x,r,p,s,tol,maxits,its,info)
c      write(*,*)'n,b,x,nel,ia,ja,a,isym,itol,tol,maxits',
c     $ n,(b(i),x(i),i=1,n),nel,(ia(i),ja(i),i=1,nelt),isym
c     $ itol,tol,maxits 
      call dslucs(n,b,x,nel,ia,ja,a,isym,itol,tol,maxits,
     $ its,err,ierr,iunit,rwork,
     $ lenw,iwork,leniw)
c      write(6,*) 'Error estimated',err
c      write(6,*) 'Error returned',ierr
c      write(6,30) its
   30 format(' number of iterations taken',i5)
c      write(6,40)(x(i),i=1,n)
   40 format('x solution',5e16.8)
 1040 FORMAT(/'  ***** SLAP Column Matrix *****'/
     $        ' Indx   ia   ja     a'/(1X,I4,1X,I4,1X,I4,1X,E16.7))
 1050 FORMAT(/'  ***** SLAP Triad Matrix *****'/
     $        ' Indx   ia   ja     a'/(1X,I4,1X,I4,1X,I4,1X,E16.7))

c      stop
      end subroutine
 
      subroutine insert(t,a,n,ia,ja,i,j,k)
c
      integer n,ia(1),ja(1),i,j,k
      integer ip1,np1
      real t,a(1)
c
c     this routine rearranges the elements in arrays a and ja
c     and updates array ia for the new element.
c
c      write(6,1000)i,j,ia(i),ia(i+1),t
1000  format('  *** from insert i,j,ia(i),ia(i+1),t ',4i5,1x,e15.8)
      l1 = k
      l2 = ia(n+1) - 1
      do 10 lb = l1,l2
         l = l2 - lb + l1
         a(l+1) = a(l)
         ja(l+1) = ja(l)
   10 continue
      a(k) = t
      ja(k) = j
      ip1 = i + 1
      np1 = n + 1
      do 20 l = ip1,np1
         ia(l) = ia(l) + 1
   20 continue
      return
      end
  
      subroutine inits(a,n,ia,ja)
c
      integer n,ia(1),ja(1)
      double precision a(1)
c
c     this routine initializes storage for the sparse
c     matrix arrays.
c     note: the matrix is initialized to have zeroes
c     on the diagonal.
c
      do 10 i = 1,n
         a(i) = 0.0e0
         ja(i) = i
         ia(i) = i
   10 continue
      ia(n+1) = n+1
      return
      end

      logical function locate(a,n,ia,ja,i,j,k)
c
      integer n,ia(1),ja(1),i,j
      real a(1)
c
c     this routine will locate the i,j-th element in the
c     sparse matrix structure.
c
      is = ia(i)
      ie = ia(i+1) - 1
c
c     row i is from is to ie in array a.
c
      do 10 l = is,ie
         if( j .gt. ja(l) ) go to 10
         if( j .ne. ja(l) ) go to 5
            locate = .true.
            k = l
            go to 20
    5    continue
         if( j .ge. ja(l) ) go to 10
            locate = .false.
            k = 0
         go to 20
   10 continue
c
c     get here if should be at the end of a row.
c
      locate = .false.
      k = 0
   20 continue
      return
      end
      subroutine put(t,a,n,nelt,ia,ja,i,j)
c
      integer nelt,n,ia(100),ja(100),i,j
      double precision t,a(100)
      integer irow,icol
c
c     this routine will insert an element into the sparse matrix
c     structure.
c
      nelt = nelt + 1
      ia(nelt) = i
      ja(nelt) = j
      a(nelt) = t
c      write(6,100)ia(nelt),ja(nelt),a(nelt)
  100 format('  *** from put ia(i),ja(j),a ',2i7,e15.8)
      return
      end
 
      subroutine scopy(n,sx,incx,sy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
