      subroutine diarfj(m,n,x,fvec,fjac,ldfjac,task,nint)
      character*(*) task
      integer n, m, ldfjac, nint
      double precision x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine diarfj
c
c     This subroutine computes the function and Jacobian matrix of the
c     data residual equations for the isomerization of alpha-pinene
c     (collocation formulation) problem.
c
c     The subroutine statement is
c
c       subroutine diarfj(m,n,x,fvec,fjac,ldfjac,task,nint,sigma)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions. m = 40.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 25*nint + 5.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a double precision array of dimension m.
c         On entry fvec need not be specified.
c         On exit fvec contains the function evaluated at x if
c            task = 'F' or 'FJ'.
c
c       fjac is a double precision array of dimension (ldfjac,n).
c         On entry fjac need not be specified.
c         On exit fjac contains the Jacobian matrix evaluated at x if
c            task = 'J' or 'FJ'.
c
c       ldfjac is an integer variable.
c          On entry ldfjac is the leading dimension of fjac.
c          On exit ldfjac is unchanged.
c
c       task is a character variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c            'F'      Evaluate the function at x.
c            'J'      Evaluate the Jacobian matrix at x.
c            'FJ'     Evaluate the function and the Jacobian at x.
c
c         On exit task is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the
c            k-stage collocation.
c         On exit nint is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer cpts, dim, maxdeg, ntimes, s
      parameter (cpts=4,maxdeg=1,ntimes=8,s=5)
      parameter (dim=maxdeg+cpts-1)
      double precision len, one, zero
      parameter (zero=0.0d0,one=1.0d0)
      parameter (len=4.0d4)

      logical feval, jeval
      integer e, i, ideg, j, k, mm, neqn, npi
      integer deg(s), eqn(s), intt(ntimes), sumdeg(0:s), var(s)
      double precision dt(ntimes,2), dw(maxdeg+1,cpts+maxdeg,s),
     +                 obs(ntimes,s), tnfhk(ntimes,0:dim,0:dim,
     +                 0:maxdeg), w(maxdeg+1,s)
      double precision h, hm, nf, tijh

      data (deg(i),i=1,s)/1, 1, 1, 1, 1/
      data (dt(i,1),i=1,ntimes)/1230.0d0, 3060.0d0, 4920.0d0, 7800.0d0,
     +     10680.0d0, 15030.0d0, 22620.0d0, 36420.0d0/
      data (obs(i,1),i=1,ntimes)/88.35d0, 76.40d0, 65.10d0, 50.40d0,
     +     37.50d0, 25.90d0, 14.00d0, 4.50d0/
      data (obs(i,2),i=1,ntimes)/7.3d0, 15.6d0, 23.1d0, 32.9d0, 42.7d0,
     +     49.1d0, 57.4d0, 63.1d0/
      data (obs(i,3),i=1,ntimes)/2.3d0, 4.5d0, 5.3d0, 6.0d0, 6.0d0,
     +     5.9d0, 5.1d0, 3.8d0/
      data (obs(i,4),i=1,ntimes)/0.4d0, 0.7d0, 1.1d0, 1.5d0, 1.9d0,
     +     2.2d0, 2.6d0, 2.9d0/
      data (obs(i,5),i=1,ntimes)/1.75d0, 2.8d0, 5.8d0, 9.3d0, 12.0d0,
     +     17.0d0, 21.0d0, 25.7d0/

c     Check input arguments for errors.

      if (m .ne. 40 .or. n .ne. 25*nint+5) then
         task = 'ERROR: M .NE. 40 OR N .NE. 25*NINT + 5 IN DIARFJ'

         return

      end if

c     Initialization.

      h = len/dble(nint)

c     Determine in which time interval each observation lies and
c     what the value of the independent variable is relative to
c     that interval.

      do 10 i = 1, ntimes
         k = int(dt(i,1)/h)
         dt(i,2) = (dt(i,1)-dble(k)*h)/h
         intt(i) = k + 1
   10 continue

c     Store all possible combinations of zeta, h, and factorial.

      hm = one
      do 50 mm = 0, maxdeg
         do 40 i = 1, ntimes
            tijh = hm
            do 30 j = 0, dim
               nf = one
               do 20 k = 0, dim
                  tnfhk(i,j,k,mm) = tijh/nf
                  nf = nf*dble(k+1)
   20          continue
               tijh = tijh*dt(i,2)
   30       continue
   40    continue
         hm = hm*h
   50 continue

      ideg = 0
      sumdeg(0) = 0
      do 60 i = 1, s
         ideg = ideg + deg(i)
         sumdeg(i) = sumdeg(i-1) + deg(i)
   60 continue
      npi = s*cpts + ideg
      neqn = nint*npi

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         feval = .true.
      else
         feval = .false.
      end if
      if (task .eq. 'J' .or. task .eq. 'FJ') then
         jeval = .true.
      else
         jeval = .false.
      end if

      if (feval) then
         do 70 i = 1, ntimes*s
            fvec(i) = zero
   70    continue
      end if
      if (jeval) then
         do 90 j = 1, neqn + s
            do 80 i = 1, ntimes*s
               fjac(i,j) = zero
   80       continue
   90    continue
      end if

c     Set up the data residual equations.

      do 170 k = 1, ntimes
         i = intt(k)
         do 130 e = 1, s
            var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
            eqn(e) = (k-1)*s + e
            do 120 mm = 1, deg(e) + 1
               w(mm,e) = zero
               do 100 j = mm, deg(e)
                  w(mm,e) = w(mm,e) + tnfhk(k,j-mm,j-mm,j-mm)*
     +                      x(var(e)+j)
                  dw(mm,j,e) = tnfhk(k,j-mm,j-mm,j-mm)
  100          continue
               do 110 j = 1, cpts
                  w(mm,e) = w(mm,e) + tnfhk(k,deg(e)+j-mm,deg(e)+j-mm,
     +                      deg(e)-mm+1)*x(var(e)+deg(e)+j)
                  dw(mm,deg(e)+j,e) = tnfhk(k,deg(e)+j-mm,deg(e)+j-mm,
     +                                deg(e)-mm+1)
  110          continue
  120       continue
  130    continue

         if (feval) then
            do 140 e = 1, s
               fvec(eqn(e)) = w(1,e) - obs(k,e)
  140       continue
         end if
         if (jeval) then
            do 160 e = 1, s
               do 150 j = 1, cpts + deg(e)
                  fjac(eqn(e),var(e)+j) = dw(1,j,e)
  150          continue
  160       continue
         end if
  170 continue

      end
