      subroutine sierfj(n,x,fvec,fjac,ldfjac,task,a,b,c,nint)
      character*(*) task
      integer n, ldfjac, nint
      real a, b, c
      real x(n), fvec(n), fjac(ldfjac,n)
c     **********
c
c     Subroutine sierfj
c
c     This subroutine computes the function and Jacobian matrix of the
c     incompressible elastic rod problem.
c
c     The subroutine statement is
c
c       subroutine sierfj(n,x,fvec,fjac,ldfjac,task,a,b,c,nint)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 15*nint+3.
c         On exit n is unchanged.
c
c       x is a real array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a real array of dimension n.
c         On entry fvec need not be specified.
c         On exit fvec contains the function evaluated at x if
c            task = 'F' or 'FJ'.
c
c       fjac is a real array of dimension (ldfjac,n).
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
c             'F'     Evaluate the function at x.
c             'J'     Evaluate the Jacobian matrix at x.
c             'FJ'    Evaluate the function and the Jacobian at x.
c             'XS'    Set x to the standard starting point xs.
c
c         On exit task is unchanged.
c
c       a is a real variable.
c         On entry a is the X-coordinate of the right end of the rod.
c         On exit a is unchanged.
c
c       b is a real variable.
c         On entry b is the Y-coordinate of the right end of the rod.
c         On exit b is unchanged.
c
c       c is a real variable.
c         On entry c is the local angle at the right end of the rod.
c         On exit c is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the
c            k-stage collocation.
c         On exit nint is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and G. L. Xue.
c
c     **********
      integer cpts, dim, maxdeg, s
      parameter (cpts=4,maxdeg=1,s=3)
      parameter (dim=maxdeg+cpts-1)
      real len, one, zero
      parameter (zero=0.0,one=1.0)
      parameter (len=one)

      logical feval, jeval
      integer e, i, ideg, ivar, j, k, l, m, npi
      integer deg(s), eqn(s), sumdeg(0:s), var(s)
      real dw(maxdeg+1,cpts+maxdeg,s), icond(s),
     +     rhnfhk(cpts,0:dim,0:dim,0:maxdeg), rho(cpts), w(maxdeg+1,s)
      real h, hm, nf, rhoijh, xt

      data (deg(i),i=1,s)/1, 1, 1/
      data (icond(i),i=1,s)/zero, zero, zero/
      data (rho(i),i=1,cpts)/0.694318413734436035E-1,
     +     0.330009490251541138, 0.669990539550781250,
     +     0.930568158626556396/

c     Check input arguments for errors.

      if (n .ne. 15*nint+3) then
         task = 'ERROR: N .NE. 15*NINT + 3 IN DIERFJ'

         return

      end if

c     Initialization.

      h = len/real(nint)

      ideg = 0
      sumdeg(0) = 0
      do 10 i = 1, s
         ideg = ideg + deg(i)
         sumdeg(i) = sumdeg(i-1) + deg(i)
   10 continue
      npi = s*cpts + ideg

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         do 20 i = 1, n
            x(i) = zero
   20    continue
         xt = zero
         ivar = 0
         do 30 i = 1, nint
            x(ivar+1) = xt
            x(ivar+2) = one
            xt = xt + h
            ivar = ivar + npi
   30    continue

         return

      end if

c     Store all possible combinations of rho, h, and n factorial.

      hm = one
      do 70 m = 0, maxdeg
         do 60 i = 1, cpts
            rhoijh = hm
            do 50 j = 0, dim
               nf = one
               do 40 k = 0, dim
                  rhnfhk(i,j,k,m) = rhoijh/nf
                  nf = nf*real(k+1)
   40          continue
               rhoijh = rhoijh*rho(i)
   50       continue
   60    continue
         hm = hm*h
   70 continue

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

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

c     Initialize arrays.

      do 90 j = 1, n
         if (feval) fvec(j) = zero
         if (jeval) then
            do 80 i = 1, n
               fjac(i,j) = zero
   80       continue
         end if
   90 continue

      do 120 k = 1, s
         do 110 j = 1, maxdeg + 1
            w(j,k) = zero
            do 100 l = 1, cpts + maxdeg
               dw(j,l,k) = zero
  100       continue
  110    continue
  120 continue

c     Satisfy initial conditions at t = 0.  yi(0) = icond(i).

      do 130 i = 1, s
         if (feval) fvec(i) = x((i-1)*cpts+sumdeg(i-1)+1) - icond(i)
         if (jeval) fjac(i,(i-1)*cpts+sumdeg(i-1)+1) = one
  130 continue

c     Set up the collocation equations.

      do 200 i = 1, nint
         do 190 k = 1, cpts
            do 170 e = 1, s
               var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
               eqn(e) = s + (i-1)*npi + (e-1)*cpts
               do 160 m = 1, deg(e) + 1
                  w(m,e) = zero
                  do 140 j = m, deg(e)
                     w(m,e) = w(m,e) + x(var(e)+j)*rhnfhk(k,j-m,j-m,j-m)
                     dw(m,j,e) = rhnfhk(k,j-m,j-m,j-m)
  140             continue
                  do 150 j = 1, cpts
                     w(m,e) = w(m,e) + x(var(e)+deg(e)+j)*
     +                        rhnfhk(k,deg(e)+j-m,deg(e)+j-m,deg(e)-m+1)
                     dw(m,deg(e)+j,e) = rhnfhk(k,deg(e)+j-m,deg(e)+j-m,
     +                                  deg(e)-m+1)
  150             continue
  160          continue
  170       continue

            if (feval) then
               fvec(eqn(1)+k) = w(2,1) - cos(w(1,3))
               fvec(eqn(2)+k) = w(2,2) - sin(w(1,3))
               fvec(eqn(3)+k) = w(2,3) - x(n-2)*w(1,1) + x(n-1)*w(1,2) -
     +                          x(n)
            end if

            if (jeval) then
               do 180 j = 1, cpts + deg(1)
                  fjac(eqn(1)+k,var(1)+j) = dw(2,j,1)
                  fjac(eqn(3)+k,var(1)+j) = -x(n-2)*dw(1,j,1)
                  fjac(eqn(2)+k,var(2)+j) = dw(2,j,2)
                  fjac(eqn(3)+k,var(2)+j) = x(n-1)*dw(1,j,2)
                  fjac(eqn(1)+k,var(3)+j) = sin(w(1,3))*dw(1,j,3)
                  fjac(eqn(2)+k,var(3)+j) = -cos(w(1,3))*dw(1,j,3)
                  fjac(eqn(3)+k,var(3)+j) = dw(2,j,3)
  180          continue
               fjac(eqn(3)+k,n-2) = -w(1,1)
               fjac(eqn(3)+k,n-1) = w(1,2)
               fjac(eqn(3)+k,n) = -one
            end if

  190    continue
  200 continue

c     Set up the continuity equations.

      do 280 i = 1, nint

         do 240 e = 1, s
            var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
            eqn(e) = s + (i-1)*npi + s*cpts + sumdeg(e-1)
            do 230 m = 1, deg(e)
               w(m,e) = zero
               do 210 j = m, deg(e)
                  w(m,e) = w(m,e) + rhnfhk(1,0,j-m,j-m)*x(var(e)+j)
                  dw(m,j,e) = rhnfhk(1,0,j-m,j-m)
  210          continue
               do 220 j = 1, cpts
                  w(m,e) = w(m,e) + x(var(e)+deg(e)+j)*
     +                     rhnfhk(1,0,deg(e)+j-m,deg(e)-m+1)
                  dw(m,deg(e)+j,e) = rhnfhk(1,0,deg(e)+j-m,deg(e)-m+1)
  220          continue
  230       continue
  240    continue
         if (i .eq. nint) go to 290

         do 270 e = 1, s
            do 260 m = 1, deg(e)
               if (feval) fvec(eqn(e)+m) = x(var(e)+npi+m) - w(m,e)
               if (jeval) then
                  fjac(eqn(e)+m,var(e)+npi+m) = one
                  do 250 j = 1, cpts + deg(e)
                     fjac(eqn(e)+m,var(e)+j) = -dw(m,j,e)
  250             continue
               end if
  260       continue
  270    continue

  280 continue
  290 continue

      if (feval) then
         fvec(n-2) = w(1,1) - a
         fvec(n-1) = w(1,2) - b
         fvec(n) = w(1,3) - c
      end if

      if (jeval) then
         var(1) = n - npi - 3
         var(2) = var(1) + 5
         var(3) = var(2) + 5
         do 300 j = 1, 5
            fjac(n-2,var(1)+j) = dw(1,j,1)
            fjac(n-1,var(2)+j) = dw(1,j,2)
            fjac(n,var(3)+j) = dw(1,j,3)
  300    continue
      end if

      end
