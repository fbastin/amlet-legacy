      subroutine ssfdfj(n,x,fvec,fjac,ldfjac,task,eps,nint)
      character*(*) task
      integer n, ldfjac, nint
      real eps
      real x(n), fvec(n), fjac(ldfjac,n)
c     **********
c
c     Subroutine ssfdfj
c
c     This subroutine computes the function and Jacobian matrix of the
c     swirling flow between disks problem.
c
c     The subroutine statement is
c
c       subroutine ssfdfj(n,x,fvec,fjac,ldfjac,task,eps,nint)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 14*nint.
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
c       eps is a real variable.
c         On entry eps is the viscosity of the fluid.
c         On exit eps is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the k-stage
c            collocation.
c         On exit nint is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer bc, cpts, dim, fdeg, gdeg, mdeg, npi
      parameter (bc=3,cpts=4,fdeg=4,gdeg=2,mdeg=4)
      parameter (dim=mdeg+cpts-1,npi=2*cpts+gdeg+fdeg)
      real omega1, omega2, one, zero
      parameter (zero=0.0,one=1.0)
      parameter (omega1=-1.0,omega2=1.0)

      logical feval, jeval
      integer eqn1, eqn2, i, j, k, m, var1, var2
      real h, hm, nf, rhoijh, xt
      real dwf(fdeg+1,cpts+fdeg), dwg(gdeg+1,cpts+gdeg),
     +     rhnfhk(cpts,0:dim,0:dim,0:mdeg), rho(cpts), wg(gdeg+1),
     +     wf(fdeg+1)

      data (rho(i),i=1,cpts)/0.694318413734436035E-1,
     +     0.330009490251541138, 0.669990539550781250,
     +     0.930568158626556396/

c     Check input arguments for errors.

      if (n .ne. 14*nint) then
         task = 'ERROR: N .NE. 14*NINT IN DSFDFJ'

         return

      end if

c     Initialization.

      h = one/real(nint)

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         do 10 i = 1, n
            x(i) = zero
   10    continue

c        The standard starting point corresponds to the solution
c        of the swirling flow problem with infinite viscosity.

         xt = zero
         do 20 i = 1, nint
            var1 = (i-1)*npi + fdeg + cpts
            x(var1+1) = omega1 + (omega2-omega1)*xt
            x(var1+2) = omega2 - omega1
            xt = xt + h
   20    continue

         return

      end if

c     Store all possible combinations of rho, h, and n factorial.

      hm = one
      do 60 m = 0, mdeg
         do 50 i = 1, cpts
            rhoijh = hm
            do 40 j = 0, dim
               nf = one
               do 30 k = 0, dim
                  rhnfhk(i,j,k,m) = rhoijh/nf
                  nf = nf*real(k+1)
   30          continue
               rhoijh = rhoijh*rho(i)
   40       continue
   50    continue
         hm = hm*h
   60 continue

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

      if (feval) then
         do 70 j = 1, n
            fvec(j) = zero
   70    continue
      end if
      if (jeval) then
         do 90 j = 1, n
            do 80 i = 1, n
               fjac(i,j) = zero
   80       continue
   90    continue
      end if
      do 110 k = 1, cpts + fdeg
         do 100 j = 1, fdeg + 1
            dwf(j,k) = zero
  100    continue
  110 continue
      do 130 k = 1, cpts + gdeg
         do 120 j = 1, gdeg + 1
            dwg(j,k) = zero
  120    continue
  130 continue

c     Set up the boundary equations at t = 0.
c     f(0) = 0, f'(0) = 0, g(0) = omega1.

      if (feval) then
         fvec(1) = x(1)
         fvec(2) = x(2)
         fvec(3) = x(cpts+fdeg+1) - omega1
      end if
      if (jeval) then
         fjac(1,1) = one
         fjac(2,2) = one
         fjac(3,cpts+fdeg+1) = one
      end if

c     Set up the collocation equations.

      do 230 i = 1, nint
         var1 = (i-1)*npi
         eqn1 = var1 + bc
         var2 = var1 + cpts + fdeg
         eqn2 = eqn1 + cpts
         do 220 k = 1, cpts
            do 160 m = 1, fdeg + 1
               wf(m) = zero
               do 140 j = m, fdeg
                  wf(m) = wf(m) + rhnfhk(k,j-m,j-m,j-m)*x(var1+j)
                  dwf(m,j) = rhnfhk(k,j-m,j-m,j-m)
  140          continue
               do 150 j = 1, cpts
                  wf(m) = wf(m) + x(var1+fdeg+j)*
     +                    rhnfhk(k,fdeg+j-m,fdeg+j-m,fdeg-m+1)
                  dwf(m,fdeg+j) = rhnfhk(k,fdeg+j-m,fdeg+j-m,fdeg-m+1)
  150          continue
  160       continue
            do 190 m = 1, gdeg + 1
               wg(m) = zero
               do 170 j = m, gdeg
                  wg(m) = wg(m) + rhnfhk(k,j-m,j-m,j-m)*x(var2+j)
                  dwg(m,j) = rhnfhk(k,j-m,j-m,j-m)
  170          continue
               do 180 j = 1, cpts
                  wg(m) = wg(m) + x(var2+gdeg+j)*
     +                    rhnfhk(k,gdeg+j-m,gdeg+j-m,gdeg-m+1)
                  dwg(m,gdeg+j) = rhnfhk(k,gdeg+j-m,gdeg+j-m,gdeg-m+1)
  180          continue
  190       continue
            if (feval) then
               fvec(eqn1+k) = eps*wf(5) + wf(4)*wf(1) + wg(2)*wg(1)
               fvec(eqn2+k) = eps*wg(3) + wf(1)*wg(2) - wf(2)*wg(1)
            end if
            if (jeval) then
               do 200 j = 1, cpts + fdeg
                  fjac(eqn1+k,var1+j) = eps*dwf(5,j) + dwf(4,j)*wf(1) +
     +                                  wf(4)*dwf(1,j)
                  fjac(eqn2+k,var1+j) = dwf(1,j)*wg(2) - dwf(2,j)*wg(1)
  200          continue
               do 210 j = 1, cpts + gdeg
                  fjac(eqn1+k,var2+j) = dwg(2,j)*wg(1) + wg(2)*dwg(1,j)
                  fjac(eqn2+k,var2+j) = eps*dwg(3,j) + wf(1)*dwg(2,j) -
     +                                  wf(2)*dwg(1,j)
  210          continue
            end if

  220    continue
  230 continue

c     Set up the continuity equations.

      do 340 i = 1, nint - 1
         var1 = (i-1)*npi
         eqn1 = var1 + bc + 2*cpts
         var2 = var1 + fdeg + cpts
         eqn2 = eqn1 + fdeg
         do 260 m = 1, fdeg
            wf(m) = zero
            do 240 j = m, fdeg
               wf(m) = wf(m) + rhnfhk(1,0,j-m,j-m)*x(var1+j)
               dwf(m,j) = rhnfhk(1,0,j-m,j-m)
  240       continue
            do 250 j = 1, cpts
               wf(m) = wf(m) + rhnfhk(1,0,fdeg+j-m,fdeg-m+1)*
     +                 x(var1+fdeg+j)
               dwf(m,fdeg+j) = rhnfhk(1,0,fdeg+j-m,fdeg-m+1)
  250       continue
  260    continue
         do 290 m = 1, gdeg
            wg(m) = zero
            do 270 j = m, gdeg
               wg(m) = wg(m) + rhnfhk(1,0,j-m,j-m)*x(var2+j)
               dwg(m,j) = rhnfhk(1,0,j-m,j-m)
  270       continue
            do 280 j = 1, cpts
               wg(m) = wg(m) + rhnfhk(1,0,gdeg+j-m,gdeg-m+1)*
     +                 x(var2+gdeg+j)
               dwg(m,gdeg+j) = rhnfhk(1,0,gdeg+j-m,gdeg-m+1)
  280       continue
  290    continue
         do 310 m = 1, fdeg
            if (feval) then
               fvec(eqn1+m) = x(var1+npi+m) - wf(m)
            end if
            if (jeval) then
               fjac(eqn1+m,var1+npi+m) = one
               do 300 j = 1, cpts + fdeg
                  fjac(eqn1+m,var1+j) = -dwf(m,j)
  300          continue
            end if
  310    continue
         do 330 m = 1, gdeg
            if (feval) fvec(eqn2+m) = x(var2+npi+m) - wg(m)
            if (jeval) then
               fjac(eqn2+m,var2+npi+m) = one
               do 320 j = 1, cpts + gdeg
                  fjac(eqn2+m,var2+j) = -dwg(m,j)
  320          continue
            end if
  330    continue
  340 continue

c     Prepare for setting up the boundary conditions at t = 1.

      var1 = n - npi
      do 370 m = 1, fdeg + 1
         wf(m) = zero
         do 350 j = m, fdeg
            wf(m) = wf(m) + rhnfhk(1,0,j-m,j-m)*x(var1+j)
            dwf(m,j) = rhnfhk(1,0,j-m,j-m)
  350    continue
         do 360 j = 1, cpts
            wf(m) = wf(m) + rhnfhk(1,0,fdeg+j-m,fdeg-m+1)*x(var1+fdeg+j)
            dwf(m,fdeg+j) = rhnfhk(1,0,fdeg+j-m,fdeg-m+1)
  360    continue
  370 continue
      var2 = var1 + fdeg + cpts
      do 400 m = 1, gdeg + 1
         wg(m) = zero
         do 380 j = m, gdeg
            wg(m) = wg(m) + rhnfhk(1,0,j-m,j-m)*x(var2+j)
            dwg(m,j) = rhnfhk(1,0,j-m,j-m)
  380    continue
         do 390 j = 1, cpts
            wg(m) = wg(m) + rhnfhk(1,0,gdeg+j-m,gdeg-m+1)*x(var2+gdeg+j)
            dwg(m,gdeg+j) = rhnfhk(1,0,gdeg+j-m,gdeg-m+1)
  390    continue
  400 continue

c     Set up the boundary equations at t = 1.
c     f(1) = 0, f'(1) = 0, g(1) = omega2.

      if (feval) then
         fvec(n-2) = wf(1)
         fvec(n-1) = wf(2)
         fvec(n) = wg(1) - omega2
      end if
      if (jeval) then
         do 410 j = 1, cpts + fdeg
            fjac(n-2,var1+j) = dwf(1,j)
            fjac(n-1,var1+j) = dwf(2,j)
  410    continue
         do 420 j = 1, cpts + gdeg
            fjac(n,var2+j) = dwg(1,j)
  420    continue
      end if

      end
