      subroutine siaofj(m,n,x,fvec,fjac,ldfjac,task,nint)
      character*(*) task
      integer n, m, ldfjac, nint
      real x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine siaofj
c
c     This subroutine computes the function and Jacobian matrix of the
c     ordinary differential constraint equations for the isomerization
c     of alpha-pinene (collocation formulation) problem.
c
c     The subroutine statement is
c
c       subroutine siaofj(m,n,x,fvec,fjac,ldfjac,task,nint)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions. m = 25*nint.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 25*nint + 5.
c         On exit n is unchanged.
c
c       x is a real array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a real array of dimension m.
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
c            'F'      Evaluate the function at x.
c            'J'      Evaluate the Jacobian matrix at x.
c            'FJ'     Evaluate the function and the Jacobian at x.
c            'XS'     Set x to the standard starting point xs.
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
      integer cpts, dim, maxdeg, s
      parameter (cpts=4,maxdeg=1,s=5)
      parameter (dim=maxdeg+cpts-1)
      real len, one, zero
      parameter (zero=0.0,one=1.0)
      parameter (len=4.0E4)

      logical feval, jeval
      integer e, i, ideg, j, k, l, mm, neqn, npi
      integer deg(s), eqn(s), sumdeg(0:s), var(s)
      real dw(maxdeg+1,cpts+maxdeg,s), icond(s),
     +     rhnfhk(cpts,0:dim,0:dim,0:maxdeg), rho(cpts), w(maxdeg+1,s)
      real h, hm, nf, rhoijh

      data (deg(i),i=1,s)/1, 1, 1, 1, 1/
      data (icond(i),i=1,s)/1.0E2, 0.0, 0.0, 0.0, 0.0/
      data (rho(i),i=1,cpts)/0.694318413734436035E-1,
     +     0.330009490251541138, 0.669990539550781250,
     +     0.930568158626556396/

c     Check input arguments for errors.

      if (m .ne. 25*nint .or. n .ne. 25*nint+5) then
         task = 'ERROR: M .NE. 25*NINT OR N .NE. 25*NINT + 5 IN DIAOFJ'

         return

      end if

c     Initialization.

      h = len/real(nint)

c     Store all possible combinations of rho, h, and n factorial.

      hm = one
      do 40 mm = 0, maxdeg
         do 30 i = 1, cpts
            rhoijh = hm
            do 20 j = 0, dim
               nf = one
               do 10 k = 0, dim
                  rhnfhk(i,j,k,mm) = rhoijh/nf
                  nf = nf*real(k+1)
   10          continue
               rhoijh = rhoijh*rho(i)
   20       continue
   30    continue
         hm = hm*h
   40 continue

      ideg = 0
      sumdeg(0) = 0
      do 50 i = 1, s
         ideg = ideg + deg(i)
         sumdeg(i) = sumdeg(i-1) + deg(i)
   50 continue
      npi = s*cpts + ideg
      neqn = nint*npi

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = icond(1)
         do 60 i = 2, neqn
            x(i) = zero
   60    continue
         x(neqn+1) = 5.84E-5
         x(neqn+2) = 2.65E-5
         x(neqn+3) = 1.63E-5
         x(neqn+4) = 2.777E-4
         x(neqn+5) = 4.61E-5

         return

      end if

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
         do 70 i = 1, neqn
            fvec(i) = zero
   70    continue
      end if
      if (jeval) then
         do 90 j = 1, neqn + s
            do 80 i = 1, neqn
               fjac(i,j) = zero
   80       continue
   90    continue
      end if

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

      do 240 i = 1, nint
         do 230 k = 1, cpts
            do 170 e = 1, s
               var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
               eqn(e) = s + (i-1)*npi + (e-1)*cpts
               do 160 mm = 1, deg(e) + 1
                  w(mm,e) = zero
                  do 140 j = mm, deg(e)
                     w(mm,e) = w(mm,e) + x(var(e)+j)*
     +                         rhnfhk(k,j-mm,j-mm,j-mm)
                     dw(mm,j,e) = rhnfhk(k,j-mm,j-mm,j-mm)
  140             continue
                  do 150 j = 1, cpts
                     w(mm,e) = w(mm,e) + x(var(e)+deg(e)+j)*
     +                         rhnfhk(k,deg(e)+j-mm,deg(e)+j-mm,
     +                         deg(e)-mm+1)
                     dw(mm,deg(e)+j,e) = rhnfhk(k,deg(e)+j-mm,
     +                                   deg(e)+j-mm,deg(e)-mm+1)
  150             continue
  160          continue
  170       continue

            if (feval) then
               fvec(eqn(1)+k) = w(2,1) + (x(neqn+1)+x(neqn+2))*w(1,1)
               fvec(eqn(2)+k) = w(2,2) - x(neqn+1)*w(1,1)
               fvec(eqn(3)+k) = w(2,3) - x(neqn+2)*w(1,1) +
     +                          (x(neqn+3)+x(neqn+4))*w(1,3) -
     +                          x(neqn+5)*w(1,5)
               fvec(eqn(4)+k) = w(2,4) - x(neqn+3)*w(1,3)
               fvec(eqn(5)+k) = w(2,5) - x(neqn+4)*w(1,3) +
     +                          x(neqn+5)*w(1,5)
            end if

            if (jeval) then
               do 180 j = 1, cpts + deg(1)
                  fjac(eqn(1)+k,var(1)+j) = dw(2,j,1) +
     +                                      (x(neqn+1)+x(neqn+2))*
     +                                      dw(1,j,1)
                  fjac(eqn(2)+k,var(1)+j) = -x(neqn+1)*dw(1,j,1)
                  fjac(eqn(3)+k,var(1)+j) = -x(neqn+2)*dw(1,j,1)
  180          continue
               do 190 j = 1, cpts + deg(2)
                  fjac(eqn(2)+k,var(2)+j) = dw(2,j,2)
  190          continue
               do 200 j = 1, cpts + deg(3)
                  fjac(eqn(3)+k,var(3)+j) = dw(2,j,3) +
     +                                      (x(neqn+3)+x(neqn+4))*
     +                                      dw(1,j,3)
                  fjac(eqn(4)+k,var(3)+j) = -x(neqn+3)*dw(1,j,3)
                  fjac(eqn(5)+k,var(3)+j) = -x(neqn+4)*dw(1,j,3)
  200          continue
               do 210 j = 1, cpts + deg(4)
                  fjac(eqn(4)+k,var(4)+j) = dw(2,j,4)
  210          continue
               do 220 j = 1, cpts + deg(5)
                  fjac(eqn(3)+k,var(5)+j) = -x(neqn+5)*dw(1,j,5)
                  fjac(eqn(5)+k,var(5)+j) = dw(2,j,5) +
     +                                      x(neqn+5)*dw(1,j,5)
  220          continue
               fjac(eqn(1)+k,neqn+1) = w(1,1)
               fjac(eqn(2)+k,neqn+1) = -w(1,1)
               fjac(eqn(1)+k,neqn+2) = w(1,1)
               fjac(eqn(3)+k,neqn+2) = -w(1,1)
               fjac(eqn(3)+k,neqn+3) = w(1,3)
               fjac(eqn(4)+k,neqn+3) = -w(1,3)
               fjac(eqn(3)+k,neqn+4) = w(1,3)
               fjac(eqn(5)+k,neqn+4) = -w(1,3)
               fjac(eqn(3)+k,neqn+5) = -w(1,5)
               fjac(eqn(5)+k,neqn+5) = w(1,5)
            end if

  230    continue
  240 continue

c     Set up the continuity equations.

      do 320 i = 1, nint - 1

         do 280 e = 1, s
            var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
            eqn(e) = s + (i-1)*npi + s*cpts + sumdeg(e-1)
            do 270 mm = 1, deg(e)
               w(mm,e) = zero
               do 250 j = mm, deg(e)
                  w(mm,e) = w(mm,e) + rhnfhk(1,0,j-mm,j-mm)*x(var(e)+j)
                  dw(mm,j,e) = rhnfhk(1,0,j-mm,j-mm)
  250          continue
               do 260 j = 1, cpts
                  w(mm,e) = w(mm,e) + x(var(e)+deg(e)+j)*
     +                      rhnfhk(1,0,deg(e)+j-mm,deg(e)-mm+1)
                  dw(mm,deg(e)+j,e) = rhnfhk(1,0,deg(e)+j-mm,
     +                                deg(e)-mm+1)
  260          continue
  270       continue
  280    continue

         do 310 e = 1, s
            do 300 mm = 1, deg(e)
               if (feval) fvec(eqn(e)+mm) = x(var(e)+npi+mm) - w(mm,e)
               if (jeval) then
                  fjac(eqn(e)+mm,var(e)+npi+mm) = one
                  do 290 j = 1, cpts + deg(e)
                     fjac(eqn(e)+mm,var(e)+j) = -dw(mm,j,e)
  290             continue
               end if
  300       continue
  310    continue

  320 continue

      end
