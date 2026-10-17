      subroutine scpffj(n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer n, ldfjac
      real x(n), fvec(n), fjac(ldfjac,n)
c     **********
c
c     Subroutine scpffj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the combustion of propane (full formulation) problem.
c
c     The subroutine statement is
c
c       subroutine scpffj(n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 11.
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
c             'XL'    Set x to the lower bound xl.
c
c         On exit task is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      real one, p, p5, rr, ten, two, zero
      parameter (zero=0.0,p5=0.5,one=1.0,two=2.0,ten=1.0E1)
      parameter (p=4.0E1,rr=ten)

      integer i, j
      real pdx, sqpdx, sqrtp, xfrac, xtau
      real k(10)

      data k/zero, zero, zero, zero, 1.930E-1, 2.597E-3, 3.448E-3,
     +     1.799E-5, 2.155E-4, 3.846E-5/

c     Check input arguments for errors.

      if (n .ne. 11) then
         task = 'ERROR: N .NE. 11 IN DCPFFJ'

         return

      end if

c     Compute a lower bound for x if task = 'XL'.

      if (task .eq. 'XL') then
         do 10 i = 1, n
            x(i) = zero
   10    continue

         return

      end if

      sqrtp = sqrt(p)

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = 5.0
         x(2) = 2.5
         x(3) = 5.0
         x(4) = 1.0E-1
         x(5) = 5.0E-2*k(5)
         x(6) = k(6)/sqrtp
         x(7) = 5.0E1*k(7)/sqrtp
         x(8) = 1.0E3*k(8)/p
         x(9) = 5.0E2*k(9)/sqrtp
         x(10) = 5.0E4*k(10)/p
         x(11) = 2.0E1

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      pdx = p/x(11)
      sqpdx = sqrt(pdx)

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         xtau = zero
         do 20 i = 1, n - 1
            xtau = xtau + x(i)
   20    continue
         fvec(1) = x(1) + x(4) - 3.0
         fvec(2) = two*x(1) + x(2) + x(4) + x(7) + x(8) + x(9) +
     +             two*x(10) - rr
         fvec(3) = two*x(2) + two*x(5) + x(6) + x(7) - 8.0
         fvec(4) = two*x(3) + x(9) - 4.0*rr
         fvec(5) = k(5)*x(2)*x(4) - x(1)*x(5)
         fvec(6) = k(6)*sqrt(x(2)*x(4)) - sqrt(x(1))*x(6)*sqpdx
         fvec(7) = k(7)*sqrt(x(1)*x(2)) - sqrt(x(4))*x(7)*sqpdx
         fvec(8) = k(8)*x(1) - x(4)*x(8)*pdx
         fvec(9) = k(9)*x(1)*sqrt(x(3)) - x(4)*x(9)*sqpdx
         fvec(10) = k(10)*x(1)**2 - (x(4)**2)*x(10)*pdx
         fvec(11) = x(11) - xtau
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         do 40 j = 1, n
            do 30 i = 1, n - 1
               fjac(i,j) = zero
   30       continue
            fjac(n,j) = -one
   40    continue
         fjac(n,n) = one

         xfrac = one/(sqrt(x(11))**3)

         fjac(1,1) = one
         fjac(2,1) = two
         fjac(5,1) = -x(5)
         fjac(6,1) = -p5*x(6)*sqpdx/sqrt(x(1))
         fjac(7,1) = p5*k(7)*sqrt(x(2))/sqrt(x(1))
         fjac(8,1) = k(8)
         fjac(9,1) = k(9)*sqrt(x(3))
         fjac(10,1) = two*k(10)*x(1)

         fjac(2,2) = one
         fjac(3,2) = two
         fjac(5,2) = k(5)*x(4)
         fjac(6,2) = p5*k(6)*sqrt(x(4))/sqrt(x(2))
         fjac(7,2) = p5*k(7)*sqrt(x(1))/sqrt(x(2))

         fjac(4,3) = two
         fjac(9,3) = p5*k(9)*x(1)/sqrt(x(3))

         fjac(1,4) = one
         fjac(2,4) = one
         fjac(5,4) = k(5)*x(2)
         fjac(6,4) = p5*k(6)*sqrt(x(2))/sqrt(x(4))
         fjac(7,4) = -p5*x(7)*sqpdx/sqrt(x(4))
         fjac(8,4) = -x(8)*pdx
         fjac(9,4) = -x(9)*sqpdx
         fjac(10,4) = -two*x(4)*x(10)*pdx

         fjac(3,5) = two
         fjac(5,5) = -x(1)

         fjac(3,6) = one
         fjac(6,6) = -sqrt(x(1))*sqpdx

         fjac(2,7) = one
         fjac(3,7) = one
         fjac(7,7) = -sqrt(x(4))*sqpdx

         fjac(2,8) = one
         fjac(8,8) = -x(4)*pdx

         fjac(2,9) = one
         fjac(4,9) = one
         fjac(9,9) = -x(4)*sqpdx

         fjac(2,10) = two
         fjac(10,10) = -(x(4)**2)*pdx

         fjac(6,11) = p5*sqrt(x(1))*x(6)*sqrtp*xfrac
         fjac(7,11) = p5*sqrt(x(4))*x(7)*sqrtp*xfrac
         fjac(8,11) = x(4)*x(8)*p/(x(11)**2)
         fjac(9,11) = p5*x(4)*x(9)*sqrtp*xfrac
         fjac(10,11) = x(4)**2*x(10)*p/(x(11)**2)
      end if

      end
