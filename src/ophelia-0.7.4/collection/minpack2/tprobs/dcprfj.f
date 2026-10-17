      subroutine dcprfj(n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer n, ldfjac
      double precision x(n), fvec(n), fjac(ldfjac,n)
c     **********
c
c     Subroutine dcprfj
c
c     This subroutine computes the function and the Jacobian matrix
c     of the combustion of propane (reduced formulation) problem.
c
c     The subroutine statement is
c
c       subroutine dcprfj(n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 5.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
c            x is set according to task.
c
c       fvec is a double precision array of dimension n.
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
      double precision eight, forty, four, one, p, p005, p05, p5, rr,
     +                 ten, three, two, zero
      parameter (zero=0.0d0,p005=5.0d-3,p05=5.0d-2,p5=0.5d0,one=1.0d0,
     +          two=2.0d0,three=3.0d0,four=4.0d0,eight=8.0d0,ten=1.0d1,
     +          forty=4.0d1)
      parameter (p=forty,rr=ten)

      integer i
      double precision k(10), r(10)
      double precision sqrtp

      data k/0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.930d-1, 2.597d-3, 3.448d-3,
     +     1.799d-5, 2.155d-4, 3.846d-5/

c     Check input arguments for errors.

      if (n .ne. 5) then
         task = 'ERROR: N .NE. 5 IN DCPRFJ'

         return

      end if

c     Compute a lower bound for x if task = 'XL'.

      if (task .eq. 'XL') then
         do 10 i = 1, n
            x(i) = zero
   10    continue

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = p005
         x(2) = p005
         x(3) = p05
         x(4) = p5
         x(5) = p05

         return

      end if

c     Initialization.

      sqrtp = sqrt(p)
      r(5) = k(5)
      r(6) = k(6)/sqrtp
      r(7) = k(7)/sqrtp
      r(8) = k(8)/p
      r(9) = k(9)/sqrtp
      r(10) = k(10)/p

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         fvec(1) = x(1)*x(2) + x(1) - three*x(5)
         fvec(2) = two*x(1)*x(2) + x(1) + two*r(10)*x(2)**2 +
     +             x(2)*x(3)**2 + r(7)*x(2)*x(3) + r(9)*x(2)*x(4) +
     +             r(8)*x(2) - rr*x(5)
         fvec(3) = two*x(2)*x(3)**2 + r(7)*x(2)*x(3) +
     +             two*r(5)*x(3)**2 + r(6)*x(3) - eight*x(5)
         fvec(4) = r(9)*x(2)*x(4) + two*x(4)**2 - four*rr*x(5)
         fvec(5) = x(1)*x(2) + x(1) + r(10)*x(2)**2 + x(2)*x(3)**2 +
     +             r(7)*x(2)*x(3) + r(9)*x(2)*x(4) + r(8)*x(2) +
     +             r(5)*x(3)**2 + r(6)*x(3) + x(4)**2 - one
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         fjac(1,1) = x(2) + one
         fjac(2,1) = two*x(2) + one
         fjac(3,1) = zero
         fjac(4,1) = zero
         fjac(5,1) = x(2) + one

         fjac(1,2) = x(1)
         fjac(2,2) = two*x(1) + four*r(10)*x(2) + x(3)**2 + r(7)*x(3) +
     +               r(9)*x(4) + r(8)
         fjac(3,2) = two*x(3)**2 + r(7)*x(3)
         fjac(4,2) = r(9)*x(4)
         fjac(5,2) = x(1) + two*r(10)*x(2) + x(3)**2 + r(7)*x(3) +
     +               r(9)*x(4) + r(8)

         fjac(1,3) = zero
         fjac(2,3) = two*x(2)*x(3) + r(7)*x(2)
         fjac(3,3) = four*x(2)*x(3) + r(7)*x(2) + four*r(5)*x(3) + r(6)
         fjac(4,3) = zero
         fjac(5,3) = two*x(2)*x(3) + r(7)*x(2) + two*r(5)*x(3) + r(6)

         fjac(1,4) = zero
         fjac(2,4) = r(9)*x(2)
         fjac(3,4) = zero
         fjac(4,4) = r(9)*x(2) + four*x(4)
         fjac(5,4) = r(9)*x(2) + two*x(4)

         fjac(1,5) = -three
         fjac(2,5) = -rr
         fjac(3,5) = -eight
         fjac(4,5) = -four*rr
         fjac(5,5) = zero
      end if

      end
