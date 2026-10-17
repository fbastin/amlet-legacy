      subroutine sedffj(m,n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer m, n, ldfjac
      real x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine sedffj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the exponential data fitting problem.
c
c     The subroutine statement is
c
c       subroutine sedffj(m,n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions. m = 33.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 5.
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
c             'F'     Evaluate the function at x.
c             'J'     Evaluate the Jacobian matrix at x.
c             'FJ'    Evaluate the function and the Jacobian at x.
c             'XS'    Set x to the standard starting point xs.
c             'XL'    Set x to the lower bound xl.
c             'XU'    Set x to the upper bound xu.
c
c         On exit task is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      real one, ten
      parameter (one=1.0,ten=1.0E1)

      integer i
      real temp, temp1, temp2
      real y(33)

      data y/8.44E-1, 9.08E-1, 9.32E-1, 9.36E-1, 9.25E-1, 9.08E-1,
     +     8.81E-1, 8.5E-1, 8.18E-1, 7.84E-1, 7.51E-1, 7.18E-1, 6.85E-1,
     +     6.58E-1, 6.28E-1, 6.03E-1, 5.8E-1, 5.58E-1, 5.38E-1, 5.22E-1,
     +     5.06E-1, 4.9E-1, 4.78E-1, 4.67E-1, 4.57E-1, 4.48E-1, 4.38E-1,
     +     4.31E-1, 4.24E-1, 4.2E-1, 4.14E-1, 4.11E-1, 4.06E-1/

c     Check input arguments for errors.

      if (m .ne. 33 .or. n .ne. 5) then
         task = 'ERROR: M .NE. 33 OR N .NE. 5 IN DEDFFJ'

         return

      end if


c     Compute a lower bound for x if task = 'XL' or an upper bound if
c     task = 'XU'.

      if (task .eq. 'XL' .or. task .eq. 'XU') then
         if (task .eq. 'XL') temp = -ten
         if (task .eq. 'XU') temp = ten
         do 10 i = 1, n
            x(i) = temp
   10    continue

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = 5.0E-1
         x(2) = 1.5
         x(3) = -1.0
         x(4) = 1.0E-2
         x(5) = 2.0E-2

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      do 20 i = 1, m
         temp = real(10*(i-1))
         temp1 = exp(-x(4)*temp)
         temp2 = exp(-x(5)*temp)
         if (task .eq. 'F' .or. task .eq. 'FJ') fvec(i) = y(i) -
     +       (x(1)+x(2)*temp1+x(3)*temp2)
         if (task .eq. 'J' .or. task .eq. 'FJ') then
            fjac(i,1) = -one
            fjac(i,2) = -temp1
            fjac(i,3) = -temp2
            fjac(i,4) = temp*x(2)*temp1
            fjac(i,5) = temp*x(3)*temp2
         end if
   20 continue

      end
