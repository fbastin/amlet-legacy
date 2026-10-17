      subroutine sgdffj(m,n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer m, n, ldfjac
      real x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine sgdffj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the Gaussian data fitting problem.
c
c     The subroutine statement is
c
c       subroutine sgdffj(m,n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions.
c            For the exponential data fitting II problem m = 65.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the exponential data fitting II problem n = 11.
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
      real ten, two, zero
      parameter (zero=0.0,two=2.0,ten=1.0E1)

      integer i
      real temp, temp1, temp2, temp3, temp4
      real y(65)

      data y/1.366, 1.191, 1.112, 1.013, 9.91E-1, 8.85E-1, 8.31E-1,
     +     8.47E-1, 7.86E-1, 7.25E-1, 7.46E-1, 6.79E-1, 6.08E-1,
     +     6.55E-1, 6.16E-1, 6.06E-1, 6.02E-1, 6.26E-1, 6.51E-1,
     +     7.24E-1, 6.49E-1, 6.49E-1, 6.94E-1, 6.44E-1, 6.24E-1,
     +     6.61E-1, 6.12E-1, 5.58E-1, 5.33E-1, 4.95E-1, 5.0E-1, 4.23E-1,
     +     3.95E-1, 3.75E-1, 3.72E-1, 3.91E-1, 3.96E-1, 4.05E-1,
     +     4.28E-1, 4.29E-1, 5.23E-1, 5.62E-1, 6.07E-1, 6.53E-1,
     +     6.72E-1, 7.08E-1, 6.33E-1, 6.68E-1, 6.45E-1, 6.32E-1,
     +     5.91E-1, 5.59E-1, 5.97E-1, 6.25E-1, 7.39E-1, 7.1E-1, 7.29E-1,
     +     7.2E-1, 6.36E-1, 5.81E-1, 4.28E-1, 2.92E-1, 1.62E-1, 9.8E-2,
     +     5.4E-2/

c     Check input arguments for errors.

      if (m .ne. 65 .or. n .ne. 11) then
         task = 'ERROR: M .NE. 65 OR N .NE. 11 IN DGDFFJ'

         return

      end if

c     Compute a lower bound for x if task = 'XL'.

      if (task .eq. 'XL') then
         do 10 i = 1, n
            x(i) = zero
   10    continue

         return

      end if

c     Compute an upper bound for x if task = 'XU'.

      if (task .eq. 'XU') then
         do 20 i = 1, n
            x(i) = ten
   20    continue

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = 1.3
         x(2) = 6.5E-1
         x(3) = 6.5E-1
         x(4) = 7.0E-1
         x(5) = 6.0E-1
         x(6) = 3.0
         x(7) = 5.0
         x(8) = 7.0
         x(9) = 2.0
         x(10) = 4.5
         x(11) = 5.5

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      do 30 i = 1, m
         temp = real(i-1)/ten
         temp1 = exp(-x(5)*temp)
         temp2 = exp(-x(6)*(temp-x(9))**2)
         temp3 = exp(-x(7)*(temp-x(10))**2)
         temp4 = exp(-x(8)*(temp-x(11))**2)
         if (task .eq. 'F' .or. task .eq. 'FJ') fvec(i) = y(i) -
     +       (x(1)*temp1+x(2)*temp2+x(3)*temp3+x(4)*temp4)
         if (task .eq. 'J' .or. task .eq. 'FJ') then
            fjac(i,1) = -temp1
            fjac(i,2) = -temp2
            fjac(i,3) = -temp3
            fjac(i,4) = -temp4
            fjac(i,5) = temp*x(1)*temp1
            fjac(i,6) = x(2)*(temp-x(9))**2*temp2
            fjac(i,7) = x(3)*(temp-x(10))**2*temp3
            fjac(i,8) = x(4)*(temp-x(11))**2*temp4
            fjac(i,9) = -two*x(2)*x(6)*(temp-x(9))*temp2
            fjac(i,10) = -two*x(3)*x(7)*(temp-x(10))*temp3
            fjac(i,11) = -two*x(4)*x(8)*(temp-x(11))*temp4
         end if
   30 continue

      end
