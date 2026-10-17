      subroutine satrfj(m,n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer m, n, ldfjac
      real x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine satrfj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the analysis of thermistor resistance problem.
c
c     The subroutine statement is
c
c       subroutine satrfj(m,n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions. m = 16.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 3.
c         On exit n is unchanged.
c
c       x is a real array of dimension n.
c         On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
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
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer i
      real temp, temp1, temp2
      real y(16)

      data y/3.478E4, 2.861E4, 2.365E4, 1.963E4, 1.637E4, 1.372E4,
     +     1.154E4, 9.744E3, 8.261E3, 7.03E3, 6.005E3, 5.147E3, 4.427E3,
     +     3.82E3, 3.307E3, 2.872E3/

c     Check input arguments for errors.

      if (m .ne. 16 .or. n .ne. 3) then
         task = 'ERROR: M .NE. 16 OR N .NE. 3 IN DATRFJ'

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = 2.0E-2
         x(2) = 4.0E3
         x(3) = 2.5E2

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      do 10 i = 1, m
         temp = real(5*i+45) + x(3)
         temp1 = x(2)/temp
         temp2 = exp(temp1)
         if (task .eq. 'F' .or. task .eq. 'FJ') fvec(i) = x(1)*temp2 -
     +       y(i)
         if (task .eq. 'J' .or. task .eq. 'FJ') then
            fjac(i,1) = temp2
            fjac(i,2) = x(1)*temp2/temp
            fjac(i,3) = -temp1*fjac(i,2)
         end if
   10 continue

      end
