      subroutine daerfj(m,n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer m, n, ldfjac
      double precision x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine daerfj
c
c     This subroutine computes the function and the Jacobian matrix
c     of the analysis of an enzyme reaction problem.
c
c     The subroutine statement is
c
c       subroutine daerfj(m,n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions. m = 11.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 4.
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
c             'F'     Evaluate the function at x.
c             'J'     Evaluate the Jacobian matrix at x.
c             'FJ'    Evaluate the function and the Jacobian at x.
c             'XS'    Set x to the standard starting point xs.
c
c         On exit task is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer i
      double precision temp1, temp2
      double precision v(11), y(11)

      data v/4.0d0, 2.0d0, 1.0d0, 5.0d-1, 2.5d-1, 1.67d-1, 1.25d-1,
     +     1.0d-1, 8.33d-2, 7.14d-2, 6.25d-2/
      data y/1.957d-1, 1.947d-1, 1.735d-1, 1.6d-1, 8.44d-2, 6.27d-2,
     +     4.56d-2, 3.42d-2, 3.23d-2, 2.35d-2, 2.46d-2/

c     Check input arguments for errors.

      if (m .ne. 11 .or. n .ne. 4) then
         task = 'ERROR: M .NE. 11 OR N .NE. 4 IN DAERFJ'

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = 2.5d-1
         x(2) = 3.9d-1
         x(3) = 4.15d-1
         x(4) = 3.9d-1

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      do 10 i = 1, m
         temp1 = v(i)*(v(i)+x(2))
         temp2 = v(i)*(v(i)+x(3)) + x(4)
         if (task .eq. 'F' .or. task .eq. 'FJ') 
     +       fvec(i) = y(i) - x(1)*temp1/temp2
         if (task .eq. 'J' .or. task .eq. 'FJ') then
            fjac(i,1) = -temp1/temp2
            fjac(i,2) = -v(i)*x(1)/temp2
            fjac(i,3) = fjac(i,1)*fjac(i,2)
            fjac(i,4) = fjac(i,3)/v(i)
         end if
   10 continue

      end
