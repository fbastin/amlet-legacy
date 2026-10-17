      subroutine diadfj(m,n,x,fvec,fjac,ldfjac,task,nh)
      character*(*) task
      integer m, n, ldfjac, nh
      double precision x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine diadfj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the isomerization of alpha-pinene (direct formulation) problem.
c
c     The subroutine statement is
c
c       subroutine diadfj(m,n,x,fvec,ldfjac,fjac,task,nh)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions. m = 40.
c         On exit m is unchanged.
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
c             'XL'    Set x to the lower bound xl.
c
c         On exit task is unchanged.
c
c       nh is an integer variable.
c           On entry nh is the number of Runge-Kutta steps taken
c              between observations. nh must be >= 1.
c           On exit nh is unchanged.
c
c     Subprograms called
c
c       MINPACK-2 ... diadde
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and the University of Minnesota.
c     Brett M. Averick and Richard G. Carter.
c
c     **********
      integer neq, neqext, nobs, nparm
      parameter (neq=5,nobs=8,nparm=5)
      parameter (neqext=neq+nparm*neq)
      double precision hundrd, six, three, two, zero
      parameter (zero=0.0d0,two=2.0d0,three=3.0d0,six=6.0d0,
     +          hundrd=1.0d2)

      integer i, iobs, j, k, ne
      double precision h, ho2, ho6, t
      double precision k1(neqext), k2(neqext), k3(neqext), k4(neqext),
     +                 obs(neq,nobs), time(nobs), work(neqext),
     +                 y(neqext)

      data obs/88.35d0, 7.3d0, 2.3d0, 0.4d0, 1.75d0, 76.40d0, 15.6d0,
     +     4.5d0, 0.7d0, 2.80d0, 65.10d0, 23.1d0, 5.3d0, 1.1d0, 5.80d0,
     +     50.40d0, 32.9d0, 6.0d0, 1.5d0, 9.30d0, 37.50d0, 42.7d0,
     +     6.0d0, 1.9d0, 12.00d0, 25.90d0, 49.1d0, 5.9d0, 2.2d0,
     +     17.00d0, 14.00d0, 57.4d0, 5.1d0, 2.6d0, 21.00d0, 4.50d0,
     +     63.1d0, 3.8d0, 2.9d0, 25.7d0/
      data time/1230.0d0, 3060.0d0, 4920.0d0, 7800.0d0, 10680.0d0,
     +     15030.0d0, 22620.0d0, 36420.0d0/

c     Check input arguments for errors.

      if (m .ne. 40 .or. n .ne. 5 .or. nh .lt. 1) then
         task = 'ERROR: M .NE. 40 OR N .NE. 5 OR NH .LT. 1 IN DIADFJ'

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
         x(1) = 5.84d-5
         x(2) = 2.65d-5
         x(3) = 1.63d-5
         x(4) = 2.777d-4
         x(5) = 4.61d-5

         return

      end if

c     Determine the number of equations to be constructed:
c     neq coupled odes for the function alone or neqext coupled odes
c     for the function and Jacobian matrix.

      ne = neq
      if (task .eq. 'J' .or. task .eq. 'FJ') ne = neqext

c     Set the initial conditions and initialize.

      y(1) = hundrd
      do 20 i = 2, ne
         y(i) = zero
   20 continue

c     Solve the ode's and evaluate the function if task = 'F', the
c     Jacobian matrix if task = 'J', or both if task = 'FJ'.

      t = zero
      do 110 iobs = 1, nobs

c        Calculate y(time(iobs)) using a fourth-order Runge-Kutta
c        scheme on nh subintervals of (time(iobs-1),time(iobs))

         h = (time(iobs)-t)/dble(nh)
         ho2 = h/two
         ho6 = h/six
         do 70 k = 1, nh

c           Compute the k1 vector.

            call diadde(ne,y,n,x,k1)

c           Compute the k2 vector.

            t = t + ho2
            do 30 i = 1, ne
               work(i) = ho2*k1(i) + y(i)
   30       continue
            call diadde(ne,work,n,x,k2)

c           Compute the k3 vector.

            do 40 i = 1, ne
               work(i) = ho2*k2(i) + y(i)
   40       continue
            call diadde(ne,work,n,x,k3)

c           Compute the k4 vector.

            t = t + ho2
            do 50 i = 1, ne
               work(i) = h*k3(i) + y(i)
   50       continue
            call diadde(ne,work,n,x,k4)

c           Compute ynew = y + h*(k1 +2*k2 + 2*k3 + k4)/6.

            do 60 i = 1, ne
               y(i) = y(i) + ho6*(k1(i)+two*(k2(i)+k3(i))+k4(i))
   60       continue

   70    continue

c        Set up the data residual equations.

         if (task .eq. 'F' .or. task .eq. 'FJ') then
            do 80 k = 1, neq
               j = k + (iobs-1)*neq
               fvec(j) = y(k) - obs(k,iobs)
   80       continue
         end if

         if (task .eq. 'J' .or. task .eq. 'FJ') then
            do 100 k = 1, neq
               j = k + (iobs-1)*neq
               do 90 i = 1, n
                  fjac(j,i) = y(k+5*i)
   90          continue
  100       continue
         end if

  110 continue

      end

      subroutine diadde(ne,y,n,x,ydot)
      integer ne, n
      double precision x(n), y(ne), ydot(ne)
c     **********
c
c     Subroutine diadde
c
c     This subroutine defines the system of ordinary differential
c     equation for the alpha-pinene problem:
c
c          ydot = h(y,x)
c
c     The subroutine statement is
c
c       subroutine diadde(ne,y,n,x,ydot)
c
c     where
c
c       ne is an integer variable.
c         On entry ne is the number of odes in the system.
c            ne = 5  denotes that the primary system is to be solved.
c            ne = 30 denotes that the extended system derived from the
c            sensitivity equations is to be solved.
c         On exit ne is unchanged.
c
c       n is an integer variable.
c         On entry n is the dimension of x. n should be 5 for the
c            alpha-pinene problem.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the array of parameters defining the ode.
c         On exit x is unchanged.
c
c       y is a double precision array of dimension ne.
c         On entry y specifies the vector valued function y evaluated at
c            the point (x).
c         On exit y is unchanged.
c
c       ydot is a double precision array of dimension ne.
c         On entry ydot need not be specified.
c         On exit ydot is set to h(y,x).
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Richard G. Carter.
c
c     **********

c     Compute components of the rhs of the ode.

      ydot(1) = -(x(1)+x(2))*y(1)
      ydot(2) = x(1)*y(1)
      ydot(3) = x(2)*y(1) - (x(3)+x(4))*y(3) + x(5)*y(5)
      ydot(4) = x(3)*y(3)
      ydot(5) = x(4)*y(3) - x(5)*y(5)

      if (ne .le. 5) return

c     Compute auxiliary equations.  y(i+5*j) is defined to be
c     the partial derivative of y(i) with respect to x(j)

      ydot(6) = -(x(1)+x(2))*y(6) - y(1)
      ydot(7) = x(1)*y(6) + y(1)
      ydot(8) = x(2)*y(6) - (x(3)+x(4))*y(8) + x(5)*y(10)
      ydot(9) = x(3)*y(8)
      ydot(10) = x(4)*y(8) - x(5)*y(10)

      ydot(11) = -(x(1)+x(2))*y(11) - y(1)
      ydot(12) = x(1)*y(11)
      ydot(13) = x(2)*y(11) + y(1) - (x(3)+x(4))*y(13) + x(5)*y(15)
      ydot(14) = x(3)*y(13)
      ydot(15) = x(4)*y(13) - x(5)*y(15)

      ydot(16) = -(x(1)+x(2))*y(16)
      ydot(17) = x(1)*y(16)
      ydot(18) = x(2)*y(16) - (x(3)+x(4))*y(18) - y(3) + x(5)*y(20)
      ydot(19) = x(3)*y(18) + y(3)
      ydot(20) = x(4)*y(18) - x(5)*y(20)

      ydot(21) = -(x(1)+x(2))*y(21)
      ydot(22) = x(1)*y(21)
      ydot(23) = x(2)*y(21) - (x(3)+x(4))*y(23) - y(3) + x(5)*y(25)
      ydot(24) = x(3)*y(23)
      ydot(25) = x(4)*y(23) - x(5)*y(25) + y(3)

      ydot(26) = -(x(1)+x(2))*y(26)
      ydot(27) = x(1)*y(26)
      ydot(28) = x(2)*y(26) - (x(3)+x(4))*y(28) + x(5)*y(30) + y(5)
      ydot(29) = x(3)*y(28)
      ydot(30) = x(4)*y(28) - x(5)*y(30) - y(5)

      end
