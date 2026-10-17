      subroutine sctsfj(m,n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer m, n, ldfjac
      real x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine sctsfj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the coating thickness standardization problem.
c
c     The subroutine statement is
c
c       subroutine sctsfj(m,n,x,fvec,fjac,ldfjac,task)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions.
c            For the coating thickness standardization problem m = 252.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables.
c            For the coating thickness standardization problem n = 134.
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
c
c         On exit task is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer mdiv2, mdiv4
      parameter (mdiv4=63,mdiv2=126)
      real one, scale1, scale2, zero
      parameter (zero=0.0,one=1.0)
      parameter (scale1=4.08,scale2=0.417)

      integer i, j
      real indvar(mdiv4,2), y(mdiv2)

      data (indvar(i,1),i=1,mdiv4)/0.7140, 0.7169, 0.7232, 0.7151,
     +     0.6848, 0.7070, 0.7177, 0.7073, 0.6734, 0.7174, 0.7125,
     +     0.6947, 0.7121, 0.7166, 0.6894, 0.6897, 0.7024, 0.7026,
     +     0.6800, 0.6957, 0.6987, 0.7111, 0.7097, 0.6809, 0.7139,
     +     0.7046, 0.6950, 0.7032, 0.7019, 0.6975, 0.6955, 0.7056,
     +     0.6965, 0.6848, 0.6995, 0.6105, 0.6027, 0.6084, 0.6081,
     +     0.6057, 0.6116, 0.6052, 0.6136, 0.6032, 0.6081, 0.6092,
     +     0.6122, 0.6157, 0.6191, 0.6169, 0.5483, 0.5371, 0.5576,
     +     0.5521, 0.5495, 0.5499, 0.4937, 0.5092, 0.5433, 0.5018,
     +     0.5363, 0.4977, 0.5296/
      data (indvar(i,2),i=1,mdiv4)/5.145, 5.241, 5.389, 5.211, 5.154,
     +     5.105, 5.191, 5.013, 5.582, 5.208, 5.142, 5.284, 5.262,
     +     6.838, 6.215, 6.817, 6.889, 6.732, 6.717, 6.468, 6.776,
     +     6.574, 6.465, 6.090, 6.350, 4.255, 4.154, 4.211, 4.287,
     +     4.104, 4.007, 4.261, 4.150, 4.040, 4.155, 5.086, 5.021,
     +     5.040, 5.247, 5.125, 5.136, 4.949, 5.253, 5.154, 5.227,
     +     5.120, 5.291, 5.294, 5.304, 5.209, 5.384, 5.490, 5.563,
     +     5.532, 5.372, 5.423, 7.237, 6.944, 6.957, 7.138, 7.009,
     +     7.074, 7.046/
      data (y(i),i=1,mdiv4)/9.3636, 9.3512, 9.4891, 9.1888, 9.3161,
     +     9.2585, 9.2913, 9.3914, 9.4524, 9.4995, 9.4179, 9.468,
     +     9.4799, 11.2917, 11.5062, 11.4579, 11.3977, 11.3688, 11.3897,
     +     11.3104, 11.3882, 11.3629, 11.3149, 11.2474, 11.2507, 8.1678,
     +     8.1017, 8.3506, 8.3651, 8.2994, 8.1514, 8.2229, 8.1027,
     +     8.3785, 8.4118, 8.0955, 8.0613, 8.0979, 8.1364, 8.1700,
     +     8.1684, 8.0885, 8.1839, 8.1478, 8.1827, 8.029, 8.1000,
     +     8.2579, 8.2248, 8.2540, 6.8518, 6.8547, 6.8831, 6.9137,
     +     6.8984, 6.8888, 8.5189, 8.5308, 8.5184, 8.5222, 8.5705,
     +     8.5353, 8.5213/
      data (y(i),i=mdiv4+1,mdiv2)/8.3158, 8.1995, 8.2283, 8.1857,
     +     8.2738, 8.2131, 8.2613, 8.2315, 8.2078, 8.2996, 8.3026,
     +     8.0995, 8.2990, 9.6753, 9.6687, 9.5704, 9.5435, 9.6780,
     +     9.7668, 9.7827, 9.7844, 9.7011, 9.8006, 9.7610, 9.7813,
     +     7.3073, 7.2572, 7.4686, 7.3659, 7.3587, 7.3132, 7.3542,
     +     7.2339, 7.4375, 7.4022, 10.7914, 10.6554, 10.7359, 10.7583,
     +     10.7735, 10.7907, 10.6465, 10.6994, 10.7756, 10.7402,
     +     10.6800, 10.7000, 10.8160, 10.6921, 10.8677, 12.3495,
     +     12.4424, 12.4303, 12.5086, 12.4513, 12.4625, 16.2290,
     +     16.2781, 16.2082, 16.2715, 16.2464, 16.1626, 16.1568/

c     Check input arguments for errors.

      if (m .ne. 252 .or. n .ne. 134) then
         task = 'ERROR: M .NE. 252 .OR. N .NE. 134 IN DCTSFJ'

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = -8.0
         x(2) = 13.0
         x(3) = 1.2
         x(4) = 0.2
         x(5) = 0.1
         x(6) = 6.0
         x(7) = 5.5
         x(8) = -5.2
         do 10 i = 1, mdiv2
            x(8+i) = zero
   10    continue

         return

      end if

c     Evaluate the function if task = 'F', the Jacobian matrix if
c     task = 'J', or both if task = 'FJ'.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         do 20 i = 1, mdiv4
            fvec(i) = x(1) + x(2)*(indvar(i,1)+x(8+i)) +
     +                x(3)*(indvar(i,2)+x(8+i+mdiv4)) +
     +                x(4)*(indvar(i,1)+x(8+i))*
     +                (indvar(i,2)+x(8+i+mdiv4)) - y(i)
            fvec(i+mdiv4) = x(5) + x(6)*(indvar(i,1)+x(8+i)) +
     +                      x(7)*(indvar(i,2)+x(8+i+mdiv4)) +
     +                      x(8)*(indvar(i,1)+x(8+i))*
     +                      (indvar(i,2)+x(8+i+mdiv4)) - y(i+mdiv4)
            fvec(i+2*mdiv4) = scale1*x(8+i)
            fvec(i+3*mdiv4) = scale2*x(8+i+mdiv4)
   20    continue
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         do 40 j = 1, 8 + 2*mdiv4
            do 30 i = 1, m
               fjac(i,j) = zero
   30       continue
   40    continue

         do 50 i = 1, mdiv4
            fjac(i,1) = one
            fjac(i,2) = indvar(i,1) + x(8+i)
            fjac(i,3) = indvar(i,2) + x(8+i+mdiv4)
            fjac(i,4) = (indvar(i,1)+x(8+i))*(indvar(i,2)+x(8+i+mdiv4))
            fjac(i,8+i) = x(2) + x(4)*(indvar(i,2)+x(8+i+mdiv4))
            fjac(i,8+i+mdiv4) = x(3) + x(4)*(indvar(i,1)+x(8+i))
            fjac(i+mdiv4,5) = one
            fjac(i+mdiv4,6) = indvar(i,1) + x(8+i)
            fjac(i+mdiv4,7) = indvar(i,2) + x(8+i+mdiv4)
            fjac(i+mdiv4,8) = (indvar(i,1)+x(8+i))*
     +                        (indvar(i,2)+x(8+i+mdiv4))
            fjac(i+mdiv4,8+i) = x(6) + x(8)*(indvar(i,2)+x(8+i+mdiv4))
            fjac(i+mdiv4,8+i+mdiv4) = x(7) + x(8)*(indvar(i,1)+x(8+i))
            fjac(i+2*mdiv4,8+i) = scale1
            fjac(i+3*mdiv4,8+i+mdiv4) = scale2
   50    continue
      end if

      end
