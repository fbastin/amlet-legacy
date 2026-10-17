      subroutine dctsfj(m,n,x,fvec,fjac,ldfjac,task)
      character*(*) task
      integer m, n, ldfjac
      double precision x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine dctsfj
c
c     This subroutine computes the function and the Jacobian matrix of
c     the coating thickness standardization problem.
c
c     The subroutine statement is
c
c       subroutine dctsfj(m,n,x,fvec,fjac,ldfjac,task)
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
      integer mdiv2, mdiv4
      parameter (mdiv4=63,mdiv2=126)
      double precision one, scale1, scale2, zero
      parameter (zero=0.0d0,one=1.0d0)
      parameter (scale1=4.08d0,scale2=0.417d0)

      integer i, j
      double precision indvar(mdiv4,2), y(mdiv2)

      data (indvar(i,1),i=1,mdiv4)/0.7140d0, 0.7169d0, 0.7232d0,
     +     0.7151d0, 0.6848d0, 0.7070d0, 0.7177d0, 0.7073d0, 0.6734d0,
     +     0.7174d0, 0.7125d0, 0.6947d0, 0.7121d0, 0.7166d0, 0.6894d0,
     +     0.6897d0, 0.7024d0, 0.7026d0, 0.6800d0, 0.6957d0, 0.6987d0,
     +     0.7111d0, 0.7097d0, 0.6809d0, 0.7139d0, 0.7046d0, 0.6950d0,
     +     0.7032d0, 0.7019d0, 0.6975d0, 0.6955d0, 0.7056d0, 0.6965d0,
     +     0.6848d0, 0.6995d0, 0.6105d0, 0.6027d0, 0.6084d0, 0.6081d0,
     +     0.6057d0, 0.6116d0, 0.6052d0, 0.6136d0, 0.6032d0, 0.6081d0,
     +     0.6092d0, 0.6122d0, 0.6157d0, 0.6191d0, 0.6169d0, 0.5483d0,
     +     0.5371d0, 0.5576d0, 0.5521d0, 0.5495d0, 0.5499d0, 0.4937d0,
     +     0.5092d0, 0.5433d0, 0.5018d0, 0.5363d0, 0.4977d0, 0.5296d0/
      data (indvar(i,2),i=1,mdiv4)/5.145d0, 5.241d0, 5.389d0, 5.211d0,
     +     5.154d0, 5.105d0, 5.191d0, 5.013d0, 5.582d0, 5.208d0,
     +     5.142d0, 5.284d0, 5.262d0, 6.838d0, 6.215d0, 6.817d0,
     +     6.889d0, 6.732d0, 6.717d0, 6.468d0, 6.776d0, 6.574d0,
     +     6.465d0, 6.090d0, 6.350d0, 4.255d0, 4.154d0, 4.211d0,
     +     4.287d0, 4.104d0, 4.007d0, 4.261d0, 4.150d0, 4.040d0,
     +     4.155d0, 5.086d0, 5.021d0, 5.040d0, 5.247d0, 5.125d0,
     +     5.136d0, 4.949d0, 5.253d0, 5.154d0, 5.227d0, 5.120d0,
     +     5.291d0, 5.294d0, 5.304d0, 5.209d0, 5.384d0, 5.490d0,
     +     5.563d0, 5.532d0, 5.372d0, 5.423d0, 7.237d0, 6.944d0,
     +     6.957d0, 7.138d0, 7.009d0, 7.074d0, 7.046d0/
      data (y(i),i=1,mdiv4)/9.3636d0, 9.3512d0, 9.4891d0, 9.1888d0,
     +     9.3161d0, 9.2585d0, 9.2913d0, 9.3914d0, 9.4524d0, 9.4995d0,
     +     9.4179d0, 9.468d0, 9.4799d0, 11.2917d0, 11.5062d0, 11.4579d0,
     +     11.3977d0, 11.3688d0, 11.3897d0, 11.3104d0, 11.3882d0,
     +     11.3629d0, 11.3149d0, 11.2474d0, 11.2507d0, 8.1678d0,
     +     8.1017d0, 8.3506d0, 8.3651d0, 8.2994d0, 8.1514d0, 8.2229d0,
     +     8.1027d0, 8.3785d0, 8.4118d0, 8.0955d0, 8.0613d0, 8.0979d0,
     +     8.1364d0, 8.1700d0, 8.1684d0, 8.0885d0, 8.1839d0, 8.1478d0,
     +     8.1827d0, 8.029d0, 8.1000d0, 8.2579d0, 8.2248d0, 8.2540d0,
     +     6.8518d0, 6.8547d0, 6.8831d0, 6.9137d0, 6.8984d0, 6.8888d0,
     +     8.5189d0, 8.5308d0, 8.5184d0, 8.5222d0, 8.5705d0, 8.5353d0,
     +     8.5213d0/
      data (y(i),i=mdiv4+1,mdiv2)/8.3158d0, 8.1995d0, 8.2283d0,
     +     8.1857d0, 8.2738d0, 8.2131d0, 8.2613d0, 8.2315d0, 8.2078d0,
     +     8.2996d0, 8.3026d0, 8.0995d0, 8.2990d0, 9.6753d0, 9.6687d0,
     +     9.5704d0, 9.5435d0, 9.6780d0, 9.7668d0, 9.7827d0, 9.7844d0,
     +     9.7011d0, 9.8006d0, 9.7610d0, 9.7813d0, 7.3073d0, 7.2572d0,
     +     7.4686d0, 7.3659d0, 7.3587d0, 7.3132d0, 7.3542d0, 7.2339d0,
     +     7.4375d0, 7.4022d0, 10.7914d0, 10.6554d0, 10.7359d0,
     +     10.7583d0, 10.7735d0, 10.7907d0, 10.6465d0, 10.6994d0,
     +     10.7756d0, 10.7402d0, 10.6800d0, 10.7000d0, 10.8160d0,
     +     10.6921d0, 10.8677d0, 12.3495d0, 12.4424d0, 12.4303d0,
     +     12.5086d0, 12.4513d0, 12.4625d0, 16.2290d0, 16.2781d0,
     +     16.2082d0, 16.2715d0, 16.2464d0, 16.1626d0, 16.1568d0/

c     Check input arguments for errors.

      if (m .ne. 252 .or. n .ne. 134) then
         task = 'ERROR: M .NE. 252 .OR. N .NE. 134 IN DCTSFJ'

         return

      end if

c     Compute the standard starting point if task = 'XS'.

      if (task .eq. 'XS') then
         x(1) = -8.0d0
         x(2) = 13.0d0
         x(3) = 1.2d0
         x(4) = 0.2d0
         x(5) = 0.1d0
         x(6) = 6.0d0
         x(7) = 5.5d0
         x(8) = -5.2d0
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
