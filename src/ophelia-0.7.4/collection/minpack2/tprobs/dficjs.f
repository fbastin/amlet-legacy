      subroutine dficjs(n,x,s,y,r,nint)
      integer n, nint
      double precision r
      double precision x(n), s(n), y(n)
c     **********
c
c     Subroutine dficjs
c
c     This subroutine computes the product f'(x)*s = y, where
c     where f'(x) is the Jacobian matrix for the flow in a channel
c     problem at the point x.
c
c     The subroutine statement is
c
c       subroutine dficjs(n,x,s,y,r,nint)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 8*nint.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x.
c         On exit x is unchanged.
c
c       s is a double precision array of dimension n.
c         On entry s contains the vector s.
c         On exit s is unchanged.
c
c       y is a double precision array of dimension n.
c         On entry out need not be specified.
c         On exit y contains the product f'(x)*s.
c
c       r is a double precision variable.
c         On entry r is the Reynolds number.
c         On exit r is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the k-stage
c            collocation.
c         On exit nint is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, R. S. Maier, G. L. Xue, R. G. Carter
c
c     **********
      integer bc, cpts, deg, dim, npi
      parameter (bc=2,cpts=4,deg=4)
      parameter (dim=deg+cpts-1,npi=cpts+deg)
      double precision one, six, three, twelve, two, zero
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0,
     +          twelve=12.0d0)

      integer eqn, i, j, k, m, var
      double precision h, hm, nf, rhoijh
      double precision dw(deg+1,cpts+deg),
     +                 rhnfhk(cpts,0:dim,0:dim,0:deg), rho(cpts),
     +                 w(deg+1)

      data (rho(i),i=1,cpts)/0.694318413734436035d-1,
     +     0.330009490251541138d0, 0.669990539550781250d0,
     +     0.930568158626556396d0/

c     Store all possible combinations of rho, h, and n factorial.

      h = one/dble(nint)
      hm = one
      do 40 m = 0, deg
         do 30 i = 1, cpts
            rhoijh = hm
            do 20 j = 0, dim
               nf = one
               do 10 k = 0, dim
                  rhnfhk(i,j,k,m) = rhoijh/nf
                  nf = nf*dble(k+1)
   10          continue
               rhoijh = rhoijh*rho(i)
   20       continue
   30    continue
         hm = hm*h
   40 continue

c     Initialize arrays.

      do 50 i = 1, n
         y(i) = zero
   50 continue
      do 70 k = 1, npi
         do 60 j = 1, deg + 1
            dw(j,k) = zero
   60    continue
   70 continue

c     Contributions from boundary equations at t = 0.
c     u(0) = 0, u'(0) = 0.

      y(1) = s(1)
      y(2) = s(2)

c    Set up the collocation equations.

      do 180 i = 1, nint
         var = (i-1)*npi
         eqn = var + bc
         do 120 k = 1, cpts
            do 100 m = 1, deg + 1
               w(m) = zero
               do 80 j = m, deg
                  w(m) = w(m) + rhnfhk(k,j-m,j-m,j-m)*x(var+j)
                  dw(m,j) = rhnfhk(k,j-m,j-m,j-m)
   80          continue
               do 90 j = 1, cpts
                  w(m) = w(m) + rhnfhk(k,deg+j-m,deg+j-m,deg-m+1)*
     +                   x(var+deg+j)
                  dw(m,deg+j) = rhnfhk(k,deg+j-m,deg+j-m,deg-m+1)
   90          continue
  100       continue

c           Contributions from collocation equations.

            do 110 j = 1, npi
               y(eqn+k) = y(eqn+k) + s(var+j)*
     +                    (dw(5,j)-r*(dw(2,j)*w(3)+w(2)*dw(3,j)-
     +                    dw(1,j)*w(4)-w(1)*dw(4,j)))
  110       continue
  120    continue

c        Set up the continuity equations.

         eqn = var + bc + cpts
         do 150 m = 1, deg
            w(m) = zero
            do 130 j = m, deg
               w(m) = w(m) + rhnfhk(1,0,j-m,j-m)*x(var+j)
               dw(m,j) = rhnfhk(1,0,j-m,j-m)
  130       continue
            do 140 j = 1, cpts
               w(m) = w(m) + rhnfhk(1,0,deg+j-m,deg-m+1)*x(var+deg+j)
               dw(m,deg+j) = rhnfhk(1,0,deg+j-m,deg-m+1)
  140       continue
  150    continue
         if (i .eq. nint) go to 190

c        Contributions from continuity equations.

         do 170 m = 1, deg
            y(eqn+m) = y(eqn+m) + s(var+cpts+deg+m)
            do 160 j = 1, npi
               y(eqn+m) = y(eqn+m) + s(var+j)*(-dw(m,j))
  160       continue
  170    continue

  180 continue

c     Contributions from boundary equations at t = 1.
c     u(1) = 1, u'(1) = 0.

  190 continue

      do 200 j = 1, npi
         y(n-1) = y(n-1) + s(var+j)*dw(1,j)
         y(n) = y(n) + s(var+j)*dw(2,j)
  200 continue

      end
