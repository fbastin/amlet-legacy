      subroutine dsfijs(nx,ny,x,s,y,lambda)
      integer nx, ny
      double precision lambda
      double precision x(nx*ny), s(nx*ny), y(nx*ny)
c     **********
c
c     Subroutine dsfijs
c
c     This subroutine computes the product f'(x)*s = y, where
c     f'(x) is the Jacobian matrix for the solid fuel ignition problem.
c
c     The subroutine statement is
c
c       subroutine dsfijs(nx,ny,x,s,y,lambda)
c
c     where
c
c       nx is an integer variable.
c         On entry nx is the number of grid points in the first
c            coordinate direction.
c         On exit nx is unchanged.
c
c       ny is an integer variable.
c         On entry ny is the number of grid points in the second
c            coordinate direction.
c         On exit ny is unchanged.
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
c       lambda is a double precision variable.
c         On entry lambda is the Frank-Kamenetski parameter.
c         On exit lambda is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     R. S. Maier and R. G. Carter.
c
c     **********
      double precision one, two, zero
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)

      integer i, j, k
      double precision hx, hxdhy, hxhy, hxhyl, hy, hydhx, sb, sc, sl,
     +                 sr, st

      hx = one/dble(nx+1)
      hy = one/dble(ny+1)
      hxdhy = hx/hy
      hydhx = hy/hx
      hxhy = hx*hy
      hxhyl = hxhy*lambda

      do 20 j = 1, ny
         do 10 i = 1, nx
            k = (j-1)*nx + i
            st = zero
            sb = zero
            sl = zero
            sr = zero
            sc = s(k)
            if (i .ne. 1) sl = s(k-1)
            if (i .ne. nx) sr = s(k+1)
            if (j .ne. 1) sb = s(k-nx)
            if (j .ne. ny) st = s(k+nx)
            y(k) = hxdhy*(-st+two*sc-sb) + hydhx*(-sr+two*sc-sl) -
     +             hxhyl*exp(x(k))*s(k)
   10    continue
   20 continue

      end
