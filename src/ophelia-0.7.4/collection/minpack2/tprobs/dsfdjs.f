      subroutine dsfdjs(n,x,s,y,eps,nint)
      integer n, nint
      double precision eps
      double precision x(n), s(n), y(n)
c     **********
c
c     Subroutine dsfdjs
c
c     This subroutine computes the product f'(x)*s = y, where
c     f'(x) is the Jacobian matrix for the swirling flow between
c     disks problem at the point x.
c
c     The subroutine statement is
c
c       subroutine dsfdjs(n,x,s,y,eps,nint)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. 14*nint.
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
c         On entry y need not be specified.
c         On exit y contains the product f'(x)*s.
c
c       eps is a double precision variable.
c         On entry eps is the viscosity of the fluid.
c         On exit eps is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the k-stage
c            collocation.
c         On exit nint is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, R. S. Maier, G. L. Xue, R. G. Carter.
c
c     **********
      integer bc, cpts, dim, fdeg, gdeg, mdeg, npi
      parameter (bc=3,cpts=4,fdeg=4,gdeg=2,mdeg=4)
      parameter (dim=mdeg+cpts-1,npi=2*cpts+gdeg+fdeg)
      double precision one, p5, two, zero
      parameter (zero=0.0d0,p5=0.5d0,one=1.0d0,two=2.0d0)
      double precision omega1, omega2
      parameter (omega1=-1.0d0,omega2=1.0d0)

      integer eqn1, eqn2, i, j, k, m, var1, var2
      double precision h, hm, nf, rhoijh, temp
      double precision dwf(fdeg+1,cpts+fdeg), dwg(gdeg+1,cpts+gdeg),
     +                 rhnfhk(cpts,0:dim,0:dim,0:mdeg), rho(cpts),
     +                 wf(fdeg+1), wg(gdeg+1)

      data (rho(i),i=1,cpts)/0.694318413734436035d-1,
     +     0.330009490251541138d0, 0.669990539550781250d0,
     +     0.930568158626556396d0/

c     Initialize.

      h = one/dble(nint)

c     Store all possible combinations of rho, h, and n factorial.

      hm = one
      do 40 m = 0, mdeg
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

      do 50 k = 1, n
         y(k) = zero
   50 continue
      do 70 k = 1, cpts + fdeg
         do 60 j = 1, fdeg + 1
            dwf(j,k) = zero
   60    continue
   70 continue
      do 90 k = 1, cpts + gdeg
         do 80 j = 1, gdeg + 1
            dwg(j,k) = zero
   80    continue
   90 continue

c     Contributions from boundary equations at t = 0.
c     f(0) = 0, f'(0) = 0, g(0) = omega1.

      y(1) = y(1) + s(1)
      y(2) = y(2) + s(2)
      y(3) = y(3) + s(cpts+fdeg+1)

c     Set up the collocation equations.

      do 190 i = 1, nint
         var1 = (i-1)*npi
         eqn1 = var1 + bc
         var2 = var1 + cpts + fdeg
         eqn2 = eqn1 + cpts
         do 180 k = 1, cpts
            do 120 m = 1, fdeg + 1
               wf(m) = zero
               do 100 j = m, fdeg
                  wf(m) = wf(m) + rhnfhk(k,j-m,j-m,j-m)*x(var1+j)
                  dwf(m,j) = rhnfhk(k,j-m,j-m,j-m)
  100          continue
               do 110 j = 1, cpts
                  wf(m) = wf(m) + x(var1+fdeg+j)*
     +                    rhnfhk(k,fdeg+j-m,fdeg+j-m,fdeg-m+1)
                  dwf(m,fdeg+j) = rhnfhk(k,fdeg+j-m,fdeg+j-m,fdeg-m+1)
  110          continue
  120       continue
            do 150 m = 1, gdeg + 1
               wg(m) = zero
               do 130 j = m, gdeg
                  wg(m) = wg(m) + rhnfhk(k,j-m,j-m,j-m)*x(var2+j)
                  dwg(m,j) = rhnfhk(k,j-m,j-m,j-m)
  130          continue
               do 140 j = 1, cpts
                  wg(m) = wg(m) + x(var2+gdeg+j)*
     +                    rhnfhk(k,gdeg+j-m,gdeg+j-m,gdeg-m+1)
                  dwg(m,gdeg+j) = rhnfhk(k,gdeg+j-m,gdeg+j-m,gdeg-m+1)
  140          continue
  150       continue

c           Contributions from collocation equations.

            do 160 j = 1, cpts + fdeg
               temp = eps*dwf(5,j) + dwf(4,j)*wf(1) + wf(4)*dwf(1,j)
               y(eqn1+k) = y(eqn1+k) + s(var1+j)*temp
               temp = dwf(1,j)*wg(2) - dwf(2,j)*wg(1)
               y(eqn2+k) = y(eqn2+k) + s(var1+j)*temp
  160       continue
            do 170 j = 1, cpts + gdeg
               temp = dwg(2,j)*wg(1) + wg(2)*dwg(1,j)
               y(eqn1+k) = y(eqn1+k) + s(var2+j)*temp
               temp = eps*dwg(3,j) + wf(1)*dwg(2,j) - wf(2)*dwg(1,j)
               y(eqn2+k) = y(eqn2+k) + s(var2+j)*temp
  170       continue
  180    continue
  190 continue

c     Set up the continuity equations.

      do 300 i = 1, nint - 1
         var1 = (i-1)*npi
         eqn1 = var1 + bc + 2*cpts
         var2 = var1 + fdeg + cpts
         eqn2 = eqn1 + fdeg
         do 220 m = 1, fdeg
            wf(m) = zero
            do 200 j = m, fdeg
               wf(m) = wf(m) + rhnfhk(1,0,j-m,j-m)*x(var1+j)
               dwf(m,j) = rhnfhk(1,0,j-m,j-m)
  200       continue
            do 210 j = 1, cpts
               wf(m) = wf(m) + rhnfhk(1,0,fdeg+j-m,fdeg-m+1)*
     +                 x(var1+fdeg+j)
               dwf(m,fdeg+j) = rhnfhk(1,0,fdeg+j-m,fdeg-m+1)
  210       continue
  220    continue
         do 250 m = 1, gdeg
            wg(m) = zero
            do 230 j = m, gdeg
               wg(m) = wg(m) + rhnfhk(1,0,j-m,j-m)*x(var2+j)
               dwg(m,j) = rhnfhk(1,0,j-m,j-m)
  230       continue
            do 240 j = 1, cpts
               wg(m) = wg(m) + rhnfhk(1,0,gdeg+j-m,gdeg-m+1)*
     +                 x(var2+gdeg+j)
               dwg(m,gdeg+j) = rhnfhk(1,0,gdeg+j-m,gdeg-m+1)
  240       continue
  250    continue

c        Contributions from continuity equations.

         do 270 m = 1, fdeg
            y(eqn1+m) = y(eqn1+m) + s(var1+npi+m)
            do 260 j = 1, cpts + fdeg
               y(eqn1+m) = y(eqn1+m) + s(var1+j)*(-dwf(m,j))
  260       continue
  270    continue
         do 290 m = 1, gdeg
            y(eqn2+m) = y(eqn2+m) + s(var2+npi+m)
            do 280 j = 1, cpts + gdeg
               y(eqn2+m) = y(eqn2+m) + s(var2+j)*(-dwg(m,j))
  280       continue
  290    continue
  300 continue

c     Prepare for setting up the boundary conditions at t = 1.

      var1 = n - npi
      var2 = var1 + fdeg + cpts
      do 330 m = 1, fdeg + 1
         wf(m) = zero
         do 310 j = m, fdeg
            wf(m) = wf(m) + rhnfhk(1,0,j-m,j-m)*x(var1+j)
            dwf(m,j) = rhnfhk(1,0,j-m,j-m)
  310    continue
         do 320 j = 1, cpts
            wf(m) = wf(m) + rhnfhk(1,0,fdeg+j-m,fdeg-m+1)*x(var1+fdeg+j)
            dwf(m,fdeg+j) = rhnfhk(1,0,fdeg+j-m,fdeg-m+1)
  320    continue
  330 continue
      do 360 m = 1, gdeg + 1
         wg(m) = zero
         do 340 j = m, gdeg
            wg(m) = wg(m) + rhnfhk(1,0,j-m,j-m)*x(var2+j)
            dwg(m,j) = rhnfhk(1,0,j-m,j-m)
  340    continue
         do 350 j = 1, cpts
            wg(m) = wg(m) + rhnfhk(1,0,gdeg+j-m,gdeg-m+1)*x(var2+gdeg+j)
            dwg(m,gdeg+j) = rhnfhk(1,0,gdeg+j-m,gdeg-m+1)
  350    continue
  360 continue

c     Contributions from the boundary equations at t = 1.
c     f(1) = 0, f'(1) = 0, g(0) = omega1.

      do 370 j = 1, cpts + fdeg
         y(n-2) = y(n-2) + s(var1+j)*dwf(1,j)
         y(n-1) = y(n-1) + s(var1+j)*dwf(2,j)
  370 continue
      do 380 j = 1, cpts + gdeg
         y(n) = y(n) + s(var2+j)*dwg(1,j)
  380 continue

      end
