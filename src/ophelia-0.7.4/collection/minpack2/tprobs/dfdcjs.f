      subroutine dfdcjs(nx,ny,x,s,y,r)
      integer nx, ny
      double precision r
      double precision x(nx*ny), s(nx*ny), y(nx*ny)
c     **********
c
c     Subroutine dfdcjs
c
c     This subroutine computes the product f'(x)*s = y, where
c     f'(x) is the Jacobian matrix for the flow in a channel
c     problem at the point x.
c
c     The subroutine statement is
c
c       subroutine dfdcjs(nx,ny,x,s,y,r)
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

c       x is a double precision array of dimension nx*ny.
c         On entry x specifies the vector x.
c         On exit x is unchanged.
c
c       s is a double precision array of dimension nx*ny.
c         On entry s contains the vector s.
c         On exit s is unchanged.
c
c       y is a double precision array of dimension nx*ny.
c         On entry out need not be specified.
c         On exit y contains the product f'(x)*s.
c
c       r is a double precision variable.
c         On entry r is the Reynolds number.
c         On exit r is unchanged.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      double precision four, one, three, two, zero
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,four=4.0d0)

      integer i, j, k, n
      double precision dpdx, dpdy, hx, hx2, hx2hy2, hx4, hy, hy2, hy4,
     +                 p, pb, pbb, pbl, pblap, pbr, pl, pll, pllap, pr,
     +                 prlap, prr, pt, ptl, ptlap, ptr, ptt

      n = nx*ny
      hx = one/dble(nx+1)
      hy = one/dble(ny+1)
      hx2 = hx*hx
      hy2 = hy*hy
      hx4 = hx2*hx2
      hy4 = hy2*hy2
      hx2hy2 = hx2*hy2

      do 20 i = 1, ny
         do 10 j = 1, nx
            k = (i-1)*nx + j
            if (i .eq. 1 .or. j .eq. 1) then
               pbl = zero
            else
               pbl = x(k-nx-1)
            end if
            if (i .eq. 1) then
               pb = zero
               pbb = x(k)
            else if (i .eq. 2) then
               pb = x(k-nx)
               pbb = zero
            else
               pb = x(k-nx)
               pbb = x(k-2*nx)
            end if
            if (i .eq. 1 .or. j .eq. nx) then
               pbr = zero
            else
               pbr = x(k-nx+1)
            end if
            if (j .eq. 1) then
               pl = zero
               pll = x(k)
            else if (j .eq. 2) then
               pl = x(k-1)
               pll = zero
            else
               pl = x(k-1)
               pll = x(k-2)
            end if
            p = x(k)
            if (j .eq. nx-1) then
               pr = x(k+1)
               prr = zero
            else if (j .eq. nx) then
               pr = zero
               prr = x(k)
            else
               pr = x(k+1)
               prr = x(k+2)
            end if
            if (i .eq. ny .or. j .eq. 1) then
               ptl = zero
            else
               ptl = x(k+nx-1)
            end if
            if (i .eq. ny-1) then
               pt = x(k+nx)
               ptt = zero
            else if (i .eq. ny) then
               pt = zero
               ptt = x(k) + two*hy
            else
               pt = x(k+nx)
               ptt = x(k+2*nx)
            end if
            if (i .eq. ny .or. j .eq. nx) then
               ptr = zero
            else
               ptr = x(k+nx+1)
            end if

            dpdy = (pt-pb)/(two*hy)
            dpdx = (pr-pl)/(two*hx)

            pblap = (pbr-two*pb+pbl)/hx2 + (p-two*pb+pbb)/hy2
            pllap = (p-two*pl+pll)/hx2 + (ptl-two*pl+pbl)/hy2
            prlap = (prr-two*pr+p)/hx2 + (ptr-two*pr+pbr)/hy2
            ptlap = (ptr-two*pt+ptl)/hx2 + (ptt-two*pt+p)/hy2

            y(k) = zero
            if (i .gt. 2) y(k) = y(k) + (one/hy4-r*dpdx/hy2/(two*hy))*
     +                           s(k-2*nx)
            if (j .gt. 2) y(k) = y(k) + (one/hx4+r*dpdy/hx2/(two*hx))*
     +                           s(k-2)
            if (j .lt. nx-1) y(k) = y(k) +
     +                              (one/hx4-r*dpdy/hx2/(two*hx))*s(k+2)
            if (i .lt. ny-1) y(k) = y(k) +
     +                              (one/hy4+r*dpdx/hy2/(two*hy))*
     +                              s(k+2*nx)
            if (i .ne. 1 .and. j .ne. 1) y(k) = y(k) +
     +          (two/hx2hy2+r*(dpdy/hy2/(two*hx)-dpdx/hx2/(two*hy)))*
     +          s(k-nx-1)
            if (i .ne. 1 .and. j .ne. nx) y(k) = y(k) +
     +          (two/hx2hy2-r*(dpdy/hy2/(two*hx)+dpdx/hx2/(two*hy)))*
     +          s(k-nx+1)
            if (i .ne. ny .and. j .ne. 1) y(k) = y(k) +
     +          (two/hx2hy2+r*(dpdy/hy2/(two*hx)+dpdx/hx2/(two*hy)))*
     +          s(k+nx-1)
            if (i .ne. ny .and. j .ne. nx) y(k) = y(k) +
     +          (two/hx2hy2-r*(dpdy/hy2/(two*hx)-dpdx/hx2/(two*hy)))*
     +          s(k+nx+1)
            if (i .ne. 1) y(k) = y(k) + (-four*(one/hx2hy2+one/hy4)+
     +                           r*((prlap-pllap)/(two*hx)/(two*hy)+
     +                           dpdx*(one/hx2+one/hy2)/hy))*s(k-nx)
            if (j .ne. 1) y(k) = y(k) + (-four*(one/hx2hy2+one/hx4)-
     +                           r*((ptlap-pblap)/(two*hx)/(two*hy)+
     +                           dpdy*(one/hx2+one/hy2)/hx))*s(k-1)
            if (j .ne. nx) y(k) = y(k) +
     +                            (-four*(one/hx2hy2+one/hx4)+r*((ptlap-
     +                            pblap)/(two*hx)/(two*hy)+
     +                            dpdy*(one/hx2+one/hy2)/hx))*s(k+1)
            if (i .ne. ny) y(k) = y(k) +
     +                            (-four*(one/hx2hy2+one/hy4)-r*((prlap-
     +                            pllap)/(two*hx)/(two*hy)+
     +                            dpdx*(one/hx2+one/hy2)/hy))*s(k+nx)
            y(k) = y(k) + (two*(four/hx2hy2+three/hx4+three/hy4))*s(k)
            if (i .eq. 1) y(k) = y(k) + (one/hy4-r*(dpdx/hy2/(two*hy)))*
     +                           s(k)
            if (j .eq. 1) y(k) = y(k) + (one/hx4+r*(dpdy/hx2/(two*hx)))*
     +                           s(k)
            if (j .eq. nx) y(k) = y(k) +
     +                            (one/hx4-r*(dpdy/hx2/(two*hx)))*s(k)
            if (i .eq. ny) y(k) = y(k) +
     +                            (one/hy4+r*(dpdx/hy2/(two*hy)))*s(k)
   10    continue
   20 continue

c     Scale the result.

      do 30 k = 1, n
         y(k) = hx2*hy2*y(k)
   30 continue

      end
