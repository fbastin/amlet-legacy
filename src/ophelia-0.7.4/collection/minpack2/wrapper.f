c     task = 1 corresponds to F
c     task = 2 corresponds to G or J
c     task = 3 corresponds to FG or FJ
c     task = 4 corresponds to XS
    
c     We use wrappers in order to cope with string conversions
c     difficulties between C and Fortran.

      subroutine wdaerfj(m, n, x, fvec, fjac, ldfjac, task)

      implicit none

      integer task
      integer m, n, ldfjac
      double precision x(n), fvec(m), fjac(ldfjac, n)

      if (task .eq. 1) then
         call daerfj(m, n, x, fvec, fjac, ldfjac, 'F');
      elseif (task .eq. 2) then
         call daerfj(m, n, x, fvec, fjac, ldfjac, 'J');
      elseif (task .eq. 3) then
         call daerfj(m, n, x, fvec, fjac, ldfjac, 'FJ');
      elseif (task .eq. 4) then
         call daerfj(m, n, x, fvec, fjac, ldfjac, 'XS')
      endif

      end

      subroutine wdmsafg(nx, ny, x, f, fgrad, task,
     +     bottom, top, left, right)

      implicit none

      integer nx, ny, task
      double precision f
      double precision x(nx*ny), fgrad(nx*ny), bottom(nx+2), top(nx+2),
     +                 left(ny+2), right(ny+2)

      if (task .eq. 1) then
         call dmsafg(nx, ny, x, f, fgrad, 'XS',
     +        bottom, top, left, right)
      endif

      end
