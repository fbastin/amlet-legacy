      subroutine diacfj(m,n,x,fvec,fjac,ldfjac,task,nint,sigma)
      character*(*) task
      integer m, n, ldfjac, nint
      double precision sigma
      double precision x(n), fvec(m), fjac(ldfjac,n)
c     **********
c
c     Subroutine diacfj
c
c     This subroutine computes the function and Jacobian matrix of the
c     isomerization of alpha-pinene (collocation formulation) problem.
c
c     The subroutine statement is
c
c       subroutine diacfj(m,n,x,fvec,fjac,ldfjac,task,nint,sigma)
c
c     where
c
c       m is an integer variable.
c         On entry m is the number of functions. m = 25*nint + 40.
c         On exit m is unchanged.
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 25*nint + 5.
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
c            'F'      Evaluate the function at x.
c            'J'      Evaluate the Jacobian matrix at x.
c            'FJ'     Evaluate the function and the Jacobian at x.
c            'XS'     Set x to the standard starting point xs.
c
c         On exit task is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the
c            k-stage collocation.
c         On exit nint is unchanged.
c
c       sigma is a double precision variable.
c         On entry sigma is the penalty constant for the constraint
c            equations.
c         On exit sigma is unchanged
c
c     Subprograms called
c
c       MINPACK-2 ... diarfj, diaofj
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer cpts, ntimes, s
      parameter (cpts=4,ntimes=8,s=5)

      integer i, ideg, j, neqn, npi
      integer deg(s)

      external diarfj, diaofj

      data (deg(i),i=1,s)/1, 1, 1, 1, 1/

c     Check input arguments for errors.

      if (m .ne. 25*nint+40 .or. n .ne. 25*nint+5) then
         task = 'ERROR: M .NE. 25*NINT + 40 OR N .NE. 25*NINT + 5 '
         task = task//'IN DIACFJ'

         return

      end if

      if (task .eq. 'XS') then
         call diaofj(m-40,n,x,fvec,fjac,ldfjac,task,nint)

         return

      end if

c     Initialize.

      ideg = 0
      do 10 i = 1, s
         ideg = ideg + deg(i)
   10 continue
      npi = s*cpts + ideg
      neqn = nint*npi

c     Evaluate the data residual equations.

      call diarfj(40,n,x,fvec,fjac,ldfjac,task,nint)

c     Move the information down.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         do 20 i = 1, s*ntimes
            fvec(neqn+i) = fvec(i)
   20    continue
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         do 40 j = 1, n
            do 30 i = 1, s*ntimes
               fjac(neqn+i,j) = fjac(i,j)
   30       continue
   40    continue
      end if

c     Evaluate the collocation constraint equations.

      call diaofj(m-40,n,x,fvec,fjac,ldfjac,task,nint)

c     Penalize the collocation constraints by sigma.

      if (task .eq. 'F' .or. task .eq. 'FJ') then
         do 50 i = 1, neqn
            fvec(i) = sigma*fvec(i)
   50    continue
      end if

      if (task .eq. 'J' .or. task .eq. 'FJ') then
         do 70 j = 1, n
            do 60 i = 1, neqn
               fjac(i,j) = sigma*fjac(i,j)
   60       continue
   70    continue
      end if

      end
