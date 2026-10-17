      subroutine sficsp(n,nint,nnz,indrow,indcol)
      integer n, nint, nnz
      integer indrow(*), indcol(*)
c     **********
c
c     Subroutine sficsp
c
c     This subroutine defines the sparsity pattern of the Jacobian
c     matrix for the flow In a channel problem.
c
c     The subroutine statement is
c
c       subroutine sficsp(n,nint,nnz,indrow,indcol)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 8*nint.
c         On exit n is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the
c            k-stage collocation.
c         On exit nint is unchanged.
c
c       nnz is an integer variable.
c         On entry nnz need not be specified
c         On exit nnz is set to the number of nonzero index pairs
c            in the sparsity structure.
c
c       indrow is an integer array of dimension at least nnz.
c         On entry indrow need not be specified.
c         On exit indrow contains the row indices of the nonzeros
c            in the sparsity structure of the Jacobian matrix.
c
c       indcol is an integer array of dimension at least nnz.
c         On entry indcol need not be specified.
c         On exit indcol contains the column indices of the nonzeros
c            in the sparsity structure of the Jacobian matrix.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, R. S. Maier, G. L. Xue, R. G. Carter
c
c     **********
      integer bc, cpts, deg, dim, npi
      parameter (bc=2,cpts=4,deg=4,dim=deg+cpts-1,npi=cpts+deg)

      integer eqn, i, j, k, m, var

c     Nonzeroes contributed by boundary equations at t = 0.

      indrow(1) = 1
      indcol(1) = 1
      indrow(2) = 2
      indcol(2) = 2
      nnz = 2

c     Nonzeroes contributed by collocation equations.

      do 30 i = 1, nint
         var = (i-1)*npi
         eqn = var + bc
         do 20 k = 1, cpts
            do 10 j = 1, npi
               nnz = nnz + 1
               indrow(nnz) = eqn + k
               indcol(nnz) = var + j
   10       continue
   20    continue
   30 continue

c     Nonzeroes contributed by the continuity equations.

      do 60 i = 1, nint - 1
         var = (i-1)*npi
         eqn = var + bc + cpts
         do 50 m = 1, deg
            nnz = nnz + 1
            indrow(nnz) = eqn + m
            indcol(nnz) = var + cpts + deg + m
            do 40 j = m, npi
               nnz = nnz + 1
               indrow(nnz) = eqn + m
               indcol(nnz) = var + j
   40       continue
   50    continue
   60 continue

c     Nonzeroes contributed by the boundary equations at t = 1.

      var = n - npi
      eqn = var + bc + cpts
      do 80 m = 1, 2
         do 70 j = m, npi
            nnz = nnz + 1
            indrow(nnz) = eqn + m
            indcol(nnz) = var + j
   70    continue
   80 continue

      end
