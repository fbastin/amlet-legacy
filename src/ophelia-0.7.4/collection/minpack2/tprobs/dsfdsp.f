      subroutine dsfdsp(n,nint,nnz,indrow,indcol)
      integer n, nint, nnz
      integer indrow(*), indcol(*)
c     **********
c
c     Subroutine dsfdsp
c
c     This subroutine defines the sparsity pattern of the Jacobian
c     matrix for the swirling flow between disks problem.
c
c     The subroutine statement is
c
c       subroutine dsfdsp(n,nint,nnz,indrow,indcol)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 14*nint.
c         On exit n is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the k-stage
c            collocation.
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
c     Brett M. Averick, R. S. Maier, G. L. Xue, R. G. Carter.
c
c     **********
      integer bc, cpts, dim, fdeg, gdeg, mdeg, npi
      parameter (bc=3,cpts=4,fdeg=4,gdeg=2,mdeg=4)
      parameter (dim=mdeg+cpts-1,npi=2*cpts+gdeg+fdeg)

      integer eqn1, eqn2, i, j, k, m, var1, var2

c     Compute the sparsity structure.

      nnz = 0

c     Nonzeroes contributed by boundary equations at t = 0.

      indrow(1) = 1
      indcol(1) = 1
      indrow(2) = 2
      indcol(2) = 2
      indrow(3) = 3
      indcol(3) = cpts + fdeg + 1
      nnz = 3

c     Nonzeroes contributed by collocation equations.

      do 40 i = 1, nint
         var1 = (i-1)*npi
         eqn1 = var1 + bc
         var2 = var1 + cpts + fdeg
         eqn2 = eqn1 + cpts
         do 30 k = 1, cpts
            do 10 j = 1, cpts + fdeg
               nnz = nnz + 1
               indrow(nnz) = eqn1 + k
               indcol(nnz) = var1 + j
               nnz = nnz + 1
               indrow(nnz) = eqn2 + k
               indcol(nnz) = var1 + j
   10       continue
            do 20 j = 1, cpts + gdeg
               nnz = nnz + 1
               indrow(nnz) = eqn1 + k
               indcol(nnz) = var2 + j
               nnz = nnz + 1
               indrow(nnz) = eqn2 + k
               indcol(nnz) = var2 + j
   20       continue
   30    continue
   40 continue

c     Nonzeroes contributed by continuity equations.

      do 90 i = 1, nint - 1
         var1 = (i-1)*npi
         eqn1 = var1 + bc + 2*cpts
         var2 = var1 + fdeg + cpts
         eqn2 = eqn1 + fdeg
         do 60 m = 1, fdeg
            nnz = nnz + 1
            indrow(nnz) = eqn1 + m
            indcol(nnz) = var1 + npi + m
            do 50 j = m, cpts + fdeg
               nnz = nnz + 1
               indrow(nnz) = eqn1 + m
               indcol(nnz) = var1 + j
   50       continue
   60    continue
         do 80 m = 1, gdeg
            nnz = nnz + 1
            indrow(nnz) = eqn2 + m
            indcol(nnz) = var2 + npi + m
            do 70 j = m, cpts + gdeg
               nnz = nnz + 1
               indrow(nnz) = eqn2 + m
               indcol(nnz) = var2 + j
   70       continue
   80    continue
   90 continue

c     Nonzeroes contributed by boundary equations at t = 1.

      var1 = n - npi
      eqn1 = var1 + bc + 2*cpts
      var2 = var1 + fdeg + cpts
      eqn2 = eqn1 + fdeg
      do 110 m = 1, 2
         do 100 j = m, cpts + fdeg
            nnz = nnz + 1
            indrow(nnz) = eqn1 + m
            indcol(nnz) = var1 + j
  100    continue
  110 continue
      do 120 j = 1, cpts + gdeg
         nnz = nnz + 1
         indrow(nnz) = n
         indcol(nnz) = var2 + j
  120 continue

      end
