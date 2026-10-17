      subroutine diersp(n,nint,nnz,indrow,indcol)
      integer n, nint, nnz
      integer indrow(*), indcol(*)
c     **********
c
c     Subroutine diersp
c
c     This subroutine computes the sparsity pattern of the Jacobian
c     matrix of the inverse elastic rod problem.
c
c     The subroutine statement is
c
c       subroutine diersp(n,nint,nnz,indrow,indcol)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables. n = 15*nint+3.
c         On exit n is unchanged.
c
c       nint is an integer variable.
c         On entry nint is the number of subintervals in the k-stage
c            collocation.
c         On exit nint is unchanged.
c
c       nnz is an integer variable.
c         On entry nnz needs not be specified.
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
c     Brett M. Averick, R. S. Maier, G. L. Xue.
c
c     **********
      integer cpts, dim, maxdeg, s
      parameter (cpts=4,maxdeg=1,s=3)
      parameter (dim=maxdeg+cpts-1)

      integer e, i, ideg, j, k, m, npi
      integer deg(s), eqn(s), sumdeg(0:s), var(s)

      data (deg(i),i=1,s)/1, 1, 1/

      nnz = 0
      ideg = 0
      sumdeg(0) = 0
      do 10 i = 1, s
         ideg = ideg + deg(i)
         sumdeg(i) = sumdeg(i-1) + deg(i)
   10 continue
      npi = s*cpts + ideg

      do 20 i = 1, s
         nnz = nnz + 1
         indrow(nnz) = i
         indcol(nnz) = (i-1)*cpts + sumdeg(i-1) + 1
   20 continue

      do 60 i = 1, nint
         do 50 k = 1, cpts
            do 30 e = 1, s
               var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
               eqn(e) = s + (i-1)*npi + (e-1)*cpts
   30       continue
            do 40 j = 1, cpts + deg(1)
               if (j .ne. 1) then
                  nnz = nnz + 1
                  indrow(nnz) = eqn(1) + k
                  indcol(nnz) = var(1) + j
               end if

               nnz = nnz + 1
               indrow(nnz) = eqn(3) + k
               indcol(nnz) = var(1) + j

               if (j .ne. 1) then
                  nnz = nnz + 1
                  indrow(nnz) = eqn(2) + k
                  indcol(nnz) = var(2) + j
               end if

               nnz = nnz + 1
               indrow(nnz) = eqn(3) + k
               indcol(nnz) = var(2) + j

               nnz = nnz + 1
               indrow(nnz) = eqn(1) + k
               indcol(nnz) = var(3) + j

               nnz = nnz + 1
               indrow(nnz) = eqn(2) + k
               indcol(nnz) = var(3) + j

               if (j .ne. 1) then
                  nnz = nnz + 1
                  indrow(nnz) = eqn(3) + k
                  indcol(nnz) = var(3) + j
               end if
   40       continue
            nnz = nnz + 1
            indrow(nnz) = eqn(3) + k
            indcol(nnz) = n - 2

            nnz = nnz + 1
            indrow(nnz) = eqn(3) + k
            indcol(nnz) = n - 1

            nnz = nnz + 1
            indrow(nnz) = eqn(3) + k
            indcol(nnz) = n

   50    continue
   60 continue

      do 110 i = 1, nint
         do 70 e = 1, s
            var(e) = (i-1)*npi + (e-1)*cpts + sumdeg(e-1)
            eqn(e) = s + (i-1)*npi + s*cpts + sumdeg(e-1)
   70    continue
         if (i .eq. nint) go to 120

         do 100 e = 1, s
            do 90 m = 1, deg(e)
               nnz = nnz + 1
               indrow(nnz) = eqn(e) + m
               indcol(nnz) = var(e) + npi + m
               do 80 j = 1, cpts + deg(e)
                  nnz = nnz + 1
                  indrow(nnz) = eqn(e) + m
                  indcol(nnz) = var(e) + j
   80          continue
   90       continue
  100    continue
  110 continue

  120 continue
      var(1) = n - npi - 3
      var(2) = var(1) + 5
      var(3) = var(2) + 5
      do 130 j = 1, 5
         nnz = nnz + 1
         indrow(nnz) = n - 2
         indcol(nnz) = var(1) + j

         nnz = nnz + 1
         indrow(nnz) = n - 1
         indcol(nnz) = var(2) + j

         nnz = nnz + 1
         indrow(nnz) = n
         indcol(nnz) = var(3) + j

  130 continue

      end
