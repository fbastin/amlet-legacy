      subroutine sfdcsp(nx,ny,nnz,indrow,indcol)
      integer nx, ny, nnz
      integer indrow(*), indcol(*)
c     **********
c
c     Subroutine sfdcsp
c
c     This subroutine defines the sparsity structure of the Jacobian
c     matrix for the flow in a driven cavity problem.
c
c     The subroutine statement is
c
c       subroutine sfdcsp(nx,ny,nnz,indrow,indcol)
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
c       nnz is an integer variable.
c         On entry nnz need not be specified.
c         On exit nnz is set to the number of nonzeros in the
c            Jacobian matrix.
c
c       indrow is an integer array of dimension at least nnz.
c         On entry indrow need not be specified.
c         On exit indrow contains the row indices of the nonzeros
c            in the Jacobian matrix.
c
c       indcol is an integer array of dimension at least nnz.
c         On entry indcol need not be specified.
c         On exit indcol contains the column indices of the nonzeros
c            in the Jacobian matrix.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer i, j, k

      nnz = 0
      do 20 i = 1, ny
         do 10 j = 1, nx
            k = (i-1)*nx + j
            nnz = nnz + 1
            indrow(nnz) = k
            indcol(nnz) = k
            if (i .gt. 2) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k - 2*nx
            end if
            if (j .gt. 2) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k - 2
            end if
            if (j .lt. nx-1) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k + 2
            end if
            if (i .lt. ny-1) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k + 2*nx
            end if
            if (i .gt. 1 .and. j .gt. 1) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k - nx - 1
            end if
            if (i .gt. 1 .and. j .lt. nx) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k - nx + 1
            end if
            if (i .lt. ny .and. j .gt. 1) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k + nx - 1
            end if
            if (i .lt. ny .and. j .lt. nx) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k + nx + 1
            end if
            if (i .ne. 1) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k - nx
            end if
            if (j .ne. 1) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k - 1
            end if
            if (j .ne. nx) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k + 1
            end if
            if (i .ne. ny) then
               nnz = nnz + 1
               indrow(nnz) = k
               indcol(nnz) = k + nx
            end if
   10    continue
   20 continue

      end
