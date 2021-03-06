module FToPanache
  use iso_c_binding
  implicit none

  type, bind(C) :: C_ShellInfo
    integer(C_INT) :: nprim
    integer(C_INT) :: am
    integer(C_INT) :: ispure
    type(C_PTR) :: exp
    type(C_PTR) :: coef
  end type C_ShellInfo 

  type, bind(C) :: C_AtomCenter
    character(C_CHAR) :: symbol(5)
    real(C_DOUBLE) :: center(3)
  end type C_AtomCenter

  ! For an array of pointer to double precision
  ! see http://stackoverflow.com/questions/8900336/arrays-of-pointers
  type dptr
    double precision, pointer :: ptr(:)
  end type dptr


  interface
    function panache_dfinit(ncenters, atoms , &
                              primary_nshellspercenter, primary_shells, &
                              aux_nshellspercenter, aux_shells, &
                              directory, optflag, bsorder, nthreads) result(res) bind(C, name="panache_dfinit")
      use iso_c_binding
      import C_ShellInfo
      import C_AtomCenter
      implicit none
      integer(C_INT), intent(in), value :: ncenters, optflag, bsorder, nthreads
      integer(C_INT), intent(in) :: primary_nshellspercenter(ncenters), &
                                    aux_nshellspercenter(ncenters)
      type(C_AtomCenter), intent(in) :: atoms(ncenters)
      type(C_ShellInfo), intent(in) :: primary_shells(*), aux_shells(*)
      character(kind=C_CHAR), intent(in) :: directory(*)
      integer(C_INT) :: res
    end function


    function panache_dfinit2(ncenters, atoms , &
                             primary_nshellspercenter, primary_shells, &
                             auxfilename, directory, optflag, bsorder, nthreads) result(res) bind(C, name="panache_dfinit2")
      use iso_c_binding
      import C_ShellInfo
      import C_AtomCenter
      implicit none
      integer(C_INT), intent(in), value :: ncenters, optflag, bsorder, nthreads
      integer(C_INT), intent(in) :: primary_nshellspercenter(ncenters)
      type(C_AtomCenter), intent(in) :: atoms(ncenters)
      type(C_ShellInfo), intent(in) :: primary_shells(*)
      character(kind=C_CHAR), intent(in) :: auxfilename(*), directory(*)
      integer(C_INT) :: res
    end function

    function panache_chinit(ncenters, atoms , &
                            primary_nshellspercenter, primary_shells, &
                            delta, directory, bsorder, nthreads) result(res) bind(C, name="panache_chinit")
      use iso_c_binding
      import C_ShellInfo
      import C_AtomCenter
      implicit none
      integer(C_INT), intent(in), value :: ncenters, bsorder, nthreads
      integer(C_INT), intent(in) :: primary_nshellspercenter(ncenters)
      type(C_AtomCenter), intent(in) :: atoms(ncenters)
      type(C_ShellInfo), intent(in) :: primary_shells(*)
      character(kind=C_CHAR), intent(in) :: directory(*)
      real(C_DOUBLE), intent(in), value :: delta
      integer(C_INT) :: res
    end function
                           

    function panache_getqbatch(handle, tensorflag, outbuf, bufsize, qstart) result(res) bind(C, name="panache_getqbatch")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, tensorflag, bufsize, qstart
      real(C_DOUBLE), intent(out) :: outbuf(bufsize)
      integer(C_INT) :: res
    end function

    function panache_getbatch(handle, tensorflag, outbuf, bufsize, ijstart) result(res) bind(C, name="panache_getbatch")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, tensorflag, bufsize, ijstart
      real(C_DOUBLE), intent(out) :: outbuf(bufsize)
      integer(C_INT) :: res
    end function

    subroutine panache_setcmatrix(handle, cmat, nmo, istrans) bind(C, name="panache_setcmatrix")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, nmo, istrans
      real(C_DOUBLE), intent(in) :: cmat(*)
    end subroutine

    function panache_tensordimensions(handle, tensorflag, naux, ndim1, ndim2) result(res) bind(C, name="panache_tensordimensions")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, tensorflag
      integer(C_INT), intent(out) :: naux, ndim1, ndim2
      integer(C_INT) :: res
    end function

    function panache_calcindex(handle, tensorflag, i, j) result(res) bind(C, name="panache_calcindex")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, tensorflag, i, j
      integer(C_INT) :: res
    end function

    function panache_ispacked(handle, tensorflag) result(res) bind(C, name="panache_ispacked")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, tensorflag
      integer(C_INT) :: res
    end function

    function panache_batchsize(handle, tensorflag) result(res) bind(C, name="panache_batchsize")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, tensorflag
      integer(C_INT) :: res
    end function

    function panache_qbatchsize(handle, tensorflag) result(res) bind(C, name="panache_qbatchsize")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, tensorflag
      integer(C_INT) :: res
    end function

    subroutine panache_genqtensors(handle, qflags, storeflags) bind(C, name="panache_genqtensors")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, qflags, storeflags
    end subroutine

    subroutine panache_delete(handle, qflags) bind(C, name="panache_delete")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, qflags
    end subroutine

    subroutine panache_printtimings(handle) bind(C, name="panache_printtimings")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle
    end subroutine

    subroutine panache_setnthread(handle, nthread, actual) bind(C, name="panache_setnthread")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, nthread
      integer(C_INT), intent(out) :: actual
    end subroutine

    subroutine panache_setnocc(handle, nocc, nfroz) bind(C, name="panache_setnocc")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, nocc, nfroz
    end subroutine

    subroutine panache_stdout() bind(C, name="panache_stdout")
      use iso_c_binding
      implicit none
    end subroutine

    subroutine panache_cleanup(handle) bind(C, name="panache_cleanup")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle
    end subroutine

    subroutine panache_cleanup_all() bind(C, name="panache_cleanup_all")
      use iso_c_binding
      implicit none
    end subroutine

  end interface

  contains
    subroutine ArraysToShellInfo(nshell, nprimpershell, am, ispure, exp, coef, shells, expptr, coefptr)
      implicit none

      integer, intent(in) :: nshell, nprimpershell(nshell), am(nshell), ispure(nshell)
      double precision, intent(in) :: exp(*), coef(*)
      type(C_ShellInfo), intent(inout) :: shells(:)
      type(dptr), intent(inout) :: expptr(:), coefptr(:)
      
      integer :: i, j, primcount, nprim

      !primcount = counts through the primitive array
      primcount = 1
      do i = 1, nshell
        nprim = nprimpershell(i)
        shells(i)%nprim     = INT(nprim, C_INT)
        shells(i)%am        = INT(am(i), C_INT)
        shells(i)%ispure    = INT(ispure(i), C_INT)
    
        allocate(expptr(i)%ptr(nprim))
        allocate(coefptr(i)%ptr(nprim))
    
        do j = 1, nprim
          expptr(i)%ptr(j) = exp(primcount)
          coefptr(i)%ptr(j) = coef(primcount)
          primcount = primcount + 1
        end do
    
        ! tricky
        ! Using the address of the first element gets around
        ! some gfortran 4.8 bugs
        shells(i)%exp  = C_LOC(expptr(i)%ptr(1))
        shells(i)%coef = C_LOC(coefptr(i)%ptr(1))
      end do
    end subroutine

    subroutine ArraysToAtoms(ncenters, symbols, xyz, atoms)
      implicit none

      integer, intent(in) :: ncenters
      double precision, intent(in) :: xyz(3,ncenters)
      character(len=*), intent(in) :: symbols(ncenters)
      type(C_AtomCenter), intent(inout) :: atoms(ncenters)

      integer :: i, j, length

      do i = 1, ncenters
        atoms(i)%center = xyz(:,i)
    
        length = min(4, len_trim(symbols(i)))
    
        do j = 1, length
          atoms(i)%symbol(j) = symbols(i)(j:j)
        end do
    
        ! safe, since symbol has a size of 5 and length
        ! is at most 4
        atoms(i)%symbol(length+1) = C_NULL_CHAR
      end do

    end subroutine

end module FToPanache


!>
!! \brief Clean up a particular calculation and free memory
!!
!! You should not attempt to use the handle afterwards
!!
!! \param [in] handle A handle (returned from an init function) for the
!!                    calculation to be cleaned up
!!
subroutine panachef_cleanup(handle)
  use FToPanache
  implicit none
  integer, intent(in) :: handle
  call panache_cleanup(INT(handle,C_INT))
end subroutine


!>
!! \brief Cleans up all calculations
!!
!! All handles are invalid after this point
!!
subroutine panachef_cleanup_all()
  use FToPanache
  implicit none
  call panache_cleanup_all()
end subroutine


!>
!! \brief Sets the text output of PANACHE to stout
!! 
subroutine panachef_stdout()
  use FToPanache
  implicit none
  call panache_stdout()
end subroutine
  

!>
!! \brief Sets the number of occupied and virtual orbitals.
!!
!! Number of virtual orbitals is taken to be the remainder after the occupied.
!! Used by Qov, etc.
!!
!! Frozen orbitals should be counted in \p nocc.
!!
!! \note You must set the C Matrix first before calling (see panache_setcmatrix())
!!
!! \param [in] handle A handle (returned from an init function) for the calculation 
!! \param [in] nocc Number of (non-frozen) occupied orbitals.
!! \param [in] nfroz Number of frozen occupied orbitals.
!!
subroutine panachef_setnocc(handle, nocc, nfroz)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, nocc, nfroz
  call panache_setnocc(INT(handle, C_INT), &
                       INT(nocc, C_INT),   &
                       INT(nfroz, C_INT))
end subroutine


!>
!! \brief Sets the maximum number of (OpenMP) threads used
!!
!! Set to zero to use the value of the environment variable OMP_NUM_THREAD (or
!! set by omp_num_threads, or the default for this machine).
!!
!! \param [in] handle A handle (returned from an init function) for the calculation
!! \param [in] nthread Max number of threads to use
!! \param [out] actual The max number of threads that will actually be used (ie if \p nthread is zero).
!! 
subroutine panachef_setnthread(handle, nthread, actual)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, nthread
  integer, intent(out) :: actual
  integer(C_INT) :: c_actual
  
  call panache_setnthread(INT(handle, C_INT),  &
                          INT(nthread, C_INT), &
                          c_actual)
  actual = INT(c_actual)
end subroutine


!>
!! \brief Prints out timing information collected so far
!!
!! All times are cumulative for all operations. The output must be set
!!  first (See Output.h)
!!
!! \param [in] handle A handle (returned from an init function) for the calculation
!!
subroutine panachef_printtimings(handle)
  use FToPanache
  implicit none
  integer, intent(in) :: handle
  call panache_printtimings(INT(handle, C_INT))
end subroutine


!>
!! \brief Delete a tensor (from memory, disk, etc)
!!
!! \param [in] handle A handle (returned from an init function) for the calculation
!! \param [in] qflags A combination of flags specifying which tensors to delete
!!
subroutine panachef_delete(handle, qflags)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, qflags
  call panache_delete(INT(handle, C_INT), INT(qflags, C_INT))
end subroutine


!>
!! \brief Generates various 3-index tensors
!!
!! For the \p qflags and \p storeflags parameters, see Flags.h. For example, to calculate the
!! Qmo and Qov tensors on disk,
!!
!! \code{.f90}
!! call panachef_genqtensors(handle, 10, 4)
!! \endcode
!!
!! Default is QSTORAGE_INMEM and not to store with QSTORAGE_BYQ
!!
!! To calculate just Qso, set \p qflags = QGEN_QSO
!!
!! \note The Qso matrix is always stored with QSTORAGE_BYQ
!! \warning Be sure to set the C-Matrix first and number of occupied orbitals first
!!          if qflags contains more than QGEN_QSO
!!
!! \note The Qso matrix is always stored with QSTORAGE_BYQ
!! \note Be sure to set the C-Matrix first!
!!
!! \param [in] handle A handle (returned from an init function) for the calculation 
!! \param [in] qflags A combination of flags specifying which tensors to generate
!! \param [in] storeflags How to store the matrix
!!
subroutine panachef_genqtensors(handle, qflags, storeflags)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, qflags, storeflags
  call panache_genqtensors(INT(handle, C_INT), &
                           INT(qflags, C_INT), &
                           INT(storeflags, C_INT))
end subroutine


!>
!! \brief Obtain the batch size for panachef_getqbatch()
!!
!! The size will be as follows
!!
!! | Tensor | Packed Size      | Unpacked Size |
!! |--------|------------------|---------------|
!! | Qso    | nso*(nso+1)/2    |               |
!! | Qmo    | nmo*(nmo+1)/2    |               |
!! | Qoo    | nocc*(nocc+1)/2  |               |
!! | Qov    |                  | nocc*nvir     |
!! | Qvv    | nvir*(nvir+1)/2  |               |
!!
!! \param [in] handle A handle (returned from an init function) for the calculation 
!! \param [in] tensorflag Which tensor to query (see Flags.h)
!! \param [out] batchsize of batches returned by panachef_getqbatch
!!
subroutine panachef_qbatchsize(handle, tensorflag, batchsize)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag
  integer, intent(out) :: batchsize
  batchsize = panache_qbatchsize(INT(handle, C_INT), INT(tensorflag, C_INT))
end subroutine


!>
!! \brief Obtain the batch size for panache_getbatch()
!!
!! The size will always be naux (number of auxiliary basis functions)
!!
!! \param [in] handle A handle (returned from an init function) for the calculation 
!! \param [in] tensorflag Which tensor to query (see Flags.h)
!! \param [out] batchsize Size of batches returned by panache_getbatch()
!!
subroutine panachef_batchsize(handle, tensorflag, batchsize)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag
  integer, intent(out) :: batchsize
  integer(C_INT) :: c_batchsize
 
  c_batchsize = panache_batchsize(INT(handle, C_INT), INT(tensorflag, C_INT))
  batchsize = INT(c_batchsize)

end subroutine


!>
!! \brief See if a particular tensor is stored packed
!!
!! \param [in] handle A handle (returned from an init function) for the calculation 
!! \param [in] tensorflag Which tensor to query (see Flags.h)
!! \param [out] ispacked Nonzero if the tensor is in packed storage
!!
subroutine panachef_ispacked(handle, tensorflag, ispacked)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag
  integer, intent(out) :: ispacked
  integer(C_INT) :: c_ispacked

  c_ispacked = panache_ispacked(INT(handle, C_INT), INT(tensorflag, C_INT))
  ispacked = INT(c_ispacked)

end subroutine


!>
!! \brief Calculate a combined orbital index
!!
!! Depends on packing
!! 
!! \note This function takes into account 1-based indexing from fortran, so you don't
!!       have to (ie this function subtracts 1 from *i and *j)
!! \param [in] handle A handle (returned from an init function) for the calculation 
!! \param [in] tensorflag Which tensor to query (see Flags.h)
!! \param [in] i First orbital index
!! \param [in] j Second orbital index
!! \param [out] ij Combined orbital index, depending on packing
!!
subroutine panachef_calcindex(handle, tensorflag, i, j, ij)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag, i, j
  integer, intent(out) :: ij
  integer(C_INT) :: c_ij

  c_ij = panache_calcindex(INT(handle, C_INT), INT(tensorflag, C_INT), &
                           INT(i-1, C_INT), INT(j-1, C_INT))
  ij = INT(c_ij)
end subroutine


!>
!! \brief Obtain the dimensions of a tensor
!!
!! \param [in] handle A handle (returned from an init function) for the calculation 
!! \param [in] tensorflag Which tensor to query (see Flags.h)
!! \param [out] naux Number of auxiliary basis functions
!! \param [out] ndim1 First dimension for a particular q
!! \param [out] ndim2 Second dimension for a particular q
!! \param [out] total Total tensor size (depends on packing)
!!
subroutine panachef_tensordimensions(handle, tensorflag, naux, ndim1, ndim2, total)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag
  integer, intent(out) :: naux, ndim1, ndim2, total
  integer(C_INT) :: c_naux, c_ndim1, c_ndim2, c_total

  c_total = panache_tensordimensions(INT(handle, C_INT), INT(tensorflag, C_INT), &
                                     c_naux, c_ndim1, c_ndim2)
  naux = INT(c_naux)
  ndim1 = INT(c_ndim1)
  ndim2 = INT(c_ndim2)
  total = INT(c_total)

end subroutine


!>
!! \brief Sets the C matrix (so-ao matrix) for use in generating Qmo, etc
!!
!! The matrix is expected be nso x nmo (MOs in the columns) in column-major order.
!! If it is nmo x nso, or the matrix is in row major order, set \p cmo_is_trans.
!!
!! \note This is different from panache_setcmatrix(), as fortran column-major order
!! matrices are handled automatically.
!!
!! The matrix is copied by the PANACHE code, so it can be safely deleted or otherwise
!! changed after calling this function.
!!
!! \param [in] handle A handle (returned from an init function) for the calculation
!! \param [in] cmat Pointer to a nso x nmo matrix representing the MO coefficients
!! \param [in] nmo Number of MOs in this C matrix
!! \param [in] istrans Set to non-zero if the matrix is the transpose (nmo x nso) or
!!                          is in row-major order.
!!
subroutine panachef_setcmatrix(handle, cmat, nmo, istrans)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, nmo, istrans
  double precision, intent(in) :: cmat(*)

  if(istrans > 0) then
    call panache_setcmatrix(INT(handle, C_INT), cmat, INT(nmo, C_INT), INT(0, C_INT))
  else
    call panache_setcmatrix(INT(handle, C_INT), cmat, INT(nmo, C_INT), INT(1, C_INT))
  end if
end subroutine


!>
!! \brief Retrieves a batch of a 3-index tensor 
!!
!! See \ref theory_page for what these tensors actually are, and memory_sec for more information
!! about memory.
!!
!! This function returns the number of batches it has stored in the buffer. The buffer
!! will contain (number of batches)*naux elements with the combined orbital index
!! as the slowest index.
!!
!! The batchsize can be obtained using BatchSize()
!!
!! Call this and process the batches, incrementing qstart by the return value,
!! until this function returns zero.
!!
!! \note Tensors are always in row-major order!
!!
!! \param [in] handle A handle (returned from an init function) for the calculation
!! \param [in] tensorflag Which tensor to get (see Flags.h)
!! \param [in] outbuf Memory location to store the tensor
!! \param [in] bufsize The size of \p outbuf (in number of doubles)
!! \param [in] qstart The starting value of q
!! \param [out] nbatch The number of batches actually stored in the buffer.
!!
subroutine panachef_getqbatch(handle, tensorflag, outbuf, bufsize, qstart, nbatch)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag, bufsize, qstart
  integer, intent(out) :: nbatch
  double precision, intent(out) :: outbuf(bufsize)
  integer(C_INT) :: c_nbatch

  c_nbatch = panache_getqbatch(INT(handle, C_INT), INT(tensorflag, C_INT), &
                               outbuf, INT(bufsize, C_INT), INT(qstart, C_INT))

  nbatch = INT(c_nbatch)
end subroutine
    

!>
!! \brief Retrieves a batch of a 3-index tensor 
!!
!! See \ref theory_page for what these tensors actually are, and memory_sec for more information
!! about memory.
!!
!! This function returns the number of batches it has stored in the buffer. The buffer
!! will contain (number of batches)*naux elements with the combined orbital index
!! as the slowest index.
!!
!! The batchsize can be obtained using BatchSize()
!!
!! Call this and process the batches, incrementing qstart by the return value,
!! until this function returns zero.
!!
!! \note Tensors are always in row-major order!
!!
!! \param [in] handle A handle (returned from an init function) for the calculation
!! \param [in] tensorflag Which tensor to get (see Flags.h)
!! \param [in] outbuf Memory location to store the tensor
!! \param [in] bufsize The size of \p outbuf (in number of doubles)
!! \param [in] ijstart The starting value of q
!! \param [out] nbatch The number of batches actually stored in the buffer.
!!
subroutine panachef_getbatch(handle, tensorflag, outbuf, bufsize, ijstart, nbatch)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag, bufsize, ijstart
  integer, intent(out) :: nbatch
  double precision, intent(out) :: outbuf(bufsize)
  integer(C_INT) :: c_nbatch

  c_nbatch = panache_getbatch(INT(handle, C_INT), INT(tensorflag, C_INT), &
                              outbuf, INT(bufsize, C_INT), INT(ijstart, C_INT))

  nbatch = INT(c_nbatch)
end subroutine



!>
!! \brief Initializes a new density-fitting calculation
!!
!! Sets up the basis set information and returns a handle that
!! is used to identify this particular calculation.
!!
!! Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
!!
!! \note Basis set coefficients should NOT be normalized
!!
!! \param [in] ncenters    The number of basis function centers
!! \param [in] xyz         Coordinates of the basis function centers. In order:
!!                         (x1, y1, z1, x2, y2, z2, ..., xN, yN, zN) (ie xyz(1:3, 1:ncenter))
!! \param [in] symbols     Atomic symbols for each center, as a set of \p ncenters strings
!! \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
!!                                       Expected to be of length ncenters.
!! \param [in] primary_am  Angular momentum of each shell (s = 0, p = 1, etc) in the primary basis. 
!!                         Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_is_pure  Whether each shell is pure/spherical or not (primary basis).
!!                              Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_nprimpershell  Number of primitives in each shell of the primary basis. 
!!                                    Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_exp  All exponents for all shells of the primary basis. 
!!                          Length should be the sum of primary_nprimpershell, with grouping
!!                          by shell.
!! \param [in] primary_coef All basis function coefficients for all shells of the primary basis. 
!!                          Length should be the sum of primary_nprimpershell, with grouping
!!                          by shell.
!!
!! \param [in] aux_nshellspercenter  Number of shells on each center for the auxiliary basis.
!!                                   Expected to be of length ncenters.
!! \param [in] aux_am  Angular momentum of each shell (s = 0, p = 1, etc) in the auxiliary basis. 
!!                     Length should be the sum of aux_nshellspercenter.
!! \param [in] aux_is_pure  Whether each shell is pure/spherical or not (auxiliary basis).
!!                          Length should be the sum of aux_nshellspercenter.
!! \param [in] aux_nprimpershell  Number of primitives in each shell of the auxiliary basis. 
!!                                Length should be the sum of aux_nshellspercenter.
!! \param [in] aux_exp  All exponents for all shells of the auxiliary basis. 
!!                      Length should be the sum of aux_nprimpershell, with grouping
!!                      by shell.
!! \param [in] aux_coef All basis function coefficients for all shells of the auxiliary basis. 
!!                      Length should be the sum of aux_nprimpershell, with grouping
!!                      by shell.
!! \param [in] directory A full path to a file to be used if storing matrices to disk.
!!                       Not referenced if the disk is not used. Should not be set to "NULL", but
!!                       may be set to an empty string if disk is not to be used.
!!                       If used, any existing files will be overwritten.
!! \param [in] optflag Flag controlling the type of metric to use
!!                     and other options. Set to zero for default (coulomb/eiginv)
!! \param [in] bsorder Basis function ordering flag
!! \param [in] nthreads Number of threads to use
!!
!! \param [out] handle A handle representing this particular density-fitting calculation.
!!
subroutine panachef_dfinit(ncenters, xyz, symbols, &
                           primary_nshellspercenter, primary_am, primary_is_pure, &
                           primary_nprimpershell, primary_exp, primary_coef, &
                           aux_nshellspercenter, aux_am, aux_is_pure, &
                           aux_nprimpershell, aux_exp, aux_coef, &
                           directory, optflag, bsorder, nthreads, handle) 
  use FToPanache

  implicit none


  integer, intent(in) :: ncenters, optflag, bsorder, nthreads, &
                         primary_nshellspercenter(ncenters), &
                         primary_am(*), primary_is_pure(*), primary_nprimpershell(*), &
                         aux_nshellspercenter(ncenters), &
                         aux_am(*), aux_is_pure(*), aux_nprimpershell(*)

  double precision, intent(in) :: xyz(3,ncenters), &
                                  primary_exp(*), primary_coef(*), &
                                  aux_exp(*), aux_coef(*)
  character(len=*), intent(in) :: symbols(ncenters), directory

  integer, intent(out) :: handle

  type(C_AtomCenter) :: atoms(ncenters)
  type(C_ShellInfo), allocatable :: pshells(:), ashells(:)
  type(dptr), allocatable :: pexpptr(:), pcoefptr(:), aexpptr(:), acoefptr(:)

  character(C_CHAR), allocatable :: directoryarr(:)

  ! for some conversion
  integer(C_INT) :: c_primary_nshellspercenter(ncenters), c_aux_nshellspercenter(ncenters)

  integer :: i, pnshells, anshells

  call ArraysToAtoms(ncenters, symbols, xyz, atoms)

  ! Count the number of shells
  ! pnshells = # of primary shells
  ! anshells = # of aux shells
  pnshells = 0
  anshells = 0
  do i = 1, ncenters
    pnshells = pnshells + primary_nshellspercenter(i)
    anshells = anshells + aux_nshellspercenter(i)

    ! convert integer types
    c_primary_nshellspercenter(i) = INT(primary_nshellspercenter(i), C_INT)
    c_aux_nshellspercenter(i) = INT(aux_nshellspercenter(i), C_INT)
  end do

  

  allocate(pshells(pnshells))
  allocate(pexpptr(pnshells))
  allocate(pcoefptr(pnshells))

  allocate(ashells(anshells))
  allocate(aexpptr(anshells))
  allocate(acoefptr(anshells))


  call ArraysToShellInfo(pnshells, primary_nprimpershell, primary_am, primary_is_pure, &
                         primary_exp, primary_coef, pshells, pexpptr, pcoefptr)
  call ArraysToShellInfo(anshells, aux_nprimpershell, aux_am, aux_is_pure, &
                         aux_exp, aux_coef, ashells, aexpptr, acoefptr)

  ! handle the filenames & paths
  allocate(directoryarr(len_trim(directory)+1))
  directoryarr(len_trim(directory)+1) = C_NULL_CHAR

  do i = 1, len_trim(directory)
    directoryarr(i) = directory(i:i)
  end do


  ! do stuff
  handle = panache_dfinit(INT(ncenters, C_INT), atoms, c_primary_nshellspercenter, &
                          pshells, c_aux_nshellspercenter, ashells, &
                          directoryarr, INT(optflag, C_INT), &
                          INT(bsorder, C_INT), INT(nthreads, C_INT)) 

  do i = 1, pnshells
    deallocate(pexpptr(i)%ptr)
    deallocate(pcoefptr(i)%ptr)
  end do
  deallocate(pshells)

  do i = 1, anshells
    deallocate(aexpptr(i)%ptr)
    deallocate(acoefptr(i)%ptr)
  end do
  deallocate(ashells)

  deallocate(directoryarr) 

  
end subroutine






!> 
!! \brief Initializes a new density-fitting calculation using an auxiliary basis set file
!!
!! Sets up the basis set information and returns a handle that
!! is used to identify this particular calculation.
!!
!! Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
!!
!! \note Basis set coefficients should NOT be normalized
!!
!! \param [in] ncenters    The number of basis function centers
!! \param [in] xyz         Coordinates of the basis function centers. In order:
!!                         (x1, y1, z1, x2, y2, z2, ..., xN, yN, zN) (ie xyz(1:3, 1:ncenter))
!! \param [in] symbols     Atomic symbols for each center, as a set of \p ncenters strings
!! \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
!!                                       Expected to be of length ncenters.
!! \param [in] primary_am  Angular momentum of each shell (s = 0, p = 1, etc) in the primary basis. 
!!                         Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_is_pure  Whether each shell is pure/spherical or not (primary basis).
!!                              Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_nprimpershell  Number of primitives in each shell of the primary basis. 
!!                                    Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_exp  All exponents for all shells of the primary basis. 
!!                          Length should be the sum of primary_nprimpershell, with grouping
!!                          by shell.
!! \param [in] primary_coef All basis function coefficients for all shells of the primary basis. 
!!                          Length should be the sum of primary_nprimpershell, with grouping
!!                          by shell.
!! \param [in] auxfilename A full path to a file containing the auxiliary basis set (in Gaussian94 format)
!!
!! \param [in] directory A full path to a file to be used if storing matrices to disk.
!!                       Not referenced if the disk is not used. Should not be set to "NULL", but
!!                       may be set to an empty string if disk is not to be used.
!! \param [in] optflag Flag controlling the type of metric to use
!!                     and other options. Set to zero for default (coulomb/eiginv)
!! \param [in] bsorder Basis function ordering flag
!! \param [in] nthreads Number of threads to use
!!
!! \param [out] handle A handle representing this particular density-fitting calculation.
!!
subroutine panachef_dfinit2(ncenters, xyz, symbols, &
                            primary_nshellspercenter, primary_am, primary_is_pure, &
                            primary_nprimpershell, primary_exp, primary_coef, &
                            auxfilename, directory, optflag, bsorder, nthreads, handle) 
  use FToPanache

  implicit none


  integer, intent(in) :: ncenters, optflag, bsorder, nthreads, &
                         primary_nshellspercenter(ncenters), &
                         primary_am(*), primary_is_pure(*), primary_nprimpershell(*)
  double precision, intent(in) :: xyz(3,ncenters), primary_exp(*), primary_coef(*)
  character(len=*), intent(in) :: symbols(ncenters), auxfilename, directory
  integer, intent(out) :: handle

  type(C_AtomCenter) :: atoms(ncenters)
  type(C_ShellInfo), allocatable :: pshells(:)
  type(dptr), allocatable :: pexpptr(:), pcoefptr(:)

  character(C_CHAR), allocatable :: auxfilearr(:), directoryarr(:)

  ! for some conversion
  integer(C_INT) :: c_primary_nshellspercenter(ncenters)

  integer :: i, pnshells

  call ArraysToAtoms(ncenters, symbols, xyz, atoms)

  ! Count the number of shells
  ! pnshells = primary # of shells
  pnshells = 0
  do i = 1, ncenters
    pnshells = pnshells + primary_nshellspercenter(i)

    ! convert integer types
    c_primary_nshellspercenter(i) = INT(primary_nshellspercenter(i), C_INT)
  end do

  allocate(pshells(pnshells))
  allocate(pexpptr(pnshells))
  allocate(pcoefptr(pnshells))


  call ArraysToShellInfo(pnshells, primary_nprimpershell, primary_am, primary_is_pure, &
                         primary_exp, primary_coef, pshells, pexpptr, pcoefptr)

  ! handle the filenames & paths
  allocate(auxfilearr(len_trim(auxfilename)+1))
  allocate(directoryarr(len_trim(directory)+1))
  auxfilearr(len_trim(auxfilename)+1) = C_NULL_CHAR
  directoryarr(len_trim(directory)+1) = C_NULL_CHAR

  do i = 1, len_trim(auxfilename)
    auxfilearr(i) = auxfilename(i:i)
  end do
  do i = 1, len_trim(directory)
    directoryarr(i) = directory(i:i)
  end do


  ! do stuff
  handle = panache_dfinit2(INT(ncenters, C_INT), atoms, c_primary_nshellspercenter, &
                           pshells, auxfilearr, directoryarr, INT(optflag, C_INT), &
                           INT(bsorder, C_INT), INT(nthreads, C_INT)) 
   

  do i = 1, pnshells
    deallocate(pexpptr(i)%ptr)
    deallocate(pcoefptr(i)%ptr)
  end do
  deallocate(pshells)

  deallocate(auxfilearr) 
  deallocate(directoryarr) 

  
end subroutine


!> 
!! \brief Initializes a new cholesky calculation
!!
!! Sets up the basis set information and returns a handle that
!! is used to identify this particular calculation.
!!
!! Information passed in is copied, so any dynamic arrays, etc, can be safely deleted afterwards
!!
!! \note Basis set coefficients should NOT be normalized
!!
!! \param [in] ncenters    The number of basis function centers
!! \param [in] xyz         Coordinates of the basis function centers. In order:
!!                         (x1, y1, z1, x2, y2, z2, ..., xN, yN, zN) (ie xyz(1:3, 1:ncenter))
!! \param [in] symbols     Atomic symbols for each center, as a set of \p ncenters strings
!! \param [in] primary_nshellspercenter  Number of shells on each center for the primary basis.
!!                                       Expected to be of length ncenters.
!! \param [in] primary_am  Angular momentum of each shell (s = 0, p = 1, etc) in the primary basis. 
!!                         Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_is_pure  Whether each shell is pure/spherical or not (primary basis).
!!                              Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_nprimpershell  Number of primitives in each shell of the primary basis. 
!!                                    Length should be the sum of primary_nshellspercenter.
!! \param [in] primary_exp  All exponents for all shells of the primary basis. 
!!                          Length should be the sum of primary_nprimpershell, with grouping
!!                          by shell.
!! \param [in] primary_coef All basis function coefficients for all shells of the primary basis. 
!!                          Length should be the sum of primary_nprimpershell, with grouping
!!                          by shell.
!!
!! \param [in] delta     Maximum error in the Cholesky procedure
!! \param [in] directory A full path to a file to be used if storing matrices to disk.
!!                       Not referenced if the disk is not used. Should not be set to "NULL", but
!!                       may be set to an empty string if disk is not to be used.
!! \param [in] bsorder Basis function ordering flag
!! \param [in] nthreads Number of threads to use
!!
!! \param [out] handle A handle representing this particular density-fitting calculation.
!!
subroutine panachef_chinit(ncenters, xyz, symbols, &
                           primary_nshellspercenter, primary_am, primary_is_pure, &
                           primary_nprimpershell, primary_exp, primary_coef, &
                           delta, directory, bsorder, nthreads, handle) 
  use FToPanache

  implicit none


  integer, intent(in) :: ncenters, bsorder, nthreads, &
                         primary_nshellspercenter(ncenters), &
                         primary_am(*), primary_is_pure(*), primary_nprimpershell(*)

  double precision, intent(in) :: xyz(3,ncenters), primary_exp(*), primary_coef(*), delta
  character(len=*), intent(in) :: symbols(ncenters), directory
  integer, intent(out) :: handle

  type(C_AtomCenter) :: atoms(ncenters)
  type(C_ShellInfo), allocatable :: pshells(:)
  type(dptr), allocatable :: pexpptr(:), pcoefptr(:)

  character(C_CHAR), allocatable :: directoryarr(:)

  ! for some conversion
  integer(C_INT) :: c_primary_nshellspercenter(ncenters)

  integer :: i, pnshells

  call ArraysToAtoms(ncenters, symbols, xyz, atoms)

  ! Count the number of shells
  ! pnshells = primary # of shells
  pnshells = 0
  do i = 1, ncenters
    pnshells = pnshells + primary_nshellspercenter(i)

    ! convert integer types
    c_primary_nshellspercenter(i) = INT(primary_nshellspercenter(i), C_INT)
  end do

  allocate(pshells(pnshells))
  allocate(pexpptr(pnshells))
  allocate(pcoefptr(pnshells))


  call ArraysToShellInfo(pnshells, primary_nprimpershell, primary_am, primary_is_pure, &
                         primary_exp, primary_coef, pshells, pexpptr, pcoefptr)

  ! handle the filenames & paths
  allocate(directoryarr(len_trim(directory)+1))
  directoryarr(len_trim(directory)+1) = C_NULL_CHAR

  do i = 1, len_trim(directory)
    directoryarr(i) = directory(i:i)
  end do


  ! do stuff
  handle = panache_chinit(INT(ncenters, C_INT), atoms, c_primary_nshellspercenter, &
                          pshells, delta, directoryarr, &
                          INT(bsorder, C_INT), INT(nthreads, C_INT)) 
   
  do i = 1, pnshells
    deallocate(pexpptr(i)%ptr)
    deallocate(pcoefptr(i)%ptr)
  end do
  deallocate(pshells)

  deallocate(directoryarr) 

  
end subroutine
