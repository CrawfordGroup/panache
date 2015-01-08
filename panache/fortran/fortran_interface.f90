module FToPanache
  use iso_c_binding, only: C_INT, C_DOUBLE, C_PTR
  implicit none

  type, bind(C) :: C_ShellInfo
    integer(C_INT) :: nprim
    integer(C_INT) :: am
    integer(C_INT) :: ispure
    type(C_PTR) :: exp
    type(C_PTR) :: coef
  end type C_ShellInfo 

  type, bind(C) :: C_AtomCenter
    type(C_PTR) :: symbol
    real(C_DOUBLE) :: center(3)
  end type C_AtomCenter


  interface
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
      integer(C_INT), intent(in), value :: handle, nthread, actual
    end subroutine

    subroutine panache_setnocc(handle, nocc, nvir) bind(C, name="panache_setnocc")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: handle, nocc, nvir
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

!    subroutine X(handle) bind(C, name="X")
!      use iso_c_binding
!      implicit none
!      integer(C_INT), intent(in), value :: handle
!    end subroutine

  end interface

     
  contains

  subroutine GetPtr(scalar_len, scalar, ptr)
    use iso_c_binding
    implicit none
  
    integer, intent(in) :: scalar_len
    character(kind=C_CHAR,len=scalar_len), intent(in), target :: scalar(1)
    character(:,kind=C_CHAR), intent(out), pointer :: ptr
    ptr => scalar(1)
  end subroutine GetPtr

  ! For converting a C null-terminated string to a fortran pointer to string
  subroutine CStringToFPtr(cstr, fstr)
    use iso_c_binding, only : C_PTR, C_CHAR, C_INT, C_F_POINTER
    implicit none

    type(C_PTR), intent(in), value :: cstr
    character(:,kind=C_CHAR), pointer, intent(out) :: fstr
    character(kind=C_CHAR), pointer :: arr(:)

    interface
      function strlen(s) BIND(C, name='strlen')
        use iso_c_binding, only: C_PTR, C_SIZE_T
        implicit none

        type(C_PTR), intent(in), value :: s
        integer(C_SIZE_T) :: strlen
      end function strlen
    end interface

    call c_f_pointer(cstr, arr, [strlen(cstr)])
    call GetPtr(size(arr), arr, fstr)
  end subroutine CStringToFPtr

end module FToPanache


!>
!! \brief Clean up a particular density-fitting calculation and free memory
!!
!! You should not attempt to use the handle afterwards
!!
!! \param [in] handle A handle (returned from an init function) for the DF
!!                    calculation to be cleaned up
!!
subroutine panachef_cleanup(handle)
  use FToPanache
  implicit none
  integer, intent(in) :: handle
  call panache_cleanup(handle)
end subroutine


!>
!! \brief Cleans up all density fitting calculations
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
!! \param [in] df_handle A handle (returned from an init function) for the DF calculation 
!! \param [in] nocc Number of (non-frozen) occupied orbitals.
!! \param [in] nfroz Number of frozen occupied orbitals.
!!
subroutine panachef_setnocc(handle, nocc, nfroz)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, nocc, nfroz
  call panache_setnocc(handle, nocc, nfroz)
end subroutine


!>
!! \brief Sets the maximum number of (OpenMP) threads used
!!
!! Set to zero to use the value of the environment variable OMP_NUM_THREAD (or
!! set by omp_num_threads, or the default for this machine).
!!
!! \param [in] df_handle A handle (returned from an init function) for this DF calculation
!! \param [in] nthread Max number of threads to use
!! \param [out] actual The max number of threads that will actually be used (ie if \p nthread is zero).
!! 
subroutine panachef_setnthread(handle, nthread, actual)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, nthread, actual
  call panache_setnthread(handle, nthread, actual)
end subroutine


!>
!! \brief Prints out timing information collected so far
!!
!! All times are cumulative for all operations. The output must be set
!!  first (See Output.h)
!!
!! \param [in] df_handle A handle (returned from an init function) for this DF calculation
!!
subroutine panachef_printtimings(handle)
  use FToPanache
  implicit none
  integer, intent(in) :: handle
  call panache_printtimings(handle)
end subroutine


!>
!! \brief Delete a tensor (from memory, disk, etc)
!!
!! \param [in] df_handle A handle (returned from an init function) for this DF calculation
!! \param [in] qflags A combination of flags specifying which tensors to delete
!!
subroutine panachef_delete(handle, qflags)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, qflags
  call panache_delete(handle, qflags)
end subroutine


!>
!! \brief Generates various 3-index tensors
!!
!! For the \p qflags and \p storeflags parameters, see Flags.h. For example, to calculate the
!! Qmo and Qov tensors on disk,
!!
!! \code{.f90}
!! call panachef_genqtensors(df_handle, 10, 4)
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
!! \param [in] df_handle A handle (returned from an init function) for the DF calculation 
!! \param [in] qflags A combination of flags specifying which tensors to generate
!! \param [in] storeflags How to store the matrix
!!
subroutine panachef_genqtensors(handle, qflags, storeflags)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, qflags, storeflags
  call panache_genqtensors(handle, qflags, storeflags)
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
!! \param [in] df_handle A handle (returned from an init function) for the DF calculation 
!! \param [in] tensorflag Which tensor to query (see Flags.h)
!! \param [out] batchsize of batches returned by panachef_getqbatch
!!
subroutine panachef_qbatchsize(handle, tensorflag, batchsize)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag
  integer, intent(out) :: batchsize
  batchsize = panache_qbatchsize(handle, tensorflag)
end subroutine


!>
!! \brief Obtain the batch size for panache_getbatch()
!!
!! The size will always be naux (number of auxiliary basis functions)
!!
!! \param [in] df_handle A handle (returned from an init function) for the DF calculation 
!! \param [in] tensorflag Which tensor to query (see Flags.h)
!! \param [out] batchsize Size of batches returned by panache_getbatch()
!!
subroutine panachef_batchsize(handle, tensorflag, batchsize)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag
  integer, intent(out) :: batchsize
  batchsize = panache_batchsize(handle, tensorflag)
end subroutine


!>
!! \brief See if a particular tensor is stored packed
!!
!! \param [in] df_handle A handle (returned from an init function) for the DF calculation 
!! \param [in] tensorflag Which tensor to query (see Flags.h)
!! \param [out] ispacked Nonzero if the tensor is in packed storage
!!
subroutine panachef_ispacked(handle, tensorflag, ispacked)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag
  integer, intent(out) :: ispacked
  ispacked = panache_ispacked(handle, tensorflag)
end subroutine


!>
!! \brief Calculate a combined orbital index
!!
!! Depends on packing
!! 
!! \note This function takes into account 1-based indexing from fortran, so you don't
!!       have to (ie this function subtracts 1 from *i and *j)
!! \param [in] df_handle A handle (returned from an init function) for the DF calculation 
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
  ij = panache_calcindex(handle, tensorflag, i-1, j-1)
end subroutine


!>
!! \brief Obtain the dimensions of a tensor
!!
!! \param [in] df_handle A handle (returned from an init function) for the DF calculation 
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
  total = panache_tensordimensions(handle, tensorflag, naux, ndim1, ndim2)
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
!! \param [in] df_handle A handle (returned from an init function) for this DF calculation
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
    call panache_setcmatrix(handle, cmat, nmo, 0)
  else
    call panache_setcmatrix(handle, cmat, nmo, 1)
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
!! \param [in] df_handle A handle (returned from an init function) for this DF calculation
!! \param [in] tensorflag Which tensor to get (see Flags.h)
!! \param [in] outbuf Memory location to store the tensor
!! \param [in] bufsize The size of \p outbuf (in number of doubles)
!! \param [in] ijstart The starting value of q
!! \param [out] nbatch The number of batches actually stored in the buffer.
!!
subroutine panachef_getqbatch(handle, tensorflag, outbuf, bufsize, qstart, nbatch)
  use FToPanache
  implicit none
  integer, intent(in) :: handle, tensorflag, bufsize, qstart
  integer, intent(out) :: nbatch
  double precision, intent(out) :: outbuf(bufsize)

  nbatch = panache_getqbatch(handle, tensorflag, outbuf, bufsize, qstart)
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
!! \param [in] df_handle A handle (returned from an init function) for this DF calculation
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

  nbatch = panache_getbatch(handle, tensorflag, outbuf, bufsize, ijstart)
end subroutine
