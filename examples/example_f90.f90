program ExampleF90
implicit none

! Counters, etc
integer :: s, i, j, a, b, ia, jb, ib, ja, q
integer :: p_nbf, a_nbf

! Filenames, etc
character(40) :: testdesc
character(256) :: testdir, fullpath



! Molecule specification
integer ncenter
double precision, allocatable :: xyz(:,:), Z(:), masses(:)
character(4), allocatable :: symbols(:)


! Basis Sets
! p_* = primary
! a_* = auxiliary
integer, allocatable :: p_nshells_per_center(:), &
                        p_nprim_per_shell(:), &
                        p_am(:), p_is_pure(:), &
                        a_nshells_per_center(:), &
                        a_nprim_per_shell(:), &
                        a_am(:), a_is_pure(:)

double precision, allocatable :: p_alpha(:), p_coef(:), &
                                 a_alpha(:), a_coef(:)

! Number of shells, primary & auxiliary
integer :: p_nshells, a_nshells

! Handle to the three index tensor stuff
integer :: handle

! Cmo matrix and orbital energies
double precision, allocatable :: cmat(:), orben(:)
integer :: nmo, nocc, nvir


! QSO buffer and info
double precision, allocatable :: ia_buf(:), jb_buf(:), ib_buf(:), ja_buf(:)
integer :: nbatch

! From tensordimensions function
integer :: matnaux, matd1, matd2, matdtotal, nbatchsize


! For the mp2 energy calculation
double precision :: e2s, e2t, iajb, ibja, denom


interface
  subroutine ReadMolecule(testdir, ncenter, symbols, xyz, Z, masses)
    implicit none
    character(*), intent(in) :: testdir
    integer, intent(out) :: ncenter
    double precision, allocatable, intent(out)  :: xyz(:,:), Z(:), masses(:)
    character(4), allocatable, intent(out) :: symbols(:)
  end subroutine ReadMolecule


  subroutine ReadBasis(testdir, typestr, nshells, nshells_per_center, nprim_per_shell, &
                       am, is_pure, alpha, coef)
    implicit none
    character(*), intent(in) :: testdir, typestr

    integer, allocatable, intent(out) :: nshells_per_center(:), &
                                         nprim_per_shell(:), &
                                         am(:), is_pure(:)
    integer, intent(out) :: nshells
    double precision, allocatable, intent(out) :: alpha(:), coef(:)
  end subroutine ReadBasis

  function CountNBF(ncenter, nshellspercenter, am, ispure) result(nbf)
    implicit none
    integer, intent(in) :: ncenter, nshellspercenter(*), am(*), ispure(*)
    integer :: nbf
  end function CountNBF

  subroutine PrintBasis(ncenter, nshellspercenter, nprimpershell, am, ispure, alpha, coef)
    implicit none
    integer, intent(in) :: ncenter, nshellspercenter(*), nprimpershell(*), am(*), ispure(*)
    double precision, intent(in) :: alpha(*), coef(*)
  end subroutine PrintBasis


  function ReadCmo(testdir, cmat) result(nmo)
    implicit none
    character(*), intent(in) :: testdir
    double precision, allocatable, intent(out) :: cmat(:)
    integer :: nmo
  end function ReadCmo

  function ReadNocc(testdir) result(nocc)
    implicit none 
    character(*), intent(in) :: testdir
    integer :: nocc
  end function ReadNOcc


  subroutine ReadOrbEn(testdir, orben)
    implicit none 
    character(*), intent(in) :: testdir
    double precision, allocatable, intent(out) :: orben(:)
  end subroutine ReadOrbEn



end interface



if (iargc() /= 1) then
  write(6,*) "Error - no path to testfiles given",iargc()
  call exit(1)
end if


call getarg(1,testdir)

testdir = trim(testdir)//"/"


! Read the description file
fullpath = trim(testdir)//'desc'
open(unit=10, file=fullpath, status='old',iostat=s)
read(10,'(A40)') testdesc
close(unit=10)

write(6,*) "Test description: ",testdesc


! Read the molecule/geometry
call ReadMolecule(testdir, ncenter, symbols, xyz, Z, masses) 

write(6,*) "NCENTER: ",ncenter
write(6,*) "SYMBOLS 1: ",symbols

call ReadBasis(testdir, "primary", p_nshells, p_nshells_per_center, &
               p_nprim_per_shell, p_am, p_is_pure, p_alpha, p_coef)

call ReadBasis(testdir, "aux", a_nshells, a_nshells_per_center, &
               a_nprim_per_shell, a_am, a_is_pure, a_alpha, a_coef)


! Print out some information
write(6,*) "Molecule:"
do i = 1,ncenter
  write(6, '(I4,A,5F15.10)') i,symbols(i),Z(i),xyz(i,1:3),masses(i)
end do

! Print out debugging information
p_nbf = CountNBF(ncenter, p_nshells_per_center, p_am, p_is_pure)
a_nbf = CountNBF(ncenter, a_nshells_per_center, a_am, a_is_pure)

write(6,*) "Primary basis"
call PrintBasis(ncenter, p_nshells_per_center, p_nprim_per_shell, &
                p_am, p_is_pure, p_alpha, p_coef)

write(6,*) "Auxiliary basis"
call PrintBasis(ncenter, a_nshells_per_center, a_nprim_per_shell, &
                a_am, a_is_pure, a_alpha, a_coef)

! Read in C matrix and occupations
! Note that in this case, nmo = nbf (no deleted orbitals, etc)
nmo = ReadCmo(testdir, cmat)
nocc = ReadNOcc(testdir)
nvir = nmo - nocc
write(6,'(A,I4)') " NBF: ", p_nbf
write(6,'(A,I4)') " NMO: ", nmo
write(6,'(A,I4)') "NOCC: ", nocc
write(6,'(A,I4)') "NVIR: ", nvir

! Read in orbital energies
call ReadOrbEn(testdir, orben)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize panache calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call panachef_stdout()
call panachef_dfinit(ncenter, xyz, symbols, &
                     p_nshells_per_center, p_am, p_is_pure, &
                     p_nprim_per_shell, p_alpha, p_coef, &
                     a_nshells_per_center, a_am, a_is_pure, &
                     a_nprim_per_shell, a_alpha, a_coef, &
                     "/tmp", 0, 0, 0, handle)

!call panachef_chinit(ncenter, xyz, symbols, &
!                     p_nshells_per_center, p_am, p_is_pure, &
!                     p_nprim_per_shell, p_alpha, p_coef, &
!                     1.0d-5, "/tmp", 0, 0, handle)

! Set the c matrix and occupations
! 1 = this is in fact a transposed matrix
!     with respect to fortran ordering
!     (it is in row major order)
call panachef_setcmatrix(handle, cmat, nmo, 1)
call panachef_setnocc(handle, nocc, 0) ! 0 = no frozen

! 8 = generate Qov
! 8 = store in memory
call panachef_genqtensors(handle, 8, 8)

! 8 = QOV
! Batch size should equal naux
call panachef_batchsize(handle, 8, nbatchsize)
call panachef_tensordimensions(handle, 8, matnaux, matd1, matd2, matdtotal)

write(6,*) "Qov Dimensions: "
write(6,*) matnaux, matd1, matd2
write(6,*) "Batchsize: ", nbatchsize

! allocate buffers for ia, jb, ib, ja
allocate(ia_buf(1:nbatchsize))
allocate(jb_buf(1:nbatchsize))
allocate(ib_buf(1:nbatchsize))
allocate(ja_buf(1:nbatchsize))

nbatch = 0
e2s = 0.0d0
e2t = 0.0d0

do i = 1, nocc
do a = 1, nvir

  ! nbatch isn't used, but should be 1
  call panachef_calcindex(handle, 8, i, a, ia)
  call panachef_getbatch(handle, 8, ia_buf, nbatchsize, ia, nbatch)

  do j = 1, nocc
  do b = 1, nvir

    call panachef_calcindex(handle, 8, j, b, jb)
    call panachef_getbatch(handle, 8, jb_buf, nbatchsize, jb, nbatch)

    call panachef_calcindex(handle, 8, i, b, ib)
    call panachef_calcindex(handle, 8, j, a, ja)
    
    call panachef_getbatch(handle, 8, ib_buf, nbatchsize, ib, nbatch)
    call panachef_getbatch(handle, 8, ja_buf, nbatchsize, ja, nbatch)
    
    iajb = 0.0d0
    ibja = 0.0d0
    
    do q = 1, nbatchsize
      iajb = iajb + ia_buf(q) * jb_buf(q)
      ibja = ibja + ib_buf(q) * ja_buf(q)
    end do
    
    denom = orben(i) + orben(j) - orben(nocc+a) - orben(nocc+b)
    
    e2s = e2s + (iajb**2)/denom
    e2t = e2t + (iajb*(iajb-ibja))/denom

  end do
  end do
end do
end do


write(6,*) "--------------------------------------------"
write(6,*) " MP2 Correlation energy"
write(6,'(A,F12.8)') "      Same spin: ", e2t
write(6,'(A,F12.8)') "  Opposite spin: ", e2s
write(6,'(A,F12.8)') "          TOTAL: ", e2s+e2t

call panachef_printtimings(handle)


call panachef_cleanup_all()



deallocate(ia_buf)
deallocate(jb_buf)
deallocate(ib_buf)
deallocate(ja_buf)
deallocate(xyz)
deallocate(symbols)
deallocate(Z)
deallocate(masses)

deallocate(cmat)
deallocate(orben)

deallocate(p_nshells_per_center)
deallocate(p_nprim_per_shell)
deallocate(p_am)
deallocate(p_is_pure)
deallocate(p_alpha)
deallocate(p_coef)
deallocate(a_nshells_per_center)
deallocate(a_nprim_per_shell)
deallocate(a_am)
deallocate(a_is_pure)
deallocate(a_alpha)
deallocate(a_coef)

end program ExampleF90

