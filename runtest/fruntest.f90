program PanacheRunTest
implicit none

! Counters, etc
integer :: s, i, j, k, prim_count, shell_count
integer :: p_nbf, a_nbf
integer :: nerrors

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
integer :: p_nshells, a_nshells

! Handle to the density fitting stuff
integer :: dfhandle

! Cmo matrix
double precision, allocatable :: cmo(:)

! QSO Matrix
integer :: matsize, expectedsize
double precision, allocatable :: qso(:), qbuf(:)

! From tensordimensions function
integer :: matd1, matd2, matdsize

! Number read from a batch and the current q value
integer :: nread, curq

! Testing a matrix
double precision :: matsum, checksum, expectedsum, expectedchecksum, values(100)
integer :: indices(100)


! Testing basis set reading
character(*) auxpath, mfilename
parameter (auxpath = "/home/ben/programming/psi4/libpanache/basis/cc-pvdz-ri.gbs")
parameter (mfilename = "Test.mat")


interface
    subroutine ReadMolecule(testdir, ncenter, symbols, xyz, Z, masses)
        character(*), intent(in) :: testdir
        integer, intent(out) :: ncenter
        double precision, allocatable, intent(out)  :: xyz(:,:), Z(:), masses(:)
        character(4), allocatable, intent(out) :: symbols(:)
    end subroutine ReadMolecule

    subroutine ReadBasis(testdir, typestr, nshells, nshells_per_center, &
                         nprim_per_shell, am, is_pure, alpha, coef)
        character(*), intent(in) :: testdir, typestr

        integer, allocatable, intent(out) :: nshells_per_center(:), &
                                             nprim_per_shell(:), &
                                             am(:), is_pure(:)

        integer, intent(out) :: nshells

        double precision, allocatable, intent(out) :: alpha(:), coef(:)
    end subroutine ReadBasis
    function TestMatrix(nelements, expected_nelements, &
                        matrix, &
                        matsum, matchecksum, &
                        expected_matsum, expected_matchecksum, &
                        indices, values)

        integer :: TestMatrix
        integer, intent(in) :: nelements, expected_nelements
        double precision, intent(in) :: matrix(*), matsum, matchecksum, expected_matsum, expected_matchecksum
        double precision, intent(in) :: values(100)
        integer, intent(in) :: indices(100)
    end function TestMatrix
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

call ReadBasis(testdir, "primary", p_nshells, p_nshells_per_center, &
               p_nprim_per_shell, p_am, p_is_pure, p_alpha, p_coef)

call ReadBasis(testdir, "aux", a_nshells, a_nshells_per_center, &
               a_nprim_per_shell, a_am, a_is_pure, a_alpha, a_coef)


write(6,*) "Molecule:"
do i = 1,ncenter
  write(6, '(I4,A,5F15.10)') i,symbols(i),Z(i),xyz(i,1:3),masses(i)
end do

! Print out debugging information
p_nbf = 0
a_nbf = 0

prim_count = 1
shell_count = 0
write(6,*)
write(6,*) "Primary basis: nshells = ",p_nshells
do i=1,ncenter
  write(6,*) "Center ",i
  write(6,*) "   nshell: ", p_nshells_per_center(i)  

  do j = 1,p_nshells_per_center(i)
    shell_count = shell_count + 1

    write(6,*) "   Shell ",shell_count
    write(6,*) "     nprim: ", p_nprim_per_shell(shell_count) 
    write(6,*) "     am = ",p_am(shell_count)
    write(6,*) "     is_pure = ",p_is_pure(shell_count)
    if (p_is_pure(shell_count) /= 0) then
      p_nbf = p_nbf + (1+2*p_am(shell_count))
    else
      p_nbf = p_nbf + ((p_am(shell_count)+2)*(p_am(shell_count)+1))/2
    end if

    do k = 1, p_nprim_per_shell(shell_count)
      write(6,'(2F15.10)') p_alpha(prim_count), p_coef(prim_count)
      prim_count = prim_count + 1
    end do

    write(6,*)
  end do
  write(6,*)
end do

prim_count = 1
shell_count = 0
write(6,*)
write(6,*) "Auxiliary basis: nshells = ",a_nshells
do i=1,ncenter
  write(6,*) "Center ",i
  write(6,*) "   nshell: ", a_nshells_per_center(i)  

  do j = 1,a_nshells_per_center(i)
    shell_count = shell_count + 1

    write(6,*) "   Shell ",shell_count
    write(6,*) "     nprim: ", a_nprim_per_shell(shell_count) 
    write(6,*) "     am = ",a_am(shell_count)
    write(6,*) "     is_pure = ",a_is_pure(shell_count)
    if (a_is_pure(shell_count) /= 0) then
      a_nbf = a_nbf + (1+2*a_am(shell_count))
    else
      a_nbf = a_nbf + ((a_am(shell_count)+2)*(a_am(shell_count)+1))/2
    end if

    do k = 1, a_nprim_per_shell(shell_count)
      write(6,'(2F15.10)') a_alpha(prim_count), a_coef(prim_count)
      prim_count = prim_count + 1
    end do

    write(6,*)
  end do
  write(6,*)
end do


matsize = a_nbf * p_nbf**2



call ReadMatrix(testdir, "qso", &
                expectedsize, expectedsum, expectedchecksum, indices, values)


allocate(qso(1:matsize))

call panachef_init(ncenter, xyz, symbols, 4, 0, &
                  p_nshells_per_center, p_am, p_is_pure, &
                  p_nprim_per_shell, p_alpha, p_coef, &
                  a_nshells_per_center, a_am, a_is_pure, &
                  a_nprim_per_shell, a_alpha, a_coef, &
                  mfilename, len(mfilename), &
                  dfhandle)

!call panachef_init2(ncenter, xyz, symbols, 4, 0, &
!                   p_nshells_per_center, p_am, p_is_pure, &
!                   p_nprim_per_shell, p_alpha, p_coef, &
!                   auxpath, len(auxpath), dfhandle)

write(6,*) "DFHANDLE: ",dfhandle

call panachef_qsodimensions(dfhandle, matd1, matd2, matdsize)

write(6,*) "Dimensions: ",matd1, matd2
write(6,*) "Three sizes: ", matsize, expectedsize, matd1*matd2
write(6,*)

allocate(cmo(1:p_nbf**2))
! TODO - read cmo

call panachef_genqso(dfhandle, 1)

!Test reading all at once
!call panachef_getbatch_qso(dfhandle, qso, matsize, nread)

! Test reading in batches
allocate(qbuf(1:5*matd2), stat=s)

call panachef_setoutputbuffer(dfhandle, qbuf, 5*matd2)

curq = 1
call panachef_getbatch_qso(dfhandle, nread)
do while (nread > 0)
  qso((curq-1)*matd2+1:(curq-1+nread)*matd2) = qbuf(:)
  curq = curq + nread
  call panachef_getbatch_qso(dfhandle, nread)
end do

deallocate(qbuf)
deallocate(cmo)

matsum = 0
checksum = 0
do i = 1,matsize
  matsum = matsum + qso(i)
  checksum = checksum + qso(i)*(i)
end do

nerrors = TestMatrix(matsize, expectedsize, &
                     qso, &
                     matsum, checksum, &
                     expectedsum, expectedchecksum, &
                     indices, values)

write(6,*)
write(6,*)
if (nerrors > 0) then
  write(6,*) "OVERALL RESULT: FAIL"
  write(6,'(A,I5,A)') "                (",nerrors," errors)"
else
  write(6,*) "OVERALL RESULT: pass"
end if
write(6,*) 
write(6,*) 

!call TestERI(dfhandle, qso, matsize, p_nshells)


call panachef_cleanup(dfhandle)

deallocate(qso)
deallocate(xyz)
deallocate(symbols)
deallocate(Z)
deallocate(masses)


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

call exit(nerrors)
 
end program PanacheRunTest









subroutine ReadMolecule(testdir, ncenter, symbols, xyz, Z, masses)
implicit none
character(*), intent(in) :: testdir
integer, intent(out) :: ncenter
double precision, allocatable, intent(out)  :: xyz(:,:), Z(:), masses(:)
character(4), allocatable, intent(out) :: symbols(:)


character(256) :: fullpath
character(10) :: symdumm
integer :: ncenterdum

integer :: s, i



fullpath = trim(testdir)//'geometry'
write(6,'(A,A)') "Reading molecule information from ",fullpath


open(unit=11, file=fullpath, status='old',iostat=s)
read(11, *) ncenter,ncenterdum,symdumm,symdumm


allocate(xyz(1:ncenter,1:3), stat=s)
allocate(symbols(1:ncenter), stat=s)
allocate(Z(1:ncenter), stat=s)
allocate(masses(1:ncenter), stat=s)

do i = 1,ncenter
    read(11, *) symbols(i), Z(i), xyz(i,1:3),masses(i)
end do

close(unit=11)

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


character(256) :: fullpath
integer :: ncenter, nprim, ndum


integer :: s, i, j, prim_count



fullpath = trim(testdir)//'basis.'//trim(typestr)
write(6,'(A,A)') "Reading basis set information from ",fullpath


open(unit=11, file=fullpath, status='old',iostat=s)
read(11, *) nshells,ndum,ndum,nprim,ncenter


allocate(nshells_per_center(1:ncenter), stat=s)
allocate(nprim_per_shell(1:nshells), stat=s)
allocate(am(1:nshells), stat=s)
allocate(is_pure(1:nshells), stat=s)
allocate(alpha(1:nprim), stat=s)
allocate(coef(1:nprim), stat=s)


read(11, *) nshells_per_center(1:ncenter)

prim_count = 1

do i = 1,nshells
  ! first number is shell index. Discard it
  read(11,*) ndum,nprim_per_shell(i),am(i),is_pure(i)

  do j = 1,nprim_per_shell(i)
    read(11,*) coef(prim_count),alpha(prim_count)
    prim_count = prim_count + 1 
  end do
end do

close(unit=11)

end subroutine ReadBasis






subroutine ReadMatrix(testdir, typestr, &
                      nelements, elementsum, checksum, indices, values)
implicit none

character(*), intent(in) :: testdir, typestr
integer, intent(out) :: nelements, indices(100)
double precision, intent(out) :: elementsum, checksum, values(100)

integer :: nrow, ncol, i, s
character(256) :: fullpath



fullpath = trim(testdir)//trim(typestr)
write(6,'(A,A)') "Reading matrix information from ",fullpath


open(unit=11, file=fullpath, status='old',iostat=s)

read(11,*) nrow, ncol, elementsum, checksum

nelements = nrow * ncol

do i = 1,100
  read(11,*) indices(i),values(i)
end do

end subroutine ReadMatrix




function TestMatrix(nelements, expected_nelements, &
                    matrix, &
                    matsum, matchecksum, &
                    expected_matsum, expected_matchecksum, &
                    indices, values)
implicit none
integer, intent(in) :: nelements, expected_nelements
double precision, intent(in) :: matrix(*), matsum, matchecksum, expected_matsum, expected_matchecksum
double precision, intent(in) :: values(100)
integer, intent(in) :: indices(100)
        
integer :: TestMatrix
character(20) :: elstr
integer :: i 

interface
    function TestInt(desc, expected, val, tolerance)
        integer :: TestInt
        character(*), intent(in) :: desc
        integer, intent(in) :: expected, val, tolerance
    end function TestInt
    function TestDouble(desc, expected, val, tolerance)
        integer :: TestDouble
        character(*), intent(in) :: desc
        double precision, intent(in) :: expected, val, tolerance
    end function TestDouble
end interface

TestMatrix = 0
 
write(6,'(A20,4A20,A10)') "Test","Value","Expected","Diff","Tolerance","Result" 
write(6,'(A)') "--------------------------------------------------------------------------------------------------------------" 
TestMatrix = TestMatrix + TestInt("# of Elements", expected_nelements, nelements, 0)
TestMatrix = TestMatrix + TestDouble("Sum", matsum, expected_matsum, 1.0d-8)
TestMatrix = TestMatrix + TestDouble("Checksum", matchecksum, expected_matchecksum, 1.0d0)

! Remeber, indices come from C/C++
! So they are zero-based
do i = 1,100
  write(elstr, '(A,I10)') "Element ",indices(i)
  TestMatrix = TestMatrix + TestDouble(elstr, matrix(indices(i)+1), values(i), 1.0d-12)
end do

return
 
end function TestMatrix


function TestInt(desc, expected, val, tolerance)
integer :: TestInt
character(*), intent(in) :: desc
integer, intent(in) :: expected, val, tolerance
integer :: diff

diff = abs(val - expected)

if (diff > tolerance) then
  write(6, 1000) desc, val, expected, diff, tolerance, "FAIL"
  TestInt = 1
else
  write(6, 1000) desc, val, expected, diff, tolerance, "pass"
  TestInt = 0
end if

1000 format(A20,4I20,A10)

return

end function TestInt




function TestDouble(desc, expected, val, tolerance)
integer :: TestDouble
character(*), intent(in) :: desc
double precision, intent(in) :: expected, val, tolerance
double precision :: diff

diff = abs(val - expected)

if (diff > tolerance) then
  write(6, 2000) desc, val, expected, diff, tolerance, "FAIL"
  TestDouble = 1
else
  write(6, 2000) desc, val, expected, diff, tolerance, "pass"
  TestDouble = 0
end if

2000 format(A20,4E20.10,A10)

return

end function TestDouble



!subroutine TestERI(dfhandle, qso, matsize, nshells)
!
!integer, intent(in)  :: matsize, nshells, dfhandle
!double precision, intent(in) :: qso(*)
!
!integer :: i,j,k,l,m,a,b,c,d,e
!integer :: neri
!integer :: s
!
!
!integer :: map(8), nmap(8), nfunc(8)
!
!
!double precision, allocatable :: eribuf(:), hondobuf(:)
!
!integer :: maxam
!
!maxam = 2
!
!
!allocate(eribuf(1:1296), stat=s)
!allocate(hondobuf(1:1296), stat=s)
!
!map(1) = 0
!map(2) = 1
!map(3) = 3
!map(4) = 5
!map(5) = 6
!map(6) = 7
!map(7) = 8
!map(8) = 9
!
!nmap(1) = 1
!nmap(2) = 2
!nmap(3) = 2
!nmap(4) = 1
!nmap(5) = 1
!nmap(6) = 1
!nmap(7) = 1
!nmap(8) = 1
!
!nfunc(1) = 1
!nfunc(2) = 4
!nfunc(3) = 4
!nfunc(4) = 6
!nfunc(5) = 1
!nfunc(6) = 1
!nfunc(7) = 1
!nfunc(8) = 1
!
!do i = 1,8
!do j = 1,8
!do k = 1,8
!do l = 1,8
!
!
!call panachef_eri_multi(dfhandle, qso, matsize, &
!                       map(i), nmap(i), &
!                       map(j), nmap(j), &
!                       map(k), nmap(k), &
!                       map(l), nmap(l), &
!                       eribuf, 1296, neri)
!
!write(6,'(A,4I3,A,I5)') "Quartet: ",i-1,j-1,k-1,l-1, " Elements = ",neri
!
!! Move to hondo buffer 
!hondobuf(:) = 0.0d0
!
!e = 0
!
!!do m = 1,neri
!!  write(6,*) eribuf(m)
!!end do
!
!write(6,*)
!
!write(6,*)
!write(6,*) "HONDO Buffer"
!do a = 1,nfunc(i)
!do b = 1,nfunc(j)
!do c = 1,nfunc(k)
!do d = 1,nfunc(l)
!  e = e + 1
!  hondobuf( (a-1)*216 + (b-1)*36 + (c-1)*6 + d) = eribuf(e)
!end do
!end do
!end do
!end do
!
!!do m = 1,1296
!!    write(6,'(F15.12)') hondobuf(m)
!!end do
!!write(6,*) 
!
!end do
!end do
!end do
!end do
!
!write(6,*)
!
!deallocate(eribuf)
!deallocate(hondobuf)
!
!end subroutine TestERI
