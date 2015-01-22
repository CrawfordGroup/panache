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
write(6,*) "HERE"

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



function CountNBF(ncenter, nshellspercenter, am, ispure) result(nbf)
implicit none
integer, intent(in) :: ncenter, nshellspercenter(*), am(*), ispure(*)
integer :: shell_count, i, j
integer :: nbf

shell_count = 1

do i=1,ncenter
  do j = 1,nshellspercenter(i)
    if (ispure(shell_count) /= 0) then
      nbf = nbf + (1+2*am(shell_count))
    else
      nbf = nbf + ((am(shell_count)+2)*(am(shell_count)+1))/2
    end if

  shell_count = shell_count + 1

  end do
end do

end function CountNBF


subroutine PrintBasis(ncenter, nshellspercenter, nprimpershell, am, ispure, alpha, coef)
implicit none
integer, intent(in) :: ncenter, nshellspercenter(*), nprimpershell(*), am(*), ispure(*)
double precision, intent(in) :: alpha(*), coef(*)

integer :: i, j, k
integer :: shell_count, prim_count, nbf

shell_count = 0
prim_count = 1

do i=1,ncenter
  write(6,*) "Center ",i
  write(6,*) "   nshell: ", nshellspercenter(i)

  do j = 1,nshellspercenter(i)
    shell_count = shell_count + 1

    write(6,*) "   Shell ",shell_count
    write(6,*) "     nprim: ", nprimpershell(shell_count)
    write(6,*) "     am = ",am(shell_count)
    write(6,*) "     is_pure = ",ispure(shell_count)

    do k = 1, nprimpershell(shell_count)
      write(6,'(2F15.10)') alpha(prim_count), coef(prim_count)
      prim_count = prim_count + 1
    end do

    write(6,*)
  end do
  write(6,*)
end do
end subroutine PrintBasis


function ReadCmo(testdir, cmat) result(nmo)
implicit none
character(*), intent(in) :: testdir
double precision, allocatable, intent(out) :: cmat(:)
character(256) :: fullpath

integer :: nao, nmo
integer :: s

fullpath = trim(testdir)//'cmat'
write(6,'(A,A)') "Reading C Matrix from ",fullpath

open(unit=11, file=fullpath, status='old',iostat=s)
read(11, *) nao, nmo

allocate(cmat(1:nao*nmo))

read(11, *) cmat(1:nao*nmo)

close(unit=11)

end function ReadCmo


subroutine ReadOrbEn(testdir, orben)
implicit none
character(*), intent(in) :: testdir
double precision, allocatable, intent(out) :: orben(:)
character(256) :: fullpath

integer :: nmo
integer :: s

fullpath = trim(testdir)//'orben'
write(6,'(A,A)') "Reading Orbital energies from ",fullpath

open(unit=11, file=fullpath, status='old',iostat=s)
read(11, *) nmo

allocate(orben(1:nmo))

read(11, *) orben(1:nmo)

close(unit=11)

end subroutine ReadOrbEn


function ReadNocc(testdir) result(nocc)
implicit none
character(*), intent(in) :: testdir
character(256) :: fullpath

integer :: nocc
integer :: s

fullpath = trim(testdir)//'nocc'
write(6,'(A,A)') "Reading nocc from ",fullpath

open(unit=11, file=fullpath, status='old',iostat=s)
read(11, *) nocc
close(unit=11)

end function ReadNOcc
  
