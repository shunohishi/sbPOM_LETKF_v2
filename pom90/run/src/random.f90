!---------------------------------------------------------------------
!
! syr: dataset start year
! eyr: dataset end year
! iyr, imon: fix random seed in a month
! iens: Ensemble member
! iyr_ens: different random year
!
!---------------------------------------------------------------------

subroutine set_ensemble_year(syr,eyr,iyr,imon,nens,iyr_ens)

  implicit none

  integer iens
  integer nseed
  integer,allocatable :: seed(:)
  real(kind = 8) noise(nens)
  
  integer,intent(in) :: syr,eyr
  integer,intent(in) :: iyr,imon
  integer,intent(in) :: nens
  integer,intent(out) :: iyr_ens(nens)



  call random_seed(size=nseed)
  allocate(seed(nseed))
  call random_seed(get=seed)
  seed(:)=iyr+imon
  call random_seed(put=seed)
  call random_number(noise)

  do iens=1,nens
     iyr_ens(iens)=syr+nint(noise(iens)*dble(eyr-syr))
  end do

end subroutine set_ensemble_year

