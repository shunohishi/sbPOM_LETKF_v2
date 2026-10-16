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
  integer iseed,nseed
  integer,allocatable :: seed(:)
  real(kind = 8) noise(nens)

  integer,intent(in) :: syr,eyr
  integer,intent(in) :: iyr,imon
  integer,intent(in) :: nens
  integer,intent(out) :: iyr_ens(nens)

  call random_seed(size=nseed)
  allocate(seed(nseed))
  do iseed=1,nseed
     seed(iseed)=iyr*100+imon+37*iseed
  end do
  call random_seed(put=seed)
  call random_number(noise)
  deallocate(seed)
  
  do iens=1,nens
     iyr_ens(iens)=syr+int(noise(iens)*dble(eyr-syr+1))
  end do

end subroutine set_ensemble_year

