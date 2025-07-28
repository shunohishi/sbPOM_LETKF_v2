!----------------------------------------------------------------
! Prepare U/V |
!----------------------------------------------------------------
!
! <GCOM-C/OLCI>
! Surface zonal and meridional velocity 
! estimated from motion vector of Chl-a
!
!-----------------------------------------------------------------
! Created by S.Ohishi 2022.12
! Modified by 
!-----------------------------------------------------------------

subroutine prepare_ssuv(iyr,imon,iday,inum,lon,lat,depth,fsm)

  use mod_rmiss
  use mod_gridinfo
  use mod_read_gcomc
  implicit none

  !---Common
  integer :: ins=31
  
  real(kind = 8) ssu(im,jm),ssv(im,jm),mask(im,jm)
  real(kind = 8) pass(im,jm),miss(im,jm)

  !GCOM-C/OLCI
  integer ifile,nfile
  integer io,no

  real(kind = 8),allocatable :: lono(:),lato(:)
  real(kind = 8),allocatable :: uo(:),vo(:)
  character(256),allocatable :: filename(:)

  !---IN
  integer,intent(in) :: iyr,imon,iday

  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(in) :: fsm(im,jm)
  
  !---INOUT
  integer,intent(inout) :: inum

  !Read GCOM-C SSUV
  write(*,'(a)') "Read motion vector"
  call ssuv_filename(iyr,imon,iday,nfile,filename)
  if(nfile == 0)then
     return
  end if

  write(*,'(a)') "Write motion vector"
  do ifile=1,nfile

     call read_gcomc_ssuv(filename(ifile),no,lono,lato,uo,vo)
     if(no == 0) cycle
     
     call write_obs_surface1d("ssu",ins,iyr,imon,iday,no,lon(1),lon(im),lono,lat(1),lat(jm),lato,uo,inum)
     call write_obs_surface1d("ssv",ins,iyr,imon,iday,no,lon(1),lon(im),lono,lat(1),lat(jm),lato,vo,inum)

     call deallocate_gcomc_ssuv(lono,lato,uo,vo)

  end do !ifile
  call deallocate_gcomc_filename(filename)

end subroutine prepare_ssuv
