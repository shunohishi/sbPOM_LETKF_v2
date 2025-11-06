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
! Modified by S.Ohishi 2025.11
!-----------------------------------------------------------------

subroutine prepare_ssuv(iyr,imon,iday,inum,lon,lat,depth,fsm)

  use mod_rmiss
  use setting, only: id_ssu, id_ssv
  use mod_gridinfo
  use mod_read_gcomc
  use mod_read_chla, im_c => im, jm_c => jm
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

  !JASMES
  real(kind = 8) lon_c(im_c),lat_c(jm_c)
  real(kind = 8) mask_c(im_c,jm_c),chl_c(im_c,jm_c)
  real(kind = 8),allocatable :: chlo(:)

  !JAMSTEM ==> GCOM-C/OLCI
  integer,allocatable :: idx(:),idy(:)
  
  !---IN
  integer,intent(in) :: iyr,imon,iday

  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(in) :: fsm(im,jm)
  
  !---INOUT
  integer,intent(inout) :: inum

  !---Read GCOM-C SSUV
  write(*,'(a)') "Read motion vector"
  call ssuv_filename(iyr,imon,iday,nfile,filename)
  if(nfile == 0)then
     return
  end if

  !---Read Chla monthly climatology
  call read_chla_clim(imon,lon_c,lat_c,mask_c,chl_c)  
  
  write(*,'(a)') "Write motion vector"
  do ifile=1,nfile

     !---Read GCOM-C SSUV
     call read_gcomc_ssuv(filename(ifile),no,lono,lato,uo,vo)
     if(no == 0) cycle

     !---Bilinear interpolation
     allocate(idx(no),idy(no),chlo(no))
     call cal_id(im_c,lon_c,no,lono,idx)
     call cal_id(jm_c,lat_c,no,lato,idy)
     call bilinear_interpolation_mask_1d(im_c,jm_c,lon_c,lat_c,mask_c,chl_c, &
          & no,lono,lato,chlo,idx,idy,rmiss)
     call remove_low_chla(no,chlo,uo)
     call remove_low_chla(no,chlo,vo)
     
     call write_obs_surface1d("ssu",ins,iyr,imon,iday, &
          & id_ssu,no,lon(1),lon(im),lono,lat(1),lat(jm),lato,uo,inum)
     call write_obs_surface1d("ssv",ins,iyr,imon,iday, &
          & id_ssv,no,lon(1),lon(im),lono,lat(1),lat(jm),lato,vo,inum)

     call deallocate_gcomc_ssuv(lono,lato,uo,vo)
     deallocate(idx,idy,chlo)
     
  end do !ifile
  call deallocate_gcomc_filename(filename)

end subroutine prepare_ssuv

!----------------------------------------------------------------------------------

subroutine remove_low_chla(n,chl,dat)

  use setting, only: chla_limit
  use mod_rmiss
  implicit none

  !---Common
  integer i

  !---IN
  integer,intent(in) :: n
  real(kind = 8),intent(in) :: chl(n)

  !---IN/OUT
  real(kind = 8),intent(inout) :: dat(n)

  do i=1,n
     if(chl(i) < chla_limit)then
        dat(i)=rmiss
     end if
  end do
  
end subroutine remove_low_chla

