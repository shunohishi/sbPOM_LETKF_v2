!-------------------------------------------------------------------
! Prepare SSS |
!-------------------------------------------------------------------
!
! 1. Read and Sum SMAP
! (2. Read and Sum SMOS)
! 3. Calculate Mean SSS
! 4. Apply fsm (nshore)
! 5. Write SSS
! 
!-------------------------------------------------------------------
! Created by S.Ohishi 2018.09
! Modified by S.Ohishi 2019.01
! Added nearshore by S.Ohishi 2020.03
! Modified by S.Ohishi 2020.04
!-------------------------------------------------------------------

subroutine prepare_sss(iyr,imon,iday,inum,lon,lat,fsm,nshore)

  use setting, only: iswitch_nearshore
  use mod_gridinfo
  use mod_read_smap,only: read_smap_filename => read_filename, read_smap,deallocate_smap,deallocate_smap_filename
  use mod_read_smos,only: read_smos_filename => read_filename, read_smos,deallocate_smos,deallocate_smos_filename
  implicit none

  !Common
  integer ihour
  integer jyr,jmon,jday
  integer ifile,nfile
  character(200),allocatable :: filename(:)

  !Model
  real(kind = 8) dat(im,jm),count(im,jm)

  !SMAP
  integer im_smap,jm_smap
  real(kind = 8),allocatable :: lon_smap(:,:),lat_smap(:,:)
  real(kind = 8),allocatable :: dat_smap(:,:)

  !SMOS
  integer ngrid_smos
  real(kind = 8),allocatable :: lon_smos(:),lat_smos(:)
  real(kind = 8),allocatable :: dat_smos(:)
  
  !---IN
  integer,intent(in) :: iyr,imon,iday

  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: fsm(im,jm),nshore(im,jm)

  !---INOUT
  integer,intent(inout) :: inum

  !Initilize
  dat(:,:)=0.d0
  count(:,:)=0.d0

  !---SMAP
  !Present day
  write(*,'(a)') "Read & Sum SMAP"
  call read_smap_filename(iyr,imon,iday,nfile,filename)
  if(nfile == 0)then
     write(*,'(a)') "Not exist SMAP data"
  else
     do ifile=1,nfile
        call read_smap(filename(ifile),jyr,jmon,jday,ihour, &
             & im_smap,jm_smap,lon_smap,lat_smap,dat_smap)
        call sum_sss(im_smap,jm_smap,lon_smap,lat_smap,dat_smap, &
             & im,jm,lon,lat,dat,count)
        call deallocate_smap(lon_smap,lat_smap,dat_smap)
     end do !ifile
     call deallocate_smap_filename(filename)
  end if
  
  !---SMOS
  !Present day
  call read_smos_filename(iyr,imon,iday,nfile,filename)
  if(nfile==0)then
     write(*,'(a)') "Not exist SMOS data"
  else
     do ifile=1,nfile
        call read_smos(filename(ifile),jyr,jmon,jday,ihour,&
             & ngrid_smos,lon_smos,lat_smos,dat_smos)
        call sum_sss(ngrid_smos,1,lon_smos,lat_smos,dat_smos, &
             & im,jm,lon,lat,dat,count)
        call deallocate_smos(lon_smos,lat_smos,dat_smos)
     end do
     call deallocate_smos_filename(filename)
  end if

  write(*,*) "Calculate SSS in model space" 
  call mean_sss(im,jm,dat,count)
  write(*,*) "Apply fsm"
  call apply_fsm(im,jm,1,dat,fsm)
  if(iswitch_nearshore == 1) call apply_fsm(im,jm,1,dat,nshore)

  write(*,*) "Write SSS"
  call write_sss(iyr,imon,iday,lon,lat,dat)
  call write_obs_surface("sss",iyr,imon,iday,inum,lon,lat,dat)

end subroutine prepare_sss

!------------------------------------------------------------------------------
! Calculate SSS |
!------------------------------------------------------------------------------

subroutine sum_sss(im_s,jm_s,lon_s,lat_s,dat_s, &
     & im,jm,lon,lat,dat,count)

  use mod_rmiss
  implicit none
  
  integer i_s,j_s
  integer i,j

  real(kind = 8) dx,dy
  
  !---IN
  !Satelllite data
  integer,intent(in) :: im_s,jm_s  
  real(kind = 8),intent(in) :: lon_s(im_s,jm_s),lat_s(im_s,jm_s),dat_s(im_s,jm_s)

  !Model
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: lon(im),lat(jm)

  !---INOUT
  real(kind = 8),intent(inout) :: dat(im,jm),count(im,jm)

  dx=lon(2)-lon(1)
  dy=lat(2)-lat(1)

  do j_s=1,jm_s
     do i_s=1,im_s

        if(lon_s(i_s,j_s) < lon(1)-0.5d0*dx &
             & .or. lon(im)+0.5d0*dx < lon_s(i_s,j_s) &
             & .or. lat_s(i_s,j_s) < lat(1)-0.5d0*dy & 
             & .or. lat(jm)+0.5d0*dy < lat_s(i_s,j_s))cycle
        if(dat_s(i_s,j_s) == rmiss)cycle

        do j=1,jm
           if(lat_s(i_s,j_s) < lat(j)-0.5d0*dy .or. lat(j)+0.5d0*dy <= lat_s(i_s,j_s))cycle
           do i=1,im
              if(lon_s(i_s,j_s) < lon(i)-0.5d0*dx .or. lon(i)+0.5d0*dx <= lon_s(i_s,j_s))cycle

              dat(i,j)=dat(i,j)+dat_s(i_s,j_s)
              count(i,j)=count(i,j)+1.d0

           end do !i
        end do !j

     end do !i_s
  end do !j_s

end subroutine sum_sss

!-------------------------------------------------------------

subroutine mean_sss(im,jm,dat,count)

  use mod_rmiss
  implicit none
  
  integer i,j

  !IN
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: count(im,jm)

  !OUT
  real(kind = 8),intent(out) :: dat(im,jm)

  do j=1,jm
     do i=1,im
        if(count(i,j) == 0.d0)then
           dat(i,j)=rmiss
        else
           dat(i,j)=dat(i,j)/count(i,j)
        end if
     end do !j
  end do !i

end subroutine mean_sss

!-----------------------------------------------------------------------
! Write Data |
!-----------------------------------------------------------------------

subroutine write_sss(iyr,imon,iday,lon,lat,dat)

  use mod_gridinfo
  use netcdf
  implicit none

  integer status,access,system
  integer ncid,varid

  character(2) mm
  character(4) yyyy
  character(100) filename

  !IN
  integer,intent(in) :: iyr,imon,iday
  real(kind = 8) lon(im),lat(jm),dat(im,jm)
  
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon

  filename="../obs/sss"//yyyy//mm//".nc"
  
  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_sss(im,jm,iyr,imon,filename)
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)

  status=nf90_inq_varid(ncid,"time",varid)
  status=nf90_put_var(ncid,varid,real(iday-0.5e0))
  
  status=nf90_inq_varid(ncid,"lon",varid)
  status=nf90_put_var(ncid,varid,real(lon))

  status=nf90_inq_varid(ncid,"lat",varid)
  status=nf90_put_var(ncid,varid,real(lat))

  status=nf90_inq_varid(ncid,"sss",varid)
  status=nf90_put_var(ncid,varid,real(dat(1:im,1:jm)),(/1,1,iday/),(/im,jm,1/))

  status=nf90_close(ncid)

end subroutine write_sss
