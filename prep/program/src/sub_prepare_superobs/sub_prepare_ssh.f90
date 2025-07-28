!------------------------------------------------------------------
! Prepare SSH |
!------------------------------------------------------------------
!
! 1. Calculate SSHA within -sshrange [day] and ssh_range-1 [day]
! 2. Calculate SSH (SSHA+MDOT)
! 3. Apply fsm
! 4. Write data
!
!-----------------------------------------------------------------
! Created by S.Ohishi 2018.09
! add when ssh_range=0 by S. Ohishi 2019.07
!------------------------------------------------------------------

subroutine prepare_ssh(iyr,imon,iday,inum,lon,lat,depth,mdot,fsm)

  use mod_rmiss
  use mod_gridinfo
  use mod_read_cmems
  implicit none

  !Common
  integer ijul
  integer jyr,jmon,jday
  integer idt,idt1,idt2
  character(200),allocatable :: filename(:)

  !Model
  real(kind = 8) ssha(im,jm),count(im,jm)
  real(kind = 8) ssh(im,jm)

  !Observation
  integer ifile,nfile
  integer ntime
  real(kind = 8),allocatable :: lon_obs(:),lat_obs(:),ssh_obs(:)

  !---IN
  integer,intent(in) :: iyr,imon,iday

  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(in) :: mdot(im,jm)
  real(kind = 8),intent(in) :: fsm(im,jm)
  
  !---INOUT
  integer,intent(inout) :: inum
  
  call ymd_julian(iyr,imon,iday,ijul)
  ssha(:,:)=0.d0
  count(:,:)=0.d0

  write(*,*) "Calculate SSHA"
  call read_filename(iyr,imon,iday,nfile,filename)

  if(nfile == 0)then
     write(*,'(a)') "Not exist satellite SSH data"
  else
     do ifile=1,nfile
        call read_cmems(filename(ifile),ntime,lon_obs,lat_obs,ssh_obs)
        call sum_ssha(ntime,lon_obs,lat_obs,ssh_obs,lon,lat,depth,ssha,count)
        call deallocate_cmems(lon_obs,lat_obs,ssh_obs)
     end do
  end if
  call deallocate_cmems_filename(filename)

  write(*,*) "Calculate SSH"
  call cal_ssh(mdot,ssha,count,ssh)

  write(*,*) "Apply fsm"
  call apply_fsm(im,jm,1,ssh,fsm)
  call apply_fsm(im,jm,1,ssha,fsm)

  write(*,*) "Write Data"
  call write_ssh(iyr,imon,iday,lon,lat,mdot,ssh,ssha)
  call write_obs_surface("ssh",iyr,imon,iday,inum,lon,lat,ssh)

end subroutine prepare_ssh

!-----------------------------------------------------------------------------
! Calculate SSH |
!-----------------------------------------------------------------------------

subroutine sum_ssha(ntime,lon_obs,lat_obs,ssh_obs,lon,lat,depth,ssha,count)

  use setting, only: ssh_depth
  use mod_rmiss
  use mod_gridinfo
  implicit none
  
  integer i,j
  integer itime
  real dx,dy

  !---IN
  !Observation
  integer,intent(in) :: ntime
  real(kind = 8),intent(in) :: lon_obs(ntime),lat_obs(ntime),ssh_obs(ntime)

  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)

  !---INOUT
  real(kind = 8),intent(inout) :: ssha(im,jm),count(im,jm)
  
  dx=lon(2)-lon(1)
  dy=lat(2)-lat(1)

  do itime=1,ntime

     if(lon_obs(itime)-0.5d0*dx < lon(1) &
          & .or. lon(im)+0.5d0*dx < lon_obs(itime) &
          & .or. lat_obs(itime)-0.5d0*dy < lat(1) &
          & .or. lat(jm)+0.5d0*dy < lat_obs(itime))cycle
     if(ssh_obs(itime) == rmiss)cycle

     do j=1,jm
        if(lat(j)-0.5d0*dy <= lat_obs(itime) &
             & .and. lat_obs(itime) < lat(j)+0.5d0*dy)then
           do i=1,im
              if(depth(i,j,km) < ssh_depth)cycle
              if(lon(i)-0.5d0*dx <= lon_obs(itime) &
                   & .and.  lon_obs(itime) < lon(i)+0.5d0*dx)then
                 ssha(i,j)=ssha(i,j)+ssh_obs(itime)
                 count(i,j)=count(i,j)+1.d0
              end if
           end do !i
        end if
     end do !j

  end do !itime

end subroutine sum_ssha

!-------------------------

subroutine cal_ssh(mdot,ssha,count,ssh)

  use mod_rmiss
  use mod_gridinfo
  implicit none

  integer i,j

  !IN
  real(kind = 8),intent(in) :: mdot(im,jm)
  real(kind = 8),intent(in) :: count(im,jm)

  !INOUT
  real(kind = 8),intent(inout) :: ssha(im,jm)
  
  !OUT
  real(kind = 8),intent(out) :: ssh(im,jm)

  do j=1,jm
     do i=1,im
        if(count(i,j) == 0.d0)then
           ssh(i,j)=rmiss
           ssha(i,j)=rmiss
        else
           ssha(i,j)=ssha(i,j)/count(i,j)
           ssh(i,j)=mdot(i,j)+ssha(i,j)
        end if
     end do
  end do

end subroutine cal_ssh

!------------------------------------------------------------------------
! Write SSH |
!------------------------------------------------------------------------

subroutine write_ssh(iyr,imon,iday,lon,lat,mdot,ssh,ssha)

  use mod_gridinfo
  use netcdf
  implicit none

  integer status,access,system
  integer ncid,varid

  !IN
  integer,intent(in) :: iyr,imon,iday
  
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: mdot(im,jm),ssh(im,jm),ssha(im,jm)

  character(2) mm,dd
  character(4) yyyy
  character(100) filename

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename="../obs/ssh"//yyyy//mm//".nc"
  
  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_ssh(im,jm,iyr,imon,filename)
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)

  status=nf90_inq_varid(ncid,"time",varid)
  status=nf90_put_var(ncid,varid,real(iday-0.5e0))
  
  status=nf90_inq_varid(ncid,"lon",varid)
  status=nf90_put_var(ncid,varid,real(lon))

  status=nf90_inq_varid(ncid,"lat",varid)
  status=nf90_put_var(ncid,varid,real(lat))

  status=nf90_inq_varid(ncid,"mdot",varid)
  status=nf90_put_var(ncid,varid,real(mdot),(/1,1,iday/),(/im,jm,1/))

  status=nf90_inq_varid(ncid,"ssh",varid)
  status=nf90_put_var(ncid,varid,real(ssh),(/1,1,iday/),(/im,jm,1/))

  status=nf90_inq_varid(ncid,"ssha",varid)
  status=nf90_put_var(ncid,varid,real(ssha),(/1,1,iday/),(/im,jm,1/))

  status=nf90_close(ncid)

end subroutine write_ssh
