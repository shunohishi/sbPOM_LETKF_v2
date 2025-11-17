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
! Modified by S. Ohishi 2019.07
!             S. Ohishi 2025.07
!------------------------------------------------------------------

subroutine prepare_ssh(iyr,imon,iday,inum,lon,lat,depth,mdot,fsm)

  use mod_rmiss
  use setting, only: id_ssh
  use mod_gridinfo
  use mod_read_cmems
  implicit none

  !---Common
  integer ifile
  integer :: ins=21
  
  character(200),allocatable :: filename(:)

  !---Model => Obs space
  integer,allocatable :: idx(:),idy(:)
  real(kind = 8),allocatable :: depth_obs(:)
  real(kind = 8),allocatable :: mdot_obs(:),ssh_obs(:)

  !---Observation
  integer nfile
  integer ntime
  real(kind = 8),allocatable :: lon_obs(:),lat_obs(:),ssha_obs(:)

  !---IN
  integer,intent(in) :: iyr,imon,iday

  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(in) :: mdot(im,jm)
  real(kind = 8),intent(in) :: fsm(im,jm)
  
  !---INOUT
  integer,intent(inout) :: inum

  !---SSHA in obs. space
  !Read filename
  write(*,'(a)') "Read SSHA filename"
  call read_filename(iyr,imon,iday,nfile,filename)

  if(nfile == 0)then
     write(*,'(a)') "No satellite SSHA data"
     return
  else
     do ifile=1,nfile
        !Read obs. SSHA
        call read_cmems(filename(ifile),ntime,lon_obs,lat_obs,ssha_obs)
        
        !Allocate
        allocate(idx(ntime),idy(ntime))
        allocate(depth_obs(ntime),mdot_obs(ntime),ssh_obs(ntime))
        
        !ID
        call cal_id(im,lon,ntime,lon_obs,idx)
        call cal_id(jm,lat,ntime,lat_obs,idy)
        
        !Simulated MDOT & Depth in obs. space
        call bilinear_interpolation_1d(im,jm,lon,lat,mdot, &
             & ntime,lon_obs,lat_obs,mdot_obs,idx,idy,rmiss)
        call bilinear_interpolation_1d(im,jm,lon,lat,depth(:,:,km), &
             & ntime,lon_obs,lat_obs,depth_obs,idx,idy,rmiss)
        
        !SSH = MDOT+SSHA
        call make_ssh_obs(ntime,depth_obs,ssha_obs,mdot_obs,ssh_obs,rmiss)

        !Write Data
        call write_obs_surface1d("ssh",ins,iyr,imon,iday,&
             & id_ssh,ntime,lon(1),lon(im),lon_obs,lat(1),lat(jm),lat_obs,ssh_obs,inum)
        
        !Deallocate
        call deallocate_cmems(lon_obs,lat_obs,ssha_obs)
        deallocate(idx,idy)
        deallocate(depth_obs,mdot_obs,ssh_obs)
        
     end do !ifile
     call deallocate_cmems_filename(filename)
  end if

end subroutine prepare_ssh

!------------------------------------------------------------------------------
! Make SSH in obs. space
!------------------------------------------------------------------------------

subroutine make_ssh_obs(ntime,depth,ssha,mdot,ssh,rmiss)

  use setting, only: ssh_depth
  implicit none

  !---Common
  integer itime
  
  !---IN
  integer,intent(in) :: ntime

  real(kind = 8),intent(in) :: depth(ntime)
  real(kind = 8),intent(in) :: ssha(ntime),mdot(ntime)
  real(kind = 8),intent(in) :: rmiss

  !---OUT
  real(kind = 8),intent(out) :: ssh(ntime)
  
  do itime=1,ntime
     if(ssha(itime) == rmiss .or. mdot(itime) == rmiss .or. abs(depth(itime)) < abs(ssh_depth))then
        ssh(itime)=rmiss
     else
        ssh(itime)=ssha(itime)+mdot(itime)
     end if
  end do
  
end subroutine make_ssh_obs
