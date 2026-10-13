!----------------------------------------------------------------
! Prepare U/V |
!----------------------------------------------------------------
!
! <GCOM-C/OLCI>
! Surface zonal and meridional velocity 
! estimated from motion vector of Chl-a
!
! 1.
! 2.
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
  call ssuv_filename(iyr,imon,iday,nfile,filename)
  if(nfile == 0)then
     ssu(:,:)=rmiss
     ssv(:,:)=rmiss
     mask(:,:)=0.d0
  else
     write(*,'(a)') "Calcaulte GCOM-C SSU/V"
     do ifile=1,nfile
        
        call read_gcomc_ssuv(filename(ifile),no,lono,lato,uo,vo)        

        call cal_ssuv(ifile,nfile,no,lono,lato,uo,vo, &
             & im,jm,lon,lat,ssu,ssv,mask,pass,miss)

        call deallocate_gcomc_ssuv(lono,lato,uo,vo)
        
     end do !ifile
     call deallocate_gcomc_filename(filename)
  end if

  !Apply
  write(*,'(a)') "Apply fsm"
  call apply_fsm(im,jm,1,ssu,fsm)
  call apply_fsm(im,jm,1,ssv,fsm)

  !Write data
  write(*,*) "Write Data"
  call write_ssuv(iyr,imon,iday,lon,lat,ssu,ssv)
  call write_obs_surface_uv("ssu",iyr,imon,iday,inum,lon,lat,ssu,mask)
  call write_obs_surface_uv("ssv",iyr,imon,iday,inum,lon,lat,ssv,mask)

end subroutine prepare_ssuv

!--------------------------------------------------------------
! Calcaulate GCOM-C SSU/V |
!--------------------------------------------------------------

subroutine cal_ssuv(ifile,nfile,no,lono,lato,uo,vo, &
     & im,jm,lon,lat,ssu,ssv,mask,pass,miss)

  use mod_rmiss
  implicit none

  !Obs.
  integer io

  integer,intent(in) :: ifile,nfile
  integer,intent(in) :: no
  real(kind = 8),intent(in) :: lono(no),lato(no)
  real(kind = 8),intent(in) :: uo(no),vo(no)

  !Model
  integer i,j
  real(kind = 8) dx,dy

  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(inout) :: ssu(im,jm),ssv(im,jm),mask(im,jm)
  real(kind = 8),intent(inout) :: pass(im,jm),miss(im,jm)

  dx=lon(2)-lon(1)
  dy=lat(2)-lat(1)

  if(ifile == 1)then
     ssu(:,:)=0.d0
     ssv(:,:)=0.d0
     pass(:,:)=0.d0
     miss(:,:)=0.d0
  end if

  do io=1,no

     if(lono(io) == rmiss .or. lato(io) == rmiss) cycle

     do j=1,jm

        if(lato(io) < lat(j)-0.5d0*dy .or. lat(j)+0.5d0*dy <= lato(io)) cycle

        do i=1,im

           if(lono(io) < lon(i)-0.5d0*dx .or. lon(i)+0.5d0*dx <= lono(io)) cycle

           ssu(i,j)=ssu(i,j)+uo(io)
           ssv(i,j)=ssv(i,j)+vo(io)
           pass(i,j)=pass(i,j)+1.d0
           
        end do !i
     end do !j
  end do !io

  if(ifile == nfile)then
     do j=1,jm
        do i=1,im
           if(pass(i,j) == 0.d0)then
              ssu(i,j)=rmiss
              ssv(i,j)=rmiss
              mask(i,j)=0.d0
           else
              ssu(i,j)=ssu(i,j)/pass(i,j)
              ssv(i,j)=ssv(i,j)/pass(i,j)
              mask(i,j)=1.d0
           end if
        end do
     end do
  end if

end subroutine cal_ssuv

!----------------------------------------------------------------------
! Write SST |
!----------------------------------------------------------------------

subroutine write_ssuv(iyr,imon,iday,lon,lat,ssu,ssv)

  use netcdf
  use mod_gridinfo
  implicit none

  integer status,access,system
  integer ncid,varid

  integer,intent(in) :: iyr,imon,iday
  
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: ssu(im,jm),ssv(im,jm)

  character(2) mm
  character(4) yyyy
  character(100) filename

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon

  filename="../obs/ssuv"//yyyy//mm//".nc"
  
  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_ssuv(im,jm,iyr,imon,filename)
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)

  status=nf90_inq_varid(ncid,"time",varid)
  status=nf90_put_var(ncid,varid,real(iday-0.5e0))
  
  status=nf90_inq_varid(ncid,"lon",varid)
  status=nf90_put_var(ncid,varid,real(lon))

  status=nf90_inq_varid(ncid,"lat",varid)
  status=nf90_put_var(ncid,varid,real(lat))

  status=nf90_inq_varid(ncid,"ssu",varid)
  status=nf90_put_var(ncid,varid,real(ssu),(/1,1,iday/),(/im,jm,1/))

  status=nf90_inq_varid(ncid,"ssv",varid)
  status=nf90_put_var(ncid,varid,real(ssv),(/1,1,iday/),(/im,jm,1/))

  status=nf90_close(ncid)

end subroutine write_ssuv

