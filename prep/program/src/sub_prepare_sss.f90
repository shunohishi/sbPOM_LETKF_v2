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
!             S.Ohishi 2025.02
!             S.Ohishi 2025.07
!-------------------------------------------------------------------

subroutine prepare_sss(iyr,imon,iday,inum,lon,lat,fsm)

  use setting, only: iswitch_nearshore
  use mod_gridinfo
  use mod_read_smap,only: read_smap_filename => read_filename, read_smap, deallocate_smap, deallocate_smap_filename
  use mod_read_smos,only: read_smos_filename => read_filename, read_smos, deallocate_smos, deallocate_smos_filename
  use mod_rmiss
  implicit none

  !---Common
  integer ihour
  integer jyr,jmon,jday
  integer ifile,nfile
  integer ins
  
  character(200),allocatable :: filename(:)

  !---SMAP
  integer im_smap,jm_smap
  
  real(kind = 8),allocatable :: lon_smap(:,:),lat_smap(:,:)
  real(kind = 8),allocatable :: dat_smap(:,:),nshore_smap(:,:)

  !---SMOS
  integer ngrid_smos
  
  real(kind = 8),allocatable :: lon_smos(:),lat_smos(:)
  real(kind = 8),allocatable :: dat_smos(:),nshore_smos(:)

  
  !---IN
  !Date
  integer,intent(in) :: iyr,imon,iday

  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: fsm(im,jm)

  !---INOUT
  integer,intent(inout) :: inum

  !---SMAP
  ins=12
  write(*,'(a)') "Read SMAP"
  call read_smap_filename(iyr,imon,iday,nfile,filename)

  if(nfile == 0)then
     write(*,'(a)') "No SMAP data"
  else

     do ifile=1,nfile
        
        !Read SMAP
        call read_smap(filename(ifile),jyr,jmon,jday,ihour, &
             & im_smap,jm_smap,lon_smap,lat_smap,dat_smap)    

        if(iyr /= jyr .or. imon /= jmon .or. iday /= jday)then
           write(*,'(a)') "***Error: Inconsistent date in Subroutine read_smap"
           write(*,'(6i6)') iyr,imon,iday,jyr,jmon,jday
           stop
        end if
        
        if(iswitch_nearshore == 1)then
           
           allocate(nshore_smap(im_smap,jm_smap))

           !Detect nearshore SSS
           call detect_nearshore_smap(lon,lat,fsm, &
                & im_smap,jm_smap,lon_smap,lat_smap,nshore_smap)

           !Replace nearshore SSS
           call replace_nearshore_smap(im_smap,jm_smap,nshore_smap,dat_smap,rmiss)

           deallocate(nshore_smap)
           
        end if
           
        !Write SSS
        call write_obs_surface2d("sss",ins,iyr,imon,iday,&
             & im_smap,jm_smap,lon(1),lon(im),lon_smap,lat(1),lat(jm),lat_smap,dat_smap,inum)

        !Deallocate
        call deallocate_smap(lon_smap,lat_smap,dat_smap)
        
     end do !ifile
     
     call deallocate_smap_filename(filename)

  end if
  
  !---SMOS
  ins=11
  write(*,'(a)') "Read SMOS"
  call read_smos_filename(iyr,imon,iday,nfile,filename)

  if(nfile==0)then
     write(*,'(a)') "No SMOS data"
  else
     do ifile=1,nfile

        !Read SMOS
        call read_smos(filename(ifile),jyr,jmon,jday,ihour,&
             & ngrid_smos,lon_smos,lat_smos,dat_smos)
        
        if(ngrid_smos == 0) cycle !For broken files

        if(iyr /= jyr .or. imon /= jmon .or. iday /= jday)then
           write(*,'(a)') "***Error: Inconsistent date in Subroutine read_smos"
           write(*,'(6i6)') iyr,imon,iday,jyr,jmon,jday
           stop
        end if

        if(iswitch_nearshore == 1)then

           allocate(nshore_smos(ngrid_smos))

           !Detect nearshore SSS
           call detect_nearshore_smos(lon,lat,fsm, &
                & ngrid_smos,lon_smos,lat_smos,nshore_smos)

           !Repalce nearshore SSS
           call replace_nearshore_smos(ngrid_smos,nshore_smos,dat_smos,rmiss)

           deallocate(nshore_smos)

        end if
           
        !Write SSS
        call write_obs_surface1d("sss",ins,iyr,imon,iday, &
             & ngrid_smos,lon(1),lon(im),lon_smos,lat(1),lat(jm),lat_smos,dat_smos,inum)

        !Deallocate
        call deallocate_smos(lon_smos,lat_smos,dat_smos)
        
     end do !ifile
     
     call deallocate_smos_filename(filename)
     
  end if

end subroutine prepare_sss

!--------------------------------------------------------------------------

subroutine replace_nearshore_smap(im,jm,nshore,dat,rmiss)

  implicit none

  !---Common
  integer i,j
  
  !---IN
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: nshore(im,jm)
  real(kind = 8),intent(in) :: rmiss
  
  !---INOUT
  real(kind = 8),intent(inout) :: dat(im,jm)

  do j=1,jm
     do i=1,im
        if(nshore(i,j) == 0.d0)then
           dat(i,j)=rmiss
        end if
     end do
  end do
  
end subroutine replace_nearshore_smap

!-------------------------------------------------------------------------

subroutine replace_nearshore_smos(im,nshore,dat,rmiss)

  implicit none

  !---Common
  integer i
  
  !---IN
  integer,intent(in) :: im

  real(kind = 8),intent(in) :: nshore(im)
  real(kind = 8),intent(in) :: rmiss
  
  !---INOUT
  real(kind = 8),intent(inout) :: dat(im)

  do i=1,im
     if(nshore(i) == 0.d0)then
        dat(i)=rmiss
     end if
  end do
  
end subroutine replace_nearshore_smos

