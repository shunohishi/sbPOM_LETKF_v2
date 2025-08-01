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

  use mod_gridinfo
  use mod_read_smap,only: read_smap_filename => read_filename, read_smap, deallocate_smap, deallocate_smap_filename
  use mod_read_smos,only: im_smos => im, jm_smos => jm, read_smos
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
  real(kind = 8),allocatable :: dat_smap(:,:)

  !---SMOS
  real(kind = 8) lon_smos(im_smos),lat_smos(jm_smos)
  real(kind = 8) dat_smos(im_smos,jm_smos)

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
  !Read SMOS (Ascending)
  call read_smos("A",iyr,imon,iday,lon_smos,lat_smos,dat_smos)
  !Write SSS
  call write_obs_surface2dg("sss",ins,iyr,imon,iday, &
       & im_smos,jm_smos,lon(1),lon(im),lon_smos,lat(1),lat(jm),lat_smos,dat_smos,inum)
  
  !Read SMOS (Descending)
  call read_smos("D",iyr,imon,iday,lon_smos,lat_smos,dat_smos)
  call write_obs_surface2dg("sss",ins,iyr,imon,iday, &
       & im_smos,jm_smos,lon(1),lon(im),lon_smos,lat(1),lat(jm),lat_smos,dat_smos,inum)

end subroutine prepare_sss
