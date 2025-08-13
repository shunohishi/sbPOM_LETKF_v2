module mod_read_jra55do

  integer,parameter :: im=640,jm=320
  
contains

  !-------------------------------------------------------------------------
  ! Read JRA55do (Tsujino et al. 2018) |
  !-------------------------------------------------------------------------
  !
  ! Horizontal resolution: ~55 km
  ! Temporal resolution: 3 hours
  !
  !--------------------------------------------------------------------------

  subroutine read_grid(long,lati,land)

    use netcdf
    implicit none

    integer status,ncid,varid
    integer i,j
    
    real(kind = 4) dland(im,jm) !Sea area percentage [%]

    real(kind = 8),intent(out) :: long(im),lati(jm)    
    real(kind = 8),intent(out) :: land(im,jm) !1: Land, 0: Ocean

    character(200) dirname,filename

    dirname="/data/R/R2402/DATA/JRA55do"
    filename=trim(dirname)//"/sftof_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr.nc"

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,long)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,lati)

    status=nf90_inq_varid(ncid,"sftof",varid)
    status=nf90_get_var(ncid,varid,dland)
    
    status=nf90_close(ncid)
    
    !$omp parallel
    !$omp do private(i,j)    
    do j=1,jm
       do i=1,im
          if(dland(i,j) == 100.e0)then
             land(i,j)=0.d0
          else
             land(i,j)=1.d0
          end if
       end do
    end do
    !$omp end do
    !$omp end parallel
    
  end subroutine read_grid

  !---------------------------------------------------------------------------
  
  subroutine read_jra55do(iyr,imon,iday,ihour,u,v,ta,qa,lw,sw,slp,prep)

    use netcdf
    use mod_rmiss
    !use julian_day
    implicit none

    real(kind = 4),parameter :: dmiss=1.e20
    integer,parameter :: nvar=8
    
    integer ijul,ijul0
    integer itime
    integer ivar
    integer status,access
    integer ncid,varid
    integer i,j

    real(kind = 4) :: tmp(im,jm)
    
    character(2) mm,dd
    character(4) yyyy,var
    character(200) dirname,filename

    !IN    
    integer,intent(in) :: iyr,imon,iday,ihour

    !OUT
    real(kind = 8),intent(out) :: u(im,jm),v(im,jm)
    real(kind = 8),intent(out) :: ta(im,jm),qa(im,jm)
    real(kind = 8),intent(out) :: lw(im,jm),sw(im,jm)
    real(kind = 8),intent(out) :: slp(im,jm),prep(im,jm)

    !ITIME
    call ymd_julian(iyr,1,1,ijul0)
    call ymd_julian(iyr,imon,iday,ijul)
    itime=(ijul-ijul0)*8+ihour/3+1

    !FILENAME
    dirname="/data/R/R2402/DATA/JRA55do"
    write(yyyy,'(i4.4)') iyr

    do ivar=1,nvar

       if(ivar == 1) var="uas" !U at 10 m [m s-1]
       if(ivar == 2) var="vas" !V at 10 m [m s-1]
       if(ivar == 3) var="tas" !Ta at 2m [K]
       if(ivar == 4) var="huss" !Qa at 2m [-]
       if(ivar == 5) var="psl"  !SLP [Pa]
       if(ivar == 6) var="rlds" !Downwelling LW [W m-2]
       if(ivar == 7) var="rsds" !Downwelling SW [W m-2]
       if(ivar == 8) var="prra" !Precpitation [kg m-2 s-1]

       if(1 <= ivar .and. ivar <= 5)then
          if(2020 <= iyr .and. iyr <= 2023)then
             filename=trim(dirname)//"/"//trim(var)//&
                  & "_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_"&
                  & //yyyy//"01010000-"//yyyy//"12312100.nc"
          else if(iyr == 2024)then
             filename=trim(dirname)//"/"//trim(var)//&
                  & "_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_"&
                  & //yyyy//"01010000-"//yyyy//"02012100.nc"
          else
             filename=trim(dirname)//"/"//trim(var)//&
                  & "_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_"&
                  & //yyyy//"01010000-"//yyyy//"12312100.nc"
          end if
       else if(6 <= ivar .and. ivar <= 8)then
          if(2020 <= iyr .and. iyr <= 2023)then
             filename=trim(dirname)//"/"//trim(var)//&
                  & "_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_"&
                  & //yyyy//"01010130-"//yyyy//"12312230.nc"
          else if(iyr == 2024)then
             filename=trim(dirname)//"/"//trim(var)//&
                  & "_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_"&
                  & //yyyy//"01010130-"//yyyy//"02012230.nc"
          else
             filename=trim(dirname)//"/"//trim(var)//&
                  & "_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_"&
                  & //yyyy//"01010130-"//yyyy//"12312230.nc"
          end if
       end if

       !Check filename
       status=access(trim(filename)," ")
       if(status == 0)then
          write(*,*) "Read: "//trim(filename)
       else
          write(*,*) "Not found:"//trim(filename)
          stop
       end if

       !Read data
       status=nf90_open(trim(filename),nf90_nowrite,ncid)
       status=nf90_inq_varid(ncid,trim(var),varid)
       status=nf90_get_var(ncid,varid,tmp,(/1,1,itime/),(/im,jm,1/))
       status=nf90_close(ncid)

       !$omp parallel
       !$omp do private(i,j)
       do j=1,jm
          do i=1,im
             if(tmp(i,j) == dmiss)then
                tmp(i,j)=real(rmiss)
             else if(ivar == 3)then !Air T [k] -> [degree C]
                tmp(i,j)=tmp(i,j)-273.15d0
             else if(ivar == 8)then !Precipitation
                ![kg m-2 s-1] / 1000.[kg/m3] = [m s-1]
                ![m s-1] * 1000 [mm/m] = ![mm s-1]
                ![m s-1] * 24*60*60*[s/day] -> [mm day -1]
                !prep(i,j)=24.d0*60.d0*60.d0*dble(tmp(i,j))
                tmp(i,j)=86400.d0*dble(tmp(i,j))
             end if
          end do
       end do
       !$omp end do
       !$omp end parallel   
       
       if(ivar == 1) u(:,:)=dble(tmp(:,:))
       if(ivar == 2) v(:,:)=dble(tmp(:,:))
       if(ivar == 3) ta(:,:)=dble(tmp(:,:))
       if(ivar == 4) qa(:,:)=dble(tmp(:,:))
       if(ivar == 5) slp(:,:)=dble(tmp(:,:))
       if(ivar == 6) lw(:,:)=dble(tmp(:,:))
       if(ivar == 7) sw(:,:)=dble(tmp(:,:))
       if(ivar == 8) prep(:,:)=dble(tmp(:,:))

    end do !ivar
    
  end subroutine read_jra55do
  
end module mod_read_jra55do
