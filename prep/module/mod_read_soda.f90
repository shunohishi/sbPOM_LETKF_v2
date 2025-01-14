module mod_read_soda

  integer,parameter :: im=720,jm=330,km=50

contains

  !------------------------------------------------------------------------
  ! Read SODA |
  !------------------------------------------------------------------------

  subroutine read_soda(iyr,imon,lont,lonu,latt,latu,deptht,depthu,t,s,u,v,ssh)

    !$use omp_lib
    use netcdf
    use mod_rmiss
    implicit none

    real(kind = 8),parameter :: dmiss=-1.d20
    real(kind = 8),parameter :: thresh=-1.d10

    integer i,j,k
    integer status,access
    integer ncid,varid

    real(kind = 4) t_tmp(im,jm,km),s_tmp(im,jm,km)
    real(kind = 4) u_tmp(im,jm,km),v_tmp(im,jm,km)
    real(kind = 4) ssh_tmp(im,jm)

    character(4) yyyy
    character(100) filename

    integer,intent(in) :: iyr,imon

    real(kind = 8),intent(out) :: lont(im),lonu(im)
    real(kind = 8),intent(out) :: latt(jm),latu(jm)
    real(kind = 8),intent(out) :: deptht(km),depthu(km)
    real(kind = 8),intent(out) :: t(im,jm,km),s(im,jm,km)
    real(kind = 8),intent(out) :: u(im,jm,km),v(im,jm,km)
    real(kind = 8),intent(out) :: ssh(im,jm)


    write(yyyy,'(i4.4)') iyr
    filename="/data/R/R2402/DATA/SODA/soda3.12.2_mn_ocean_reg_"//yyyy//".nc"

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not Found: "//trim(filename)
       stop
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !Longitude    
    status=nf90_inq_varid(ncid,"xt_ocean",varid)
    status=nf90_get_var(ncid,varid,lont)

    status=nf90_inq_varid(ncid,"xu_ocean",varid)
    status=nf90_get_var(ncid,varid,lonu)

    !Latitude
    status=nf90_inq_varid(ncid,"yt_ocean",varid)
    status=nf90_get_var(ncid,varid,latt)

    status=nf90_inq_varid(ncid,"yu_ocean",varid)
    status=nf90_get_var(ncid,varid,latu)

    !Depth
    status=nf90_inq_varid(ncid,"st_ocean",varid)
    status=nf90_get_var(ncid,varid,deptht)

    status=nf90_inq_varid(ncid,"sw_ocean",varid)
    status=nf90_get_var(ncid,varid,depthu)

    !Temp
    status=nf90_inq_varid(ncid,"temp",varid)
    status=nf90_get_var(ncid,varid,t_tmp,(/1,1,1,imon/),(/im,jm,km,1/))

    !Salinity
    status=nf90_inq_varid(ncid,"salt",varid)
    status=nf90_get_var(ncid,varid,s_tmp,(/1,1,1,imon/),(/im,jm,km,1/))

    !Zonal Velocity
    status=nf90_inq_varid(ncid,"u",varid)
    status=nf90_get_var(ncid,varid,u_tmp,(/1,1,1,imon/),(/im,jm,km,1/))

    !Meridional velocity
    status=nf90_inq_varid(ncid,"v",varid)
    status=nf90_get_var(ncid,varid,v_tmp,(/1,1,1,imon/),(/im,jm,km,1/))

    !SSH
    status=nf90_inq_varid(ncid,"ssh",varid)
    status=nf90_get_var(ncid,varid,ssh_tmp,(/1,1,imon/),(/im,jm,1/))

    status=nf90_close(ncid)

    !Pre-process
    !$omp parallel
    !$omp do private(i,j)  
    do j=1,jm
       do i=1,im
          if(nint(ssh_tmp(i,j)) == nint(dmiss) .or. ssh_tmp(i,j) < thresh)then
             ssh(i,j)=rmiss
          else
             ssh(i,j)=dble(ssh_tmp(i,j))
          end if
       end do
    end do
    !$omp end do

    !$omp do private(i,j,k)    
    do k=1,km
       do j=1,jm
          do i=1,im

             if(nint(t_tmp(i,j,k)) == nint(dmiss) .or. t_tmp(i,j,k) < thresh)then
                t(i,j,k)=rmiss
             else
                t(i,j,k)=dble(t_tmp(i,j,k))
             end if

             if(nint(s_tmp(i,j,k)) == nint(dmiss) .or. s_tmp(i,j,k) < thresh)then
                s(i,j,k)=rmiss
             else
                s(i,j,k)=dble(s_tmp(i,j,k))
             end if

             if(nint(u_tmp(i,j,k)) == nint(dmiss) .or. u_tmp(i,j,k) < thresh)then
                u(i,j,k)=rmiss
             else
                u(i,j,k)=dble(u_tmp(i,j,k))
             end if

             if(nint(v_tmp(i,j,k)) == nint(dmiss) .or. v_tmp(i,j,k) < thresh)then
                v(i,j,k)=rmiss
             else
                v(i,j,k)=dble(v_tmp(i,j,k))
             end if

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine read_soda

  !----------------------------------------------------------------------------------

    subroutine read_soda_info(iyr,lont,lonu,latt,latu,deptht,depthu)

    !$use omp_lib
    use netcdf
    implicit none

    integer status,access
    integer ncid,varid

    character(4) yyyy
    character(100) filename

    integer,intent(in) :: iyr

    real(kind = 8),intent(out) :: lont(im),lonu(im)
    real(kind = 8),intent(out) :: latt(jm),latu(jm)
    real(kind = 8),intent(out) :: deptht(km),depthu(km)


    write(yyyy,'(i4.4)') iyr
    filename="/data/R/R2402/DATA/SODA/soda3.12.2_mn_ocean_reg_"//yyyy//".nc"

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not Found: "//trim(filename)
       stop
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !Longitude    
    status=nf90_inq_varid(ncid,"xt_ocean",varid)
    status=nf90_get_var(ncid,varid,lont)

    status=nf90_inq_varid(ncid,"xu_ocean",varid)
    status=nf90_get_var(ncid,varid,lonu)

    !Latitude
    status=nf90_inq_varid(ncid,"yt_ocean",varid)
    status=nf90_get_var(ncid,varid,latt)

    status=nf90_inq_varid(ncid,"yu_ocean",varid)
    status=nf90_get_var(ncid,varid,latu)

    !Depth
    status=nf90_inq_varid(ncid,"st_ocean",varid)
    status=nf90_get_var(ncid,varid,deptht)

    status=nf90_inq_varid(ncid,"sw_ocean",varid)
    status=nf90_get_var(ncid,varid,depthu)

  end subroutine read_soda_info
  
  !----------------------------------------------------------------------------------

  subroutine read_soda_2d(varname,iyr,imon,dat)

    !$use omp_lib
    use netcdf
    use mod_rmiss
    implicit none

    real(kind = 8),parameter :: dmiss=-1.d20
    real(kind = 8),parameter :: thresh=-1.d10

    integer i,j
    integer status,access
    integer ncid,varid

    real(kind = 4) tmp(im,jm)

    character(4) yyyy
    character(100) filename

    integer,intent(in) :: iyr,imon
    character(3),intent(in) :: varname
    
    real(kind = 8),intent(out) :: dat(im,jm)

    write(yyyy,'(i4.4)') iyr
    filename="/data/R/R2402/DATA/SODA/soda3.12.2_mn_ocean_reg_"//yyyy//".nc"

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not Found: "//trim(filename)
       stop
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !SSH
    if(varname == "ssh")then
       status=nf90_inq_varid(ncid,"ssh",varid)
    else
       write(*,*) "***Error: Not found "//varname
       stop
    end if
       
    status=nf90_get_var(ncid,varid,tmp,(/1,1,imon/),(/im,jm,1/))

    status=nf90_close(ncid)

    !Pre-process
    !$omp parallel
    !$omp do private(i,j)  
    do j=1,jm
       do i=1,im
          if(nint(tmp(i,j)) == nint(dmiss) .or. tmp(i,j) < thresh)then
             dat(i,j)=rmiss
          else
             dat(i,j)=dble(tmp(i,j))
          end if
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine read_soda_2d

  !---------------------------------------------------------------------------------
  
  subroutine read_soda_3d(varname,iyr,imon,dat)

    !$use omp_lib
    use netcdf
    use mod_rmiss
    implicit none

    real(kind = 8),parameter :: dmiss=-1.d20
    real(kind = 8),parameter :: thresh=-1.d10

    integer i,j,k
    integer status,access
    integer ncid,varid

    real(kind = 4) tmp(im,jm,km)

    character(4) yyyy
    character(100) filename

    integer,intent(in) :: iyr,imon
    character(1),intent(in) :: varname
    
    real(kind = 8),intent(out) :: dat(im,jm,km)

    write(yyyy,'(i4.4)') iyr
    filename="/data/R/R2402/DATA/SODA/soda3.12.2_mn_ocean_reg_"//yyyy//".nc"

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not Found: "//trim(filename)
       stop
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    if(varname == "t")then
       status=nf90_inq_varid(ncid,"temp",varid)
    else if(varname == "s")then
       status=nf90_inq_varid(ncid,"salt",varid)
    else if(varname == "u")then
       status=nf90_inq_varid(ncid,"u",varid)
    else if(varname == "v")then
       status=nf90_inq_varid(ncid,"v",varid)
    else
       write(*,*) "***Error: Not found "//varname
       stop
    end if

    status=nf90_get_var(ncid,varid,tmp,(/1,1,1,imon/),(/im,jm,km,1/))

    status=nf90_close(ncid)

    !Pre-process
    !$omp parallel
    !$omp do private(i,j,k)    
    do k=1,km
       do j=1,jm
          do i=1,im

             if(nint(tmp(i,j,k)) == nint(dmiss) .or. tmp(i,j,k) < thresh)then
                dat(i,j,k)=rmiss
             else
                dat(i,j,k)=dble(tmp(i,j,k))
             end if

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine read_soda_3d
   
end module mod_read_soda
