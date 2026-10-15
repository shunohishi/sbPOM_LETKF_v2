module mod_io_common

contains

  !-----------------------------------------------------------------------
  ! Read Argument |
  !-----------------------------------------------------------------------

  subroutine read_argument(syr,smon,sday,eyr,emon,eday)

    implicit none

    !---Common
    integer i,length,status

    character(:),allocatable :: arg

    intrinsic :: command_argument_count, get_command_argument

    !---Out
    integer,intent(out) :: syr,smon,sday
    integer,intent(out) :: eyr,emon,eday

    if(command_argument_count() /= 6)then
       write(*,*) "***Error: Argument size => ",command_argument_count()
       stop
    endif

    do i=1,command_argument_count()

       call get_command_argument(i,length=length,status=status)

       if(status /= 0)then
          write(*,*) "Error: arugument ",status
       else

          allocate(character(length) :: arg)

          call get_command_argument(i,arg,status=status)

          if(i == 1)then
             read(arg,'(i4)') syr
          else if(i == 2)then
             read(arg,'(i2)') smon
          else if(i == 3)then
             read(arg,'(i2)') sday
          else if(i == 4)then
             read(arg,'(i4)') eyr
          else if(i == 5)then
             read(arg,'(i2)') emon
          else if(i == 6)then
             read(arg,'(i2)') eday
          end if

          deallocate(arg)

       end if

    end do

  end subroutine read_argument

  !---------------------------------------------------------------------------------
  ! Read data in obs. space |
  !---------------------------------------------------------------------------------

  subroutine read_hdata(buoyname,varname,idat_a,iyr,imon,iday, &
       & km_o,lon_o,lat_o,dep_o,dat_o,hmean_a,hsprd_a)

    use setting, only: datname
    use mod_rmiss
    use mod_check_netcdf
    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,dimid,varid

    character(100) filename
    character(4) yyyy
    character(2) mm

    !---IN
    integer,intent(in) :: idat_a
    integer,intent(in) :: iyr,imon,iday

    character(10),intent(in) :: buoyname
    character(1),intent(in) :: varname

    !---OUT    
    integer,intent(out) :: km_o

    !---IN/OUT (*for no file case)
    real(kind = 8),intent(inout) :: lon_o,lat_o
    real(kind = 8),allocatable,intent(inout) :: dep_o(:),dat_o(:)
    real(kind = 8),allocatable,intent(inout) :: hmean_a(:),hsprd_a(:)

    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    filename="dat/"//trim(buoyname)//"/"//trim(varname)//"/"//trim(datname(idat_a))//"."//yyyy//mm//".nc"

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read "//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       km_o=0
       return
    end if

    !---Open file
    status=nf90_open(trim(filename),nf90_nowrite,ncid)
    call check_netcdf(status)

    !---Get km_o
    status=nf90_inq_dimid(ncid,"z",dimid)
    call check_netcdf(status)
    status=nf90_inquire_dimension(ncid,dimid,len = km_o)
    call check_netcdf(status)

    !---Allocate
    if(allocated(dep_o)) deallocate(dep_o)
    if(allocated(dat_o)) deallocate(dat_o)
    if(allocated(hmean_a)) deallocate(hmean_a)
    if(allocated(hsprd_a)) deallocate(hsprd_a)
    allocate(dep_o(km_o),dat_o(km_o))
    allocate(hmean_a(km_o),hsprd_a(km_o))

    !---Read data
    status=nf90_inq_varid(ncid,"lon_o",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,lon_o)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat_o",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,lat_o)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"dep_o",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,dep_o)
    call check_netcdf(status)

    !status=nf90_inq_varid(ncid,"pres_o",varid)
    !status=nf90_get_var(ncid,varid,pres_o,(/1,iday/),(/km_o,1/))

    status=nf90_inq_varid(ncid,"h"//trim(varname)//"mean_a",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,hmean_a,(/1,iday/),(/km_o,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"h"//trim(varname)//"sprd_a",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,hsprd_a,(/1,iday/),(/km_o,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,trim(varname)//"_o",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,dat_o,(/1,iday/),(/km_o,1/))
    call check_netcdf(status)

    status=nf90_close(ncid)
    call check_netcdf(status)

  end subroutine read_hdata

  !------------------------------------

  subroutine deallocate_hdata(dep_o,dat_o,hmean_a,hsprd_a)

    implicit none

    real(kind = 8),allocatable,intent(inout) :: dep_o(:),dat_o(:)
    real(kind = 8),allocatable,intent(inout) :: hmean_a(:),hsprd_a(:)

    if(allocated(dep_o)) deallocate(dep_o)
    if(allocated(dat_o)) deallocate(dat_o)
    if(allocated(hmean_a)) deallocate(hmean_a)
    if(allocated(hsprd_a)) deallocate(hsprd_a)

  end subroutine deallocate_hdata

  !---------------------------------------------------------------------------------
  ! Write data in obs. space |
  !---------------------------------------------------------------------------------

  subroutine write_hdata(buoyname,varname,idat_a,iyr,imon,iday, &
       & km_o,lon_o,lat_o,dep_o,pres_o,dat_o,hmean_a,hsprd_a)

    use setting, only: datname
    use mod_make_ncfile
    use mod_check_netcdf
    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,varid

    character(100) filename
    character(4) yyyy
    character(2) mm

    !---IN
    integer,intent(in) :: idat_a
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: km_o

    real(kind = 8),intent(in) :: lon_o,lat_o,dep_o(km_o),pres_o(km_o)
    real(kind = 8),intent(in) :: dat_o(km_o)
    real(kind = 8),intent(in) :: hmean_a(km_o),hsprd_a(km_o)

    character(10),intent(in) :: buoyname
    character(1),intent(in) :: varname

    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    filename="dat/"//trim(buoyname)//"/"//trim(varname)//"/"//trim(datname(idat_a))//"."//yyyy//mm//".nc"

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Output to "//trim(filename)
    else
       write(*,*) "Make & Output to "//trim(filename)
       call make_ncfile(km_o,varname,filename)
    end if

    !---Write data
    status=nf90_open(trim(filename),nf90_write,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lon_o",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,lon_o)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat_o",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,lat_o)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"dep_o",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,dep_o)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"pres_o",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,pres_o,(/1,iday/),(/km_o,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"h"//trim(varname)//"mean_a",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,hmean_a,(/1,iday/),(/km_o,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"h"//trim(varname)//"sprd_a",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,hsprd_a,(/1,iday/),(/km_o,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,trim(varname)//"_o",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,dat_o,(/1,iday/),(/km_o,1/))
    call check_netcdf(status)

    status=nf90_close(ncid)
    call check_netcdf(status)

  end subroutine write_hdata

  !---------------------------------------------------------------------------------
  ! Write data in obs. space |
  !---------------------------------------------------------------------------------

  subroutine write_obs(buoyname,varname,datname,iyr,imon,iday,km_o,dep_o,hmean_a,dat_o)

    implicit none

    !---Common
    integer k

    character(100) format
    character(10) yyyymmdd
    character(4) yyyy
    character(2) mm,dd

    !---IN
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: km_o

    real(kind = 8),intent(in) :: dep_o(km_o)
    real(kind = 8),intent(in) :: hmean_a(km_o),dat_o(km_o)

    character(10),intent(in) :: buoyname,datname
    character(1),intent(in) :: varname    

    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday

    yyyymmdd=yyyy//"-"//mm//"-"//dd

    write(format,'(a)') "(a,3f12.5)" 

    open(11,file="dat/"//trim(buoyname)//"/"//trim(varname)//"/"//trim(datname)//".dat",access="append")
    do k=1,km_o       
       write(11,trim(format)) yyyymmdd,dep_o(k),hmean_a(k),dat_o(k)
    end do
    close(11)

  end subroutine write_obs

  !---------------------------------------------------------------------------------
  ! Write Monthly statistics |
  !---------------------------------------------------------------------------------

  subroutine write_mave(buoyname,varname,syr,eyr,ndat_a,km_o,dep_o, &
       & num_stat_mave,num_sprd_mave,bias_mave,rmsd_mave,sprd_mave)

    implicit none

    !---Common
    integer iyr,imon
    integer k

    character(100) format
    character(10) yyyymmdd

    !---IN
    integer,intent(in) :: syr,eyr
    integer,intent(in) :: ndat_a
    integer,intent(in) :: km_o
    integer,intent(in) :: num_stat_mave(km_o,ndat_a,12,syr:eyr)
    integer,intent(in) :: num_sprd_mave(km_o,ndat_a,12,syr:eyr)

    real(kind = 8),intent(in) :: dep_o(km_o)
    real(kind = 8),intent(in) :: bias_mave(km_o,ndat_a,12,syr:eyr)
    real(kind = 8),intent(in) :: rmsd_mave(km_o,ndat_a,12,syr:eyr)
    real(kind = 8),intent(in) :: sprd_mave(km_o,ndat_a,12,syr:eyr)

    character(10),intent(in) :: buoyname
    character(1),intent(in) :: varname

    write(format,'(a,I0,a,I0,a)') "(a,f12.5,",ndat_a,"i10,",ndat_a,"f12.5)"

    open(11,file="dat/"//trim(buoyname)//"/"//trim(varname)//"bias_mave.dat",status="replace")
    open(12,file="dat/"//trim(buoyname)//"/"//trim(varname)//"rmsd_mave.dat",status="replace")
    open(13,file="dat/"//trim(buoyname)//"/"//trim(varname)//"sprd_mave.dat",status="replace")
    do iyr=syr,eyr
       do imon=1,12

          write(yyyymmdd,'(i4.4,a,i2.2,a)') iyr,"-",imon,"-15"

          do k=1,km_o

             write(11,trim(format)) yyyymmdd,dep_o(k),num_stat_mave(k,:,imon,iyr),bias_mave(k,:,imon,iyr)
             write(12,trim(format)) yyyymmdd,dep_o(k),num_stat_mave(k,:,imon,iyr),rmsd_mave(k,:,imon,iyr)
             write(13,trim(format)) yyyymmdd,dep_o(k),num_sprd_mave(k,:,imon,iyr),sprd_mave(k,:,imon,iyr)

          end do
       end do
    end do
    close(11)
    close(12)
    close(13)             

  end subroutine write_mave

  !---------------------------------------------------------------------------------
  ! Write statistics for whole analysis period |
  !---------------------------------------------------------------------------------

  subroutine write_ave(buoyname,varname,sjul,ejul,ndat_a,km_o,dep_o, &
       & num_stat_ave,num_sprd_ave, &
       & bias_ave,abias_dif_low,abias_dif_ave,abias_dif_upp, &
       & rmsd_ave,rmsd_dif_low,rmsd_dif_ave,rmsd_dif_upp, &
       & sprd_ave,cor_ave)

    implicit none

    !---Common
    integer k
    integer idat_a

    integer num_min(km_o,ndat_a)

    character(100) format

    !---IN
    integer,intent(in) :: sjul,ejul
    integer,intent(in) :: ndat_a
    integer,intent(in) :: km_o
    integer,intent(in) :: num_stat_ave(km_o,ndat_a)
    integer,intent(in) :: num_sprd_ave(km_o,ndat_a)

    real(kind = 8),intent(in) :: dep_o(km_o)
    real(kind = 8),intent(in) :: bias_ave(km_o,ndat_a)
    real(kind = 8),intent(in) :: abias_dif_low(km_o,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: abias_dif_ave(km_o,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: abias_dif_upp(km_o,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_ave(km_o,ndat_a)
    real(kind = 8),intent(in) :: rmsd_dif_low(km_o,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_dif_ave(km_o,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_dif_upp(km_o,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: sprd_ave(km_o,ndat_a)
    real(kind = 8),intent(in) :: cor_ave(km_o,ndat_a)

    character(10),intent(in) :: buoyname
    character(1),intent(in) :: varname

    write(format,'(a,I0,a)') "(f12.5,",ndat_a*4,"f12.5)"

    open(11,file="dat/"//trim(buoyname)//"/"//trim(varname)//"bias_ave.dat",status="replace")
    open(12,file="dat/"//trim(buoyname)//"/"//trim(varname)//"rmsd_ave.dat",status="replace")
    do k=1,km_o

       write(11,trim(format)) dep_o(k),num_stat_ave(k,:)*100.d0/dble(ejul-sjul+1),bias_ave(k,:),abias_dif_low(k,1,:),abias_dif_upp(k,1,:)
       write(12,trim(format)) dep_o(k),num_stat_ave(k,:)*100.d0/dble(ejul-sjul+1),rmsd_ave(k,:),rmsd_dif_low(k,1,:),rmsd_dif_upp(k,1,:)

    end do
    close(11)
    close(12)

    write(format,'(a,I0,a)') "(f12.5,",2*ndat_a,"f12.5)"    

    open(13,file="dat/"//trim(buoyname)//"/"//trim(varname)//"sprd_ave.dat",status="replace")              
    do k=1,km_o
       write(13,trim(format)) dep_o(k),num_sprd_ave(k,:)*100.d0/dble(ejul-sjul+1),sprd_ave(k,:)
    end do
    close(13)             

    write(format,'(a,I0,a)') "(f12.5,",2*ndat_a,"f12.5)"    

    do idat_a=1,ndat_a
       do k=1,km_o
          num_min(k,idat_a)=min(num_stat_ave(k,idat_a),num_sprd_ave(k,idat_a))
       end do
    end do

    open(14,file="dat/"//trim(buoyname)//"/"//trim(varname)//"cor_ave.dat",status="replace")              
    do k=1,km_o
       write(14,trim(format)) dep_o(k),num_min(k,:)*100.d0/dble(ejul-sjul+1),cor_ave(k,:)
    end do
    close(14)             

  end subroutine write_ave

  !---------------------------------------------------------------------------------
  ! Write statistics (bias, RMSD, spread difference) |
  !---------------------------------------------------------------------------------  

  subroutine write_dave(buoyname,varname,ndat_a, &
       & bias_dave, &
       & rmsd_dave,rmsd_dif_dlow,rmsd_dif_dave,rmsd_dif_dupp, &
       & sprd_dave,cor_dave)

    implicit none

    !---Common
    integer idat_a

    character(100) format

    !---IN
    integer,intent(in) :: ndat_a

    real(kind = 8),intent(in) :: bias_dave(ndat_a)
    real(kind = 8),intent(in) :: rmsd_dave(ndat_a)
    real(kind = 8),intent(in) :: rmsd_dif_dlow(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_dif_dave(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_dif_dupp(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: sprd_dave(ndat_a)
    real(kind = 8),intent(in) :: cor_dave(ndat_a)

    character(10),intent(in) :: buoyname
    character(1),intent(in) :: varname

    write(format,'(a,I0,a)') "(",3*ndat_a,"f12.5)"

    open(11,file="dat/"//trim(buoyname)//"/"//trim(varname)//"rmsd_dave.dat",status="replace")
    write(11,trim(format)) rmsd_dave(:),rmsd_dif_dlow(1,:),rmsd_dif_dupp(1,:)
    close(11)

    write(format,'(a,I0,a)') "(",ndat_a,"f12.5)"

    open(11,file="dat/"//trim(buoyname)//"/"//trim(varname)//"bias_dave.dat",status="replace")
    open(12,file="dat/"//trim(buoyname)//"/"//trim(varname)//"sprd_dave.dat",status="replace")
    open(13,file="dat/"//trim(buoyname)//"/"//trim(varname)//"cor_dave.dat",status="replace")
    write(11,trim(format)) bias_dave(:)
    write(12,trim(format)) sprd_dave(:)
    write(13,trim(format)) cor_dave(:)
    close(11)
    close(12)             
    close(13)             

  end subroutine write_dave

end module mod_io_common
