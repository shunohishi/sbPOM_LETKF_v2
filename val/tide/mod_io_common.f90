module mod_io_common

contains

  !--------------------------------------------------------------------------
  ! Read argument
  !--------------------------------------------------------------------------  
  
  subroutine read_argument(syr,smon,sday,eyr,emon,eday)

    implicit none

    !---Common
    integer i,length,status

    character(:),allocatable :: arg

    intrinsic :: command_argument_count, get_command_argument

    !---Out
    integer,intent(out) :: syr,smon,sday
    integer,intent(out) :: eyr,emon,eday

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

  !--------------------------------------------------------------------------
  ! Write data in obs. space
  !--------------------------------------------------------------------------  

  subroutine write_hdat(idat_a,nst,ijul,lon_a,lat_a,hdat_a,hsprd_a, &
       & lon_o,lat_o,dat_o,dist)

    use mod_check_netcdf
    use mod_io, only: get_datname
    use mod_julian
    use mod_make_ncfile
    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,varid
    integer iyr,imon,iday

    character(100) filename
    character(10) datname
    character(4) yyyy
    character(2) mm

    !---IN
    integer,intent(in) :: idat_a
    integer,intent(in) :: nst
    integer,intent(in) :: ijul

    real(kind = 8),intent(in) :: lon_a(nst),lat_a(nst)
    real(kind = 8),intent(in) :: hdat_a(nst),hsprd_a(nst)

    real(kind = 8),intent(in) :: lon_o(nst),lat_o(nst)
    real(kind = 8),intent(in) :: dat_o(nst)

    real(kind = 8),intent(in) :: dist(nst)

    call get_datname(idat_a,datname)

    !---Julian day => iyr,imon,iday
    call julian_ymd(ijul,iyr,imon,iday)

    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    !---Filename
    filename="dat/"//trim(datname)//"."//yyyy//mm//".nc"

    !---Make file
    status=access(trim(filename)," ")
    if(status /= 0)then
       call make_obsfile(nst,filename)
    end if

    status=nf90_open(trim(filename),nf90_write,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lon_a",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,lon_a(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat_a",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,lat_a(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"hmean_a",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,hdat_a(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"hsprd_a",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,hsprd_a(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lon_o",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,lon_o(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat_o",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,lat_o(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"h_o",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,dat_o(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"dist",varid)
    call check_netcdf(status)
    status=nf90_put_var(ncid,varid,dist(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_close(ncid)    
    call check_netcdf(status)

  end subroutine write_hdat

  !---------------------------------------------------------------------
  ! Read data in obs. space |
  !---------------------------------------------------------------------

  subroutine read_hdat(idat_a,nst,ijul,lon_a,lat_a,hdat_a,hsprd_a, &
       & lon_o,lat_o,dat_o,dist)

    use mod_check_netcdf
    use mod_io, only: get_datname
    use mod_julian
    use mod_rmiss
    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,varid
    integer iyr,imon,iday

    character(100) filename
    character(10) datname
    character(4) yyyy
    character(2) mm

    !---IN
    integer,intent(in) :: idat_a
    integer,intent(in) :: nst
    integer,intent(in) :: ijul

    !---OUT
    real(kind = 8),intent(out) :: lon_a(nst),lat_a(nst)
    real(kind = 8),intent(out) :: hdat_a(nst),hsprd_a(nst)

    real(kind = 8),intent(out) :: lon_o(nst),lat_o(nst)
    real(kind = 8),intent(out) :: dat_o(nst)

    real(kind = 8),intent(out) :: dist(nst)

    call get_datname(idat_a,datname)

    !---Julian day => iyr,imon,iday
    call julian_ymd(ijul,iyr,imon,iday)

    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    !---Filename
    filename="dat/"//trim(datname)//"."//yyyy//mm//".nc"

    !---Make file
    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "Not found "//trim(filename)
       lon_a(:)=rmiss
       lat_a(:)=rmiss
       hdat_a(:)=rmiss
       hsprd_a(:)=rmiss
       lon_o(:)=rmiss
       lat_o(:)=rmiss
       dat_o(:)=rmiss
       dist(:)=rmiss
       return
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lon_a",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,lon_a(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat_a",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,lat_a(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"hmean_a",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,hdat_a(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"hsprd_a",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,hsprd_a(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lon_o",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,lon_o(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat_o",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,lat_o(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"h_o",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,dat_o(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"dist",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,dist(1:nst),(/1,iday/),(/nst,1/))
    call check_netcdf(status)

    status=nf90_close(ncid)    
    call check_netcdf(status)

  end subroutine read_hdat

  !---------------------------------------------------------------------
  ! Write time serires |
  !---------------------------------------------------------------------

  subroutine write_ts(ndat_a,nst,sjul,ejul,dat_a,sprd_a,dat_o)

    use mod_julian
    use mod_rmiss
    implicit none

    !---Common
    integer ist
    integer ijul
    integer iyr,imon,iday
    integer pass(nst)

    character(100) format
    character(10) yyyymmdd
    character(4) yyyy
    character(3) nnn
    character(2) mm,dd

    !---IN
    integer,intent(in) :: ndat_a
    integer,intent(in) :: nst
    integer,intent(in) :: sjul,ejul

    real(kind = 8),intent(in) :: dat_a(nst,sjul:ejul,ndat_a),sprd_a(nst,sjul:ejul,ndat_a)
    real(kind = 8),intent(in) :: dat_o(nst,sjul:ejul)

    !--Format
    write(format,'(a,I0,a)') "(a,",2*ndat_a+1,"f12.5)"

    !---Count
    pass(:)=0

    do ist=1,nst
       do ijul=sjul,ejul

          if(dat_o(ist,ijul) == rmiss) cycle
          pass(ist)=pass(ist)+1

       end do
    end do

    !---Write
    do ist=1,nst

       if(pass(ist) == 0) cycle

       write(nnn,'(i3.3)') ist

       open(1,file="dat/"//nnn//".dat",status="replace")
       do ijul=sjul,ejul
          call julian_ymd(ijul,iyr,imon,iday)
          write(yyyy,'(i4.4)') iyr
          write(mm,'(i2.2)') imon
          write(dd,'(i2.2)') iday
          yyyymmdd=yyyy//"-"//mm//"-"//dd
          write(1,trim(format)) yyyymmdd,dat_a(ist,ijul,1:ndat_a),sprd_a(ist,ijul,1:ndat_a),dat_o(ist,ijul)
       end do
       close(1)

    end do !ist

  end subroutine write_ts

  !---------------------------------------------------------------------
  ! Write Statistics (Monthly) |
  !---------------------------------------------------------------------

  subroutine write_stat_mave(nst,ndat_a,syr,eyr,lon,lat,num_stat,num_sprd,bias,rmsd,sprd)

    implicit none

    !---Common
    integer ist
    integer iyr,imon

    character(10) yyyymmdd
    character(4) yyyy
    character(2) mm
    character(100) format

    !---IN
    integer,intent(in) :: nst
    integer,intent(in) :: ndat_a
    integer,intent(in) :: syr,eyr

    integer,intent(in) :: num_stat(nst,ndat_a,12,syr:eyr)
    integer,intent(in) :: num_sprd(nst,ndat_a,12,syr:eyr)

    real(kind = 8),intent(in) :: lon(nst),lat(nst)
    real(kind = 8),intent(in) :: bias(nst,ndat_a,12,syr:eyr)
    real(kind = 8),intent(in) :: rmsd(nst,ndat_a,12,syr:eyr)
    real(kind = 8),intent(in) :: sprd(nst,ndat_a,12,syr:eyr)

    write(format,'(a,I0,a,I0,a)') "(i6,2f12.5,x,a,",ndat_a,"i10,",ndat_a,"f12.5)"

    open(1,file="dat/bias_mave.dat",status="replace")
    open(2,file="dat/rmsd_mave.dat",status="replace")
    open(3,file="dat/sprd_mave.dat",status="replace")
    do ist=1,nst
       do iyr=syr,eyr

          write(yyyy,'(i4.4)') iyr          

          do imon=1,12

             write(mm,'(i2.2)') imon
             yyyymmdd=yyyy//"-"//mm//"-15"

             write(1,trim(format)) ist,lon(ist),lat(ist),yyyymmdd,num_stat(ist,:,imon,iyr),bias(ist,:,imon,iyr)
             write(2,trim(format)) ist,lon(ist),lat(ist),yyyymmdd,num_stat(ist,:,imon,iyr),rmsd(ist,:,imon,iyr)
             write(3,trim(format)) ist,lon(ist),lat(ist),yyyymmdd,num_sprd(ist,:,imon,iyr),sprd(ist,:,imon,iyr)
          end do
       end do
    end do
    close(1)
    close(2)
    close(3)

  end subroutine write_stat_mave

  !---------------------------------------------------------------------
  ! Write Statistics (ALL period) |
  !---------------------------------------------------------------------

  subroutine write_stat_ave(nst,ndat_a,lon,lat,num_stat,num_sprd,bias,rmsd,sprd,cor, &
       & abias_dif_low,abias_dif_ave,abias_dif_upp, &
       & rmsd_dif_low,rmsd_dif_ave,rmsd_dif_upp)

    implicit none

    !---Common
    integer ist
    integer idat_a

    integer num_min(nst,ndat_a)

    character(100) format1,format2

    !---IN
    integer,intent(in) :: nst
    integer,intent(in) :: ndat_a

    integer,intent(in) :: num_stat(nst,ndat_a)
    integer,intent(in) :: num_sprd(nst,ndat_a)

    real(kind = 8),intent(in) :: lon(nst),lat(nst)
    real(kind = 8),intent(in) :: bias(nst,ndat_a)
    real(kind = 8),intent(in) :: rmsd(nst,ndat_a)    
    real(kind = 8),intent(in) :: sprd(nst,ndat_a)
    real(kind = 8),intent(in) :: cor(nst,ndat_a)

    real(kind = 8),intent(in) :: abias_dif_low(nst,ndat_a,ndat_a),abias_dif_ave(nst,ndat_a,ndat_a),abias_dif_upp(nst,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_dif_low(nst,ndat_a,ndat_a),rmsd_dif_ave(nst,ndat_a,ndat_a),rmsd_dif_upp(nst,ndat_a,ndat_a)

    do idat_a=1,ndat_a
       do ist=1,nst
          num_min(ist,idat_a)=min(num_stat(ist,idat_a),num_sprd(ist,idat_a))
       end do
    end do

    write(format1,'(a,I0,a,I0,a)') "(i6,2f12.5,",ndat_a,"i10,",ndat_a*3,"f12.5)"
    write(format2,'(a,I0,a,I0,a)') "(i6,2f12.5,",ndat_a,"i10,",ndat_a,"f12.5)"

    open(1,file="dat/bias_ave.dat",status="replace")
    open(2,file="dat/rmsd_ave.dat",status="replace")
    open(3,file="dat/sprd_ave.dat",status="replace")
    open(4,file="dat/cor_ave.dat",status="replace")
    do ist=1,nst
       write(1,trim(format1)) ist,lon(ist),lat(ist),num_stat(ist,:),bias(ist,:),abias_dif_low(ist,1,:),abias_dif_upp(ist,1,:)
       write(2,trim(format1)) ist,lon(ist),lat(ist),num_stat(ist,:),rmsd(ist,:),rmsd_dif_low(ist,1,:),rmsd_dif_upp(ist,1,:)
       write(3,trim(format2)) ist,lon(ist),lat(ist),num_sprd(ist,:),sprd(ist,:)
       write(4,trim(format2)) ist,lon(ist),lat(ist),num_min(ist,:),cor(ist,:)       
    end do
    close(1)
    close(2)
    close(3)
    close(4)

  end subroutine write_stat_ave

  !---------------------------------------------------------------------
  ! Write Statistics (Station & Time average) |
  !---------------------------------------------------------------------

  subroutine write_stat_ave_all(ndat_a,bias,rmsd,sprd,cor1,cor2, &
       & abias_dif_low,abias_dif_ave,abias_dif_upp, &
       & rmsd_dif_low,rmsd_dif_ave,rmsd_dif_upp)

    implicit none

    !---Common
    character(100) format1,format2,format3

    !---IN
    integer,intent(in) :: ndat_a

    real(kind = 8),intent(in) :: bias(ndat_a)
    real(kind = 8),intent(in) :: rmsd(ndat_a)    
    real(kind = 8),intent(in) :: sprd(ndat_a)
    real(kind = 8),intent(in) :: cor1(ndat_a),cor2(ndat_a)

    real(kind = 8),intent(in) :: abias_dif_low(ndat_a,ndat_a),abias_dif_ave(ndat_a,ndat_a),abias_dif_upp(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_dif_low(ndat_a,ndat_a),rmsd_dif_ave(ndat_a,ndat_a),rmsd_dif_upp(ndat_a,ndat_a)

    write(format1,'(a,I0,a)') "(",ndat_a*3,"f12.5)"
    write(format2,'(a,I0,a)') "(",ndat_a,"f12.5)"
    write(format3,'(a,I0,a)') "(",2*ndat_a,"f12.5)"

    open(1,file="dat/bias_ave_all.dat",status="replace")
    open(2,file="dat/rmsd_ave_all.dat",status="replace")
    open(3,file="dat/sprd_ave_all.dat",status="replace")
    open(4,file="dat/cor_ave_all.dat",status="replace")
    write(1,trim(format1)) bias(:),abias_dif_low(1,:),abias_dif_upp(1,:)
    write(2,trim(format1)) rmsd(:),rmsd_dif_low(1,:),rmsd_dif_upp(1,:)
    write(3,trim(format2)) sprd(:)
    write(4,trim(format3)) cor1(:),cor2(:)
    close(1)
    close(2)
    close(3)
    close(4)

  end subroutine write_stat_ave_all
  
end module mod_io_common
