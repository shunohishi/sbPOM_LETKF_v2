module mod_io

contains

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

  !---------------------------------------------------------------------------------

  subroutine get_datname(idat,datname)

    implicit none

    !---IN
    integer,intent(in) :: idat

    !---OUT
    character(10),intent(out) :: datname

    if(idat == 1) datname="lora"
    if(idat == 2) datname="glorys"
    if(idat == 3) datname="oras5"
    if(idat == 4) datname="cglors"
    
  end subroutine get_datname
  
  !---------------------------------------------------------------------------------

  subroutine read_grid(idat,im,jm,km,lon,lat,mask)

    use mod_read_lora, only: read_grid_lora => read_grid
    use mod_read_glorys025, only: read_glorys025
    implicit none

    !---Common
    integer i,j
    
    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) :: dir="QGLOBAL"

    !GLORYS
    character(10) datname
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: im,jm,km

    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm),mask(im,jm)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))
    allocate(tmp3d(im,jm,km))

    call get_datname(idat,datname)
    
    if(idat == 1)then
       dir="QGLOBAL"
       call read_grid_lora(dir,lon,tmp1dx,tmp1dx, &
            & lat,tmp1dy,tmp1dy, &
            & tmp3d,tmp3d,tmp3d, &
            & mask,tmp2d,tmp2d)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       varname="t"
       call read_glorys025(datname,varname,2003,1,1,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,tmp3d)
       lon(:)=tmp1dx(:)
       lat(:)=tmp1dy(:)
       mask(:,:)=tmp2d(:,:)
    end if

    !open(1,file="mask.dat",status="replace")
    !do j=1,jm
    !   do i=1,im             
    !      write(1,'(3f12.5)') lon(i),lat(j), mask(i,j)
    !   end do
    !end do
    !close(1)
    !stop

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)
    deallocate(tmp3d)

  end subroutine read_grid

  !-------------------------------------------------

  subroutine read_ssh_anal(idat,iyr,imon,iday,im,jm,km,mask,dat,sprd)

    use mod_read_lora, only: read_anal
    use mod_read_glorys025, only: read_glorys025
    use mod_rmiss
    implicit none

    !---Common
    integer k

    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)

    !LORA
    integer imem !Dummy
    character(10) :: dir="QGLOBAL"
    character(10) :: letkf="letkf"
    character(10) :: region="qglobal"
    character(10) :: ms

    !GLORYS
    character(10) datname
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: im,jm,km

    real(kind = 8),intent(in) :: mask(im,jm)

    !---OUT
    real(kind = 8),intent(out) :: dat(im,jm)
    real(kind = 8),intent(out) :: sprd(im,jm)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))

    call get_datname(idat,datname)    

    if(idat == 1)then
       imem=0
       k=1
       ms="mean"
       call read_anal(dir,letkf,region,ms,imem,"el",iyr,imon,iday,im,jm,k,mask,dat)
       ms="sprd"
       call read_anal(dir,letkf,region,ms,imem,"el",iyr,imon,iday,im,jm,k,mask,sprd)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       k=1
       varname="h"
       call read_glorys025(datname,varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,dat)
       sprd=rmiss
    end if    
    
    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine read_ssh_anal

  !----------------------------------------------------

  subroutine write_hdat(idat_a,nst,ijul,lon_a,lat_a,hdat_a,hsprd_a, &
       & lon_o,lat_o,dat_o,dist)

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

    status=nf90_inq_varid(ncid,"lon_a",varid)
    status=nf90_put_var(ncid,varid,lon_a(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"lat_a",varid)
    status=nf90_put_var(ncid,varid,lat_a(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"hmean_a",varid)
    status=nf90_put_var(ncid,varid,hdat_a(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"hsprd_a",varid)
    status=nf90_put_var(ncid,varid,hsprd_a(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"lon_o",varid)
    status=nf90_put_var(ncid,varid,lon_o(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"lat_o",varid)
    status=nf90_put_var(ncid,varid,lat_o(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"h_o",varid)
    status=nf90_put_var(ncid,varid,dat_o(1:nst),(/1,iday/),(/nst,1/))
    
    status=nf90_inq_varid(ncid,"dist",varid)
    status=nf90_put_var(ncid,varid,dist(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_close(ncid)    
    
  end subroutine write_hdat

  !----------------------------------------------------

  subroutine read_hdat(idat_a,nst,ijul,lon_a,lat_a,hdat_a,hsprd_a, &
       & lon_o,lat_o,dat_o,dist)

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
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"lon_a",varid)
    status=nf90_get_var(ncid,varid,lon_a(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"lat_a",varid)
    status=nf90_get_var(ncid,varid,lat_a(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"hmean_a",varid)
    status=nf90_get_var(ncid,varid,hdat_a(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"hsprd_a",varid)
    status=nf90_get_var(ncid,varid,hsprd_a(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"lon_o",varid)
    status=nf90_get_var(ncid,varid,lon_o(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"lat_o",varid)
    status=nf90_get_var(ncid,varid,lat_o(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_inq_varid(ncid,"h_o",varid)
    status=nf90_get_var(ncid,varid,dat_o(1:nst),(/1,iday/),(/nst,1/))
    
    status=nf90_inq_varid(ncid,"dist",varid)
    status=nf90_get_var(ncid,varid,dist(1:nst),(/1,iday/),(/nst,1/))

    status=nf90_close(ncid)    
    
  end subroutine read_hdat

  !-----------------------------------------------

  subroutine write_ts(ndat_a,nst,sjul,ejul,dat_a,dat_o)

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

    real(kind = 8),intent(in) :: dat_a(nst,sjul:ejul,ndat_a)
    real(kind = 8),intent(in) :: dat_o(nst,sjul:ejul)

    !--Format
    write(format,'(a,I0,a)') "(a,",ndat_a+1,"f12.5)"

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
          write(1,trim(format)) yyyymmdd,dat_a(ist,ijul,1:ndat_a),dat_o(ist,ijul)
       end do
       close(1)

    end do !ist
    
  end subroutine write_ts

  !----------------------------------------------

  subroutine write_stat_mave(nst,ndat_a,syr,eyr,lon,lat,num,bias,rmsd,sprd)

    implicit none

    !---Common
    integer ist
    integer iyr,imon
    
    character(100) format
    
    !---IN
    integer,intent(in) :: nst
    integer,intent(in) :: ndat_a
    integer,intent(in) :: syr,eyr

    integer,intent(in) :: num(nst,ndat_a,12,syr:eyr)

    real(kind = 8),intent(in) :: lon(nst),lat(nst)
    real(kind = 8),intent(in) :: bias(nst,ndat_a,12,syr:eyr)
    real(kind = 8),intent(in) :: rmsd(nst,ndat_a,12,syr:eyr)
    real(kind = 8),intent(in) :: sprd(nst,ndat_a,12,syr:eyr)

    write(format,'(a,I0,a,I0,a)') "(i6,2f12.5,2i6,",ndat_a,"i10,",ndat_a,"f12.5)"
    
    open(1,file="dat/bias_mave.dat",status="replace")
    open(2,file="dat/rmsd_mave.dat",status="replace")
    open(3,file="dat/sprd_mave.dat",status="replace")
    do ist=1,nst
       do iyr=syr,eyr
          do imon=1,12
             write(1,trim(format)) ist,lon(ist),lat(ist),iyr,imon,num(ist,:,imon,iyr),bias(ist,:,imon,iyr)
             write(2,trim(format)) ist,lon(ist),lat(ist),iyr,imon,num(ist,:,imon,iyr),rmsd(ist,:,imon,iyr)
             write(3,trim(format)) ist,lon(ist),lat(ist),iyr,imon,num(ist,:,imon,iyr),sprd(ist,:,imon,iyr)
          end do
       end do
    end do
    close(1)
    close(2)
    close(3)

  end subroutine write_stat_mave

  !-----------------------------------------------
  
  subroutine write_stat_ave(nst,ndat_a,lon,lat,num,bias,rmsd,sprd, &
       & bias_dof,bias_tcrit,bias_tval, &
       & rmsd_dof,rmsd_tcrit,rmsd_tval)

    implicit none

    !---Common
    integer ist
    integer idat_a

    character(100) format1,format2

    !---IN
    integer,intent(in) :: nst
    integer,intent(in) :: ndat_a

    integer,intent(in) :: num(nst,ndat_a)
    integer,intent(in) :: bias_dof(nst,ndat_a,ndat_a)
    integer,intent(in) :: rmsd_dof(nst,ndat_a,ndat_a)

    
    real(kind = 8),intent(in) :: lon(nst),lat(nst)
    real(kind = 8),intent(in) :: bias(nst,ndat_a)
    real(kind = 8),intent(in) :: rmsd(nst,ndat_a)    
    real(kind = 8),intent(in) :: sprd(nst,ndat_a)

    real(kind = 8),intent(in) :: bias_tcrit(nst,ndat_a,ndat_a),bias_tval(nst,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_tcrit(nst,ndat_a,ndat_a),rmsd_tval(nst,ndat_a,ndat_a)

    write(format1,'(a,I0,a,I0,a)') "(i6,2f12.5,",ndat_a,"i10,",ndat_a*3,"f12.5)"
    write(format2,'(a,I0,a,I0,a)') "(i6,2f12.5,",ndat_a,"i10,",ndat_a,"f12.5)"
    
    open(1,file="dat/bias_ave.dat",status="replace")
    open(2,file="dat/rmsd_ave.dat",status="replace")
    open(3,file="dat/sprd_ave.dat",status="replace")
    do ist=1,nst
       write(1,trim(format1)) ist,lon(ist),lat(ist),num(ist,:),bias(ist,:),bias_tcrit(ist,1,:),bias_tval(ist,1,:)
       write(2,trim(format1)) ist,lon(ist),lat(ist),num(ist,:),rmsd(ist,:),rmsd_tcrit(ist,1,:),rmsd_tval(ist,1,:)
       write(3,trim(format2)) ist,lon(ist),lat(ist),num(ist,:),sprd(ist,:)
    end do
    close(1)
    close(2)
    close(3)

  end subroutine write_stat_ave

  !----------------------------------------------

  subroutine write_stat_ave_all(ndat_a,bias,rmsd,sprd, &
       & bias_dof,bias_tcrit,bias_tval, &
       & rmsd_dof,rmsd_tcrit,rmsd_tval)

    implicit none

    !---Common
    character(100) format1,format2

    !---IN
    integer,intent(in) :: ndat_a

    integer,intent(in) :: bias_dof(ndat_a,ndat_a)
    integer,intent(in) :: rmsd_dof(ndat_a,ndat_a)
    
    real(kind = 8),intent(in) :: bias(ndat_a)
    real(kind = 8),intent(in) :: rmsd(ndat_a)    
    real(kind = 8),intent(in) :: sprd(ndat_a)

    real(kind = 8),intent(in) :: bias_tcrit(ndat_a,ndat_a),bias_tval(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_tcrit(ndat_a,ndat_a),rmsd_tval(ndat_a,ndat_a)

    write(format1,'(a,I0,a)') "(",ndat_a*3,"f12.5)"
    write(format2,'(a,I0,a)') "(",ndat_a,"f12.5)"
    
    open(1,file="dat/bias_ave_all.dat",status="replace")
    open(2,file="dat/rmsd_ave_all.dat",status="replace")
    open(3,file="dat/sprd_ave_all.dat",status="replace")
    write(1,trim(format1)) bias(:),bias_tcrit(1,:),bias_tval(1,:)
    write(2,trim(format1)) rmsd(:),rmsd_tcrit(1,:),rmsd_tval(1,:)
    write(3,trim(format2)) sprd(:)
    close(1)
    close(2)
    close(3)

  end subroutine write_stat_ave_all
    
end module mod_io
