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

  subroutine read_grid(idat,im,jm,km,lont,lonu,lonv,latt,latu,latv,dept,depu,depv,maskt,masku,maskv)

    use setting, only: datname
    use mod_read_lora, only: read_grid_lora => read_grid
    use mod_read_glorys025, only: extract_glorys025
    implicit none

    !---Common
    integer i,j
    
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) :: dir="QGLOBAL"

    !GLORYS
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: im,jm,km

    !---OUT
    real(kind = 8),intent(out) :: lont(im),lonu(im),lonv(im)
    real(kind = 8),intent(out) :: latt(jm),latu(jm),latv(jm)
    real(kind = 8),intent(out) :: dept(im,jm,km),depu(im,jm,km),depv(im,jm,km)
    real(kind = 8),intent(out) :: maskt(im,jm),masku(im,jm),maskv(im,jm)

    allocate(tmp3d(im,jm,km))

    if(idat == 1)then
       dir="QGLOBAL"
       call read_grid_lora(dir,lont,lonu,lonv, &
            & latt,latu,latv, &
            & dept,depu,depv, &
            & maskt,masku,maskv)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       varname="t"
       call extract_glorys025(datname(idat),varname,2003,1,1,1,im,1,jm,1,km,lont,latt,dept(1,1,:),maskt,tmp3d)
       lonu(:)=lont(:)
       lonv(:)=lont(:)
       latu(:)=latt(:)
       latv(:)=latt(:)
       do j=1,jm
          do i=1,im
             dept(i,j,:)=dept(1,1,:)
             depu(i,j,:)=dept(1,1,:)
             depv(i,j,:)=dept(1,1,:)
          end do
       end do
       masku(:,:)=maskt(:,:)
       maskv(:,:)=maskt(:,:)
    end if

    deallocate(tmp3d)

  end subroutine read_grid

  !-------------------------------------------------

  subroutine extract_data(varname,idat,iyr,imon,iday,is,im,js,jm,ks,km,mean,sprd)

    use setting, only: datname
    use mod_read_lora, only: extract_anal
    use mod_read_glorys025, only: extract_glorys025
    use mod_rmiss
    implicit none

    !---Common
    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)

    !LORA
    integer imem !Dummy
    character(10) :: dir="QGLOBAL"
    character(10) :: letkf="letkf"
    character(10) :: region="qglobal"
    character(10) :: ms

    !GLORYS
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: is,im,js,jm,ks,km

    !---OUT
    real(kind = 8),intent(out) :: mean(im,jm,km),sprd(im,jm,km)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))

    if(idat == 1)then
       imem=0
       ms="mean"
       call extract_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,is,im,js,jm,ks,km,mean)
       ms="sprd"
       call extract_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,is,im,js,jm,ks,km,sprd)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       call extract_glorys025(datname(idat),varname, &
            & iyr,imon,iday,is,im,js,jm,ks,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)
       sprd=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine extract_data

  !------------------------------------

  subroutine read_hdata(buoyname,varname,idat_a,iyr,imon,iday, &
       & km_o,lon_o,lat_o,dep_o,dat_o,hmean_a,hsprd_a)

    use setting, only: datname
    use mod_rmiss
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

    !---Get km_o
    status=nf90_inq_dimid(ncid,"z",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = km_o)
    
    !---Allocate
    allocate(dep_o(km_o),dat_o(km_o))
    allocate(hmean_a(km_o),hsprd_a(km_o))
    
    !---Read data
    status=nf90_inq_varid(ncid,"lon_o",varid)
    status=nf90_get_var(ncid,varid,lon_o)

    status=nf90_inq_varid(ncid,"lat_o",varid)
    status=nf90_get_var(ncid,varid,lat_o)

    status=nf90_inq_varid(ncid,"dep_o",varid)
    status=nf90_get_var(ncid,varid,dep_o)

    status=nf90_inq_varid(ncid,"h"//trim(varname)//"mean_a",varid)
    status=nf90_get_var(ncid,varid,hmean_a,(/1,iday/),(/km_o,1/))

    status=nf90_inq_varid(ncid,"h"//trim(varname)//"sprd_a",varid)
    status=nf90_get_var(ncid,varid,hsprd_a,(/1,iday/),(/km_o,1/))

    status=nf90_inq_varid(ncid,trim(varname)//"_o",varid)
    status=nf90_get_var(ncid,varid,dat_o,(/1,iday/),(/km_o,1/))
        
    status=nf90_close(ncid)
        
  end subroutine read_hdata

  !------------------------------------

  subroutine deallcate_hdata(dep_o,dat_o,hmean_a,hsprd_a)

    implicit none

    real(kind = 8),allocatable,intent(inout) :: dep_o(:),dat_o(:)
    real(kind = 8),allocatable,intent(inout) :: hmean_a(:),hsprd_a(:)

    if(allocated(dep_o)) deallocate(dep_o)
    if(allocated(dat_o)) deallocate(dat_o)
    if(allocated(hmean_a)) deallocate(hmean_a)
    if(allocated(hsprd_a)) deallocate(hsprd_a)
    
  end subroutine deallcate_hdata    
  
  !------------------------------------
  
  subroutine write_hdata(buoyname,varname,idat_a,iyr,imon,iday, &
       & km_o,lon_o,lat_o,dep_o,dat_o,hmean_a,hsprd_a)

    use setting, only: datname
    use mod_make_ncfile
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

    real(kind = 8),intent(in) :: lon_o,lat_o,dep_o(km_o)
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

    status=nf90_inq_varid(ncid,"lon_o",varid)
    status=nf90_put_var(ncid,varid,lon_o)

    status=nf90_inq_varid(ncid,"lat_o",varid)
    status=nf90_put_var(ncid,varid,lat_o)

    status=nf90_inq_varid(ncid,"dep_o",varid)
    status=nf90_put_var(ncid,varid,dep_o)

    status=nf90_inq_varid(ncid,"h"//trim(varname)//"mean_a",varid)
    status=nf90_put_var(ncid,varid,hmean_a,(/1,iday/),(/km_o,1/))

    status=nf90_inq_varid(ncid,"h"//trim(varname)//"sprd_a",varid)
    status=nf90_put_var(ncid,varid,hsprd_a,(/1,iday/),(/km_o,1/))

    status=nf90_inq_varid(ncid,trim(varname)//"_o",varid)
    status=nf90_put_var(ncid,varid,dat_o,(/1,iday/),(/km_o,1/))
        
    status=nf90_close(ncid)
        
  end subroutine write_hdata

  !------------------------------------------------------------------

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
  
  !------------------------------------------------------------------

  subroutine write_mave(buoyname,varname,syr,eyr,ndat_a,km_o,dep_o,num_mave,bias_mave,rmsd_mave,sprd_mave)

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
    integer,intent(in) :: num_mave(km_o,ndat_a,12,syr:eyr)
    
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
             
             write(11,trim(format)) yyyymmdd,dep_o(k),num_mave(k,:,imon,iyr),bias_mave(k,:,imon,iyr)
             write(12,trim(format)) yyyymmdd,dep_o(k),num_mave(k,:,imon,iyr),rmsd_mave(k,:,imon,iyr)
             write(13,trim(format)) yyyymmdd,dep_o(k),num_mave(k,:,imon,iyr),sprd_mave(k,:,imon,iyr)

          end do
       end do
    end do
    close(11)
    close(12)
    close(13)             

  end subroutine write_mave

  !---------------------------------------------------------

  subroutine write_ave(buoyname,varname,sjul,ejul,ndat_a,km_o,dep_o,num_ave, &
       & bias_ave,bias_dof_ave,bias_tcrit_ave,bias_tval_ave, &
       & rmsd_ave,rmsd_dof_ave,rmsd_tcrit_ave,rmsd_tval_ave, &
       & sprd_ave)

    implicit none

    !---Common
    integer k
    integer idat_a
    
    character(100) format

    !---IN
    integer,intent(in) :: sjul,ejul
    integer,intent(in) :: ndat_a
    integer,intent(in) :: km_o
    integer,intent(in) :: num_ave(km_o,ndat_a)
    integer,intent(in) :: bias_dof_ave(km_o,ndat_a,ndat_a)
    integer,intent(in) :: rmsd_dof_ave(km_o,ndat_a,ndat_a)
    
    real(kind = 8),intent(in) :: dep_o(km_o)
    real(kind = 8),intent(in) :: bias_ave(km_o,ndat_a)
    real(kind = 8),intent(in) :: bias_tcrit_ave(km_o,ndat_a,ndat_a),bias_tval_ave(km_o,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: rmsd_ave(km_o,ndat_a)
    real(kind = 8),intent(in) :: rmsd_tcrit_ave(km_o,ndat_a,ndat_a),rmsd_tval_ave(km_o,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: sprd_ave(km_o,ndat_a)

    character(10),intent(in) :: buoyname
    character(1),intent(in) :: varname
    
    write(format,'(a,I0,a)') "(f12.5,",ndat_a*4,"f12.5)"
    
    open(11,file="dat/"//trim(buoyname)//"/"//trim(varname)//"bias_ave.dat",status="replace")
    open(12,file="dat/"//trim(buoyname)//"/"//trim(varname)//"rmsd_ave.dat",status="replace")
    do k=1,km_o

       write(11,trim(format)) dep_o(k),num_ave(k,:)*100.d0/dble(ejul-sjul+1),bias_ave(k,:),bias_tcrit_ave(k,1,:),bias_tval_ave(k,1,:)
       write(12,trim(format)) dep_o(k),num_ave(k,:)*100.d0/dble(ejul-sjul+1),rmsd_ave(k,:),rmsd_tcrit_ave(k,1,:),rmsd_tval_ave(k,1,:)

    end do
    close(11)
    close(12)

    write(format,'(a,I0,a)') "(f12.5,",2*ndat_a,"f12.5)"    

    open(13,file="dat/"//trim(buoyname)//"/"//trim(varname)//"sprd_ave.dat",status="replace")              
    do k=1,km_o
       write(13,trim(format)) dep_o(k),num_ave(k,:)*100.d0/dble(ejul-sjul+1),sprd_ave(k,:)
    end do
    close(13)             
    
    !write(format,'(a,I0,a,I0,a)') "(f12.5,i6,",ndat_a,"i10,",2*ndat_a,"f12.5)"
    
    !open(11,file="dat/"//trim(buoyname)//"/"//trim(varname)//"bias_tval_ave.dat",status="replace")
    !open(12,file="dat/"//trim(buoyname)//"/"//trim(varname)//"rmsd_tval_ave.dat",status="replace")
    !do k=1,km_o
    !   do idat_a=1,ndat_a

    !      write(11,trim(format)) dep_o(k),idat_a, &
    !           & bias_dof_ave(k,idat_a,:),bias_tcrit_ave(k,idat_a,:),bias_tval_ave(k,idat_a,:)
    !      write(12,trim(format)) dep_o(k),idat_a, &
    !           & rmsd_dof_ave(k,idat_a,:),rmsd_tcrit_ave(k,idat_a,:),rmsd_tval_ave(k,idat_a,:)
          
    !   end do
    !end do
    !close(11)
    !close(12)
    
  end subroutine write_ave
  
end module mod_io
