module mod_read_jcope_fgo

  integer,parameter :: im=3602,jm=1502,km=44

  character(100),parameter :: jcope_fgo_dir="/lvs0/rccs-dart/ohishi/DATA/JCOPE-FGO"
  character(100),parameter :: grid_file="/lvs0/rccs-dart/ohishi/DATA/JCOPE-FGO/basic.dat"

contains

  !---------------------------------------------------------------------------
  ! Read JCOPE-FGO |
  !---------------------------------------------------------------------------

  subroutine get_jcope_fgo_info(var_in,varname)

    implicit none

    !---IN
    character(1),intent(in) :: var_in

    !---OUT
    character(2),intent(out) :: varname

    if(var_in == "h")then
       varname="EL"
    else if(var_in == "t")then
       varname="T"       
    else if(var_in == "s")then
       varname="S"
    else if(var_in == "u")then
       varname="U"
    else if(var_in == "v")then
       varname="V"
    else
       write(*,*) "***Error: Incorrect varname => "//trim(var_in)
       stop
    end if

  end subroutine get_jcope_fgo_info

  !----------------------------------------------------------------------------

  subroutine read_grid(lont,lonu,lonv,latt,latu,latv,maskt,masku,maskv,dep)

    implicit none

    !---Parameter
    real(kind = 4),parameter :: dmiss=1.e22
    real(kind = 8),parameter :: res=0.1d0    
    character(10),parameter :: yyyymmddhh="2003010112"
    
    !---Common
    integer i,j,k
    integer iyy,imm,idd,ihh
    integer ivar
    integer status,access
    
    real(kind = 4),allocatable :: z(:,:,:),zz(:,:,:),dz(:,:,:)
    real(kind = 4),allocatable :: tmp(:,:,:)

    real(kind = 8) mask(im,jm)
    
    character(200) filename
    character(4) cnull,param4

    !---OUT
    real(kind = 8),intent(out) :: lont(im),lonu(im),lonv(im)
    real(kind = 8),intent(out) :: latt(jm),latu(jm),latv(jm)
    real(kind = 8),intent(out) :: maskt(im,jm),masku(im,jm),maskv(im,jm)
    real(kind = 8),intent(out) :: dep(im,jm,km)

    !---Longitude & Latitude
    do i=1,im
       lont(i)=-0.1d0+dble(i-1)*res
       lonu(i)=-0.05d0+dble(i-1)*res
       lonv(i)=lont(i)
    end do

    do j=1,jm
       latt(j)=-75.05d0+dble(j-1)*res
       latu(j)=latt(j)
       latv(j)=-75.d0+dble(j-1)*res
    end do

    !---Check file
    status=access(trim(grid_file)," ")
    if(status == 0)then
       write(*,*) "Read :"//trim(grid_file)
    else
       write(*,*) "***Error: Not found "//trim(grid_file)
       stop
    end if

    !---Read Depth
    allocate(z(im,jm,km),zz(im,jm,km),dz(im,jm,km))

    open(1,file=trim(grid_file),status="old",access="stream",form="unformatted", &
         action="read",convert="big_endian")
    read(1,pos=5) cnull
    read(1) z
    read(1) zz
    read(1) dz
    close(1)

    !$omp parallel    
    !$omp do private(i,j,k) collapse(3)
    do k=1,km
       do j=1,jm
          do i=1,im
             if(abs(zz(i,j,km)) <= 1.e0)then
                dep(i,j,k)=0.d0
             else
                dep(i,j,k)=dble(abs(zz(i,j,k)))
             end if
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(z,zz,dz)

    !---Mask    
    allocate(tmp(im,jm,km))

    do ivar=1,3

       if(ivar == 1)then
          filename=trim(jcope_fgo_dir)//"/T_"//yyyymmddhh
       else if(ivar == 2)then
          filename=trim(jcope_fgo_dir)//"/U_"//yyyymmddhh
       else if(ivar == 3)then
          filename=trim(jcope_fgo_dir)//"/V_"//yyyymmddhh
       end if
          
       open(1,file=trim(filename),status="old",access="sequential",form="unformatted", &
            action="read",convert="big_endian")
       read(1) param4,iyy,imm,idd,ihh,tmp
       close(1)
    
       !$omp parallel
       !$omp do private(i,j) collapse(2)
       do j=1,jm
          do i=1,im
             if(tmp(i,j,1) == dmiss)then
                mask(i,j)=0.d0
             else
                mask(i,j)=1.d0
             end if
          end do
       end do
       !$omp end do
       !$omp end parallel

       if(ivar == 1)then
          maskt(:,:)=mask(:,:)
       else if(ivar == 2)then
          masku(:,:)=mask(:,:)
       else if(ivar == 3)then
          maskv(:,:)=mask(:,:)
       end if
       
    end do !ivar

  end subroutine read_grid

  !----------------------------------------------------------------------------

  subroutine read_jcope_fgo(var_in,iyr,imon,iday,km_in,mask,dat)

    !$use omp_lib    
    use mod_rmiss
    implicit none

    !---Parameter
    real(kind = 4),parameter :: dmiss=1.e22

    !---Common
    integer i,j,k
    integer iyy,imm,idd,ihh
    integer status,access

    real(kind = 4),allocatable :: tmp(:,:,:)

    character(200) filename
    character(4) param4
    character(4) yyyy
    character(2) mm,dd,hh
    character(2) varname    

    !---IN
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: km_in

    real(kind = 8),intent(in) :: mask(im,jm)

    character(1),intent(in) :: var_in

    !---OUT
    real(kind = 8),intent(out) :: dat(im,jm,km_in)

    !---Get information
    call get_jcope_fgo_info(var_in,varname)

    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') 12

    filename=trim(jcope_fgo_dir)//"/"//trim(varname)//"_"//yyyy//mm//dd//hh

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read :"//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if

    !---Read data
    allocate(tmp(im,jm,km_in))

    open(1,file=trim(filename),status="old",access="sequential",form="unformatted", &
         action="read",convert="big_endian")
    read(1) param4,iyy,imm,idd,ihh,tmp
    close(1)

    !---Post process
    !$omp parallel
    !$omp do private(i,j,k) collapse(3)
    do k=1,km_in
       do j=1,jm
          do i=1,im
             if(mask(i,j) == 0.d0 .or. tmp(i,j,k) == dmiss)then
                dat(i,j,k)=rmiss
             else
                dat(i,j,k)=dble(tmp(i,j,k))
             end if
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
    
    deallocate(tmp)

  end subroutine read_jcope_fgo

  !----------------------------------------------------------------------------

  subroutine extract_jcope_fgo(var_in,iyr,imon,iday,is,im_in,js,jm_in,ks,km_in,dat)

    !$use omp_lib    
    use mod_rmiss
    implicit none

    !---Parameter
    real(kind = 4),parameter :: dmiss=1.e22

    !---Common
    integer i,j,k
    integer iyy,imm,idd,ihh
    integer status,access

    real(kind = 4),allocatable :: tmp(:,:,:)

    character(200) filename
    character(4) param4
    character(4) yyyy
    character(2) mm,dd,hh
    character(2) varname    

    !---IN
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: is,im_in
    integer,intent(in) :: js,jm_in
    integer,intent(in) :: ks,km_in

    character(1),intent(in) :: var_in

    !---OUT
    real(kind = 8),intent(out) :: dat(im_in,jm_in,km_in)

    !---Get information
    call get_jcope_fgo_info(var_in,varname)

    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') 12

    filename=trim(jcope_fgo_dir)//"/"//trim(varname)//"_"//yyyy//mm//dd//hh

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read :"//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if

    !---Read data
    allocate(tmp(im,jm,km))

    open(1,file=trim(filename),status="old",access="sequential",form="unformatted", &
         action="read",convert="big_endian")
    read(1) param4,iyy,imm,idd,ihh,tmp
    close(1)
    
    !---Post process    
    !$omp parallel
    !$omp do private(i,j,k) collapse(3)
    do k=1,km_in
       do j=1,jm_in
          do i=1,im_in
             if(tmp(is+i-1,js+j-1,ks+k-1) == dmiss .or. 1.e10 <= abs(tmp(is+i-1,js+j-1,ks+k-1)))then
                dat(i,j,k)=rmiss
             else
                dat(i,j,k)=dble(tmp(is+i-1,js+j-1,ks+k-1))
             end if
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
    
    deallocate(tmp)

  end subroutine extract_jcope_fgo
  
end module mod_read_jcope_fgo
