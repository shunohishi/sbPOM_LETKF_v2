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

  subroutine read_grid(idat,im,jm,km,lont,lonu,lonv,latt,latu,latv,maskt,masku,maskv)

    use mod_read_lora, only: read_grid_lora => read_grid
    use mod_read_glorys025, only: read_glorys025
    implicit none

    !---Common
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
    real(kind = 8),intent(out) :: lont(im),lonu(im),lonv(im)
    real(kind = 8),intent(out) :: latt(jm),latu(jm),latv(jm)
    real(kind = 8),intent(out) :: maskt(im,jm),masku(im,jm),maskv(im,jm)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))
    allocate(tmp3d(im,jm,km))

    if(idat == 1)then
       dir="QGLOBAL"
       call read_grid_lora(dir,lont,lonu,lonv, &
            & latt,latu,latv, &
            & tmp3d,tmp3d,tmp3d, &
            & maskt,masku,maskv)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       if(idat == 2) datname="glorys"
       if(idat == 3) datname="oras5"
       if(idat == 4) datname="cglors"
       varname="t"
       call read_glorys025(datname,varname,2003,1,1,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,tmp3d)
       lont(:)=tmp1dx(:)
       lonu(:)=tmp1dx(:)
       lonv(:)=tmp1dx(:)
       latt(:)=tmp1dy(:)
       latu(:)=tmp1dy(:)
       latv(:)=tmp1dy(:)
       maskt(:,:)=tmp2d(:,:)
       masku(:,:)=tmp2d(:,:)
       maskv(:,:)=tmp2d(:,:)
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)
    deallocate(tmp3d)

  end subroutine read_grid

  !-------------------------------------------------

  subroutine read_surface_data(idat,iyr,imon,iday,im,jm,km,maskt,masku,maskv,t,u,v,tsprd,usprd,vsprd)

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

    real(kind = 8),intent(in) :: maskt(im,jm),masku(im,jm),maskv(im,jm)

    !---OUT
    real(kind = 8),intent(out) :: t(im,jm),u(im,jm),v(im,jm)
    real(kind = 8),intent(out) :: tsprd(im,jm),usprd(im,jm),vsprd(im,jm)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))

    if(idat == 1)then
       imem=0
       k=1
       ms="mean"
       call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im,jm,k,maskt,t)
       call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im,jm,k,masku,u)
       call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im,jm,k,maskv,v)
       ms="sprd"
       call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im,jm,k,maskt,tsprd)
       call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im,jm,k,masku,usprd)
       call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im,jm,k,maskv,vsprd)     
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       if(idat == 2) datname="glorys"
       if(idat == 3) datname="oras5"
       if(idat == 4) datname="cglors"
       k=1
       varname="t"
       call read_glorys025(datname,varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,t)
       varname="u"
       call read_glorys025(datname,varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,u)
       varname="v"
       call read_glorys025(datname,varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine read_surface_data

  !-------------------------------------------------

  subroutine read_obs(idat_a,iyr,imon,iday,nobs,lon_o,lat_o, &
       & ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o,u_o,v_o)

    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,dimid,varid

    character(100) filename
    character(10) datname
    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd

    !---IN
    integer,intent(in) :: idat_a
    integer,intent(in) :: iyr,imon,iday

    !---OUT
    integer,intent(out) :: nobs

    real(kind = 8),intent(out),allocatable :: lon_o(:),lat_o(:)
    real(kind = 8),intent(out),allocatable :: ht_a(:),hu_a(:),hv_a(:)
    real(kind = 8),intent(out),allocatable :: htsprd_a(:),husprd_a(:),hvsprd_a(:)
    real(kind = 8),intent(out),allocatable :: t_o(:),u_o(:),v_o(:)

    write(*,*) "Start"

    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd

    if(idat_a == 1)then
       datname="lora"
    else if(idat_a == 2)then
       datname="glorys"
    else if(idat_a == 3)then
       datname="oras5"
    else if(idat_a == 4)then
       datname="cglors"
    else
       write(*,*) "***Error: Incorrect datname"
       stop
    end if

    filename="dat/"//yyyy//mm//"/"//trim(datname)//"."//yyyymmdd//".nc"    

    write(*,*) trim(filename)

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read: "//trim(filename)
    else
       write(*,*) "***Error: Not found"//trim(filename)
       stop
    end if

    !---Read
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !nobs
    status=nf90_inq_dimid(ncid,"nobs",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = nobs)

    if(nobs == 0)then
       status=nf90_close(ncid)
       return
    end if

    !Allocate
    allocate(lon_o(nobs),lat_o(nobs))
    allocate(ht_a(nobs),hu_a(nobs),hv_a(nobs))
    allocate(htsprd_a(nobs),husprd_a(nobs),hvsprd_a(nobs))
    allocate(t_o(nobs),u_o(nobs),v_o(nobs))

    status=nf90_inq_varid(ncid,"lon_o",varid)
    status=nf90_get_var(ncid,varid,lon_o)

    status=nf90_inq_varid(ncid,"lat_o",varid)
    status=nf90_get_var(ncid,varid,lat_o)

    status=nf90_inq_varid(ncid,"ht_a",varid)
    status=nf90_get_var(ncid,varid,ht_a)

    status=nf90_inq_varid(ncid,"hu_a",varid)
    status=nf90_get_var(ncid,varid,hu_a)

    status=nf90_inq_varid(ncid,"hv_a",varid)
    status=nf90_get_var(ncid,varid,hv_a)

    status=nf90_inq_varid(ncid,"htsprd_a",varid)
    status=nf90_get_var(ncid,varid,htsprd_a)

    status=nf90_inq_varid(ncid,"husprd_a",varid)
    status=nf90_get_var(ncid,varid,husprd_a)

    status=nf90_inq_varid(ncid,"hvsprd_a",varid)
    status=nf90_get_var(ncid,varid,hvsprd_a)

    status=nf90_inq_varid(ncid,"t_o",varid)
    status=nf90_get_var(ncid,varid,t_o)

    status=nf90_inq_varid(ncid,"u_o",varid)
    status=nf90_get_var(ncid,varid,u_o)

    status=nf90_inq_varid(ncid,"v_o",varid)
    status=nf90_get_var(ncid,varid,v_o)

    status=nf90_close(ncid)

  end subroutine read_obs

  !-------------------------------------------------

  subroutine deallocate_obs(lon_o,lat_o, &
       & ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o,u_o,v_o)

    implicit none

    real(kind = 8),intent(inout),allocatable :: lon_o(:),lat_o(:)
    real(kind = 8),intent(inout),allocatable :: ht_a(:),hu_a(:),hv_a(:)
    real(kind = 8),intent(inout),allocatable :: htsprd_a(:),husprd_a(:),hvsprd_a(:)
    real(kind = 8),intent(inout),allocatable :: t_o(:),u_o(:),v_o(:)

    if(allocated(lon_o)) deallocate(lon_o)
    if(allocated(lat_o)) deallocate(lat_o)
    if(allocated(ht_a)) deallocate(ht_a)
    if(allocated(hu_a)) deallocate(hu_a)
    if(allocated(hv_a)) deallocate(hv_a)
    if(allocated(htsprd_a)) deallocate(htsprd_a)
    if(allocated(husprd_a)) deallocate(husprd_a)
    if(allocated(hvsprd_a)) deallocate(hvsprd_a)
    if(allocated(t_o)) deallocate(t_o)
    if(allocated(u_o)) deallocate(u_o)
    if(allocated(v_o)) deallocate(v_o)

  end subroutine deallocate_obs

  !-------------------------------------------------

  subroutine write_obs(idat_a,ijul,nobs,ijul_o,lon_o,lat_o, &
       & ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o,u_o,v_o, &
       & ncid,inum)

    use netcdf
    use mod_julian
    use mod_rmiss
    use mod_make_ncfile
    implicit none

    !---Common
    integer status,access,system
    integer varid
    integer iyr,imon,iday
    integer n

    real(kind = 8),allocatable :: lon_o_tmp(:),lat_o_tmp(:)
    real(kind = 8),allocatable :: ht_a_tmp(:),hu_a_tmp(:),hv_a_tmp(:)
    real(kind = 8),allocatable :: htsprd_a_tmp(:),husprd_a_tmp(:),hvsprd_a_tmp(:)
    real(kind = 8),allocatable :: t_o_tmp(:),u_o_tmp(:),v_o_tmp(:)
    
    character(100) filename
    character(10) datname
    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd

    logical, allocatable :: mask(:)
    
    !---IN
    integer,intent(in) :: idat_a
    integer,intent(in) :: ijul
    integer,intent(in) :: nobs
    integer,intent(in) :: ijul_o(nobs)

    real(kind = 8),intent(in) :: lon_o(nobs),lat_o(nobs)
    real(kind = 8),intent(in) :: ht_a(nobs),hu_a(nobs),hv_a(nobs)
    real(kind = 8),intent(in) :: htsprd_a(nobs),husprd_a(nobs),hvsprd_a(nobs)
    real(kind = 8),intent(in) :: t_o(nobs),u_o(nobs),v_o(nobs)

    !---INOUT
    integer,intent(inout) :: ncid
    integer,intent(inout) :: inum

    !---Filename
    call julian_ymd(ijul,iyr,imon,iday)
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd

    if(idat_a == 1)then
       datname="lora"
    else if(idat_a == 2)then
       datname="glorys"
    else if(idat_a == 3)then
       datname="oras5"
    else if(idat_a == 4)then
       datname="cglors"
    else
       write(*,*) "***Error: Incorrect datname"
       stop
    end if

    filename="dat/"//yyyy//mm//"/"//trim(datname)//"."//yyyymmdd//".nc"

    !---Make obs. file
    if(inum == 0)then

       status=access(trim(filename)," ")
       if(status == 0)then
          write(*,*) "Output to existing file:"//trim(filename)
          !status=system("rm -f "//trim(filename))
       else
          call make_obsfile(filename)
       end if
          
    end if

    !---Extract
    allocate(mask(nobs))
    mask = (ijul == ijul_o) .and. .not.(t_o == rmiss .and. u_o == rmiss .and. v_o == rmiss)

    lon_o_tmp    = pack(lon_o,    mask)
    lat_o_tmp    = pack(lat_o,    mask)
    ht_a_tmp     = pack(ht_a,     mask)
    hu_a_tmp     = pack(hu_a,     mask)
    hv_a_tmp     = pack(hv_a,     mask)
    htsprd_a_tmp = pack(htsprd_a, mask)
    husprd_a_tmp = pack(husprd_a, mask)
    hvsprd_a_tmp = pack(hvsprd_a, mask)
    t_o_tmp      = pack(t_o,      mask)
    u_o_tmp      = pack(u_o,      mask)
    v_o_tmp      = pack(v_o,      mask)

    n=size(lon_o_tmp)

    !---Write    
    if(n /= 0)then
       
       if(ncid == 0)then
          status=nf90_open(trim(filename),nf90_write,ncid)
       end if
          
       status=nf90_inq_varid(ncid,"lon_o",varid)
       status=nf90_put_var(ncid,varid,lon_o_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"lat_o",varid)
       status=nf90_put_var(ncid,varid,lat_o_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"ht_a",varid)
       status=nf90_put_var(ncid,varid,ht_a_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"hu_a",varid)
       status=nf90_put_var(ncid,varid,hu_a_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"hv_a",varid)
       status=nf90_put_var(ncid,varid,hv_a_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"htsprd_a",varid)
       status=nf90_put_var(ncid,varid,htsprd_a_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"husprd_a",varid)
       status=nf90_put_var(ncid,varid,husprd_a_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"hvsprd_a",varid)
       status=nf90_put_var(ncid,varid,hvsprd_a_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"t_o",varid)
       status=nf90_put_var(ncid,varid,t_o_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"u_o",varid)
       status=nf90_put_var(ncid,varid,u_o_tmp,(/inum+1/),(/n/))

       status=nf90_inq_varid(ncid,"v_o",varid)
       status=nf90_put_var(ncid,varid,v_o_tmp,(/inum+1/),(/n/))
       
    end if

    !---Count up
    inum=inum+n
    
    !---Deallocate
    if(allocated(mask))      deallocate(mask)
    if(allocated(lon_o_tmp)) deallocate(lon_o_tmp)
    if(allocated(lat_o_tmp)) deallocate(lat_o_tmp)
    if(allocated(ht_a_tmp))  deallocate(ht_a_tmp)
    if(allocated(hu_a_tmp))  deallocate(hu_a_tmp)
    if(allocated(hv_a_tmp))  deallocate(hv_a_tmp)
    if(allocated(htsprd_a_tmp))  deallocate(htsprd_a_tmp)
    if(allocated(husprd_a_tmp))  deallocate(husprd_a_tmp)
    if(allocated(hvsprd_a_tmp))  deallocate(hvsprd_a_tmp)
    if(allocated(t_o_tmp))  deallocate(t_o_tmp)
    if(allocated(u_o_tmp))  deallocate(u_o_tmp)
    if(allocated(v_o_tmp))  deallocate(v_o_tmp)
    
  end subroutine write_obs

  !-------------------------------------------------

  subroutine write_bin(im_bin,jm_bin,ndat_a,dx_bin,dy_bin,lon_bin,lat_bin, &
       & unum_bin,ubias_bin,urmsd_bin,usprd_bin,       &
       & ubias_dof_bin,ubias_tcrit_bin,ubias_tval_bin, &
       & urmsd_dof_bin,urmsd_tcrit_bin,urmsd_tval_bin, &
       & vnum_bin,vbias_bin,vrmsd_bin,vsprd_bin,       &
       & vbias_dof_bin,vbias_tcrit_bin,vbias_tval_bin, &
       & vrmsd_dof_bin,vrmsd_tcrit_bin,vrmsd_tval_bin, &
       & tnum_bin,tbias_bin,trmsd_bin,tsprd_bin,       &
       & tbias_dof_bin,tbias_tcrit_bin,tbias_tval_bin, &
       & trmsd_dof_bin,trmsd_tcrit_bin,trmsd_tval_bin)
    
    implicit none

    !---Common
    integer i_bin,j_bin
    integer idat_a

    character(100) format

    !---IN
    integer,intent(in) :: im_bin,jm_bin,ndat_a
    real(kind = 8),intent(in) :: dx_bin,dy_bin

    integer,intent(in) :: unum_bin(im_bin,jm_bin,ndat_a)
    integer,intent(in) :: ubias_dof_bin(im_bin,jm_bin,ndat_a,ndat_a),urmsd_dof_bin(im_bin,jm_bin,ndat_a,ndat_a)
    integer,intent(in) :: vnum_bin(im_bin,jm_bin,ndat_a)
    integer,intent(in) :: vbias_dof_bin(im_bin,jm_bin,ndat_a,ndat_a),vrmsd_dof_bin(im_bin,jm_bin,ndat_a,ndat_a)
    integer,intent(in) :: tnum_bin(im_bin,jm_bin,ndat_a)
    integer,intent(in) :: tbias_dof_bin(im_bin,jm_bin,ndat_a,ndat_a),trmsd_dof_bin(im_bin,jm_bin,ndat_a,ndat_a)

    real(kind = 8),intent(in) :: lon_bin(im_bin),lat_bin(jm_bin)    
    real(kind = 8),intent(in) :: ubias_bin(im_bin,jm_bin,ndat_a),urmsd_bin(im_bin,jm_bin,ndat_a),usprd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: ubias_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),ubias_tval_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: urmsd_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),urmsd_tval_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: vbias_bin(im_bin,jm_bin,ndat_a),vrmsd_bin(im_bin,jm_bin,ndat_a),vsprd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: vbias_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),vbias_tval_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: vrmsd_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),vrmsd_tval_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: tbias_bin(im_bin,jm_bin,ndat_a),trmsd_bin(im_bin,jm_bin,ndat_a),tsprd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: tbias_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),tbias_tval_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: trmsd_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),trmsd_tval_bin(im_bin,jm_bin,ndat_a,ndat_a)
    
    write(format,'(a,I0,a,I0,a)') '(2f12.5,',ndat_a,'i10,',ndat_a,'f12.5)'

    open(11,file="dat/ubias_bin.dat",status="replace")
    open(12,file="dat/vbias_bin.dat",status="replace")
    open(13,file="dat/tbias_bin.dat",status="replace")
    open(14,file="dat/urmsd_bin.dat",status="replace")
    open(15,file="dat/vrmsd_bin.dat",status="replace")
    open(16,file="dat/trmsd_bin.dat",status="replace")
    open(17,file="dat/usprd_bin.dat",status="replace")
    open(18,file="dat/vsprd_bin.dat",status="replace")
    open(19,file="dat/tsprd_bin.dat",status="replace")
    do j_bin=1,jm_bin-1
       do i_bin=1,im_bin-1
          
          write(11,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & unum_bin(i_bin,j_bin,:),ubias_bin(i_bin,j_bin,:)
          write(12,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & vnum_bin(i_bin,j_bin,:),vbias_bin(i_bin,j_bin,:)
          write(13,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & tnum_bin(i_bin,j_bin,:),tbias_bin(i_bin,j_bin,:)

          write(14,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & unum_bin(i_bin,j_bin,:),urmsd_bin(i_bin,j_bin,:)
          write(15,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & vnum_bin(i_bin,j_bin,:),vrmsd_bin(i_bin,j_bin,:)
          write(16,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & tnum_bin(i_bin,j_bin,:),trmsd_bin(i_bin,j_bin,:)

          write(17,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & unum_bin(i_bin,j_bin,:),usprd_bin(i_bin,j_bin,:)
          write(18,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & vnum_bin(i_bin,j_bin,:),vsprd_bin(i_bin,j_bin,:)
          write(19,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & tnum_bin(i_bin,j_bin,:),tsprd_bin(i_bin,j_bin,:)          
          
       end do !i_bin
    end do    !j_bin
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)

    write(format,'(a,I0,a,I0,a)') '(2f12.5,i6,',ndat_a,'i6,',2*ndat_a,'f12.5)'

    open(11,file="dat/ubias_tval_bin.dat",status="replace")
    open(12,file="dat/vbias_tval_bin.dat",status="replace")
    open(13,file="dat/tbias_tval_bin.dat",status="replace")
    open(14,file="dat/urmsd_tval_bin.dat",status="replace")
    open(15,file="dat/vrmsd_tval_bin.dat",status="replace")
    open(16,file="dat/trmsd_tval_bin.dat",status="replace")
    do j_bin=1,jm_bin-1
       do i_bin=1,im_bin-1
          do idat_a=1,ndat_a
             
             write(11,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,ubias_dof_bin(i_bin,j_bin,idat_a,:),ubias_tcrit_bin(i_bin,j_bin,idat_a,:),ubias_tval_bin(i_bin,j_bin,idat_a,:)
             write(12,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,vbias_dof_bin(i_bin,j_bin,idat_a,:),vbias_tcrit_bin(i_bin,j_bin,idat_a,:),vbias_tval_bin(i_bin,j_bin,idat_a,:)
             write(13,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,tbias_dof_bin(i_bin,j_bin,idat_a,:),tbias_tcrit_bin(i_bin,j_bin,idat_a,:),tbias_tval_bin(i_bin,j_bin,idat_a,:)

             write(14,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,urmsd_dof_bin(i_bin,j_bin,idat_a,:),urmsd_tcrit_bin(i_bin,j_bin,idat_a,:),urmsd_tval_bin(i_bin,j_bin,idat_a,:)
             write(15,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,vrmsd_dof_bin(i_bin,j_bin,idat_a,:),vrmsd_tcrit_bin(i_bin,j_bin,idat_a,:),vrmsd_tval_bin(i_bin,j_bin,idat_a,:)
             write(16,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,trmsd_dof_bin(i_bin,j_bin,idat_a,:),trmsd_tcrit_bin(i_bin,j_bin,idat_a,:),trmsd_tval_bin(i_bin,j_bin,idat_a,:)

          end do
       end do
    end do
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
          
  end subroutine write_bin

  !---------------------------------------------------

  subroutine write_ave(syr,eyr,ndat_a, &
     & unum_mave,ubias_mave,urmsd_mave,usprd_mave, &
     & vnum_mave,vbias_mave,vrmsd_mave,vsprd_mave, &
     & tnum_mave,tbias_mave,trmsd_mave,tsprd_mave, &
     & unum_yave,ubias_yave,urmsd_yave,usprd_yave, &
     & vnum_yave,vbias_yave,vrmsd_yave,vsprd_yave, &
     & tnum_yave,tbias_yave,trmsd_yave,tsprd_yave, &
     & unum_ave,ubias_ave,urmsd_ave,usprd_ave,       &
     & ubias_dof_ave,ubias_tcrit_ave,ubias_tval_ave, &
     & urmsd_dof_ave,urmsd_tcrit_ave,urmsd_tval_ave, &
     & vnum_ave,vbias_ave,vrmsd_ave,vsprd_ave,       &
     & vbias_dof_ave,vbias_tcrit_ave,vbias_tval_ave, &
     & vrmsd_dof_ave,vrmsd_tcrit_ave,vrmsd_tval_ave, &
     & tnum_ave,tbias_ave,trmsd_ave,tsprd_ave,       &
     & tbias_dof_ave,tbias_tcrit_ave,tbias_tval_ave, &
     & trmsd_dof_ave,trmsd_tcrit_ave,trmsd_tval_ave)

    
    implicit none

    !---Common
    integer iyr,imon
    integer idat_a

    character(100) format
    character(10) yyyymmdd

    !---IN  
    integer,intent(in) :: syr,eyr,ndat_a

    !Monthly
    integer,intent(in) :: unum_mave(ndat_a,12,syr:eyr),vnum_mave(ndat_a,12,syr:eyr),tnum_mave(ndat_a,12,syr:eyr)

    !Yearly
    integer,intent(in) :: unum_yave(ndat_a,syr:eyr),vnum_yave(ndat_a,syr:eyr),tnum_yave(ndat_a,syr:eyr)

    !ALL
    integer,intent(in) :: unum_ave(ndat_a),vnum_ave(ndat_a),tnum_ave(ndat_a)
    integer,intent(in) :: ubias_dof_ave(ndat_a,ndat_a),vbias_dof_ave(ndat_a,ndat_a),tbias_dof_ave(ndat_a,ndat_a)
    integer,intent(in) :: urmsd_dof_ave(ndat_a,ndat_a),vrmsd_dof_ave(ndat_a,ndat_a),trmsd_dof_ave(ndat_a,ndat_a)
    
    !Monthly
    real(kind = 8),intent(in) :: ubias_mave(ndat_a,12,syr:eyr),vbias_mave(ndat_a,12,syr:eyr),tbias_mave(ndat_a,12,syr:eyr)
    real(kind = 8),intent(in) :: urmsd_mave(ndat_a,12,syr:eyr),vrmsd_mave(ndat_a,12,syr:eyr),trmsd_mave(ndat_a,12,syr:eyr)
    real(kind = 8),intent(in) :: usprd_mave(ndat_a,12,syr:eyr),vsprd_mave(ndat_a,12,syr:eyr),tsprd_mave(ndat_a,12,syr:eyr)

    !Yearly
    real(kind = 8),intent(in) :: ubias_yave(ndat_a,syr:eyr),vbias_yave(ndat_a,syr:eyr),tbias_yave(ndat_a,syr:eyr)
    real(kind = 8),intent(in) :: urmsd_yave(ndat_a,syr:eyr),vrmsd_yave(ndat_a,syr:eyr),trmsd_yave(ndat_a,syr:eyr)
    real(kind = 8),intent(in) :: usprd_yave(ndat_a,syr:eyr),vsprd_yave(ndat_a,syr:eyr),tsprd_yave(ndat_a,syr:eyr)

    !ALL
    real(kind = 8),intent(in) :: ubias_ave(ndat_a),vbias_ave(ndat_a),tbias_ave(ndat_a)
    real(kind = 8),intent(in) :: urmsd_ave(ndat_a),vrmsd_ave(ndat_a),trmsd_ave(ndat_a)
    real(kind = 8),intent(in) :: usprd_ave(ndat_a),vsprd_ave(ndat_a),tsprd_ave(ndat_a)

    real(kind = 8),intent(in) :: ubias_tcrit_ave(ndat_a,ndat_a),vbias_tcrit_ave(ndat_a,ndat_a),tbias_tcrit_ave(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: ubias_tval_ave(ndat_a,ndat_a),vbias_tval_ave(ndat_a,ndat_a),tbias_tval_ave(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: urmsd_tcrit_ave(ndat_a,ndat_a),vrmsd_tcrit_ave(ndat_a,ndat_a),trmsd_tcrit_ave(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: urmsd_tval_ave(ndat_a,ndat_a),vrmsd_tval_ave(ndat_a,ndat_a),trmsd_tval_ave(ndat_a,ndat_a)
    
    !---Monthly
    write(format,'(a,I0,a,I0,a)') "(a,",ndat_a,"i10,",ndat_a,"f12.5)"

    open(11,file="dat/ubias_mave.dat",status="replace")
    open(12,file="dat/vbias_mave.dat",status="replace")
    open(13,file="dat/tbias_mave.dat",status="replace")
    open(14,file="dat/urmsd_mave.dat",status="replace")
    open(15,file="dat/vrmsd_mave.dat",status="replace")
    open(16,file="dat/trmsd_mave.dat",status="replace")
    open(17,file="dat/usprd_mave.dat",status="replace")
    open(18,file="dat/vsprd_mave.dat",status="replace")
    open(19,file="dat/tsprd_mave.dat",status="replace")
    do iyr=syr,eyr
       do imon=1,12

          write(yyyymmdd,'(i4.4,a,i2.2,a)') iyr,"-",imon,"-15"
          
          write(11,trim(format)) yyyymmdd,unum_mave(:,imon,iyr),ubias_mave(:,imon,iyr)
          write(12,trim(format)) yyyymmdd,vnum_mave(:,imon,iyr),vbias_mave(:,imon,iyr)
          write(13,trim(format)) yyyymmdd,tnum_mave(:,imon,iyr),tbias_mave(:,imon,iyr)

          write(14,trim(format)) yyyymmdd,unum_mave(:,imon,iyr),urmsd_mave(:,imon,iyr)
          write(15,trim(format)) yyyymmdd,vnum_mave(:,imon,iyr),vrmsd_mave(:,imon,iyr)
          write(16,trim(format)) yyyymmdd,tnum_mave(:,imon,iyr),trmsd_mave(:,imon,iyr)

          write(17,trim(format)) yyyymmdd,unum_mave(:,imon,iyr),usprd_mave(:,imon,iyr)
          write(18,trim(format)) yyyymmdd,vnum_mave(:,imon,iyr),vsprd_mave(:,imon,iyr)
          write(19,trim(format)) yyyymmdd,tnum_mave(:,imon,iyr),tsprd_mave(:,imon,iyr)
          
       end do
    end do
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)

    !---Yearly
    write(format,'(a,I0,a,I0,a)') "(i6,",ndat_a,"i10,",ndat_a,"f12.5)"

    open(11,file="dat/ubias_yave.dat",status="replace")
    open(12,file="dat/vbias_yave.dat",status="replace")
    open(13,file="dat/tbias_yave.dat",status="replace")
    open(14,file="dat/urmsd_yave.dat",status="replace")
    open(15,file="dat/vrmsd_yave.dat",status="replace")
    open(16,file="dat/trmsd_yave.dat",status="replace")
    open(17,file="dat/usprd_yave.dat",status="replace")
    open(18,file="dat/vsprd_yave.dat",status="replace")
    open(19,file="dat/tsprd_yave.dat",status="replace")
    do iyr=syr,eyr
       
       write(11,trim(format)) iyr,unum_yave(:,iyr),ubias_yave(:,iyr)
       write(12,trim(format)) iyr,vnum_yave(:,iyr),vbias_yave(:,iyr)
       write(13,trim(format)) iyr,tnum_yave(:,iyr),tbias_yave(:,iyr)

       write(14,trim(format)) iyr,unum_yave(:,iyr),urmsd_yave(:,iyr)
       write(15,trim(format)) iyr,vnum_yave(:,iyr),vrmsd_yave(:,iyr)
       write(16,trim(format)) iyr,tnum_yave(:,iyr),trmsd_yave(:,iyr)

       write(17,trim(format)) iyr,unum_yave(:,iyr),usprd_yave(:,iyr)
       write(18,trim(format)) iyr,vnum_yave(:,iyr),vsprd_yave(:,iyr)
       write(19,trim(format)) iyr,tnum_yave(:,iyr),tsprd_yave(:,iyr)
 
    end do
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)

    !---ALL
    write(format,'(a,I0,a,I0,a)') "(",ndat_a,"i10,",ndat_a,"f12.5)"

    open(11,file="dat/ubias_ave.dat",status="replace")
    open(12,file="dat/vbias_ave.dat",status="replace")
    open(13,file="dat/tbias_ave.dat",status="replace")
    open(14,file="dat/urmsd_ave.dat",status="replace")
    open(15,file="dat/vrmsd_ave.dat",status="replace")
    open(16,file="dat/trmsd_ave.dat",status="replace")
    open(17,file="dat/usprd_ave.dat",status="replace")
    open(18,file="dat/vsprd_ave.dat",status="replace")
    open(19,file="dat/tsprd_ave.dat",status="replace")

    write(11,trim(format)) unum_ave(:),ubias_ave(:)
    write(12,trim(format)) vnum_ave(:),vbias_ave(:)
    write(13,trim(format)) tnum_ave(:),tbias_ave(:)

    write(14,trim(format)) unum_ave(:),urmsd_ave(:)
    write(15,trim(format)) vnum_ave(:),vrmsd_ave(:)
    write(16,trim(format)) tnum_ave(:),trmsd_ave(:)

    write(17,trim(format)) unum_ave(:),usprd_ave(:)
    write(18,trim(format)) vnum_ave(:),vsprd_ave(:)
    write(19,trim(format)) tnum_ave(:),tsprd_ave(:)
    
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)

    write(format,'(a,I0,a,I0,a,I0,a)') "(i6,",ndat_a,"i6,",ndat_a,"f12.5,",ndat_a,"f12.5)"

    open(11,file="dat/ubias_tval_ave.dat",status="replace")
    open(12,file="dat/vbias_tval_ave.dat",status="replace")
    open(13,file="dat/tbias_tval_ave.dat",status="replace")
    open(14,file="dat/urmsd_tval_ave.dat",status="replace")
    open(15,file="dat/vrmsd_tval_ave.dat",status="replace")
    open(16,file="dat/trmsd_tval_ave.dat",status="replace")
    do idat_a=1,ndat_a
       
       write(11,trim(format)) idat_a,ubias_dof_ave(idat_a,:),ubias_tcrit_ave(idat_a,:),ubias_tval_ave(idat_a,:)
       write(12,trim(format)) idat_a,vbias_dof_ave(idat_a,:),vbias_tcrit_ave(idat_a,:),vbias_tval_ave(idat_a,:)
       write(13,trim(format)) idat_a,tbias_dof_ave(idat_a,:),tbias_tcrit_ave(idat_a,:),tbias_tval_ave(idat_a,:)

       write(14,trim(format)) idat_a,urmsd_dof_ave(idat_a,:),urmsd_tcrit_ave(idat_a,:),urmsd_tval_ave(idat_a,:)
       write(15,trim(format)) idat_a,vrmsd_dof_ave(idat_a,:),vrmsd_tcrit_ave(idat_a,:),vrmsd_tval_ave(idat_a,:)
       write(16,trim(format)) idat_a,trmsd_dof_ave(idat_a,:),trmsd_tcrit_ave(idat_a,:),trmsd_tval_ave(idat_a,:)
       
    end do
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    
  end subroutine write_ave

  
end module mod_io
