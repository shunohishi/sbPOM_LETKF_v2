module mod_io

contains

  !----------------------------------------------------------------------------
  ! Read Argument |
  !----------------------------------------------------------------------------

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
  ! Grid size (Analysis data) |
  !---------------------------------------------------------------------------------

  ! *** To be modified ***
  subroutine get_grid_size(idat,im,jm,km)

    use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
    use mod_read_bran2020,   only: im_bran => im, jm_bran => jm, km_bran => km
    use mod_read_glorys12v1, only: im_g010 => im, jm_g010 => jm, km_g010 => km
    use mod_read_glorys025,  only: im_g025 => im, jm_g025 => jm, km_g025 => km
    use mod_read_jcope_fgo,  only: im_jcope => im, jm_jcope => jm, km_jcope => km   
    implicit none

    !---IN
    integer,intent(in) :: idat

    !---OUT
    integer,intent(out) :: im,jm,km

     if(idat == 1)then !---LORA-NP
        im=im_lora
        jm=jm_lora
        km=km_lora
     else if(idat == 2)then !---BRAN2020
        im=im_bran
        jm=jm_bran
        km=km_bran
     else if(idat == 3)then !---GLORYS010
        im=im_g010
        jm=jm_g010
        km=km_g010
     else if(idat == 4)then !---JCOPE-FGO
        im=im_jcope
        jm=jm_jcope
        km=km_jcope
     else
        write(*,*) "***Error: Incorecot idat_a => ",idat
        stop
     end if        

  end subroutine get_grid_size
  
  !---------------------------------------------------------------------------------
  ! Read Grid data |
  !---------------------------------------------------------------------------------

  ! ***To be modified ***
  subroutine read_grid(idat,im,jm,km,lont,lonu,lonv,latt,latu,latv,maskt,masku,maskv)

    use setting, only: datname
    use mod_read_lora,       only: read_grid_lora => read_grid
    use mod_read_bran2020,   only: read_bran2020
    use mod_read_glorys12v1, only: read_glorys12v1
    use mod_read_jcope_fgo,  only: read_grid_jcope => read_grid
    implicit none

    !---Common
    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) dir

    !Others
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

    if(idat == 1)then !--- LORA-NP
       dir="NP"
       call read_grid_lora(dir,lont,lonu,lonv, &
            & latt,latu,latv, &
            & tmp3d,tmp3d,tmp3d,tmp3d, &
            & maskt,masku,maskv)
    else if(idat == 2)then !--- BRAN
       varname="t"
       call read_bran2020(varname,2003,1,1,km,lont,latt,tmp1dz,maskt,tmp3d)
       varname="u"
       call read_bran2020(varname,2003,1,1,km,lonu,latu,tmp1dz,masku,tmp3d)
       varname="v"
       call read_bran2020(varname,2003,1,1,km,lonv,latv,tmp1dz,maskv,tmp3d)
    else if(idat == 3)then !---GLORYS12V1 (Probably Regridded)
       varname="t"
       call read_glorys12v1(varname,2003,1,1,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,tmp3d)
       lont(:)=tmp1dx(:)
       lonu(:)=tmp1dx(:)
       lonv(:)=tmp1dx(:)
       latt(:)=tmp1dy(:)
       latu(:)=tmp1dy(:)
       latv(:)=tmp1dy(:)
       maskt(:,:)=tmp2d(:,:)
       masku(:,:)=tmp2d(:,:)
       maskv(:,:)=tmp2d(:,:)
    else if(idat == 4)then !---JCOPE-FGO
       call read_grid_jcope(lont,lonu,lonv,latt,latu,latv,tmp2d,tmp3d)
       maskt(:,:)=tmp2d(:,:)
       masku(:,:)=tmp2d(:,:)
       maskv(:,:)=tmp2d(:,:)
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)
    deallocate(tmp3d)

  end subroutine read_grid

  !---------------------------------------------------------------------------------
  ! Read sea-surface analysis data |
  !---------------------------------------------------------------------------------

  !*** To be modified ***
  subroutine read_surface_data(idat,iyr,imon,iday,im,jm,km,maskt,masku,maskv,t,u,v,tsprd,usprd,vsprd)

    use setting, only: datname
    use mod_read_lora,       only: read_anal
    use mod_read_bran2020,   only: read_bran2020
    use mod_read_glorys12v1, only: read_glorys12v1
    use mod_read_jcope_fgo,  only: read_jcope_fgo
    use mod_rmiss
    implicit none

    !---Common
    integer k

    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    integer imem !Dummy
    character(10) dir
    character(10) letkf
    character(10) region
    character(10) ms

    !Others
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
    allocate(tmp3d(im,jm,km))

    if(idat == 1)then !---LORA-NP
       imem=0
       k=1
       dir="NP"
       letkf="letkf"
       region="np"       
       ms="mean"
       call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im,jm,k,maskt,t)
       call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im,jm,k,masku,u)
       call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im,jm,k,maskv,v)
       ms="sprd"
       call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im,jm,k,maskt,tsprd)
       call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im,jm,k,masku,usprd)
       call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im,jm,k,maskv,vsprd)
    else if(idat == 2)then !---BRAN
       k=1
       varname="t"       
       call read_bran2020(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,t)
       varname="u"       
       call read_bran2020(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,u)
       varname="v"       
       call read_bran2020(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss
    else if(idat == 3)then !---GLORYS12V1
       k=1
       varname="t"
       call read_glorys12v1(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,t)
       varname="u"
       call read_glorys12v1(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,u)
       varname="v"
       call read_glorys12v1(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss
    else if(idat == 4)then !---JCOPE-FGO
       k=1
       varname="t"
       call read_jcope_fgo(varname,iyr,imon,iday,k,maskt,t)
       varname="u"
       call read_jcope_fgo(varname,iyr,imon,iday,k,masku,u)
       varname="v"
       call read_jcope_fgo(varname,iyr,imon,iday,k,maskv,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)
    deallocate(tmp3d)

  end subroutine read_surface_data

  !---------------------------------------------------------------------------------
  ! Read data in observation space |
  !---------------------------------------------------------------------------------

  subroutine read_obs(idat_a,iyr,imon,iday,nobs,lon_o,lat_o, &
       & ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o,u_o,v_o)

    use netcdf
    use setting, only: ndat_a, datname
    implicit none

    !---Common
    integer status,access
    integer ncid,dimid,varid

    character(100) filename
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

    if(idat_a < 1 .or. ndat_a < idat_a)then
       write(*,*) "***Error: Incorrect datname"
       stop
    end if

    filename="dat/"//yyyy//mm//"/"//trim(datname(idat_a))//"."//yyyymmdd//".nc"    

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

  !---------------------------------------------------------------------------------
  ! Write data in observation space |
  !---------------------------------------------------------------------------------

  subroutine write_obs(idat_a,ijul,nobs,ijul_o,lon_o,lat_o, &
       & ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o,u_o,v_o, &
       & ncid,inum)

    use netcdf
    use mod_julian
    use mod_rmiss
    use setting, only: ndat_a, datname
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

    if(idat_a < 1 .or. ndat_a < idat_a)then
       write(*,*) "***Error: Incorrect datname"
       stop
    end if

    filename="dat/"//yyyy//mm//"/"//trim(datname(idat_a))//"."//yyyymmdd//".nc"

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

  !---------------------------------------------------------------------------------
  ! Write bin data |
  !---------------------------------------------------------------------------------

  subroutine write_bin(im_bin,jm_bin,ndat_a,dx_bin,dy_bin,lon_bin,lat_bin, &
       & unum_stat_bin,unum_sprd_bin,ubias_bin,urmsd_bin,usprd_bin,ucor_bin, &
       & uabias_dif_low_bin,uabias_dif_ave_bin,uabias_dif_upp_bin, &
       & urmsd_dif_low_bin,urmsd_dif_ave_bin,urmsd_dif_upp_bin, &
       & vnum_stat_bin,vnum_sprd_bin,vbias_bin,vrmsd_bin,vsprd_bin,vcor_bin, &
       & vabias_dif_low_bin,vabias_dif_ave_bin,vabias_dif_upp_bin, &
       & vrmsd_dif_low_bin,vrmsd_dif_ave_bin,vrmsd_dif_upp_bin, &
       & tnum_stat_bin,tnum_sprd_bin,tbias_bin,trmsd_bin,tsprd_bin,tcor_bin, &
       & tabias_dif_low_bin,tabias_dif_ave_bin,tabias_dif_upp_bin, &
       & trmsd_dif_low_bin,trmsd_dif_ave_bin,trmsd_dif_upp_bin)

    implicit none

    !---Common
    integer i_bin,j_bin
    integer idat_a

    integer unum_min_bin(im_bin,jm_bin,ndat_a)
    integer vnum_min_bin(im_bin,jm_bin,ndat_a)
    integer tnum_min_bin(im_bin,jm_bin,ndat_a)
    
    character(100) format

    !---IN
    integer,intent(in) :: im_bin,jm_bin,ndat_a
    real(kind = 8),intent(in) :: dx_bin,dy_bin

    integer,intent(in) :: unum_stat_bin(im_bin,jm_bin,ndat_a)
    integer,intent(in) :: vnum_stat_bin(im_bin,jm_bin,ndat_a)
    integer,intent(in) :: tnum_stat_bin(im_bin,jm_bin,ndat_a)

    integer,intent(in) :: unum_sprd_bin(im_bin,jm_bin,ndat_a)
    integer,intent(in) :: vnum_sprd_bin(im_bin,jm_bin,ndat_a)
    integer,intent(in) :: tnum_sprd_bin(im_bin,jm_bin,ndat_a)

    real(kind = 8),intent(in) :: lon_bin(im_bin),lat_bin(jm_bin)    

    real(kind = 8),intent(in) :: ubias_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: urmsd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: usprd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: ucor_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: uabias_dif_low_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: uabias_dif_ave_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: uabias_dif_upp_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: urmsd_dif_low_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: urmsd_dif_ave_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: urmsd_dif_upp_bin(im_bin,jm_bin,ndat_a,ndat_a)

    real(kind = 8),intent(in) :: vbias_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: vrmsd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: vsprd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: vcor_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: vabias_dif_low_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: vabias_dif_ave_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: vabias_dif_upp_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: vrmsd_dif_low_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: vrmsd_dif_ave_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: vrmsd_dif_upp_bin(im_bin,jm_bin,ndat_a,ndat_a)

    real(kind = 8),intent(in) :: tbias_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: trmsd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: tsprd_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: tcor_bin(im_bin,jm_bin,ndat_a)
    real(kind = 8),intent(in) :: tabias_dif_low_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: tabias_dif_ave_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: tabias_dif_upp_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: trmsd_dif_low_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: trmsd_dif_ave_bin(im_bin,jm_bin,ndat_a,ndat_a)
    real(kind = 8),intent(in) :: trmsd_dif_upp_bin(im_bin,jm_bin,ndat_a,ndat_a)

    do idat_a=1,ndat_a
       do j_bin=1,jm_bin
          do i_bin=1,im_bin
             unum_min_bin(i_bin,j_bin,idat_a)=min(unum_stat_bin(i_bin,j_bin,idat_a),unum_sprd_bin(i_bin,j_bin,idat_a))
             vnum_min_bin(i_bin,j_bin,idat_a)=min(vnum_stat_bin(i_bin,j_bin,idat_a),vnum_sprd_bin(i_bin,j_bin,idat_a))
             tnum_min_bin(i_bin,j_bin,idat_a)=min(tnum_stat_bin(i_bin,j_bin,idat_a),tnum_sprd_bin(i_bin,j_bin,idat_a))
          end do
       end do
    end do
    
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
    open(20,file="dat/ucor_bin.dat",status="replace")
    open(21,file="dat/vcor_bin.dat",status="replace")
    open(22,file="dat/tcor_bin.dat",status="replace")

    do j_bin=1,jm_bin-1
       do i_bin=1,im_bin-1

          write(11,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & unum_stat_bin(i_bin,j_bin,:),ubias_bin(i_bin,j_bin,:)
          write(12,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & vnum_stat_bin(i_bin,j_bin,:),vbias_bin(i_bin,j_bin,:)
          write(13,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & tnum_stat_bin(i_bin,j_bin,:),tbias_bin(i_bin,j_bin,:)

          write(14,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & unum_stat_bin(i_bin,j_bin,:),urmsd_bin(i_bin,j_bin,:)
          write(15,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & vnum_stat_bin(i_bin,j_bin,:),vrmsd_bin(i_bin,j_bin,:)
          write(16,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & tnum_stat_bin(i_bin,j_bin,:),trmsd_bin(i_bin,j_bin,:)

          write(17,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & unum_sprd_bin(i_bin,j_bin,:),usprd_bin(i_bin,j_bin,:)
          write(18,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & vnum_sprd_bin(i_bin,j_bin,:),vsprd_bin(i_bin,j_bin,:)
          write(19,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & tnum_sprd_bin(i_bin,j_bin,:),tsprd_bin(i_bin,j_bin,:)          

          write(20,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & unum_min_bin(i_bin,j_bin,:),ucor_bin(i_bin,j_bin,:)
          write(21,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & vnum_min_bin(i_bin,j_bin,:),vcor_bin(i_bin,j_bin,:)
          write(22,trim(format)) &
               & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
               & tnum_min_bin(i_bin,j_bin,:),tcor_bin(i_bin,j_bin,:)                    
          
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
    close(20)
    close(21)
    close(22)

    write(format,'(a,I0,a,I0,a)') '(2f12.5,i6,',ndat_a,'i10,',2*ndat_a,'f12.5)'

    open(11,file="dat/ubias_dif_bin.dat",status="replace")
    open(12,file="dat/vbias_dif_bin.dat",status="replace")
    open(13,file="dat/tbias_dif_bin.dat",status="replace")
    open(14,file="dat/urmsd_dif_bin.dat",status="replace")
    open(15,file="dat/vrmsd_dif_bin.dat",status="replace")
    open(16,file="dat/trmsd_dif_bin.dat",status="replace")
    do j_bin=1,jm_bin-1
       do i_bin=1,im_bin-1
          do idat_a=1,ndat_a

             write(11,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,unum_stat_bin(i_bin,j_bin,:),uabias_dif_low_bin(i_bin,j_bin,idat_a,:),uabias_dif_upp_bin(i_bin,j_bin,idat_a,:)
             write(12,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,vnum_stat_bin(i_bin,j_bin,:),vabias_dif_low_bin(i_bin,j_bin,idat_a,:),vabias_dif_upp_bin(i_bin,j_bin,idat_a,:)
             write(13,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,tnum_stat_bin(i_bin,j_bin,:),tabias_dif_low_bin(i_bin,j_bin,idat_a,:),tabias_dif_upp_bin(i_bin,j_bin,idat_a,:)

             write(14,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,unum_stat_bin(i_bin,j_bin,:),urmsd_dif_low_bin(i_bin,j_bin,idat_a,:),urmsd_dif_upp_bin(i_bin,j_bin,idat_a,:)
             write(15,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,vnum_stat_bin(i_bin,j_bin,:),vrmsd_dif_low_bin(i_bin,j_bin,idat_a,:),vrmsd_dif_upp_bin(i_bin,j_bin,idat_a,:)
             write(16,trim(format)) &
                  & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
                  & idat_a,tnum_stat_bin(i_bin,j_bin,:),trmsd_dif_low_bin(i_bin,j_bin,idat_a,:),trmsd_dif_upp_bin(i_bin,j_bin,idat_a,:)

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

  !---------------------------------------------------------------------------------
  ! Write static data |
  !---------------------------------------------------------------------------------

  subroutine write_ave(syr,eyr,ndat_a, &
       & unum_stat_mave,unum_sprd_mave,ubias_mave,urmsd_mave,usprd_mave, &
       & vnum_stat_mave,vnum_sprd_mave,vbias_mave,vrmsd_mave,vsprd_mave, &
       & tnum_stat_mave,tnum_sprd_mave,tbias_mave,trmsd_mave,tsprd_mave, &
       & unum_stat_yave,unum_sprd_yave,ubias_yave,urmsd_yave,usprd_yave, &
       & vnum_stat_yave,vnum_sprd_yave,vbias_yave,vrmsd_yave,vsprd_yave, &
       & tnum_stat_yave,tnum_sprd_yave,tbias_yave,trmsd_yave,tsprd_yave, &
       & unum_stat_ave,unum_sprd_ave,ubias_ave,urmsd_ave,usprd_ave,ucor_ave, &
       & uabias_dif_low_ave,uabias_dif_ave_ave,uabias_dif_upp_ave, &
       & urmsd_dif_low_ave,urmsd_dif_ave_ave,urmsd_dif_upp_ave, &
       & vnum_stat_ave,vnum_sprd_ave,vbias_ave,vrmsd_ave,vsprd_ave,vcor_ave, &
       & vabias_dif_low_ave,vabias_dif_ave_ave,vabias_dif_upp_ave, &
       & vrmsd_dif_low_ave,vrmsd_dif_ave_ave,vrmsd_dif_upp_ave, &
       & tnum_stat_ave,tnum_sprd_ave,tbias_ave,trmsd_ave,tsprd_ave,tcor_ave, &
       & tabias_dif_low_ave,tabias_dif_ave_ave,tabias_dif_upp_ave, &
       & trmsd_dif_low_ave,trmsd_dif_ave_ave,trmsd_dif_upp_ave)

    implicit none

    !---Common
    integer iyr,imon
    integer idat_a

    integer unum_min_ave(ndat_a),vnum_min_ave(ndat_a),tnum_min_ave(ndat_a)    
    
    character(100) format
    character(10) yyyymmdd

    !---IN  
    integer,intent(in) :: syr,eyr,ndat_a

    !Monthly
    integer,intent(in) :: unum_stat_mave(ndat_a,12,syr:eyr),vnum_stat_mave(ndat_a,12,syr:eyr),tnum_stat_mave(ndat_a,12,syr:eyr)
    integer,intent(in) :: unum_sprd_mave(ndat_a,12,syr:eyr),vnum_sprd_mave(ndat_a,12,syr:eyr),tnum_sprd_mave(ndat_a,12,syr:eyr)

    !Yearly
    integer,intent(in) :: unum_stat_yave(ndat_a,syr:eyr),vnum_stat_yave(ndat_a,syr:eyr),tnum_stat_yave(ndat_a,syr:eyr)
    integer,intent(in) :: unum_sprd_yave(ndat_a,syr:eyr),vnum_sprd_yave(ndat_a,syr:eyr),tnum_sprd_yave(ndat_a,syr:eyr)

    !ALL
    integer,intent(in) :: unum_stat_ave(ndat_a),vnum_stat_ave(ndat_a),tnum_stat_ave(ndat_a)
    integer,intent(in) :: unum_sprd_ave(ndat_a),vnum_sprd_ave(ndat_a),tnum_sprd_ave(ndat_a)

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
    real(kind = 8),intent(in) :: ucor_ave(ndat_a),vcor_ave(ndat_a),tcor_ave(ndat_a)

    real(kind = 8),intent(in) :: uabias_dif_low_ave(ndat_a,ndat_a),vabias_dif_low_ave(ndat_a,ndat_a),tabias_dif_low_ave(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: uabias_dif_ave_ave(ndat_a,ndat_a),vabias_dif_ave_ave(ndat_a,ndat_a),tabias_dif_ave_ave(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: uabias_dif_upp_ave(ndat_a,ndat_a),vabias_dif_upp_ave(ndat_a,ndat_a),tabias_dif_upp_ave(ndat_a,ndat_a)

    real(kind = 8),intent(in) :: urmsd_dif_low_ave(ndat_a,ndat_a),vrmsd_dif_low_ave(ndat_a,ndat_a),trmsd_dif_low_ave(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: urmsd_dif_ave_ave(ndat_a,ndat_a),vrmsd_dif_ave_ave(ndat_a,ndat_a),trmsd_dif_ave_ave(ndat_a,ndat_a)
    real(kind = 8),intent(in) :: urmsd_dif_upp_ave(ndat_a,ndat_a),vrmsd_dif_upp_ave(ndat_a,ndat_a),trmsd_dif_upp_ave(ndat_a,ndat_a)

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

          write(11,trim(format)) yyyymmdd,unum_stat_mave(:,imon,iyr),ubias_mave(:,imon,iyr)
          write(12,trim(format)) yyyymmdd,vnum_stat_mave(:,imon,iyr),vbias_mave(:,imon,iyr)
          write(13,trim(format)) yyyymmdd,tnum_stat_mave(:,imon,iyr),tbias_mave(:,imon,iyr)

          write(14,trim(format)) yyyymmdd,unum_stat_mave(:,imon,iyr),urmsd_mave(:,imon,iyr)
          write(15,trim(format)) yyyymmdd,vnum_stat_mave(:,imon,iyr),vrmsd_mave(:,imon,iyr)
          write(16,trim(format)) yyyymmdd,tnum_stat_mave(:,imon,iyr),trmsd_mave(:,imon,iyr)

          write(17,trim(format)) yyyymmdd,unum_sprd_mave(:,imon,iyr),usprd_mave(:,imon,iyr)
          write(18,trim(format)) yyyymmdd,vnum_sprd_mave(:,imon,iyr),vsprd_mave(:,imon,iyr)
          write(19,trim(format)) yyyymmdd,tnum_sprd_mave(:,imon,iyr),tsprd_mave(:,imon,iyr)

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

       write(11,trim(format)) iyr,unum_stat_yave(:,iyr),ubias_yave(:,iyr)
       write(12,trim(format)) iyr,vnum_stat_yave(:,iyr),vbias_yave(:,iyr)
       write(13,trim(format)) iyr,tnum_stat_yave(:,iyr),tbias_yave(:,iyr)

       write(14,trim(format)) iyr,unum_stat_yave(:,iyr),urmsd_yave(:,iyr)
       write(15,trim(format)) iyr,vnum_stat_yave(:,iyr),vrmsd_yave(:,iyr)
       write(16,trim(format)) iyr,tnum_stat_yave(:,iyr),trmsd_yave(:,iyr)

       write(17,trim(format)) iyr,unum_sprd_yave(:,iyr),usprd_yave(:,iyr)
       write(18,trim(format)) iyr,vnum_sprd_yave(:,iyr),vsprd_yave(:,iyr)
       write(19,trim(format)) iyr,tnum_sprd_yave(:,iyr),tsprd_yave(:,iyr)

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
    do idat_a=1,ndat_a
       unum_min_ave(idat_a)=min(unum_stat_ave(idat_a),unum_sprd_ave(idat_a))
       vnum_min_ave(idat_a)=min(vnum_stat_ave(idat_a),vnum_sprd_ave(idat_a))
       tnum_min_ave(idat_a)=min(tnum_stat_ave(idat_a),tnum_sprd_ave(idat_a))
    end do
    
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
    open(20,file="dat/ucor_ave.dat",status="replace")
    open(21,file="dat/vcor_ave.dat",status="replace")
    open(22,file="dat/tcor_ave.dat",status="replace")

    write(11,trim(format)) unum_stat_ave(:),ubias_ave(:)
    write(12,trim(format)) vnum_stat_ave(:),vbias_ave(:)
    write(13,trim(format)) tnum_stat_ave(:),tbias_ave(:)

    write(14,trim(format)) unum_stat_ave(:),urmsd_ave(:)
    write(15,trim(format)) vnum_stat_ave(:),vrmsd_ave(:)
    write(16,trim(format)) tnum_stat_ave(:),trmsd_ave(:)

    write(17,trim(format)) unum_sprd_ave(:),usprd_ave(:)
    write(18,trim(format)) vnum_sprd_ave(:),vsprd_ave(:)
    write(19,trim(format)) tnum_sprd_ave(:),tsprd_ave(:)

    write(20,trim(format)) unum_min_ave(:),ucor_ave(:)
    write(21,trim(format)) vnum_min_ave(:),vcor_ave(:)
    write(22,trim(format)) tnum_min_ave(:),tcor_ave(:)
    
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)

    write(format,'(a,I0,a)') "(i6,",2*ndat_a,"f12.5)"

    open(11,file="dat/ubias_dif_ave.dat",status="replace")
    open(12,file="dat/vbias_dif_ave.dat",status="replace")
    open(13,file="dat/tbias_dif_ave.dat",status="replace")
    open(14,file="dat/urmsd_dif_ave.dat",status="replace")
    open(15,file="dat/vrmsd_dif_ave.dat",status="replace")
    open(16,file="dat/trmsd_dif_ave.dat",status="replace")
    do idat_a=1,ndat_a

       write(11,trim(format)) idat_a,uabias_dif_low_ave(idat_a,:),uabias_dif_upp_ave(idat_a,:)
       write(12,trim(format)) idat_a,vabias_dif_low_ave(idat_a,:),vabias_dif_upp_ave(idat_a,:)
       write(13,trim(format)) idat_a,tabias_dif_low_ave(idat_a,:),tabias_dif_upp_ave(idat_a,:)

       write(14,trim(format)) idat_a,urmsd_dif_low_ave(idat_a,:),urmsd_dif_upp_ave(idat_a,:)
       write(15,trim(format)) idat_a,vrmsd_dif_low_ave(idat_a,:),vrmsd_dif_upp_ave(idat_a,:)
       write(16,trim(format)) idat_a,trmsd_dif_low_ave(idat_a,:),trmsd_dif_upp_ave(idat_a,:)

    end do
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)

  end subroutine write_ave

end module mod_io
