module mod_read_gtspp

  !Modified S.Ohishi 2020.05
  !Modified S.Ohishi 2025.02
  !Modified S.Ohishi 2025.07

  character(100),parameter :: gtspp_dir="/data/R/R2402/DATA/GTSPP"

contains

  !-----------------------------------------------------------
  ! Make & Read filename |
  !-----------------------------------------------------------
  
  subroutine make_filename(iyr,imon)

    implicit none

    integer jyr,jmon,jday
    integer status,access,system
    integer ifile,nfile

    real(kind = 8) long,lati

    character(100) filename
    character(8) yyyymmdd
    character(6) yyyymm
    character(4) yyyy
    character(2) mm,dd

    !---IN
    integer,intent(in) :: iyr,imon
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    yyyymm=yyyy//mm

    status=access(trim(gtspp_dir)//"/filename/"//yyyymm//".txt"," ")
    if(status == 0)then
       write(*,'(a)') "Exist "//trim(gtspp_dir)//"/filename/"//yyyymm//".txt"
       return
    end if

    write(*,'(a)') "Start: Make monthly filename list at "//yyyymm
    status=system("ls "//trim(gtspp_dir)//"/"//yyyymm// &
         & " > "//trim(gtspp_dir)//"/filename/"//yyyymm//".txt")
    write(*,'(a)') "End: Make monthly filename list"//yyyymm

    nfile=0
    open(1,file=trim(gtspp_dir)//"/filename/"//yyyymm//".txt",status="old")
    do 
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,'(a,i10)') "Total number of file:",nfile

    write(*,'(a)') "Start: Make daily filename list at "//yyyymm
    open(1,file=trim(gtspp_dir)//"/filename/"//yyyymm//".txt",status="old")
    do ifile=1,nfile

       if(mod(ifile,100) == 1) write(*,'(i10,a,i10)') ifile,"/",nfile
       read(1,'(a)') filename
       call read_info(filename,iyr,imon,jyr,jmon,jday,long,lati)

       if(iyr == jyr .and. imon == jmon)then
          write(yyyy,'(i4.4)') jyr
          write(mm,'(i2.2)') jmon
          write(dd,'(i2.2)') jday
          yyyymmdd=yyyy//mm//dd
          open(11, &
               & file=trim(gtspp_dir)//"/filename/"//yyyymmdd//".txt", &
               & access="append")
          write(11,'(2f12.5,x,a)') long,lati,trim(filename)
          close(11)
       else
          write(*,'(a)') "***Error: Inconsistent date"
          write(*,'(3i6)') iyr,imon
          write(*,'(3i6)') jyr,jmon,jday
          stop
       end if
    end do
    close(1)

    write(*,'(a)') "End: Make daily filename list at "//yyyymm
    
  end subroutine make_filename

  !-------------------------

  subroutine read_filename(iyr,imon,iday,nfile,long,lati,filename)

    implicit none

    integer ifile
    integer status,access

    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd

    !---IN
    integer,intent(in) :: iyr,imon,iday

    !---OUT
    integer,intent(out) :: nfile
    real(kind = 8),allocatable,intent(out) :: long(:),lati(:)
    character(100),allocatable,intent(out) :: filename(:)


    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd

    status=access(trim(gtspp_dir)//"/filename/"//yyyymmdd//".txt"," ")
    if(status /= 0)then
       write(*,'(a)') "***Error: Not found "//trim(gtspp_dir)//"/filename/"&
            &//yyyymmdd//".txt"
       stop
    end if

    nfile=0
    open(1,file=trim(gtspp_dir)//"/filename/"//yyyymmdd//".txt",status="old")
    do 
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,'(a,i10)') "The number of file:",nfile

    if(nfile == 0)then
       return
    else
       
       allocate(long(nfile),lati(nfile),filename(nfile))

       open(1,file=trim(gtspp_dir)//"/filename/"//yyyymmdd//".txt",status="old")
       do ifile=1,nfile
          read(1,'(2f12.5,x,a)') long(ifile),lati(ifile),filename(ifile)
       end do
       close(1)

    end if

  end subroutine read_filename

  !-----------------------------------------------------------------
  ! Read information(date, position)
  !-----------------------------------------------------------------

  subroutine read_info(filename,iyr,imon,jyr,jmon,jday,long,lati)

    use netcdf
    implicit none

    integer status,access
    integer ncid,varid
    integer ijul

    real(kind = 4) rlon,rlat
    real(kind = 8) :: dtime

    character(200) fullfilename
    character(4) yyyy
    character(2) mm

    !IN
    integer,intent(in) :: iyr,imon
    character(100),intent(in) :: filename

    !OUT
    integer,intent(out) :: jyr,jmon,jday
    real(kind = 8),intent(out) :: long,lati

    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    fullfilename=trim(gtspp_dir)//"/"//yyyy//mm//"/"//trim(filename)

    status=access(trim(fullfilename)," ")
    if(status /= 0)then
       write(*,'(a)') "***Error: Not found "//trim(fullfilename)
       stop
    end if

    status=nf90_open(trim(fullfilename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"time",varid)
    status=nf90_get_var(ncid,varid,dtime)

    status=nf90_inq_varid(ncid,"longitude",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,rlat)

    status=nf90_close(ncid)

    long=dble(rlon)
    lati=dble(rlat)
    
    if(long < 0.d0)then
       long=long+360.d0
    end if

    call ymd_julian(1900,1,1,ijul)
    call julian_ymd(int(dtime)+ijul,jyr,jmon,jday)

  end subroutine read_info

  !----------------------------------------------------------------
  ! Read data |
  !------------
  !
  ! <flag>
  ! 0: No QC done, 1: Good data, 2:Probably good data
  ! 3: Probably bad data, 4: bad data
  ! 5: Changed, 6-8 Reserved, 9: Element missing
  !
  !----------------------------------------------------------------
  
  subroutine read_gtspp(iyr,imon,iday,filename,km,long,lati,depth,t,s)

    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: iflag=1    
    real(kind = 4),parameter :: dmiss=99999.e0

    !---Common
    integer status,access
    integer ncid,dimid,varid
    integer jyr,jmon,jday
    integer ifile,nfile
    integer k

    integer,allocatable :: dflag(:),tflag(:),sflag(:)

    real(kind = 4),allocatable :: rdep(:),rt(:),rs(:)
    
    character(200) fullfilename
    character(4) yyyy
    character(2) mm

    !---IN
    integer,intent(in) :: iyr,imon,iday
    character(100),intent(in) :: filename

    !---OUT
    integer,intent(out) :: km
    real(kind = 8),allocatable,intent(out) :: depth(:),t(:),s(:)

    !---INOUT
    real(kind = 8),intent(inout) :: long,lati
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    fullfilename=trim(gtspp_dir)//"/"//yyyy//mm//"/"//trim(filename)

    status=access(trim(fullfilename)," ")
    if(status /= 0)then
       write(*,'(a)') "***Error: Not Found "//trim(fullfilename)
       stop
    end if

    call read_info(filename,iyr,imon,jyr,jmon,jday,long,lati)

    status=nf90_open(trim(fullfilename),nf90_nowrite,ncid)

    !Allocate
    status=nf90_inq_dimid(ncid,"z",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = km)
    allocate(rdep(km),rt(km),rs(km))
    allocate(depth(km),t(km),s(km))
    allocate(dflag(km),tflag(km),sflag(km))

    !Read data
    !depth
    status=nf90_inq_varid(ncid,"z",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,rdep)
       depth(:)=dble(rdep(:))
    else
       depth(:)=rmiss
    end if

    status=nf90_inq_varid(ncid,"z_variable_quality_flag",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,dflag)
    else
       dflag(:)=rmiss
    end if

    !temperature
    status=nf90_inq_varid(ncid,"temperature",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,rt,(/1,1,1,1/),(/1,1,km,1/))
       t(:)=dble(rt(:))
    else
       t(:)=rmiss
    end if

    status=nf90_inq_varid(ncid,"temperature_quality_flag",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,tflag)
    else
       tflag(:)=rmiss
    end if

    !salinity
    status=nf90_inq_varid(ncid,"salinity",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,rs,(/1,1,1,1/),(/1,1,km,1/))
       s(:)=dble(rs(:))
    else
       s(:)=rmiss
    end if

    status=nf90_inq_varid(ncid,"salinity_quality_flag",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,sflag)
    else
       sflag(:)=rmiss
    end if

    status=nf90_close(ncid)

    !QC
    do k=1,km
       if(tflag(k) /= iflag) t(k)=rmiss
       if(sflag(k) /= iflag) s(k)=rmiss
    end do

    deallocate(rdep,rt,rs)
    deallocate(dflag,tflag,sflag)

  end subroutine read_gtspp

  !-----------------------------------------------------------
  ! Observation type |
  !-------------------
  !
  ! ins:
  !
  ! -1=XBT/MBT
  ! -2=XCTD
  ! -3=CTD
  ! -4=Ship (Not used in this subroutine)
  ! -5=Profiling float (Floating buoy)
  ! -6=Drifter Buoy
  ! -7=Mooring Buoy
  ! -8=Animal
  ! 0=Others
  !
  !-----------------------------------------------------------

  subroutine read_gtspp_obs_type(iyr,imon,filename,ins)

    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: nxbt=39
    integer,parameter :: nxctd=10
    integer,parameter :: nctd=14
    integer,parameter :: npf=52
    integer,parameter :: nb=1
    integer,parameter :: na=2
    integer,parameter :: no=4

    integer,parameter :: ndb=9
    integer,parameter :: nfb=13
    integer,parameter :: nmb=17
    
    !---Common
    integer i,j
    integer status
    integer ncid,dimid,varid
    integer ns
    integer itype

    character(200) fullfilename
    character(4) yyyy
    character(2) mm

    character(10) gticode
    character(10),allocatable,dimension(:) :: cparm
    character(4),allocatable,dimension(:) :: pcode

    
    !---PEQ$
    integer :: xbt(nxbt)=(/1, 2, 11, 21, 31, 32, 41, 42, 51, 52, 61, 71, 81, &
         & 201, 202, 211, 212, 221, 222, 231, 241, 251, 252, 261, &
         & 401, 411, 421, 431, 441, 451, 461, 462, 471, 481, 491, 501, 510, 800, 900/)
    integer :: xctd(nxctd)=(/700, 710, 720, 730, 741, 742, 743, 744, 745, 751/)
    integer :: ctd(nctd)=(/780, 781, &
         & 810, 830, &
         & 901, 902, 903, 904, 905, 906, 907, 908, 909, 910/)
    integer :: pf(npf)=(/831, 834, 835, 836, 837, 838, 839, 840, &
              & 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, &
              & 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, &
              & 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, &
              & 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, &
              & 881, 882, 884, 886/)
    integer :: buoy(nb)=(/820/)
    integer :: animal(na)=(/995,996/)
    integer :: other(no)=(/760, 825, 999, 1023/)

    !---BTP$
    integer :: db(ndb)=(/0, 1, 2, 3, 4, 5, 6, 38, 39/)
    integer :: fb(nfb)=(/8, 9, 10, 11, 12, 13, 14, 15, &
         & 26, 27, 28, 29, 30/)
    integer :: mb(nmb)=(/16, 17, 18, 19, 20,  &
         & 21, 22, 23, 24, 25, &
         & 31, 32, 33, 34, 35, 36, 37/)
    
    !---IN
    integer,intent(in) :: iyr,imon
    character(100),intent(in) :: filename

    !---OUT
    integer,intent(out) :: ins
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    fullfilename=trim(gtspp_dir)//"/"//yyyy//mm//"/"//trim(filename)
    
    !---Read pcode and cparm
    status=nf90_open(trim(fullfilename),nf90_nowrite,ncid)

    status=nf90_inq_dimid(ncid,"num_surf",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = ns)

    allocate(cparm(ns),pcode(ns))
    
    status=nf90_inq_varid(ncid,"surfacecodes_pcode",varid)
    status=nf90_get_var(ncid,varid,pcode)

    status=nf90_inq_varid(ncid,"surfacecodes_cparm",varid)
    status=nf90_get_var(ncid,varid,cparm)
    
    status=nf90_close(ncid)
    
    !---Detect observation
    !Initial
    ins=0

    !Phase1: BTP$
    do i=1,ns
       
       if(trim(pcode(i)) == "BTP$")then
          
          read(cparm(i),*) itype

          !Drifter buoy
          do j=1,ndb
             if(itype == db(j))then
                ins=-6
                exit
             end if             
          end do

          !Floating buoy
          do j=1,nfb
             if(itype == fb(j))then
                ins=-5
                exit
             end if             
          end do

          
          !Mooring buoy
          do j=1,nmb
             if(itype == mb(j))then
                ins=-7
                exit
             end if             
          end do
       end if
       
    end do !i

    if(ins /= 0)then
       deallocate(cparm,pcode)
       return
    end if

    !Phase2: BUOY
    do i=1,ns
       
       if(trim(pcode(i)) == "BUOY")THEN

          read(cparm(i),*) itype
          
          if(itype == 0)then !Drifting buoy
             ins=-6
          else if(itype == 1)then !Mooring buoy
             ins=-7
          else if(itype == 2)then !Floating buoy
             ins=-5
          end if             
       end if
    end do
    
    if(ins /= 0)then
       deallocate(cparm,pcode)
       return
    end if
    
    !Phase3: PEQ$
    do i=1,ns
       
       if(trim(pcode(i)) == "PEQ$")then

          read(cparm(i),*) itype

          !XBT
          do j=1,nxbt
             if(itype == xbt(j))then
                ins=-1
                exit
             end if            
          end do

          !XCTD
          do j=1,nxctd
             if(itype == xctd(j))then
                ins=-2
                exit
             end if
          end do

          !CTD
          do j=1,nctd
             if(itype == ctd(j))then
                ins=-3
                exit
             end if
          end do

          !Profiling float
          do j=1,npf
             if(itype == pf(j))then
                ins=-5
                exit
             end if
          end do

          !Buoy
          do j=1,nb
             if(itype == buoy(j))then
                ins=-7
                exit
             end if
          end do

          !Animal
          do j=1,na
             if(itype == animal(j))then
                ins=-8
                exit
             end if
          end do
          
          !Others
          do j=1,no
             if(itype == other(j))then
                ins=0
                exit
             end if
          end do
          
       end if
       
    end do !i

    if(ins /= 0)then
       deallocate(cparm,pcode)
       return
    end if

    !Phase4
    do i=1,ns
       if(trim(pcode(i)) == "DPC$" .or. trim(pcode(i)) == "LTN#" .or. trim(pcode(i)) == "MFD#" &
            & .or. trim(pcode(i)) == "PFR$" .or. trim(pcode(i)) == "PGSS" .or. trim(pcode(i)) == "PRT$" &
            & .or. trim(pcode(i)) == "RCT$" .or. trim(pcode(i)) == "RCT#" .or. trim(pcode(i)) == "XDIN" &
            & .or. trim(pcode(i)) == "XDMT" .or. trim(pcode(i)) == "XTRS")then
          ins=-1 !XBT
          exit
       else if(trim(pcode(i)) == "SCT$")then
          ins=-3 !CTD
          exit
       end if
    end do
       
    deallocate(cparm,pcode)
    
  end subroutine read_gtspp_obs_type
  
  !-------------------------------

  subroutine deallocate_gtspp(depth,t,s)

    implicit none
    real(kind = 8),allocatable :: depth(:),t(:),s(:)

    deallocate(depth,t,s)

  end subroutine deallocate_gtspp

  !-------------------------------

  subroutine deallocate_gtspp_filename(long,lati,filename)

    implicit none

    real(kind = 8),allocatable :: long(:),lati(:)
    character(100),allocatable :: filename(:)

    deallocate(long,lati)
    deallocate(filename)

  end subroutine deallocate_gtspp_filename

end module mod_read_gtspp
