module mod_read_en4

  character(100),parameter :: en4_dir="/data/R/R2402/DATA/EN4"

contains

  !-------------------------------------------------------
  ! Read EN4 dataset |
  !-------------------------------------------------------
  !
  ! Website: https://www.metoffice.gov.uk/hadobs/en4/download-en4-2-2.html
  ! Provider: Met Office
  ! *2-3 months delay from realtime
  !
  !-------------------------------------------------------
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
  !-------------------------------------------------------
  !
  ! Created by S.Ohishi 2025.07
  ! Modified by S.Ohishi 2026.10: Read QC flag
  ! 
  !--------------------------------------------------------
  
  subroutine read_en4(iyr,imon,iday,np,km,long,lati,depth,t,s,ins)

    use mod_rmiss
    use mod_check_netcdf
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

    real(kind = 8),parameter :: dmiss=99999.d0
    
    !---Common
    integer status,access
    integer ncid,dimid,varid
    integer ip,j,k
    integer itype
    integer ijul,sjul
    integer jyr,jmon

    real(kind = 8),allocatable :: jul(:)
    
    real(kind = 4),allocatable :: rdep(:,:)
    real(kind = 4),allocatable :: rt(:,:),rs(:,:)
    
    character(200) filename
    character(4) yyyy
    character(2) mm
    character(len=:),allocatable :: cbuf

    character(4),allocatable,dimension(:) :: cins    
    character(1),allocatable,dimension(:) :: pos_qc
    character(1),allocatable,dimension(:,:) :: t_qc,s_qc
    
    !---WMO instrument code (Table 1770)
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
    
    !---IN
    integer,intent(in) :: iyr,imon

    !---OUT
    integer,intent(out) :: np !Number of profile
    integer,intent(out) :: km !Number of level

    integer,allocatable,intent(out) :: iday(:)
    integer,allocatable,intent(out) :: ins(:)
    
    real(kind = 8),allocatable,intent(out) :: long(:),lati(:),depth(:,:)
    real(kind = 8),allocatable,intent(out) :: t(:,:),s(:,:)
    
    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    filename=trim(en4_dir)//"/EN.4.2.2.f.profiles.g10."//yyyy//mm//".nc "

    !---Check Access
    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,'(a)') "Read "//trim(filename)
    else
       write(*,'(a)') "Skip to Read ==> Reason: Not found "//trim(filename)
       np=0
       km=0
       return
    end if

    !---Read np and km
    status=nf90_open(trim(filename),nf90_nowrite,ncid)
    call check_netcdf(status)

    !Number of Profile
    status=nf90_inq_dimid(ncid,"N_PROF",dimid)
    call check_netcdf(status)
    status=nf90_inquire_dimension(ncid,dimid,len = np)
    call check_netcdf(status)

    !Number of Layer
    status=nf90_inq_dimid(ncid,"N_LEVELS",dimid)
    call check_netcdf(status)
    status=nf90_inquire_dimension(ncid,dimid,len = km)
    call check_netcdf(status)

    !---Allocate
    allocate(jul(np))
    allocate(pos_qc(np),t_qc(km,np),s_qc(km,np))
    allocate(rdep(km,np))
    allocate(rt(km,np),rs(km,np))
    allocate(cins(np))
    
    allocate(iday(np))
    allocate(long(np),lati(np),depth(km,np))
    allocate(t(km,np),s(km,np))
    allocate(ins(np))
    
    !---Read data
    !JULD since 1950-01-01
    status=nf90_inq_varid(ncid,"JULD",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,jul)
    call check_netcdf(status)
    
    !LONGITUDE
    status=nf90_inq_varid(ncid,"LONGITUDE",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,long)
    call check_netcdf(status)
    
    !LATITUDE
    status=nf90_inq_varid(ncid,"LATITUDE",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,lati)
    call check_netcdf(status)
    
    !DEPH_CORRECTED (plus) => Depth
    status=nf90_inq_varid(ncid,"DEPH_CORRECTED",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,rdep)
    call check_netcdf(status)
    
    !POSITION_QC (4: reject)
    allocate(character(len=np) :: cbuf)
    status=nf90_inq_varid(ncid,"POSITION_QC",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,cbuf)
    call check_netcdf(status)
    do ip=1,np
       pos_qc(ip)=cbuf(ip:ip)
    end do
    deallocate(cbuf)
    
    !POTM_CORRECTED => Potemtial temperature
    status=nf90_inq_varid(ncid,"POTM_CORRECTED",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,rt)
    call check_netcdf(status)
    
    !POTM_CORRECTED_QC (4: reject)
    allocate(character(len=km*np) :: cbuf)
    status=nf90_inq_varid(ncid,"POTM_CORRECTED_QC",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,cbuf)
    call check_netcdf(status)
    do ip=1,np
       do k=1,km
          t_qc(k,ip)=cbuf((ip-1)*km+k:(ip-1)*km+k)
       end do
    end do
        
    !PSAL_CORRECTED => Practical salinity
    status=nf90_inq_varid(ncid,"PSAL_CORRECTED",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,rs)
    call check_netcdf(status)

    !PSAL_CORRECTED_QC (4:reject)
    status=nf90_inq_varid(ncid,"PSAL_CORRECTED_QC",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,cbuf)
    call check_netcdf(status)
    do ip=1,np
       do k=1,km
          s_qc(k,ip)=cbuf((ip-1)*km+k:(ip-1)*km+k)
       end do
    end do
    deallocate(cbuf)
    
    !WMO_INST_TYPE
    status=nf90_inq_varid(ncid,"WMO_INST_TYPE",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,cins)
    call check_netcdf(status)
    
    !INST_REFERENCE
    !status=nf90_inq_varid(ncid,"INST_REFERENCE",varid)
    !status=nf90_get_var(ncid,varid,cins)
    
    status=nf90_close(ncid)
    call check_netcdf(status)

    !---sjul
    call ymd_julian(1950,1,1,sjul)

    !---Post Process
    do ip=1,np

       !Date
       call julian_ymd(sjul+int(jul(ip)),jyr,jmon,iday(ip))       
       if(iyr /= jyr .or. imon /= jmon)then
          write(*,*) "***Error: Inconsistent date"
          write(*,*) iyr,imon
          write(*,*) jyr,jmon,iday(ip)
          stop
       end if
       
       !Longitude & Latitude
       if(pos_qc(ip) == "4" .or. long(ip) == dmiss .or. lati(ip) == dmiss)then
          long(ip)=rmiss
          lati(ip)=rmiss
          depth(:,ip)=rmiss
          t(:,ip)=rmiss
          s(:,ip)=rmiss
          ins(ip)=nint(rmiss)
          cycle
       end if

       !Longitude
       if(long(ip) < 0.d0)then
          long(ip)=long(ip)+360.d0
       end if

       !Instrument
       ins(ip)=0
       read(cins(ip),*,iostat=status) itype
       if(status /= 0) itype=-999
          
       !XBT
       do j=1,nxbt
          if(itype == xbt(j))then
             ins(ip)=-1
          end if
       end do

       !XCTD
       do j=1,nxctd
          if(itype == xctd(j))then
             ins(ip)=-2
          end if
       end do
       
       !CTD
       do j=1,nctd
          if(itype == ctd(j))then
             ins(ip)=-3
          end if
       end do
          
       !Profiling float
       do j=1,npf
          if(itype == pf(j))then
             ins(ip)=-5
          end if
       end do
       
       !Buoy
       do j=1,nb
          if(itype == buoy(j))then
             ins(ip)=-7
          end if
       end do

       !Animal
       do j=1,na
          if(itype == animal(j))then
             ins(ip)=-8
          end if
       end do
       
       !Others
       do j=1,no
          if(itype == other(j))then
             ins(ip)=0
          end if
       end do

       do k=1,km

          !Depth
          if(rdep(k,ip) == real(dmiss))then
             depth(k,ip)=rmiss
             t(k,ip)=rmiss
             s(k,ip)=rmiss
             cycle
          else
             depth(k,ip)=dble(rdep(k,ip))
          end if

          !Temperature
          if(rt(k,ip) == real(dmiss) .or. t_qc(k,ip) == "0" .or. t_qc(k,ip) == "4")then
             t(k,ip)=rmiss
          else
             t(k,ip)=dble(rt(k,ip))
          end if

          !Salinity
          if(rs(k,ip) == real(dmiss) .or. s_qc(k,ip) == "0" .or. s_qc(k,ip) == "4")then
             s(k,ip)=rmiss
          else
             s(k,ip)=dble(rs(k,ip))
          end if
          
       end do !k
    end do    !ip
    
    !---Deallocate
    deallocate(jul)
    deallocate(pos_qc,t_qc,s_qc)
    deallocate(rdep,rt,rs)
    deallocate(cins)
    
  end subroutine read_en4

  !--------------------------------------------------------

  subroutine deallocate_en4(iday,long,lati,depth,t,s,ins)

    implicit none
    
    integer,allocatable,intent(inout) :: iday(:),ins(:)
    
    real(kind = 8),allocatable,intent(inout) :: long(:),lati(:),depth(:,:)
    real(kind = 8),allocatable,intent(inout) :: t(:,:),s(:,:)

    if(allocated(iday))  deallocate(iday)
    if(allocated(ins))   deallocate(ins)
    if(allocated(long))  deallocate(long)
    if(allocated(lati))  deallocate(lati)
    if(allocated(depth)) deallocate(depth)
    if(allocated(t))     deallocate(t)
    if(allocated(s))     deallocate(s)
    
  end subroutine deallocate_en4

  !---------------------------------------------------------
  
end module mod_read_en4
