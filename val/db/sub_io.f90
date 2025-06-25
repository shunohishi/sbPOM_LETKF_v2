subroutine write_bin(im_bin,jm_bin,ndat_a,dx_bin,dy_bin,lon_bin,lat_bin, &
     & unum_bin,ubias_bin,urmsd_bin, &
     & vnum_bin,vbias_bin,vrmsd_bin, &
     & tnum_bin,tbias_bin,trmsd_bin)

  implicit none

  !---Common
  integer i_bin,j_bin

  character(20) format
  
  !---IN
  integer,intent(in) :: im_bin,jm_bin,ndat_a
  real(kind = 8),intent(in) :: dx_bin,dy_bin
    
  integer,intent(in) :: unum_bin(im_bin,jm_bin,ndat_a)
  integer,intent(in) :: vnum_bin(im_bin,jm_bin,ndat_a)
  integer,intent(in) :: tnum_bin(im_bin,jm_bin,ndat_a)

  real(kind = 8),intent(in) :: lon_bin(im_bin),lat_bin(jm_bin)  
  real(kind = 8),intent(in) :: ubias_bin(im_bin,jm_bin,ndat_a),urmsd_bin(im_bin,jm_bin,ndat_a)
  real(kind = 8),intent(in) :: vbias_bin(im_bin,jm_bin,ndat_a),vrmsd_bin(im_bin,jm_bin,ndat_a)
  real(kind = 8),intent(in) :: tbias_bin(im_bin,jm_bin,ndat_a),trmsd_bin(im_bin,jm_bin,ndat_a)

  write(format,'(a,I0,a,I0,a)') '(2f12.5,',ndat_a,'i6,',2*ndat_a,'f12.5)'
  
  open(11,file="dat/u_bin.dat",status="replace")
  open(12,file="dat/v_bin.dat",status="replace")
  open(13,file="dat/t_bin.dat",status="replace")
  do j_bin=1,jm_bin-1
     do i_bin=1,im_bin-1
        write(11,trim(format)) &
             & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
             & unum_bin(i_bin,j_bin,:),ubias_bin(i_bin,j_bin,:),urmsd_bin(i_bin,j_bin,:)
        write(12,trim(format)) &
             & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
             & vnum_bin(i_bin,j_bin,:),vbias_bin(i_bin,j_bin,:),vrmsd_bin(i_bin,j_bin,:)
        write(13,trim(format)) &
             & lon_bin(i_bin)+0.5d0*dx_bin,lat_bin(j_bin)+0.5d0*dy_bin, &
             & tnum_bin(i_bin,j_bin,:),tbias_bin(i_bin,j_bin,:),trmsd_bin(i_bin,j_bin,:)
     end do !i_bin
  end do    !j_bin
  close(11)
  close(12)
  close(13)
  
end subroutine write_bin

!---------------------------------------------------

subroutine write_ave(syr,eyr,ndat_a, &
     & unum_mave,ubias_mave,urmsd_mave, &
     & vnum_mave,vbias_mave,vrmsd_mave, &
     & tnum_mave,tbias_mave,trmsd_mave, &
     & unum_yave,ubias_yave,urmsd_yave, &
     & vnum_yave,vbias_yave,vrmsd_yave, &
     & tnum_yave,tbias_yave,trmsd_yave, &
     & unum_ave,ubias_ave,urmsd_ave, &
     & vnum_ave,vbias_ave,vrmsd_ave, &
     & tnum_ave,tbias_ave,trmsd_ave)

  implicit none

  !---Common
  integer iyr,imon

  character(20) format
  character(10) yyyymmdd
  
  !---IN  
  integer,intent(in) :: syr,eyr,ndat_a

  !Monthly
  integer,intent(in) :: unum_mave(ndat_a,12,syr:eyr)
  integer,intent(in) :: vnum_mave(ndat_a,12,syr:eyr)
  integer,intent(in) :: tnum_mave(ndat_a,12,syr:eyr)

  !Yearly
  integer,intent(in) :: unum_yave(ndat_a,syr:eyr)
  integer,intent(in) :: vnum_yave(ndat_a,syr:eyr)
  integer,intent(in) :: tnum_yave(ndat_a,syr:eyr)

  !ALL
  integer,intent(in) :: unum_ave(ndat_a)
  integer,intent(in) :: vnum_ave(ndat_a)
  integer,intent(in) :: tnum_ave(ndat_a)
  
  !Monthly
  real(kind = 8),intent(in) :: ubias_mave(ndat_a,12,syr:eyr),vbias_mave(ndat_a,12,syr:eyr),tbias_mave(ndat_a,12,syr:eyr)
  real(kind = 8),intent(in) :: urmsd_mave(ndat_a,12,syr:eyr),vrmsd_mave(ndat_a,12,syr:eyr),trmsd_mave(ndat_a,12,syr:eyr)

  !Yearly
  real(kind = 8),intent(in) :: ubias_yave(ndat_a,syr:eyr),vbias_yave(ndat_a,syr:eyr),tbias_yave(ndat_a,syr:eyr)
  real(kind = 8),intent(in) :: urmsd_yave(ndat_a,syr:eyr),vrmsd_yave(ndat_a,syr:eyr),trmsd_yave(ndat_a,syr:eyr)

  !ALL
  real(kind = 8),intent(in) :: ubias_ave(ndat_a),vbias_ave(ndat_a),tbias_ave(ndat_a)
  real(kind = 8),intent(in) :: urmsd_ave(ndat_a),vrmsd_ave(ndat_a),trmsd_ave(ndat_a)
  
  !---Monthly
  write(format,'(a,I0,a,I0,a)') "(a,",ndat_a,"i6,",2*ndat_a,"f12.5)"
  
  open(11,file="dat/u_mave.dat",status="replace")
  open(12,file="dat/v_mave.dat",status="replace")
  open(13,file="dat/t_mave.dat",status="replace")
  do iyr=syr,eyr
     do imon=1,12
        write(yyyymmdd,'(i4.4,a,i2.2,a)') iyr,"-",imon,"-15"
        write(11,trim(format)) &
             & yyyymmdd,unum_mave(:,imon,iyr),ubias_mave(:,imon,iyr),urmsd_mave(:,imon,iyr)
        write(12,trim(format)) &
             & yyyymmdd,vnum_mave(:,imon,iyr),vbias_mave(:,imon,iyr),vrmsd_mave(:,imon,iyr)
        write(13,trim(format)) &
             & yyyymmdd,tnum_mave(:,imon,iyr),tbias_mave(:,imon,iyr),trmsd_mave(:,imon,iyr)
     end do
  end do
  close(11)
  close(12)
  close(13)

  !---Yearly
  write(format,'(a,I0,a,I0,a)') "(i6,",ndat_a,"i6,",2*ndat_a,"f12.5)"
  
  open(11,file="dat/u_yave.dat",status="replace")
  open(12,file="dat/v_yave.dat",status="replace")
  open(13,file="dat/t_yave.dat",status="replace")
  do iyr=syr,eyr
     write(11,trim(format)) &
          & iyr,unum_yave(:,iyr),ubias_yave(:,iyr),urmsd_yave(:,iyr)
     write(12,trim(format)) &
          & iyr,vnum_yave(:,iyr),vbias_yave(:,iyr),vrmsd_yave(:,iyr)
     write(13,trim(format)) &
          & iyr,tnum_yave(:,iyr),tbias_yave(:,iyr),trmsd_yave(:,iyr)
  end do
  close(11)
  close(12)
  close(13)

  !---ALL
  write(format,'(a,I0,a,I0,a)') "(",ndat_a,"i6,",2*ndat_a,"f12.5)"
  
  open(11,file="dat/u_ave.dat",status="replace")
  open(12,file="dat/v_ave.dat",status="replace")
  open(13,file="dat/t_ave.dat",status="replace")
  write(11,trim(format)) &
       & unum_ave(:),ubias_ave(:),urmsd_ave(:)
  write(12,trim(format)) &
       & vnum_ave(:),vbias_ave(:),vrmsd_ave(:)
  write(13,trim(format)) &
       & tnum_ave(:),tbias_ave(:),trmsd_ave(:)
  close(11)
  close(12)
  close(13)
  
end subroutine write_ave
