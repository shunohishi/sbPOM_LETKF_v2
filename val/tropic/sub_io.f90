subroutine write_month(idata,iyr,imon,im,jm,km,lon,lat,depth,mmean,manom)

  implicit none

  !---Common
  integer i,j,k

  character(4) yyyy
  character(2) mm,dd
  
  !---IN
  integer,intent(in) :: idata
  integer,intent(in) :: iyr,imon
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(in) :: mmean(im,jm,km),manom(im,jm,km)

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') idata
  
  open(1,file="dat/mmean"//yyyy//mm//"-d"//dd//".dat",status="replace")
  do k=1,km
     do j=1,jm
        do i=1,im
           write(1,'(5f12.5)') lon(i),lat(j),depth(i,j,k),mmean(i,j,k),manom(i,j,k) 
        end do
     end do
  end do
  close(1)
  
end subroutine write_month

!-----------------------------------------------------------------------------------

subroutine write_mclim(idata,imon,im,jm,km,lon,lat,depth,mclim)

  implicit none

  !---Common
  integer i,j,k
  
  character(2) mm,dd
  
  !---IN
  integer,intent(in) :: idata
  integer,intent(in) :: imon
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(in) :: mclim(im,jm,km)

  write(dd,'(i2.2)') idata
  write(mm,'(i2.2)') imon
  
  open(1,file="dat/mclim"//mm//"-d"//dd//".dat",status="replace")
  do k=1,km
     do j=1,jm
        do i=1,im
           write(1,'(4f12.5)') lon(i),lat(j),depth(i,j,k),mclim(i,j,k) 
        end do
     end do
  end do
  close(1)
  
end subroutine write_mclim

!--------------------------------------------------------------------------------------

subroutine write_index(nino,nino_norm,iod,iod_norm,anino,anino_norm)

  use setting
  implicit none

  !---Common
  integer iyr,imon

  character(10) yyyymmdd
  character(4) yyyy
  character(2) mm
  
  !---IN
  real(kind = 8),intent(in) :: nino(syr:eyr,12,ndata),nino_norm(syr:eyr,12,ndata)
  real(kind = 8),intent(in) :: iod(syr:eyr,12,ndata),iod_norm(syr:eyr,12,ndata)
  real(kind = 8),intent(in) :: anino(syr:eyr,12,ndata),anino_norm(syr:eyr,12,ndata)

  open(1,file="dat/nino34.dat",status="replace")
  open(2,file="dat/iod.dat",status="replace")
  open(3,file="dat/anino.dat",status="replace")
  do iyr=syr,eyr
     do imon=1,12

        write(yyyy,'(i4.4)') iyr
        write(mm,'(i2.2)') imon
        yyyymmdd=yyyy//"-"//mm//"-"//"15"
        
        write(1,'(a,8f12.5)') yyyymmdd,nino(iyr,imon,:),nino_norm(iyr,imon,:)
        write(2,'(a,8f12.5)') yyyymmdd,iod(iyr,imon,:),iod_norm(iyr,imon,:)
        write(3,'(a,8f12.5)') yyyymmdd,anino(iyr,imon,:),anino_norm(iyr,imon,:)
        
     end do
  end do
  close(1)
  close(2)
  close(3)
  
end subroutine write_index

