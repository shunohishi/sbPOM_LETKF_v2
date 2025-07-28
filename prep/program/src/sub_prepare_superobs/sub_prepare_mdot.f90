!----------------------------------------------------------------------
! Prepare MDOT |
!----------------------------------------------------------------------
! <Model SSH>
! 1. Calculate sshmean from model outputs averaging in syr-eyr
!
!----------------------------------------------------------------------
! Created by S.Ohishi 2018.09
! Modified by S.Ohishi 2019.10 for ensemble pom
!----------------------------------------------------------------------

subroutine prepare_mdot(lon,lat,mdot)

  use mod_rmiss
  use setting, only : syr => syr_ssh, eyr => eyr_ssh
  use mod_gridinfo
  implicit none
  
  !Common
  integer iyr,imon,iday,nday

  real(kind = 8) ssh(im,jm)
  real(kind = 8) sshmean(im,jm),sshpass(im,jm),sshmiss(im,jm)
  
  !IN
  real(kind = 8),intent(in) :: lon(im),lat(jm)

  !OUT
  real(kind = 8),intent(out) :: mdot(im,jm)
  

  write(*,*) "Calculate model ssh mean"
  do iyr=syr,eyr
     do imon=1,12

        write(*,'(i4.4,i2.2)') iyr,imon
        call number_of_day(iyr,imon,nday)

        do iday=1,nday
           call read_ssh(iyr,imon,iday,ssh)
           call cal_clim(iyr,syr,eyr,imon,iday, &
                & im,jm,1,ssh,sshmean,sshpass,sshmiss)
        end do !iday

     end do !imon
  end do !iyr

  mdot(:,:)=sshmean(:,:)

end subroutine prepare_mdot

!----------------------------------------------------------------------
! Calculate Mean |
!----------------------------------------------------------------------

subroutine cal_clim(iyr,syr,eyr,imon,iday,im,jm,km,dat,mean,pass,miss)

  use mod_rmiss
  implicit none

  integer iyr,syr,eyr,imon,iday
  integer i,j,k,im,jm,km

  real(kind = 8) dat(im,jm,km)
  real(kind = 8) mean(im,jm,km),pass(im,jm,km),miss(im,jm,km)

  if(iyr == syr .and. imon == 1 .and. iday == 1)then
     mean(:,:,:)=0.d0
     pass(:,:,:)=0.d0
     miss(:,:,:)=0.d0
  end if

  do k=1,km
     do j=1,jm
        do i=1,im
           if(dat(i,j,k) == rmiss)then
              miss(i,j,k)=miss(i,j,k)+1.d0
           else
              mean(i,j,k)=mean(i,j,k)+dat(i,j,k)
              pass(i,j,k)=pass(i,j,k)+1.d0
           end if
        end do
     end do
  end do

  if(iyr == eyr .and. imon == 12 .and. iday == 31)then
     do k=1,km
        do j=1,jm
           do i=1,im
              if(pass(i,j,k) == 0.d0)then
                 mean(i,j,k)=rmiss
              else
                 mean(i,j,k)=mean(i,j,k)/pass(i,j,k)
              end if
           end do
        end do
     end do
  end if

end subroutine cal_clim
