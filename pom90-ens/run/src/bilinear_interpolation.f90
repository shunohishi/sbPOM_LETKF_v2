!---------------------------------------------------------------
! Calculate ID (longitude)|
!--------------------------
!
! 1  ----i2----2
! i1 ----i2----i1+1
! im1----i2----1
! 
! * id denotes i1, when dat2 is at i2
!
!---------------------------------------------------------------

subroutine set_idlon(im1,lon1,im2,lon2,id)

  !$use omp_lib
  use common_pom_var, only: r_size
  implicit none

  !---Common
  integer i1,i2
  integer n1,n2

  !---IN
  integer,intent(in) :: im1,im2
  real(kind = r_size),intent(in) :: lon1(im1),lon2(im2)

  !---OUT
  integer,intent(out) :: id(im2)

  id(:)=0

  !$omp parallel
  !$omp do private(i1,i2,n1,n2)
  do i2=1,im2

     n1=floor(lon1(1)/360.d0)
     n2=floor(lon2(i2)/360.d0)

     !ex. -0.5 < lon2 < 0.5
     if(lon1(im1)-(n1+1)*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0  <= lon1(1)-n1*360.d0)then
        id(i2)=im1
        !ex. 359.5 < lon2 < 360.5
     else if(lon1(im1)-n1*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0  <= lon1(1)-(n1-1)*360.d0)then
        id(i2)=im1
     else

        do i1=1,im1-1
           n1=floor(lon1(i1)/360.d0)
           if(lon1(i1)-n1*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0 <= lon1(i1+1)-n1*360.d0)then
              id(i2)=i1
           end if
        end do !i1

     end if
  end do !i2
  !$omp end do
  !$omp end parallel

  do i2=1,im2
     if(id(i2) == 0)then
        write(*,*) "***Error: idx is missed at ",i2
        stop
     end if
  end do

end subroutine set_idlon

!------------------------

subroutine set_idlat(jm1,lat1,jm2,lat2,id)

  !$use omp_lib
  use common_pom_var, only: r_size
  implicit none

  !---Common
  integer j1,j2

  !---IN
  integer,intent(in) :: jm1,jm2
  real(kind = r_size),intent(in) :: lat1(jm1),lat2(jm2)

  !---OUT
  integer,intent(out) :: id(jm2)

  id(:)=0

  !$omp parallel
  !$omp do private(j1,j2)  
  do j2=1,jm2
     do j1=1,jm1-1

        if(lat1(j1) <= lat2(j2) .and. lat2(j2) <= lat1(j1+1))then
           id(j2)=j1
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel

  do j2=1,jm2
     if(id(j2) == 0)then
        write(*,*) "***Error: idy is missed at ",j2
        stop
     end if
  end do

end subroutine set_idlat

!------------------------------------------------------------
! Coast line |
!------------------------------------------------------------
!Coast line gird:1, other: 0
!------------------------------------------------------------

subroutine coast_line(im,jm,fsm,cline)

  use common_pom_var, only: r_size  
  !$use omp_lib  
  implicit none

  !---Common
  integer i,j,i1,i2,j1,j2
  integer pass(im,jm),miss(im,jm)

  !---IN
  integer,intent(in) :: im,jm
  real(kind = r_size),intent(in) :: fsm(im,jm)

  !---OUT
  integer,intent(out) :: cline(im,jm)

  cline(:,:)=0

  !$omp parallel
  !$omp do private(i1,j1,i2,j2)
  do j1=1,jm
     do i1=1,im

        pass(i1,j1)=0
        miss(i1,j1)=0

        do j2=j1-1,j1+1
           if(j2 < 1 .or. jm < j2)cycle
           do i2=i1-1,i1+1
              if(i2 < 1 .or. im < i2)cycle

              if(fsm(i2,j2) == REAL(0.d0,r_size))then
                 miss(i1,j1)=miss(i1,j1)+1
              else
                 pass(i1,j1)=pass(i1,j1)+1
              end if

           end do
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j)  
  do i=1,im
     do j=1,jm
        !Make cline (miss = 0.)
        if(fsm(i,j) == REAL(1.d0,r_size) .and. miss(i,j) /= 0)then
           cline(i,j)=1
        endif
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine coast_line
!----------------------------------------------------
! Distance between (x1,y1) and (x2,y2)
!----------------------------------------------------
!
! x1,x2: longitude (degree)
! y1,y2: latitude (degree)
!---------------------------------------------------

subroutine distance(x1,y1,x2,y2,r)

  use common_pom_var, only: r_size,pi,earth
  implicit none

  !---IN  
  real(kind = r_size),intent(in) :: x1,y1
  real(kind = r_size),intent(in) :: x2,y2

  !---OUT
  real(kind = r_size),intent(out) :: r

  if(x1 == x2 .and. y1 == y2)then
     r=0.d0
  else
     r=earth*acos( &
          & cos(pi*y1/180.d0)*cos(pi*y2/180.d0)*cos(pi*x1/180.d0)*cos(pi*x2/180.d0) &
          & +cos(pi*y1/180.d0)*cos(pi*y2/180.d0)*sin(pi*x1/180.d0)*sin(pi*x2/180.d0) &
          & +sin(pi*y1/180.d0)*sin(pi*y2/180.d0) &
          & )
  end if

end subroutine distance

!--------------------------------------------------------------
! Calculate CaMa-Flood ID at model coast line |
!--------------------------------------------------------------
!
!--- INPUT
! im,jm,lon,lat: model grid
! cline: model coast line
! 
! *_cama: CaMa-Flood grid
!
!--- OUTPUT
! idx,idy: ID
! dnum: duplication number
!
!--------------------------------------------------------------

subroutine set_id_cama(im,jm,lon,lat,cline, &
     & im_cama,jm_cama,lon_cama,lat_cama,cline_cama, &
     & idx,idy,dnum)

  use common_pom_var, only: r_size
  implicit none

  real(kind = r_size),parameter :: r_max=50.d3

  integer i1,i2,di,j1,j2,dj
  real(kind = r_size) r,r0 !distance between 2 point

  !---IN
  integer,intent(in) :: im,jm  
  integer,intent(in) :: cline(im,jm)

  real(kind = r_size),intent(in) :: lon(im),lat(jm)

  integer,intent(in) :: im_cama,jm_cama
  real(kind = r_size),intent(in) :: lon_cama(im_cama),lat_cama(jm_cama)
  integer,intent(in) :: cline_cama(im_cama,jm_cama)

  !---OUT
  integer,intent(out) :: idx(im,jm),idy(im,jm) !river position
  integer,intent(out) :: dnum(im,jm)           !duplication number

  idx(:,:)=0
  idy(:,:)=0
  dnum(:,:)=0

  !Grid id with minimum distance
  !$omp parallel
  !$omp do private(i1,j1,i2,j2,r0)
  do j1=1,jm
     do i1=1,im

        if(cline(i1,j1) == 0) cycle
        r0=r_max
        
        do j2=1,jm_cama
           do i2=1,im_cama


              if(cline_cama(i2,j2) == 0) cycle

              call distance(lon(i1),lat(j1),lon_cama(i2),lat_cama(j2),r)

              if(r < r0)then
                 r0 = r
                 idx(i1,j1)=i2
                 idy(i1,j1)=j2
              end if

           end do
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i1,j1,di,dj,i2,j2)  
  do j1=1,jm
     do i1=1,im

        if(cline(i1,j1) == 0) cycle

        do j2=1,jm
           do i2=1,im        
              if(idx(i1,j1) /= idx(i2,j2) .or. idy(i1,j1) /= idy(i2,j2))cycle
              dnum(i1,j1)=dnum(i1,j1)+1
           end do
        end do

     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine set_id_cama

!-----------------------------------------------------------------------
! Fill value in horizontal direction |
!----------------------------------------
!
! ncount: the number of times in fill value
!
!-----------------------------------------------------------------------

subroutine fillvalue_2d(ncount,im,jm,dat,fsm)

  use common_pom_var, only: r_size
  !$use omp_lib
  implicit none

  !---Common
  integer i,j
  integer di,dj
  integer count,pass

  real(kind = r_size) tmp(im,jm)

  !---IN
  integer,intent(in) :: ncount
  integer,intent(in) :: im,jm

  real(kind = r_size),intent(inout) :: dat(im,jm),fsm(im,jm)

  count=1

  do while(count <= ncount)

     !     write(*,*) "The number of times filling value:",count,"/",ncount 

     !$omp parallel
     !$omp do private(i,j,di,dj,pass)
     do j=1,jm
        do i=1,im

           if(dat(i,j) == 0.d0 .and. fsm(i,j) == REAL(0.d0,r_size))then

              tmp(i,j)=0.d0
              pass=0

              do dj=-1,1
                 do di=-1,1
                    if(i+di < 1 .or. im < i+di)cycle
                    if(j+dj < 1 .or. jm < j+dj)cycle
                    if(fsm(i+di,j+dj) == REAL(1.d0,r_size))then
                       tmp(i,j)=tmp(i,j)+dat(i+di,j+dj)
                       pass=pass+1
                    end if
                 end do
              end do

              if(pass == 0)then
                 tmp(i,j)=REAL(0.d0,r_size)
              else
                 tmp(i,j)=tmp(i,j)/dble(pass)
              end if

           else
              tmp(i,j)=dat(i,j)
           end if
        end do
     end do
     !$omp end do
     !$omp end parallel

     dat(:,:)=tmp(:,:)
     count=count+1

  end do

end subroutine fillvalue_2d

!-----------------------------------------------------------
! Bilinear Interpolation |
!-----------------------------------------------------------
!
! f01(x0,y1) ___________ f11(x1,y1)
!  |                       |
! y|-----------fxy(x,y)    |
!  |            |          |
!  |            |          |
! f00(x0,y0) ___|_______ f10(x1,y0)
!               x
!-----------------------------------------------------------

subroutine bilinear_interpolation_2d &
     & (idlon,idlat,im1,jm1,lon1,lat1,dat1,fsm1,im2,jm2,lon2,lat2,dat2,fsm2)

  !$use omp_lib
  use common_pom_var, only: r_size  
  implicit none

  !---Common
  integer i10,i11,j10,j11
  integer i2,j2
  integer n1,n10,n11,n2

  real(kind = r_size) x0,x1,y0,y1,x,y
  real(kind = r_size) f00,f10,f01,f11

  !---IN
  integer,intent(in) :: idlon(im2),idlat(jm2)
  integer,intent(in) :: im1,jm1
  integer,intent(in) :: im2,jm2

  real(kind = r_size),intent(in) :: lon1(im1),lat1(jm1),dat1(im1,jm1),fsm1(im1,jm1)
  real(kind = r_size),intent(in) :: lon2(im2),lat2(jm2),fsm2(im2,jm2)

  !---OUT
  real(kind = r_size),intent(out) :: dat2(im2,jm2)

  !$omp parallel
  !$omp do private(i2,j2,n2,i10,i11,n1,n10,n11,j10,j11,x0,x1,y0,y1,x,y,f00,f10,f01,f11)    
  do j2=1,jm2
     do i2=1,im2

        if(fsm2(i2,j2) == REAL(0.d0,r_size))then
           dat2(i2,j2)=0.d0
           cycle
        end if
        
        n2=floor(lon2(i2)/360.d0)

        if(idlon(i2) == im1)then

           i10=im1
           i11=1
           n1=floor(lon1(1)/360.d0)

           !ex. -0.5 < lon2 < 0.5
           if(lon1(im1)-(n1+1)*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0 <= lon1(1)-n1*360.d0)then
              n10=floor(lon1(im1)/360.d0)-1
              n11=floor(lon1(1)/360.d0)
              !ex.359.5 < lon2 < 360.5
           else if(lon1(im1)-n1*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0 <= lon1(1)-(n1-1)*360.d0)then
              n10=floor(lon1(im1)/360.d0)
              n11=floor(lon1(1)/360.d0)+1
           end if

        else

           i10=idlon(i2)
           i11=idlon(i2)+1
           n10=floor(lon1(i10)/360.d0)
           n11=floor(lon1(i11)/360.d0)

        end if

        j10=idlat(j2)
        j11=idlat(j2)+1

        x0=REAL(lon1(i10)-n10*360.d0,r_size)
        x1=REAL(lon1(i11)-n11*360.d0,r_size)
        y0=lat1(j10)
        y1=lat1(j11)
        x=REAL(lon2(i2)-n2*360.d0,r_size)
        y=lat2(j2)
        f00=dat1(i10,j10)
        f10=dat1(i11,j10)
        f01=dat1(i10,j11)
        f11=dat1(i11,j11)
        
        dat2(i2,j2)=(y1-y)/(y1-y0)*((x1-x)/(x1-x0)*f00+(x-x0)/(x1-x0)*f10) &
             & +(y-y0)/(y1-y0)*((x1-x)/(x1-x0)*f01+(x-x0)/(x1-x0)*f11)
        
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine bilinear_interpolation_2d

!-----------------------------------------------------------------------------
! disrtribution River flux 
!-----------------------------------------------------------------------------

subroutine distribute_river(im_cama,jm_cama,dat_cama,im,jm,dat,idx,idy,dnum)

  use common_pom_var, only: r_size  
  implicit none

  !---Common
  integer i,j

  !---IN
  integer,intent(in) :: im_cama,jm_cama  
  integer,intent(in) :: im,jm
  integer,intent(in) :: idx(im,jm),idy(im,jm),dnum(im,jm)

  real(kind = r_size),intent(in) :: dat_cama(im_cama,jm_cama)

  !---OUT
  real(kind = r_size),intent(out) :: dat(im,jm)

  dat(:,:)=0.d0

  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im

        if(idx(i,j) == 0 .or. idy(i,j) == 0 .or. dnum(i,j) == 0)cycle

        dat(i,j)=dat_cama(idx(i,j),idy(i,j))/REAL(dnum(i,j),r_size)

        if(dat(i,j) < 0.d0) dat(i,j)=0.d0

     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine distribute_river
