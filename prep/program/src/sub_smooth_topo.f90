!--------------------------------------------------------------------
! Smoothing topography |
!--------------------------------------------------------------------
! Mellor et al. (1994) JAOT
!
! To reduce errors caused by steep slope,
! repeat smoothing procedure whenever (H2-H1)/(H2+H1) < dh is satisfied
!
! Created by S.Ohishi 2018.07
! Modified by S.Ohishi 2023.03
! Modified by S.Ohishi 2023.12
!
!--------------------------------------------------------------------

subroutine smooth_topo(im,jm,lon,lat,topo,topo_dif,ratio,fsm)

  use setting, only: dh,escale,delta
  implicit none
  
  integer,parameter :: nmax=200
  
  integer n

  integer count(0:nmax)
  
  real(kind = 8) max(0:nmax)
  real(kind = 8) sigma(im,jm) !e-folding scale in Gaussian filter
  real(kind = 8) tmp(im,jm)

  !---IN
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: lon(im),lat(jm),fsm(im,jm)

  !---IN/OUT
  real(kind = 8),intent(inout) :: topo(im,jm)

  !---OUT
  real(kind = 8),intent(out) :: topo_dif(im,jm),ratio(im,jm)
  
  n=0
  count(:)=0
  max(:)=1.d0
  tmp(:,:)=topo(:,:)

  !---Gaussian filter
  if(escale /= 0.d0)then
     sigma(:,:)=escale
     call gauss_filter(im,jm,lon,lat,fsm,sigma,topo)
  end if
  
  !---Smoothing to satisfy dh condition
  if(dh /= 0.d0)then
     do while(dh < max(n))

        n=n+1

        write(*,'(i6,a)') n," Times"
        write(*,'(a)') "-----Before Smoothing Filter-----"
        call calculate_ratio(im,jm,topo,fsm,count(n),max(n),ratio)  
        write(*,'(a)') "--------------------------------"

        call smooth_filter(im,jm,max(n),topo,fsm)

        write(*,'(a)') "-----After Smoothing Filter-----"
        call calculate_ratio(im,jm,topo,fsm,count(n),max(n),ratio)
        write(*,'(a)') "-------------------------------"

        if(max(n) <= dh+delta .or. n == nmax) exit

     end do !do while
  end if
  
  !---Summary Information
  write(*,'(a)') "---Summary information on smoothing topography---"
  call calculate_ratio(im,jm,topo,fsm,count(0),max(0),ratio)
  write(*,'(a)') "------------------------------------------------"

  topo_dif(:,:)=tmp(:,:)-topo(:,:)
  
end subroutine smooth_topo

!------------------------------------------------------------------

subroutine calculate_ratio(im,jm,topo,fsm,count,max,ratio)

  use setting, only: lons,lone,dh,delta
  implicit none

  integer i,j,i0,i1,i2,j0,j1,j2,di,dj
  integer i_max,j_max
  integer im_tmp
  
  real(kind = 8) xratio,xratio1,xratio2
  real(kind = 8) yratio,yratio1,yratio2

  logical lcyc
  
  !---IN
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: topo(im,jm),fsm(im,jm)

  !---OUT
  integer,intent(out) :: count
  
  real(kind = 8),intent(out) :: ratio(im,jm)

  !---IN/OUT
  real(kind = 8),intent(inout) :: max

  count=0
  max=0.d0
  ratio(:,:)=0.d0

  if(lons == lone-360.d0)then
     lcyc=.true.
  else
     lcyc=.false.
  end if
  
  if(lcyc)then
     im_tmp=im-2
  else
     im_tmp=im
  end if
  
  !calculate count (the number of grid exceeding dh)
  do j=1,jm

     j0=j-1
     j1=j
     j2=j+1

     do i=1,im_tmp
           
        i0=i-1
        i1=i
        i2=i+1

        if(lcyc)then

           if(i0 < 1)then
              i0=im_tmp
           end if

           if(im_tmp < i2)then
              i2=1
           end if
           
        end if        

        if(fsm(i,j) == 0.d0)then
           ratio(i,j)=0.d0
           cycle
        end if
        
        if(i0 < 1 .or. im_tmp < i0 .or. i1 < 1 .or. im_tmp < i1)then
           xratio1=0.d0
        else if(fsm(i0,j) == 0.d0 .or. fsm(i1,j) == 0.d0)then
           xratio1=0.d0
        else
           xratio1=abs(topo(i1,j)-topo(i0,j))/abs(topo(i1,j)+topo(i0,j))
        end if

        if(i1 < 1 .or. im_tmp < i1 .or. i2 < 1 .or. im_tmp < i2)then
           xratio2=0.d0
        else if(fsm(i1,j) == 0.d0 .or. fsm(i2,j) == 0.d0)then
           xratio2=0.d0
        else
           xratio2=abs(topo(i2,j)-topo(i1,j))/abs(topo(i2,j)+topo(i1,j))
        end if

        if(xratio1 < xratio2)then
           xratio=xratio2
        else
           xratio=xratio1
        end if
        
        if(j0 < 1 .or. jm < j0 .or. j1 < 1 .or. jm < j1)then
           yratio1=0.d0
        else if(fsm(i,j0) == 0.d0 .or. fsm(i,j1) == 0.d0)then
           yratio1=0.d0
        else
           yratio1=abs(topo(i,j1)-topo(i,j0))/abs(topo(i,j1)+topo(i,j0))
        end if

        if(j1 < 1 .or. jm < j1 .or. j2 < 1 .or. jm < j2)then
           yratio2=0.d0
        else if(fsm(i,j1) == 0.d0 .or. fsm(i,j2) == 0.d0)then
           yratio2=0.d0
        else
           yratio2=abs(topo(i,j2)-topo(i,j1))/abs(topo(i,j2)+topo(i,j1))
        end if

        if(yratio1 < yratio2)then
           yratio=yratio2
        else
           yratio=yratio1
        end if
        
        if(yratio < xratio)then
           ratio(i,j)=xratio
        else
           ratio(i,j)=yratio
        end if
        
        if(ratio(i,j) > max)then
           i_max=i
           j_max=j
           max=ratio(i,j)
        end if
           
        if(ratio(i,j) > dh+delta)then
           count=count+1
        end if

     end do
  end do

  write(*,'(a,i10,a,i10)') "#Topography to be smoothed:",count,"/",im*jm
  write(*,'(a,f12.5)') "Max of (H2-H1)/(H2+H1):",max
    
end subroutine calculate_ratio

!------------------------------------------------------------------------------

subroutine smooth_filter(im,jm,dh_max,topo,fsm)

  use setting, only: lons,lone,dh,delta
  implicit none

  integer itime  
  integer i,si,ei,di,i1,i2
  integer j,sj,ej,dj,j1,j2
  integer im1,im2
  integer count

  real(kind = 8) dh_tmp,dh_cor
  real(kind = 8) ratio(im,jm)

  logical lcyc

  !---IN
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: fsm(im,jm)

  !---OUT
  real(kind = 8),intent(out) :: dh_max

  !---IN/OUT
  real(kind = 8),intent(inout) :: topo(im,jm)

  if(lons == lone-360.d0)then
     lcyc=.true.
  else
     lcyc=.false.
  end if

  if(lcyc)then
     im1=im-2
     im2=im-2
  else
     im1=im
     im2=im-1
  end if
     
  call calculate_ratio(im,jm,topo,fsm,count,dh_max,ratio)

  do itime=1,4
     
     if(itime == 1)then !Eastward
        si=1
        ei=im2
        di=1
        sj=1
        ej=jm
        dj=1
     else if(itime == 2)then !Westward
        si=im2
        ei=1
        di=-1
        sj=1
        ej=jm
        dj=1
     else if(itime == 3)then !Northward
        si=1
        ei=im1
        di=1
        sj=1
        ej=jm-1
        dj=1
     else if(itime == 4)then !Southward
        si=1
        ei=im1
        di=1
        sj=jm-1
        ej=1
        dj=-1
     end if
     
     do j=sj,ej,dj
        do i=si,ei,di

           if(itime == 1 .or. itime == 2)then
              i1=i
              i2=i+1
              if(lcyc .and. im1 < i2)then
                 i2=1
              end if              
              j1=j
              j2=j
           else if(itime == 3 .or. itime == 4)then
              i1=i
              i2=i
              j1=j
              j2=j+1
           end if
           
           if(fsm(i1,j1) == 0.d0 .or. fsm(i2,j2) == 0.d0) cycle

           dh_tmp=abs(topo(i2,j2)-topo(i1,j1))/abs(topo(i2,j2)+topo(i1,j1))
           if(dh_tmp <= dh) cycle
           
           !Modification
           dh_cor=0.5d0*(dh_tmp-dh)*abs(topo(i2,j2)+topo(i1,j1))
           
           if(topo(i1,j1) < topo(i2,j2))then
              topo(i2,j2)=topo(i2,j2)-dh_cor
              topo(i1,j1)=topo(i1,j1)+dh_cor
           else
              topo(i2,j2)=topo(i2,j2)+dh_cor
              topo(i1,j1)=topo(i1,j1)-dh_cor
           end if

        end do !i
     end do !j

  end do !itime

  if(lcyc)then
     topo(im-1:im,:)=topo(1:2,:)
  end if
  
  call calculate_ratio(im,jm,topo,fsm,count,dh_max,ratio)
     
end subroutine smooth_filter
