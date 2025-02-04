subroutine gradient(itype,im,jm,lon,lat,dat,fsm,dx,dy,ddat)
    
  use mod_parameter
  use mod_rmiss
  implicit none

  integer itype !1:global, 2:regional
  integer i,j
  integer i1,i2,j1,j2

  real(kind = 8) dlon,dlat
  
  !---IN
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: dat(im,jm),fsm(im,jm)
  
  !---OUT  
  real(kind = 8) dx(im,jm),dy(im,jm),ddat(im,jm)

  if(itype /= 1 .and. itype /= 2)then
     write(*,*) "***Error: ",itype
     stop
  end if
  
  do j=1,jm

     if(j==1)then
        j1=j
        j2=j+1
     elseif(j==jm)then
        j1=j-1
        j2=j
     else
        j1=j-1
        j2=j+1
     end if

     dlat=lat(j2)-lat(j1)
       
     do i=1,im

        if(itype == 1)then !global
           if(i==1)then
              i1=im
              i2=i+1
              dlon=lon(i2)-lon(i1)+360.d0
           elseif(i==im)then
              i1=im-1
              i2=1
              dlon=lon(i2)-lon(i1)+360.d0
           else
              i1=i-1
              i2=i+1
              dlon=lon(i2)-lon(i1)
           end if
        elseif(itype == 2)then
           if(i==1)then
              i1=i
              i2=i+1
              dlon=lon(i2)-lon(i1)
           elseif(i==im)then
              i1=i-1
              i2=i
              dlon=lon(i2)-lon(i1)             
           else
              i1=i-1
              i2=i+1
              dlon=lon(i2)-lon(i1)
           end if
        end if
             
        if(fsm(i1,j)==0.d0 .or. fsm(i2,j)==0.d0)then
           dx(i,j)=rmiss
        else
           dx(i,j)=180.d5*(dat(i2,j)-dat(i1,j))/(dlon*pi*earth*cos(pi*lat(j)/180.d0))
        end if

        if(fsm(i,j1)==0.d0 .or. fsm(i,j2)==0.d0)then
           dy(i,j)=rmiss
        else
           dy(i,j)=180.d5*(dat(i,j2)-dat(i,j1))/(dlat*pi*earth)
        end if

        if((fsm(i1,j)==0.d0 .or. fsm(i1,j) == 0.d0) &
             & .and. (fsm(i,j1)==0.d0 .or. fsm(i,j2)==0.d0))then
           ddat(i,j)=rmiss
        else if(fsm(i1,j)==0.d0 .or. fsm(i1,j) == 0.d0)then
           ddat(i,j)=abs(dy(i,j))
        else if(fsm(i,j1)==0.d0 .or. fsm(i,j2)==0.d0)then
           ddat(i,j)=abs(dx(i,j))
        else
           ddat(i,j)=sqrt(dx(i,j)*dx(i,j)+dy(i,j)*dy(i,j))
        end if
          
     end do
  end do
  
end subroutine gradient
