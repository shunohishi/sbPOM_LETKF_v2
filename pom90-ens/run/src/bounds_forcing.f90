!_______________________________________________________________________
subroutine bcond(idx)

  ! apply open boundary conditions by nesting
  ! closed boundary conditions are automatically enabled through
  ! specification of the masks, dum, dvm and fsm, in which case the open
  ! boundary conditions, included below, will be overwritten

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer,intent(in) :: idx

  integer i,j,k
  real(kind = r_size) u1
  real(kind = r_size) hmax

  if(idx == 1) then

     ! External (2-D) elevation boundary conditions
     if(n_west == -1)then
        elf(1,1:jm)=elw(1:jm)
     end if

     if(n_east == -1)then
        elf(im,1:jm)=ele(1:jm)
     end if

     if(n_south == -1)then
        elf(1:im,1)=els(1:im)
     end if

     if(n_north == -1)then
        elf(1:im,jm)=eln(1:im)
     end if

     !$omp parallel
     !$omp do private(i,j)
     do j=1,jm
        do i=1,im
           elf(i,j)=elf(i,j)*fsm(i,j)
        end do
     end do
     !$omp end do
     !$omp end parallel

     return

  else if(idx == 2)then

     ! external (2-D) velocity boundary conditions
     ! west
     !$omp parallel
     if(n_west == -1)then

        !$omp do private(j)
        do j=2,jmm1
           uaf(2,j)=uabw(j)-rfw*sqrt(grav/h(2,j))*(el(2,j)-elw(j))
           uaf(2,j)=ramp*uaf(2,j)
           uaf(1,j)=uaf(2,j)
           vaf(1,j)=vabw(j)
        end do
        !$omp end do
        
     end if

     ! east
     if(n_east == -1)then

        !$omp do private(j)        
        do j=2,jmm1
           uaf(im,j)=uabe(j)+rfe*sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j))
           uaf(im,j)=ramp*uaf(im,j)
           vaf(im,j)=vabe(j)
        end do
        !$omp end do
        
     end if

     ! south
     if(n_south == -1) then

        !$omp do private(i)                
        do i=2,imm1
           vaf(i,2)=vabs(i)-rfs*sqrt(grav/h(i,2))*(el(i,2)-els(i))
           vaf(i,2)=ramp*vaf(i,2)
           vaf(i,1)=vaf(i,2)
           uaf(i,1)=uabs(i)
        end do
        !$omp end do
        
     end if

     ! north
     if(n_north == -1) then

        !$omp do private(i)                        
        do i=2,imm1
           vaf(i,jm)=vabn(i)+rfn*sqrt(grav/h(i,jmm1))*(el(i,jmm1)-eln(i))
           vaf(i,jm)=ramp*vaf(i,jm)
           uaf(i,jm)=uabn(i)
        end do
        !$omp end do

     end if

     !$omp do private(i,j)                             
     do j=1,jm
        do i=1,im
           uaf(i,j)=uaf(i,j)*dum(i,j)
           vaf(i,j)=vaf(i,j)*dvm(i,j)
        end do
     end do
     !$omp end do
     !$omp end parallel

     return

  else if(idx == 3) then
     
     ! internal (3-D) velocity boundary conditions
     ! radiation conditions
     ! smoothing is used in the direction tangential to the boundaries
     hmax=4500.d0

     ! east
     if(n_east == -1)then
        uf(im,2:jmm1,1:kbm1)=ube(2:jmm1,1:kbm1)
        vf(im,2:jmm1,1:kbm1)=vbe(2:jmm1,1:kbm1)
     end if

     ! west
     if(n_west == -1) then
        uf(2,2:jmm1,1:kbm1)=ubw(2:jmm1,1:kbm1)
        uf(1,2:jmm1,1:kbm1)=uf(2,2:jmm1,1:kbm1)
        vf(1,2:jmm1,1:kbm1)=vbw(2:jmm1,1:kbm1)
     end if
     
     ! south
     if(n_south == -1) then
        vf(2:imm1,2,1:kbm1)=vbs(2:imm1,1:kbm1)
        vf(2:imm1,1,1:kbm1)=vf(2:imm1,2,1:kbm1)
        uf(2:imm1,1,1:kbm1)=ubs(2:imm1,1:kbm1)
     end if

     ! north
     if(n_north == -1) then
        vf(2:imm1,jm,1:kbm1)=vbn(2:imm1,1:kbm1)
        uf(2:imm1,jm,1:kbm1)=ubn(2:imm1,1:kbm1)
     end if

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,kbm1
        do j=1,jm
           do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

     return
     
  else if(idx == 4)then

     ! temperature and salinity boundary conditions (using uf and vf,
     ! respectively)
     
     ! east
     if(n_east == -1)then
        uf(im,1:jm,1:kbm1)=tbe(1:jm,1:kbm1) 
        vf(im,1:jm,1:kbm1)=sbe(1:jm,1:kbm1)
     end if

     ! west
     if(n_west == -1)then
        uf(1,1:jm,1:kbm1)=tbw(1:jm,1:kbm1)
        vf(1,1:jm,1:kbm1)=sbw(1:jm,1:kbm1)
     end if
     
     ! south
     if(n_south == -1)then
        uf(1:im,1,1:kbm1)=tbs(1:im,1:kbm1)
        vf(1:im,1,1:kbm1)=sbs(1:im,1:kbm1)
     end if

     ! north
     if(n_north == -1) then
        uf(1:im,jm,1:kbm1)=tbn(1:im,1:kbm1)
        vf(1:im,jm,1:kbm1)=sbn(1:im,1:kbm1)
     end if

     !$omp parallel
     !$omp do private(i,j,k)          
     do k=1,kbm1
        do j=1,jm
           do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

     return
     
  else if(idx == 5) then
     
     ! vertical velocity boundary conditions
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,kbm1
        do j=1,jm
           do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

     return
     
  else if(idx == 6) then

     ! q2 and q2l boundary conditions
     ! west

     !$omp parallel
     if(n_west == -1)then
        !$omp do private(j,k,u1)
        do k=1,kb
           do j=1,jm

              u1=2.d0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1 >= 0.d0) then
                 uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
                 vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
              else
                 uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
                 vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
              end if
              
           end do
        end do
        !$omp end do
     end if

     ! east
     if(n_east == -1)then
        !$omp do private(j,k,u1)
        do k=1,kb
           do j=1,jm

              u1=2.d0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1 <= 0.d0) then
                 uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
                 vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
              else
                 uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
                 vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
              end if
           end do
        end do
        !$omp end do
     end if

     ! south
     if(n_south == -1) then
        !$omp do private(i,k,u1)
        do k=1,kb
           do i=1,im
              u1=2.d0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1 >= 0.d0) then
                 uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
                 vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
              else
                 uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
                 vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
              end if
           end do
        end do
        !$omp end do
     end if
           
     ! north
     if(n_north == -1)then
        !$omp do private(i,k,u1)
        do k=1,kb
           do i=1,im
              u1=2.d0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1 <= 0.d0) then
                 uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
                 vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
              else
                 uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
                 vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
              end if
           end do
        end do
        !$omp end do
     end if

     !$omp do private(i,j,k)     
     do k=1,kb
        do j=1,jm
           do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.d-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.d-10
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

  endif

end subroutine bcond
