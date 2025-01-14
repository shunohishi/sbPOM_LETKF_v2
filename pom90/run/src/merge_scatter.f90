!  merge2d() merges a 2-D array into a 2-D array in rank0 process.

subroutine merge2d(dest,src)

  !$use omp_lib  
  use mpi
  use common_pom_var
  implicit none

  !Common
  integer nblk
  integer ijloc,iloc,jloc
  integer iglb,jglb
  integer m,n,ierr
  
  real(kind = r_size) work(im*jm*n_proc)

  !IN/OUT
  real(kind = r_size),intent(in)  :: src(im,jm)
  real(kind = r_sngl),intent(out) :: dest(im_global,jm_global)

  nblk=im*jm

  if(r_size == kind(0.0d0))then
     call MPI_Gather(src,nblk,MPI_DOUBLE_PRECISION, &
          & work,nblk,MPI_DOUBLE_PRECISION, &
          & master_task,pom_comm,ierr)
  else if(r_size == kind(0.0e0))then
     call MPI_Gather(src,nblk,MPI_REAL, &
          & work,nblk,MPI_REAL, &
          & master_task,pom_comm,ierr)     
  end if
  
  if(my_task == master_task)then

     !$omp parallel
     !$omp do private(m,n,ijloc,iloc,jloc,iglb,jglb)
     do m=1,nblk*n_proc
        
        n=int(m/nblk)
        
        ijloc=mod(m,nblk)
        if(ijloc == 0)then
           n=n-1
           ijloc=nblk
        end if
            
        jloc=int(ijloc/im)+1
        iloc=mod(ijloc,im)
        if(iloc == 0)then
           jloc=jloc-1
           iloc=im
        end if
            
        iglb=iloc+mod(n,nproc_x)*(im-2)
        jglb=jloc+(n/nproc_x)*(jm-2)
            
        dest(iglb,jglb)=work(m)
            
     end do
     !$omp end do
     !$omp end parallel
     
  endif
      
  call MPI_Barrier(pom_comm,ierr)

end subroutine merge2d

!____________________________________________________________________________
!  merge3d() merges a 3D array (src) in each process into
!  a array (dest) in rank0 process.

subroutine merge3d(dest,src)

  !$use omp_lib  
  use mpi
  use common_pom_var
  
  implicit none

  !Common
  integer nblk
  integer kloc,ijkloc,ijloc,iloc,jloc
  integer iglb,jglb
  integer m,n
  integer ierr
  
  real(kind = r_size) work(im*jm*kb*n_proc)

  !IN/OUT
  real(kind = r_size),intent(in)  :: src(im,jm,kb)
  real(kind = r_sngl),intent(out) :: dest(im_global,jm_global,kb)


  nblk=im*jm*kb

  if(r_size == kind(0.0d0))then
     call MPI_Gather(src,nblk,MPI_DOUBLE_PRECISION, &
          & work,nblk,MPI_DOUBLE_PRECISION, &
          & master_task,pom_comm,ierr)
  else if(r_size == kind(0.0e0))then
     call MPI_Gather(src,nblk,MPI_REAL, &
          & work,nblk,MPI_REAL, &
          & master_task,pom_comm,ierr)
  end if
     
  if(my_task == master_task)then

     !$omp parallel
     !$omp do private(m,n,ijkloc,ijloc,iloc,jloc,kloc,iglb,jglb)
     do m=1,nblk*n_proc

        n=int(m/nblk)
        
        ijkloc=mod(m,nblk)
        if(ijkloc == 0)then
           n=n-1
           ijkloc=nblk
        end if
        
        kloc=int(ijkloc/(im*jm))+1
        ijloc=mod(ijkloc,im*jm)
        
        if(ijloc == 0)then
           kloc=kloc-1
           ijloc=im*jm
        end if
        
        jloc=int(ijloc/im)+1
        
        iloc=mod(ijloc,im)
        if(iloc == 0)then
           jloc=jloc-1
           iloc=im
        end if

        iglb=iloc+mod(n,nproc_x)*(im-2)
        jglb=jloc+(n/nproc_x)*(jm-2)

            
        dest(iglb,jglb,kloc)=work(m)
     end do
     !$omp end do
     !$omp end parallel
     
  endif
      
  call MPI_Barrier(pom_comm,ierr)

end subroutine merge3d

!_______________________________________________________________________________________
!  scatter2d() scatter a global 2-D array into a local 2-D array

subroutine scatter2d(dest,src)

  !$use omp_lib  
  use mpi
  use common_pom_var
  
  implicit none

  !Common
  integer iloc,jloc

  !IN/OUT
  real(kind = r_dble),intent(in)  :: src(im_global,jm_global)
  real(kind = r_size),intent(out) :: dest(im,jm)

  !$omp parallel
  !$omp do private(iloc,jloc)
  do jloc=1,jm
     do iloc=1,im
        dest(iloc,jloc)=src(i_global(iloc),j_global(jloc))
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine scatter2d

!________________________________________________________________________________________
!  scatter3d() scatter a global 3-D array into a local 3-D array

subroutine scatter3d(dest,src)

  !$use omp_lib  
  use mpi
  use common_pom_var
  implicit none

  !Common
  integer iloc,jloc,k

  !IN/OUT
  real(kind = r_dble),intent(in)  :: src(im_global,jm_global,kb)
  real(kind = r_size),intent(out) :: dest(im,jm,kb)

  !$omp parallel
  !$omp do private(iloc,jloc,k)  
  do k=1,kb
     do jloc=1,jm
        do iloc=1,im
           dest(iloc,jloc,k)=src(i_global(iloc),j_global(jloc),k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine scatter3d

!_____________________________________________________________________________
!  scatter_bnd_ew() 
!    scatter a global 2-D/3-D array into 
!    a eastern/western boundary of local 2-D array

subroutine scatter_bnd_ew(dest,src,kmax)

  !$use omp_lib  
  use mpi
  use common_pom_var
  
  implicit none

  !Common
  integer kmax
  integer jloc,k

  !IN/OUT
  real(kind = r_dble),intent(in)  :: src(jm_global,kmax)
  real(kind = r_size),intent(out) :: dest(jm,kmax)

  if(n_west == -1 .or. n_east == -1)then

     !$omp parallel
     !$omp do private(jloc,k)       
     do k=1,kmax
        do jloc=1,jm
           dest(jloc,k)=src(j_global(jloc),k)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if
  
end subroutine scatter_bnd_ew

!_______________________________________________________________________________
!  scatter_bnd_ns() 
!    scatter a global 2-D/3-D array into 
!    a northern/southern boundary of local 2-D array

subroutine scatter_bnd_ns(dest,src,kmax)

  !$use omp_lib  
  use mpi
  use common_pom_var
  
  implicit none

  !Common
  integer kmax
  integer iloc,k

  !IN/OUT
  real(kind = r_dble),intent(in)  :: src(im_global,kmax)
  real(kind = r_size),intent(out) :: dest(im,kmax)

  if(n_north == -1 .or. n_south == -1)then

     !$omp parallel
     !$omp do private(iloc,k)            
     do k=1,kmax
        do iloc=1,im
           dest(iloc,k)=src(i_global(iloc),k)
        end do
     end do
     !$omp end do
     !$omp end parallel

  end if
  
end subroutine scatter_bnd_ns
