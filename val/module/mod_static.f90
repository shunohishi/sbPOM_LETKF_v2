module mod_static

contains

  !----------------------------------------------------------------------------------
  ! Bin Static |
  !----------------------------------------------------------------------------------

  subroutine static_bin_ini( &
       & im_bin,jm_bin, &
       & num_bin,bias_bin,rmsd_bin)

    implicit none

    !---IN
    integer,intent(in) :: im_bin,jm_bin

    !---OUT
    integer,intent(out) :: num_bin(im_bin,jm_bin)
    real(kind = 8),intent(out) :: bias_bin(im_bin,jm_bin),rmsd_bin(im_bin,jm_bin)

    num_bin(:,:)=0
    bias_bin(:,:)=0.d0
    rmsd_bin(:,:)=0.d0
    
  end subroutine static_bin_ini

  !------------------------
  
  subroutine static_bin_add(ijul,n_o,ijul_o,lon_o,lat_o,hdat_a,dat_o, &
       & im_bin,jm_bin,lon_bin,lat_bin,num_bin,bias_bin,rmsd_bin)

    use mod_rmiss
    implicit none

    !---Common
    integer i_o
    integer i_bin,j_bin

    !---IN
    integer,intent(in) :: ijul
    integer,intent(in) :: n_o
    integer,intent(in) :: im_bin,jm_bin

    integer,intent(in) :: ijul_o(n_o)
    
    real(kind = 8),intent(in) :: lon_o(n_o),lat_o(n_o)
    real(kind = 8),intent(in) :: hdat_a(n_o),dat_o(n_o)

    real(kind = 8),intent(in) :: lon_bin(im_bin),lat_bin(jm_bin)
    
    !---INOUT
    integer,intent(inout) :: num_bin(im_bin,jm_bin)

    real(kind = 8),intent(inout) :: bias_bin(im_bin,jm_bin)
    real(kind = 8),intent(inout) :: rmsd_bin(im_bin,jm_bin)

    do i_o=1,n_o

       if(ijul /= ijul_o(i_o)) cycle
       if(lon_o(i_o) < lon_bin(1) .or. lon_bin(im_bin) < lon_o(i_o)) cycle
       if(lat_o(i_o) < lat_bin(1) .or. lat_bin(jm_bin) < lat_o(i_o)) cycle
       if(hdat_a(i_o) == rmiss .or. dat_o(i_o) == rmiss) cycle
       
       do j_bin=1,jm_bin-1
          
          if(lat_o(i_o) < lat_bin(j_bin) .or. lat_bin(j_bin+1) <= lat_o(i_o)) cycle
          
          do i_bin=1,im_bin-1
             
             if(lon_o(i_o) < lon_bin(i_bin) .or. lon_bin(i_bin+1) <= lon_o(i_o)) cycle
             
             num_bin(i_bin,j_bin)=num_bin(i_bin,j_bin)+1
             bias_bin(i_bin,j_bin)=bias_bin(i_bin,j_bin)+(hdat_a(i_o)-dat_o(i_o))
             rmsd_bin(i_bin,j_bin)=rmsd_bin(i_bin,j_bin)+(hdat_a(i_o)-dat_o(i_o))**2
             
          end do !i_bin
       end do !j_bin
       
    end do !i_o
    
  end subroutine static_bin_add

  !------------------------

  subroutine static_bin_end(im_bin,jm_bin,num_bin,bias_bin,rmsd_bin)

    use mod_rmiss
    implicit none

    !---Common
    integer i_bin,j_bin
    
    !---IN
    integer,intent(in) :: im_bin,jm_bin
    integer,intent(in) :: num_bin(im_bin,jm_bin)

    !---INOUT
    real(kind = 8),intent(inout) :: bias_bin(im_bin,jm_bin)
    real(kind = 8),intent(inout) :: rmsd_bin(im_bin,jm_bin)
    
    do j_bin=1,jm_bin
       do i_bin=1,im_bin

          if(num_bin(i_bin,j_bin) == 0)then
             bias_bin(i_bin,j_bin)=rmiss
             rmsd_bin(i_bin,j_bin)=rmiss
          else
             bias_bin(i_bin,j_bin)=bias_bin(i_bin,j_bin)/dble(num_bin(i_bin,j_bin))
             rmsd_bin(i_bin,j_bin)=sqrt(rmsd_bin(i_bin,j_bin)/dble(num_bin(i_bin,j_bin)))
          end if
          
       end do !i_bin
    end do !j_bin
                 
  end subroutine static_bin_end
  
  !-----------------------------------------------------------------------------
  ! Spatial & Temporal Average
  !-----------------------------------------------------------------------------

  subroutine static_ave_ini(num_ave,bias_ave,rmsd_ave)

    implicit none

    !---OUT
    integer,intent(out) :: num_ave

    real(kind = 8),intent(out) :: bias_ave,rmsd_ave

    num_ave=0
    bias_ave=0.d0
    rmsd_ave=0.d0
    
  end subroutine static_ave_ini

  !---------------------------

  subroutine static_ave_add(ijul,n_o,ijul_o,hdat_a,dat_o, &
       & num_ave,bias_ave,rmsd_ave)

    use mod_rmiss
    implicit none

    !---Common
    integer i_o
    
    !---IN
    integer,intent(in) :: ijul
    integer,intent(in) :: n_o
    integer,intent(in) :: ijul_o(n_o)
    
    real(kind = 8),intent(inout) :: hdat_a(n_o),dat_o(n_o)
    
    !---INOUT
    integer,intent(inout) :: num_ave

    real(kind = 8),intent(inout) :: bias_ave,rmsd_ave

    
    do i_o=1,n_o

       if(ijul /= ijul_o(i_o)) cycle
       if(hdat_a(i_o) == rmiss .or. dat_o(i_o) == rmiss) cycle
       
       num_ave=num_ave+1
       bias_ave=bias_ave+hdat_a(i_o)-dat_o(i_o)
       rmsd_ave=rmsd_ave+(hdat_a(i_o)-dat_o(i_o))**2

    end do
    
  end subroutine static_ave_add

  !---------------------------

  subroutine static_ave_end(num_ave,bias_ave,rmsd_ave)

    use mod_rmiss
    implicit none

    !---IN
    integer,intent(in) :: num_ave

    !---INOUT
    real(kind = 8),intent(inout) :: bias_ave,rmsd_ave

    if(num_ave == 0)then
       bias_ave=rmiss
       rmsd_ave=rmiss
    else
       bias_ave=bias_ave/dble(num_ave)
       rmsd_ave=sqrt(rmsd_ave/dble(num_ave))
    end if
    
  end subroutine static_ave_end
  
end module mod_static
