module mod_stat

contains

  !----------------------------------------------------------------------------------
  ! Bin Stat |
  !----------------------------------------------------------------------------------

  subroutine stat_bin_ini( &
       & im_bin,jm_bin, &
       & num_bin,bias_bin,rmsd_bin,sprd_bin)

    implicit none

    !---IN
    integer,intent(in) :: im_bin,jm_bin

    !---OUT
    integer,intent(out) :: num_bin(im_bin,jm_bin)
    real(kind = 8),intent(out) :: bias_bin(im_bin,jm_bin),rmsd_bin(im_bin,jm_bin),sprd_bin(im_bin,jm_bin)

    num_bin(:,:)=0
    bias_bin(:,:)=0.d0
    rmsd_bin(:,:)=0.d0
    sprd_bin(:,:)=0.d0

  end subroutine stat_bin_ini

  !------------------------

  subroutine stat_bin_add(n_o,lon_o,lat_o,hdat_a,hsprd_a,dat_o, &
       & im_bin,jm_bin,lon_bin,lat_bin, &
       & num_bin,bias_bin,rmsd_bin,sprd_bin)

    use mod_rmiss
    implicit none

    !---Common
    integer i_o
    integer i_bin,j_bin

    !---IN
    integer,intent(in) :: n_o
    integer,intent(in) :: im_bin,jm_bin

    real(kind = 8),intent(in) :: lon_o(n_o),lat_o(n_o)
    real(kind = 8),intent(in) :: hdat_a(n_o),hsprd_a(n_o),dat_o(n_o)

    real(kind = 8),intent(in) :: lon_bin(im_bin),lat_bin(jm_bin)

    !---INOUT
    integer,intent(inout) :: num_bin(im_bin,jm_bin)

    real(kind = 8),intent(inout) :: bias_bin(im_bin,jm_bin)
    real(kind = 8),intent(inout) :: rmsd_bin(im_bin,jm_bin)
    real(kind = 8),intent(inout) :: sprd_bin(im_bin,jm_bin)

    do i_o=1,n_o

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
             sprd_bin(i_bin,j_bin)=sprd_bin(i_bin,j_bin)+hsprd_a(i_o)

          end do !i_bin
       end do !j_bin

    end do !i_o

  end subroutine stat_bin_add

  !------------------------

  subroutine stat_bin_end(im_bin,jm_bin,num_bin,bias_bin,rmsd_bin,sprd_bin)

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
    real(kind = 8),intent(inout) :: sprd_bin(im_bin,jm_bin)

    do j_bin=1,jm_bin
       do i_bin=1,im_bin

          if(num_bin(i_bin,j_bin) == 0)then
             bias_bin(i_bin,j_bin)=rmiss
             rmsd_bin(i_bin,j_bin)=rmiss
             sprd_bin(i_bin,j_bin)=rmiss
          else
             bias_bin(i_bin,j_bin)=bias_bin(i_bin,j_bin)/dble(num_bin(i_bin,j_bin))
             rmsd_bin(i_bin,j_bin)=sqrt(rmsd_bin(i_bin,j_bin)/dble(num_bin(i_bin,j_bin)))
             sprd_bin(i_bin,j_bin)=sprd_bin(i_bin,j_bin)/dble(num_bin(i_bin,j_bin))
          end if

       end do !i_bin
    end do !j_bin

  end subroutine stat_bin_end

  !-----------------------------------------------------------------------------
  ! Spatial & Temporal Average
  !-----------------------------------------------------------------------------

  subroutine stat_ave_ini(num_ave,bias_ave,rmsd_ave,sprd_ave)

    implicit none

    !---OUT
    integer,intent(out) :: num_ave

    real(kind = 8),intent(out) :: bias_ave,rmsd_ave,sprd_ave

    num_ave=0
    bias_ave=0.d0
    rmsd_ave=0.d0
    sprd_ave=0.d0

  end subroutine stat_ave_ini

  !---------------------------

  subroutine stat_ave_add(n_o,hdat_a,hsprd_a,dat_o, &
       & num_ave,bias_ave,rmsd_ave,sprd_ave)

    use mod_rmiss
    implicit none

    !---Common
    integer i_o

    !---IN
    integer,intent(in) :: n_o

    real(kind = 8),intent(inout) :: hdat_a(n_o),hsprd_a(n_o),dat_o(n_o)

    !---INOUT
    integer,intent(inout) :: num_ave

    real(kind = 8),intent(inout) :: bias_ave,rmsd_ave,sprd_ave


    do i_o=1,n_o

       if(hdat_a(i_o) == rmiss .or. dat_o(i_o) == rmiss) cycle

       num_ave=num_ave+1
       bias_ave=bias_ave+hdat_a(i_o)-dat_o(i_o)
       rmsd_ave=rmsd_ave+(hdat_a(i_o)-dat_o(i_o))**2
       sprd_ave=sprd_ave+hsprd_a(i_o)

    end do

  end subroutine stat_ave_add

  !---------------------------

  subroutine stat_ave_end(num_ave,bias_ave,rmsd_ave,sprd_ave)

    use mod_rmiss
    implicit none

    !---IN
    integer,intent(in) :: num_ave

    !---INOUT
    real(kind = 8),intent(inout) :: bias_ave,rmsd_ave,sprd_ave

    if(num_ave == 0)then
       bias_ave=rmiss
       rmsd_ave=rmiss
       sprd_ave=rmiss
    else
       bias_ave=bias_ave/dble(num_ave)
       rmsd_ave=sqrt(rmsd_ave/dble(num_ave))
       sprd_ave=sprd_ave/dble(num_ave)
    end if

  end subroutine stat_ave_end

  !-----------------------------------------------------------------------------
  ! Spatial & Temporal 1D Average
  !-----------------------------------------------------------------------------

  subroutine stat_1d_ini(im,num_ave,bias_ave,rmsd_ave,sprd_ave)

    implicit none

    !---IN
    integer,intent(in) :: im

    !---OUT
    integer,intent(out) :: num_ave(im)

    real(kind = 8),intent(out) :: bias_ave(im),rmsd_ave(im),sprd_ave(im)

    num_ave(:)=0
    bias_ave(:)=0.d0
    rmsd_ave(:)=0.d0
    sprd_ave(:)=0.d0

  end subroutine stat_1d_ini

  !---------------------------

  subroutine stat_1d_add(im,hdat_a,hsprd_a,dat_o, &
       & num_ave,bias_ave,rmsd_ave,sprd_ave)

    use mod_rmiss
    implicit none

    !---Common
    integer i

    !---IN
    integer,intent(in) :: im

    real(kind = 8),intent(inout) :: hdat_a(im),hsprd_a(im),dat_o(im)

    !---INOUT
    integer,intent(inout) :: num_ave(im)

    real(kind = 8),intent(inout) :: bias_ave(im),rmsd_ave(im),sprd_ave(im)


    do i=1,im

       if(hdat_a(i) == rmiss .or. dat_o(i) == rmiss) cycle

       num_ave(i)=num_ave(i)+1
       bias_ave(i)=bias_ave(i)+hdat_a(i)-dat_o(i)
       rmsd_ave(i)=rmsd_ave(i)+(hdat_a(i)-dat_o(i))**2
       sprd_ave(i)=sprd_ave(i)+hsprd_a(i)

    end do

  end subroutine stat_1d_add

  !---------------------------

  subroutine stat_1d_end(im,num_ave,bias_ave,rmsd_ave,sprd_ave)

    use mod_rmiss
    implicit none

    !---Common
    integer i

    !---IN
    integer,intent(in) :: im
    integer,intent(in) :: num_ave(im)

    !---INOUT
    real(kind = 8),intent(inout) :: bias_ave(im),rmsd_ave(im),sprd_ave(im)

    do i=1,im
       if(num_ave(i) == 0)then
          bias_ave(i)=rmiss
          rmsd_ave(i)=rmiss
          sprd_ave(i)=rmiss
       else
          bias_ave(i)=bias_ave(i)/dble(num_ave(i))
          rmsd_ave(i)=sqrt(rmsd_ave(i)/dble(num_ave(i)))
          sprd_ave(i)=sprd_ave(i)/dble(num_ave(i))
       end if
    end do

  end subroutine stat_1d_end

  !-----------------------------------------------------------------------
  ! Basic 
  !-----------------------------------------------------------------------

  subroutine average(n,dat,ave)

    implicit none

    !---Common
    integer i

    !---IN
    integer,intent(in) :: n

    real(kind = 8),intent(in) :: dat(n)

    !---OUT
    real(kind = 8),intent(out) :: ave

    ave=0.d0
    do i=1,n
       ave=ave+dat(i)
    end do

    ave=ave/dble(n)

  end subroutine average

  !-------------------------

  subroutine standard_deviation(n,dat,ave,std)

    implicit none

    !---Common
    integer i

    !---IN
    integer,intent(in) :: n

    real(kind = 8),intent(in) :: dat(n),ave

    !---OUT
    real(kind = 8),intent(out) :: std

    std=0.d0
    do i=1,n
       std=std+(dat(i)-ave)*(dat(i)-ave)
    end do

    std=sqrt(std/dble(n-1))

  end subroutine standard_deviation

  !-------------------------------------------------------------

  subroutine correlation(n,dat1,dat2,cor)

    implicit none

    !---Common
    integer i

    real(kind = 8) ave1,ave2
    real(kind = 8) std1,std2
    real(kind = 8) cov

    !---IN
    integer,intent(in) :: n

    real(kind = 8),intent(in) :: dat1(n),dat2(n)

    !---OUT
    real(kind = 8) cor

    !Average
    call average(n,dat1,ave1)
    call average(n,dat2,ave2)

    !Standard deviation
    call standard_deviation(n,dat1,ave1,std1)
    call standard_deviation(n,dat2,ave2,std2)

    !Covariance
    cov=0.d0
    do i=1,n
       cov=cov+(dat1(i)-ave1)*(dat2(i)-ave2)
    end do

    cov=cov/dble(n-1)

    !Correlation
    cor=cov/(std1*std2)

  end subroutine correlation

  !-----------------------------------------------------------

  subroutine average_rmiss(n,dat,ave)

    use mod_rmiss
    implicit none

    !---Common
    integer i
    integer np

    !---IN
    integer,intent(in) :: n

    real(kind = 8),intent(in) :: dat(n)

    !---OUT
    real(kind = 8),intent(out) :: ave

    np=0    
    ave=0.d0

    do i=1,n
       if(dat(i) == rmiss) cycle
       ave=ave+dat(i)
       np=np+1
    end do

    if(np == 0)then
       ave=rmiss
    else
       ave=ave/dble(np)
    end if
       
  end subroutine average_rmiss

  !-------------------------

  subroutine standard_deviation_rmiss(n,dat,ave,std)

    use mod_rmiss
    implicit none

    !---Common
    integer i
    integer np

    !---IN
    integer,intent(in) :: n

    real(kind = 8),intent(in) :: dat(n),ave

    !---OUT
    real(kind = 8),intent(out) :: std

    np=0
    std=0.d0

    do i=1,n
       if(dat(i) == rmiss) cycle
       std=std+(dat(i)-ave)*(dat(i)-ave)
       np=np+1
    end do

    if(np <= 1)then
       std=rmiss
    else
       std=sqrt(std/dble(np-1))
    end if

  end subroutine standard_deviation_rmiss

  !--------------------------

  subroutine correlation_rmiss(n,dat1,dat2,cor)

    use mod_rmiss
    implicit none

    !---Common
    integer i
    integer np

    real(kind = 8) ave1,ave2
    real(kind = 8) std1,std2
    real(kind = 8) cov

    !---IN
    integer,intent(in) :: n

    real(kind = 8),intent(in) :: dat1(n),dat2(n)

    !---OUT
    real(kind = 8) cor

    !Average
    call average(n,dat1,ave1)
    call average(n,dat2,ave2)

    !Standard deviation
    call standard_deviation(n,dat1,ave1,std1)
    call standard_deviation(n,dat2,ave2,std2)

    if(ave1 == rmiss .or. ave2 == rmiss .or. std1 == rmiss .or. std2 == rmiss)then

       cov=rmiss

    else

       !Covariance
       np=0
       cov=0.d0
       do i=1,n
          if(dat1(i) == rmiss) cycle
          cov=cov+(dat1(i)-ave1)*(dat2(i)-ave2)
          np=np+1
       end do

       cov=cov/dble(np-1)

       !Correlation
       cor=cov/(std1*std2)

    end if

  end subroutine correlation_rmiss
  
  !---------------------------------------------------------------------
  ! Paired t-test
  !---------------------------------------------------------------------

  subroutine paired_t_test(nall,dat1,dat2,dof,tcrit,tval)

    use mod_rmiss
    implicit none

    !---Parameter
    integer,parameter :: ndof=100
    real(kind=8), parameter :: tcrit_table(ndof) = (/ &
         63.657d0, 9.925d0, 5.841d0, 4.604d0, 4.032d0, &  ! df=1-5
         3.707d0, 3.499d0, 3.355d0, 3.250d0, 3.169d0, &  ! df=6-10
         3.106d0, 3.055d0, 3.012d0, 2.977d0, 2.947d0, &  ! df=11-15
         2.921d0, 2.898d0, 2.878d0, 2.861d0, 2.845d0, &  ! df=16-20
         2.831d0, 2.819d0, 2.807d0, 2.797d0, 2.787d0, &  ! df=21-25
         2.779d0, 2.771d0, 2.763d0, 2.756d0, 2.750d0, &  ! df=26-30
         2.744d0, 2.738d0, 2.733d0, 2.728d0, 2.724d0, &  ! df=31-35
         2.719d0, 2.715d0, 2.712d0, 2.708d0, 2.704d0, &  ! df=36-40
         2.701d0, 2.698d0, 2.695d0, 2.692d0, 2.690d0, &  ! df=41-45
         2.687d0, 2.685d0, 2.682d0, 2.680d0, 2.678d0, &  ! df=46-50
         2.676d0, 2.674d0, 2.672d0, 2.670d0, 2.668d0, &  ! df=51-55
         2.667d0, 2.665d0, 2.663d0, 2.662d0, 2.660d0, &  ! df=56-60
         2.659d0, 2.657d0, 2.656d0, 2.655d0, 2.653d0, &  ! df=61-65
         2.652d0, 2.651d0, 2.650d0, 2.648d0, 2.647d0, &  ! df=66-70
         2.646d0, 2.645d0, 2.644d0, 2.643d0, 2.642d0, &  ! df=71-75
         2.641d0, 2.640d0, 2.639d0, 2.638d0, 2.638d0, &  ! df=76-80
         2.637d0, 2.636d0, 2.635d0, 2.634d0, 2.634d0, &  ! df=81-85
         2.633d0, 2.632d0, 2.632d0, 2.631d0, 2.630d0, &  ! df=86-90
         2.630d0, 2.629d0, 2.629d0, 2.628d0, 2.627d0, &  ! df=91-95
         2.627d0, 2.626d0, 2.626d0, 2.626d0, 2.626d0  /)  ! df=96-100 (99%)

    !    real(kind = 8),parameter ::  tcrit_table(ndof) = (/   &  !Paired t-test, two-side test
    !         & 12.706d0, 4.303d0, 3.182d0, 2.776d0, 2.571d0,  &  ! df=1-5
    !         &  2.447d0, 2.365d0, 2.306d0, 2.262d0, 2.228d0,  &  ! df=6-10
    !         &  2.201d0, 2.179d0, 2.160d0, 2.145d0, 2.131d0,  &  ! df=11-15
    !         &  2.120d0, 2.110d0, 2.101d0, 2.093d0, 2.086d0,  &  ! df=16-20
    !         &  2.080d0, 2.074d0, 2.069d0, 2.064d0, 2.060d0,  &  ! df=21-25
    !         &  2.056d0, 2.052d0, 2.048d0, 2.045d0, 2.042d0,  &  ! df=26-30
    !         &  2.040d0, 2.037d0, 2.035d0, 2.032d0, 2.030d0,  &  ! df=31-35
    !         &  2.028d0, 2.026d0, 2.024d0, 2.023d0, 2.021d0,  &  ! df=36-40
    !         &  2.020d0, 2.018d0, 2.017d0, 2.015d0, 2.014d0,  &  ! df=41-45
    !         &  2.013d0, 2.012d0, 2.011d0, 2.010d0, 2.009d0 /)   ! df=46-50

    !---Common
    integer i,n

    real(kind = 8) ave,std,cor
    real(kind = 8),allocatable :: dif(:)

    !---IN
    integer,intent(in) :: nall

    real(kind = 8),intent(in) :: dat1(nall),dat2(nall)

    !---OUT
    integer,intent(out) :: dof
    real(kind = 8),intent(out) :: tcrit,tval

    !Initial
    dof=int(rmiss)
    tcrit=rmiss
    tval=rmiss

    !Check missing value
    n=0
    do i=1,nall
       if(dat1(i) == rmiss .or. dat2(i) == rmiss)then
          cycle
       else
          n=n+1
       end if
    end do

    if(n <= 3)then
       return
    end if

    !Difference
    allocate(dif(n))

    n=0
    do i=1,nall
       if(dat1(i) == rmiss .or. dat2(i) == rmiss)then
          cycle
       else
          n=n+1
          dif(n)=dat2(i)-dat1(i)
       end if
    end do

    call average(n,dif,ave)
    call standard_deviation(n,dif,ave,std)

    !Lag correlation
    call correlation(n-1,dif(1:n-1),dif(2:n),cor)

    !Effective degree of freedom (Bretherton et al. 1999)
    dof=int(n*(1.d0-cor)/(1.d0+cor))
    dof=max(dof,2) !Minimum DOF = 2 becasue of t value

    !T critical value
    if(dof <= 0)then
       return
    else if(ndof < dof)then
       tcrit=1.96d0
    else
       tcrit=tcrit_table(dof)
    end if

    !T-value
    tval=sqrt(dof-1.d0)*ave/std

    write(*,'(a,f12.5,x,a,f12.5)') "Ave:",ave,"Std:",std
    write(*,'(a,f12.5)') "Lag-correlaction:",cor
    write(*,'(a,i6)') "Effective degree of freedom:",dof
    write(*,'(a,f12.5)') "T value:",tval

    deallocate(dif)

  end subroutine paired_t_test

  !---------------------------------------------------------------------
  ! Bootstrap |
  !---------------------------------------------------------------------
  ! call init_random() in main program
  !---------------------------------------------------------------------

  subroutine bootstrap(nall,dat1,dat2,dat_low,dat_ave,dat_upp)

    use mod_rmiss
    implicit none

    !---Parameter
    integer,parameter :: nboot=1000 !number of sample
    real,parameter :: sig_lev=0.01 !two-sided significance level [-]
    
    !---Common
    integer i,n
    integer iboot
    integer id(nboot),id_s,id_e

    real(kind = 8),allocatable :: dif(:)
    real(kind = 8),allocatable :: random(:),sample(:)
    real(kind = 8),allocatable :: sample_ave(:)
    
    !---IN
    integer,intent(in) :: nall

    real(kind = 8),intent(in) :: dat1(nall),dat2(nall)

    !---OUT
    real(kind = 8),intent(out) :: dat_low,dat_ave,dat_upp

    !---Check missing value
    n=0
    do i=1,nall
       if(dat1(i) == rmiss .or. dat2(i) == rmiss)then
          cycle
       else
          n=n+1
       end if
    end do

    if(n == 0)then
       dat_low=rmiss
       dat_ave=rmiss
       dat_upp=rmiss       
       return
    end if

    !---Difference
    allocate(dif(n))

    n=0
    do i=1,nall
       if(dat1(i) == rmiss .or. dat2(i) == rmiss)then
          cycle
       else
          n=n+1
          dif(n)=dat2(i)-dat1(i)
       end if
    end do

    !---Average
    call average(n,dif,dat_ave)

    !---Sampling
    allocate(random(n),sample(n),sample_ave(nboot))
    !call init_random() ==> Call Main Program
    do iboot=1,nboot

       call random_number(random) !0<= random < 1

       do i=1,n
          sample(i)=dif( 1+ int(random(i)*n) )
       end do

       id(iboot)=iboot
       call average(n,sample,sample_ave(iboot))

    end do

    !---Sort
    id_s=1
    id_e=nboot
    call quick_sort_asnd(nboot,sample_ave,id,id_s,id_e)

    !---Low & Upper
    i=max(1, int( (0.5d0*sig_lev)*nboot))
    dat_low=sample_ave(i)

    i=min(int( (1.d0-0.5d0*sig_lev)*nboot), nboot)
    dat_upp=sample_ave(i)
    
    deallocate(dif)
    deallocate(random,sample)
    deallocate(sample_ave)
    
  end subroutine bootstrap

  !---------------------------

  subroutine init_random()

    implicit none

    !---Common
    integer i
    integer seedsize
    integer,allocatable :: seed(:)

    call random_seed(size=seedsize)
    allocate(seed(seedsize))

    do i=1,seedsize
       seed(i)=12345+1000*i
    end do
    
    call random_seed(put=seed)
        
    deallocate(seed)
    
  end subroutine init_random

  !--------------------------

  RECURSIVE SUBROUTINE quick_sort_asnd(n,dat,idx,idx_s,idx_e)

    IMPLICIT NONE

    !---IN
    INTEGER,INTENT(IN) :: n
    INTEGER,INTENT(IN) :: idx_s,idx_e

    !---INOUT
    INTEGER,INTENT(INOUT) :: idx(n)

    REAL(kind = 8),INTENT(INOUT) :: dat(n)

    !---WORK
    INTEGER i,j,it
    REAL(kind = 8) x,t

    x = dat((idx_s+idx_e)/2)
    i = idx_s
    j = idx_e

    do

       do while (dat(i) < x)
          i=i+1
       end do

       do while (x < dat(j))
          j=j-1
       end do

       if(j <= i) exit

       t=dat(i)
       dat(i)=dat(j)
       dat(j)=t

       it=idx(i)
       idx(i)=idx(j)
       idx(j)=it

       i=i+1
       j=j-1

    end do

    if(idx_s < i-1) call quick_sort_asnd(n,dat,idx,idx_s,i-1)
    if(j+1 < idx_e) call quick_sort_asnd(n,dat,idx,j+1,idx_e)

  END SUBROUTINE quick_sort_asnd
  
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  subroutine substitute_rmiss(ndat,dof,tcrit,tval)

    use mod_rmiss
    implicit none

    !---IN
    integer,intent(in) :: ndat

    !---OUT
    integer,intent(out) :: dof(ndat,ndat)
    real(kind = 8),intent(out) :: tcrit(ndat,ndat),tval(ndat,ndat)

    dof(:,:)=int(rmiss)
    tcrit(:,:)=rmiss
    tval(:,:)=rmiss

  end subroutine substitute_rmiss

  !---------------------------------------------------

  subroutine convert_absolute_bias(n,bias,abias)

    use mod_rmiss
    implicit none

    !---Common
    integer i

    !---IN
    integer,intent(in) :: n
    real(kind = 8),intent(in) :: bias(n)

    !---OUT
    real(kind = 8),intent(out) :: abias(n)

    do i=1,n
       if(bias(i) == rmiss)then
          abias(i)=rmiss
       else
          abias(i)=abs(bias(i))
       end if
    end do

  end subroutine convert_absolute_bias

end module mod_stat
