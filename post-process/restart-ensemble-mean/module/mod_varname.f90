module mod_varname

contains

  subroutine varname2d(ivar,varname)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname

    select case(ivar)
    case (1)
       varname = "wubot"
    case (2)
       varname = "wvbot"
    case (3)
       varname = "aam2d"
    case (4)
       varname = "ua"
    case (5)
       varname = "uab"
    case (6)
       varname = "va"
    case (7)
       varname = "vab"
    case (8)
       varname = "el"
    case (9)
       varname = "elb"
    case (10)
       varname = "et"
    case (11)
       varname = "etb"
    case (12)
       varname = "egb"
    case (13)
       varname = "utb"
    case (14)
       varname = "vtb"
    case (15)
       varname = "adx2d"
    case (16)
       varname = "ady2d"
    case (17)
       varname = "advua"
    case (18)
       varname = "advva"
    case default
       varname = ""
       write(*,*) "***Error: invalid ivar in varname2d =", ivar
    end select

  end subroutine varname2d

  !----------------------------------------

  subroutine varname3d(ivar,varname)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname

    select case(ivar)
    case (1)
       varname = "u"
    case (2)
       varname = "ub"
    case (3)
       varname = "v"
    case (4)
       varname = "vb"
    case (5)
       varname = "w"
    case (6)
       varname = "t"
    case (7)
       varname = "tb"
    case (8)
       varname = "s"
    case (9)
       varname = "sb"
    case (10)
       varname = "rho"
    case (11)
       varname = "km"
    case (12)
       varname = "kh"
    case (13)
       varname = "kq"
    case (14)
       varname = "l"
    case (15)
       varname = "q2"
    case (16)
       varname = "q2b"
    case (17)
       varname = "q2l"
    case (18)
       varname = "q2lb"
    case (19)
       varname = "aam"
    case default
       varname = ""
       write(*,*) "***Error: invalid ivar in varname3d =", ivar
    end select

  end subroutine varname3d

end module mod_varname
