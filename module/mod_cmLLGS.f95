module mod_cmLLGS !module of LLG+spin-transfer torque for macrospin
  use constant
  use mod_mLLGS
  use mod_rfpmLLGS
  implicit none
  type cmLLGS
    type(rfpmLLGS) r
    logical lout
    character(80) filename,filename_mh
    double precision Hmin,Hmax,dH
  contains
    procedure::init=>cmLLGS_init
    procedure::initconf=>cmLLGS_initconf
    procedure::readparam=>cmLLGS_readparam
    procedure::calc=>cmLLGS_calc
    procedure::calcmh=>cmLLGS_calcmh
  end type
contains

  subroutine cmLLGS_init(s)
    use constant
    implicit none
    class(cmLLGS) s

    s%filename='mLLGS.txt'
    s%filename_mh='mh_mLLGS.txt'

    call type_mLLGS%init()

  end subroutine

  subroutine cmLLGS_readparam(s)
    implicit none
    class(cmLLGS) s
    integer i,j

      call s%r%init()
      call s%r%read()
      call s%r%ana()
      call s%r%print()
      call s%r%close()

      type_mLLGS%nl=s%r%nl
      
      call type_mLLGS%alloc()

      type_mLLGS%dt = s%r%dt
      type_mLLGS%N = s%r%n
      type_mLLGS%Neq = s%r%neq
      type_mLLGS%Hext = s%r%Hext
      type_mLLGS%theta_h = s%r%theta_h
      type_mLLGS%phi_h = s%r%phi_h
      type_mLLGS%alpha = s%r%alpha
      type_mLLGS%Ms = s%r%Ms
      type_mLLGS%Ku = s%r%Ku
      type_mLLGS%Kuni = s%r%Kuni
      type_mLLGS%dF = s%r%dF
      type_mLLGS%Jex = s%r%Jex
      type_mLLGS%Jex2 = s%r%Jex2
      type_mLLGS%DMI = s%r%DMI

      s%Hmin=s%r%Hmin
      s%Hmax=s%r%Hmax
      s%dH=s%r%dH

      type_mLLGS%dir_D =s%r%dir_D
      type_mLLGS%dir_uni = s%r%dir_uni

      call type_mLLGS%calcparam()
      
      s%filename = s%r%cfilename
      s%filename_mh = s%r%cfilename_mh

      type_mLLGS%mx = s%r%mx
      type_mLLGS%my = s%r%my
      type_mLLGS%mz = s%r%mz

  end subroutine

  subroutine cmLLGS_initconf(s)
    implicit none
    class(cmLLGS) s

      type_mLLGS%mx = s%r%mx
      type_mLLGS%my = s%r%my
      type_mLLGS%mz = s%r%mz

  end subroutine

  subroutine cmLLGS_calcmh(s)
    implicit none
    class(cmLLGS) s
    integer i,j,n
    double precision Hext

    type_mLLGS%alpha=s%r%alpha_mh
    n=(s%Hmax-s%Hmin)/s%dH
    open(20,file='data/'//trim(s%filename_mh),status='replace')
    if(s%r%print_mh_all .eqv. .true.)then
      open(25,file='data/all_'//trim(s%filename_mh),status='replace')
    else 
    end if
    do i=1,n
      Hext=s%Hmax-dble(i-1)*s%dH
      type_mLLGS%Hext=Hext
      if(s%r%nint_mh==1)then
        call s%initconf()
      end if
      call type_mLLGS%calcequiv()
      call type_mLLGS%calcsum()
      if(s%r%print_mh_all .eqv. .true.)then
        write(20,'(100e16.5)')Hext,(type_mLLGS%mvec(j),j=1,3),(type_mLLGS%mx(j),type_mLLGS%my(j),type_mLLGS%mz(j),j=1,type_mLLGS%nl)
        write(25,'(100e16.5)')(0.0d0,0.0d0,0.0d0+dble(j),&
        type_mLLGS%mx(j),type_mLLGS%my(j),type_mLLGS%mz(j)+dble(j),j=1,type_mLLGS%nl)
      else
        write(20,'(4e16.5)')Hext,(type_mLLGS%mvec(j),j=1,3)
      end if
    end do

    do i=1,n
      Hext=s%Hmin+dble(i-1)*s%dH
      type_mLLGS%Hext=Hext
      if(s%r%nint_mh==1)then
        call s%initconf()
      end if
      call type_mLLGS%calcequiv()
      call type_mLLGS%calcsum()
      if(s%r%print_mh_all .eqv. .true.)then
        write(20,'(100e16.5)')Hext,(type_mLLGS%mvec(j),j=1,3),(type_mLLGS%mx(j),type_mLLGS%my(j),type_mLLGS%mz(j),j=1,type_mLLGS%nl)
        write(25,'(100e16.5)')(0.0d0,0.0d0,0.0d0+dble(j),&
        type_mLLGS%mx(j),type_mLLGS%my(j),type_mLLGS%mz(j)+dble(j),j=1,type_mLLGS%nl)
      else
        write(20,'(4e16.5)')Hext,(type_mLLGS%mvec(j),j=1,3)
      end if
    end do
    close(20)
    if(s%r%print_mh_all .eqv. .true.)then
      close(25)
    else 
    end if

    return
  end subroutine

  subroutine cmLLGS_calc(s)
    implicit none
    class(cmLLGS) s
    integer i,j,k

    if(s%lout .eqv. .true.)then
      do k=1,type_mLLGS%nl
        write(*,'(A,i2,4e12.4)')'initial magnetization = ',k,type_mLLGS%mx(k),type_mLLGS%my(k),type_mLLGS%mz(k),type_mLLGS%norm(k)
      end do
    end if

    type_mLLGS%hpx=0.0d0
    type_mLLGS%hpy=0.0d0
    type_mLLGS%hpz=0.0d0
    !write(*,'(6e12.4)')(s%mx(j),s%my(j),s%mz(j),j=1,2)
    open(10,file=trim(s%filename),status='replace')
    do i=1,type_mLLGS%N

      type_mLLGS%time=type_mLLGS%dt*dble(i)

      call type_mLLGS%calctimestep()
      !write(*,'(6e12.4)')(s%mx(j),s%my(j),s%mz(j),j=1,2)

      if(isnan(type_mLLGS%mx(1)) .eqv. .true.)then
      write(*,'(A,7e12.4)')'Nan is found',type_mLLGS%time,(type_mLLGS%mx(j),type_mLLGS%my(j),type_mLLGS%mz(j),j=1,2)
      exit
      end if

      do k=1,type_mLLGS%nl
        write(*,'(A,i2,4e16.5)')'magnetization = ',k,type_mLLGS%mx(k),type_mLLGS%my(k),type_mLLGS%mz(k),type_mLLGS%norm(k)
      end do

    end do
    close(10)

    if(s%lout .eqv. .true.)then
      do k=1,type_mLLGS%nl
        write(*,'(A,i2,4e12.4)')'magnetization = ',k,type_mLLGS%mx(k),type_mLLGS%my(k),type_mLLGS%mz(k),type_mLLGS%norm(k)
      end do
    end if

    return
  end subroutine

end module


  