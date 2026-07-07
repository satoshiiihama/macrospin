module mod_twoLLGSdynamics !module of LLG+spin-transfer torque for macrospin
  use constant
  use mod_twoLLGS
  use mod_rfptwoLLGS
  implicit none
  type twoLLGSdynamics
    type(rfptwoLLGS) r
    double precision,allocatable::time_twoLLGS(:),mvec_twoLLGS(:,:,:),mvec(:,:),nvec(:,:)
    double precision time_exc,jc,width,torque
    logical lout
    character(80) filename
  contains
    procedure::init=>twoLLGSdynamics_init
    procedure::readparam=>twoLLGSdynamics_readparam
    procedure::calcparam=>twoLLGSdynamics_calcparam
    procedure::setini=>twoLLGSdynamics_setini
    procedure::alloc=>twoLLGSdynamics_alloc
    procedure::dealloc=>twoLLGSdynamics_dealloc
    procedure::calc=>twoLLGSdynamics_calc
  end type
  type(twoLLGSdynamics) type_twoLLGSdynamics
contains

  subroutine twoLLGSdynamics_init(s)
    use constant
    implicit none
    class(twoLLGSdynamics) s

    s%filename='twoLLGS.txt'

    call type_twoLLGS%init()

    type_twoLLGS%dt=1.0d0
    type_twoLLGS%N=1000

    type_twoLLGS%Hext=0.2d0
    type_twoLLGS%theta_h=90.0d0*pi/180.0d0
    type_twoLLGS%phi_h=90.0d0*pi/180.0d0

    type_twoLLGS%alpha=0.01
    type_twoLLGS%gfac=2.1d0
    type_twoLLGS%gamma=type_twoLLGS%gfac*myuB/hbar*1.0d-12
    type_twoLLGS%gamma_p=type_twoLLGS%gamma/(1.0d0+type_twoLLGS%alpha**2.0d0)
    
    type_twoLLGS%Ms=1.1d6
    type_twoLLGS%Hkeff=-1.0d0*myu0*type_twoLLGS%Ms
    type_twoLLGS%dF=5.0d-9
    type_twoLLGS%Jex=-1.7d-2
    type_twoLLGS%sigma1(1)=0.0d0
    type_twoLLGS%sigma1(2)=0.0d0
    type_twoLLGS%sigma1(3)=1.0d0

    type_twoLLGS%mx(1)=1.0d0
    type_twoLLGS%my(1)=0.0d0
    type_twoLLGS%mz(1)=0.0d0
    type_twoLLGS%mx(2)=-1.0d0
    type_twoLLGS%my(2)=0.0d0
    type_twoLLGS%mz(2)=0.0d0


    s%width=20.0d-9
    s%jc=(1000.0d-6)/(pi*s%width*s%width)
    s%torque=0.5d0*hbar/e/type_twoLLGS%Ms(1)/type_twoLLGS%dF(1)*s%jc
    s%time_exc=10.0d0

    return
  end subroutine

  subroutine twoLLGSdynamics_readparam(s)
    implicit none
    class(twoLLGSdynamics) s

      call s%r%init()
      call s%r%read()
      call s%r%ana()
      call s%r%print()
      call s%r%close()

      type_twoLLGS%dt = s%r%dt
      type_twoLLGS%N = s%r%n
      type_twoLLGS%Neq = s%r%neq
      type_twoLLGS%Hext = s%r%Hext
      type_twoLLGS%theta_h = s%r%theta_h
      type_twoLLGS%phi_h = s%r%phi_h
      type_twoLLGS%alpha = s%r%alpha
      type_twoLLGS%Ms = s%r%Ms
      type_twoLLGS%Ku = s%r%Ku
      type_twoLLGS%dF = s%r%dF
      type_twoLLGS%Jex = s%r%Jex
      
      type_twoLLGS%sigma1 = s%r%sigma
      s%jc = s%r%jc
      s%time_exc = s%r%time_exc 
      s%filename = s%r%cfilename

    return
  end subroutine

  subroutine twoLLGSdynamics_calcparam(s)
    implicit none
    class(twoLLGSdynamics) s

    s%torque=0.5d0*hbar/e/type_twoLLGS%Ms(1)/type_twoLLGS%dF(1)*s%jc
    call type_twoLLGS%calcparam()

    return
  end subroutine

  subroutine twoLLGSdynamics_setini(s)
    implicit none
    class(twoLLGSdynamics) s

      type_twoLLGS%mx = s%r%mx
      type_twoLLGS%my = s%r%my
      type_twoLLGS%mz = s%r%mz

    return
  end subroutine

  subroutine twoLLGSdynamics_alloc(s)
    implicit none
    class(twoLLGSdynamics) s

    allocate(s%time_twoLLGS(type_twoLLGS%N),s%mvec_twoLLGS(2,3,type_twoLLGS%N),s%mvec(3,type_twoLLGS%N),s%nvec(3,type_twoLLGS%N))

    return
  end subroutine

  subroutine twoLLGSdynamics_dealloc(s)
    implicit none
    class(twoLLGSdynamics) s

    deallocate(s%time_twoLLGS,s%mvec_twoLLGS,s%mvec,s%nvec)

    return
  end subroutine

  subroutine twoLLGSdynamics_calc(s)
    implicit none
    class(twoLLGSdynamics) s
    integer i,j

    if(s%lout .eqv. .true.)then
      write(*,'(A,4e12.4)')'initial magnetization 1 = ',type_twoLLGS%mx(1),type_twoLLGS%my(1),type_twoLLGS%mz(1)
      write(*,'(A,4e12.4)')'initial magnetization 2 = ',type_twoLLGS%mx(2),type_twoLLGS%my(2),type_twoLLGS%mz(2)
    end if

    type_twoLLGS%hpx=0.0d0
    type_twoLLGS%hpy=0.0d0
    type_twoLLGS%hpz=0.0d0
    !write(*,'(6e12.4)')(s%mx(j),s%my(j),s%mz(j),j=1,2)
    open(10,file=trim(s%filename),status='replace')
    do i=1,type_twoLLGS%N

      type_twoLLGS%time=type_twoLLGS%dt*dble(i)

      call type_twoLLGS%calctimestep()
      !write(*,'(6e12.4)')(s%mx(j),s%my(j),s%mz(j),j=1,2)

      if(isnan(type_twoLLGS%mx(1)) .eqv. .true.)then
      write(*,'(A,7e12.4)')'Nan is found',type_twoLLGS%time,(type_twoLLGS%mx(j),type_twoLLGS%my(j),type_twoLLGS%mz(j),j=1,2)
      exit
      end if

      if(type_twoLLGS%time>s%time_exc)then
        type_twoLLGS%torque1(1)=s%torque
      else
        type_twoLLGS%torque1(1)=0.0d0
      end if

      write(10,'(8e16.5)')type_twoLLGS%time,type_twoLLGS%torque1(1),(type_twoLLGS%mx(j),type_twoLLGS%my(j),type_twoLLGS%mz(j),j=1,2)

    end do
    close(10)

    if(s%lout .eqv. .true.)then
      write(*,'(A,4e12.4)')'equilibrium magnetization 1 = ',type_twoLLGS%mx(1),type_twoLLGS%my(1),type_twoLLGS%mz(1)
      write(*,'(A,4e12.4)')'equilibrium magnetization 2 = ',type_twoLLGS%mx(2),type_twoLLGS%my(2),type_twoLLGS%mz(2)
    end if

    return
  end subroutine

end module


  