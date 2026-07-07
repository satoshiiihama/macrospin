module mod_mLLGS !module of LLG+spin-transfer torque for multilayer macrospins coupled by Jex and DMI
  implicit none
  type mLLGS
    integer iout,N,Neq,nl
    double precision time,dt,Hext,theta_h,phi_h
    double precision,allocatable::hpx(:),hpy(:),hpz(:)
    double precision,allocatable::Hk(:),Hkeff(:),Ms(:),dF(:),Ku(:),Jex(:),Jex2(:),DMI(:)
    double precision,allocatable::Hex(:),Hex2(:),Hdmi(:),Kuni(:),Huni(:)
    double precision,allocatable::dir_D(:,:),dir_uni(:,:)
    double precision alpha,gamma,gamma_p,gfac
    double precision,allocatable::mx(:),my(:),mz(:),Hx(:),Hy(:),Hz(:),dmx(:),dmy(:),dmz(:),dm_norm(:),norm(:),mvec(:)
    double precision,allocatable::torque(:),sigma(:,:),taux(:),tauy(:),tauz(:)
    logical lout
  contains
    procedure::init=>mLLGS_init
    procedure::alloc=>mLLGS_alloc
    procedure::dealloc=>mLLGS_dealloc
    procedure::calcparam=>mLLGS_calcparam
    procedure::calcequiv=>mLLGS_calcequiv
    procedure::calcsum=>mLLGS_calcsum
    procedure::calctimestep=>mLLGS_calctimestep
    procedure::calcRK=>mLLGS_calcRK
    procedure::calceachstep=>mLLGS_calceachstep
    procedure::effectivefield=>mLLGS_effectivefield
  end type
  type(mLLGS) type_mLLGS

  contains

    subroutine mLLGS_init(s)
    use constant
    implicit none
    class(mLLGS) s

    s%nl=2

    s%dt=1.0d0
    s%N=50000
    s%Neq=50000

    s%iout=50
    s%gfac=2.1d0
    s%gamma=s%gfac*myuB/hbar*1.0d-12
    s%alpha=0.01
    s%gamma_p=s%gamma/(1.0d0+s%alpha**2.0d0)

    s%Hext=0.2d0 !T
    s%theta_h=0.5d0*pi !rad
    s%phi_h=0.0d0 !rad

    s%lout=.true.

  end subroutine

  subroutine mLLGS_alloc(s)
    use constant
    implicit none
    class(mLLGS) s

    allocate(s%hpx(s%nl),s%hpy(s%nl),s%hpz(s%nl))
    allocate(s%Hk(s%nl),s%Hkeff(s%nl),s%Huni(s%nl),s%Ms(s%nl),s%dF(s%nl),s%Ku(s%nl),s%Kuni(s%nl),s%dir_uni(3,s%nl),s%Jex(s%nl-1),&
    s%Jex2(s%nl-1),s%DMI(s%nl-1),s%dir_D(3,s%nl-1),s%Hex(s%nl-1),s%Hex2(s%nl-1),s%Hdmi(s%nl-1))
    allocate(s%mx(s%nl),s%my(s%nl),s%mz(s%nl),s%Hx(s%nl),s%Hy(s%nl),s%Hz(s%nl),&
    s%dmx(s%nl),s%dmy(s%nl),s%dmz(s%nl),s%dm_norm(s%nl),s%norm(s%nl),s%mvec(3))
    allocate(s%torque(s%nl),s%sigma(3,s%nl),s%taux(s%nl),s%tauy(s%nl),s%tauz(s%nl))

    s%sigma=1.0d0
    s%sigma=0.0d0
    s%sigma=0.0d0
    s%torque=0.0d0

    s%Ms=1.0d6
    s%Ku=0.0d0 !J/m2
    s%Hkeff=2.0d0*s%Ku/s%Ms-myu0*s%Ms

    s%dF=1.0d-9
    s%Jex=-0.9d-3
    s%Jex2=0.0d0
    s%DMI=0.1d-3

    s%mx=0.0d0
    s%my=0.0d0
    s%mz=1.0d0

  end subroutine

  subroutine mLLGS_dealloc(s)
    implicit none
    class(mLLGS) s

    deallocate(s%hpx,s%hpy,s%hpz)
    deallocate(s%Hk,s%Hkeff,s%Huni,s%Ms,s%dF,s%Ku,s%Kuni,s%dir_uni,s%Jex,s%Jex2,s%DMI,s%dir_D,s%Hex,s%Hex2,s%Hdmi)
    deallocate(s%mx,s%my,s%mz,s%Hx,s%Hy,s%Hz,s%dmx,s%dmy,s%dmz,s%dm_norm,s%norm,s%mvec)
    deallocate(s%torque,s%sigma,s%taux,s%tauy,s%tauz)

  end subroutine

  subroutine mLLGS_calcparam(s)
    use constant
    implicit none
    class(mLLGS) s

    s%gamma_p=s%gamma/(1.0d0+s%alpha**2.0d0)
    s%Hkeff=2.0d0*s%Ku/s%Ms-myu0*s%Ms
    s%Huni=2.0d0*s%Kuni/s%Ms
    s%Hex=s%Jex/s%Ms/s%dF
    s%Hex2=2.0d0*s%Jex2/s%Ms/s%dF
    s%Hdmi=s%DMI/s%Ms/s%dF
    
    return
  end subroutine

  subroutine mLLGS_calcequiv(s)
    implicit none
    class(mLLGS) s
    integer i,j,k

    if(s%lout .eqv. .true.)then
      do k=1,s%nl
        write(*,'(A,i2,4e12.4)')'initial magnetization = ',k,s%mx(k),s%my(k),s%mz(k),s%norm(k)
      end do
    end if

    s%hpx=0.0d0
    s%hpy=0.0d0
    s%hpz=0.0d0
    !write(*,'(6e12.4)')(s%mx(j),s%my(j),s%mz(j),j=1,2)
    do i=1,s%Neq

      s%time=s%dt*dble(i)

      call s%calctimestep()
      !write(*,'(6e12.4)')(s%mx(j),s%my(j),s%mz(j),j=1,2)

      if(isnan(s%mx(1)) .eqv. .true.)then
      write(*,'(A,100e12.4)')'Nan is found',s%time,(s%mx(j),s%my(j),s%mz(j),j=1,s%nl)
      exit
      end if
               
      if(maxval(s%dm_norm)<1.0d-8)then
        if(s%lout .eqv. .true.)then
          write(*,*)'calculation is converged when i = ',i
        end if
        exit
      end if

      if(i==s%Neq)then
        write(*,*)'calculation is not converged'
      end if

    end do

    if(s%lout .eqv. .true.)then
      do k=1,s%nl
        write(*,'(A,i2,4e12.4)')'equilibrium magnetization 1 = ',k,s%mx(k),s%my(k),s%mz(k),s%norm(k)
      end do
    end if

    return
  end subroutine

  subroutine mLLGS_calcsum(s)
    implicit none
    class(mLLGS) s
    integer i
    
    s%mvec=0.0d0
    do i=1,s%nl
      s%mvec(1)=s%mvec(1)+s%mx(i)
      s%mvec(2)=s%mvec(2)+s%my(i)
      s%mvec(3)=s%mvec(3)+s%mz(i)  
    end do
    s%mvec=s%mvec/dble(s%nl)
    
  end subroutine

  subroutine mLLGS_calctimestep(s)
    implicit none
    class(mLLGS) s
    integer i
    double precision dm1(s%nl),dm2(s%nl),dm3(s%nl)
    
    call s%calcRK(s%mx,s%my,s%mz,dm1,dm2,dm3)
    s%mx=s%mx+dm1
    s%my=s%my+dm2
    s%mz=s%mz+dm3
    
    do i=1,s%nl
      s%norm(i)=sqrt(s%mx(i)*s%mx(i)+s%my(i)*s%my(i)+s%mz(i)*s%mz(i))
    end do

    s%mx=s%mx/s%norm
    s%my=s%my/s%norm
    s%mz=s%mz/s%norm

    s%dm_norm=sqrt(dm1**2.0d0+dm2**2.0d0+dm3**2.0d0)

  end subroutine

  subroutine mLLGS_calcRK(s,m1,m2,m3,dm1,dm2,dm3)
    implicit none
    class(mLLGS) s
    integer i
    double precision m1(s%nl),m2(s%nl),m3(s%nl),dm1(s%nl),dm2(s%nl),dm3(s%nl),mp1(s%nl),mp2(s%nl),mp3(s%nl)
    double precision l1x(s%nl),l1y(s%nl),l1z(s%nl)
    double precision l2x(s%nl),l2y(s%nl),l2z(s%nl)
    double precision l3x(s%nl),l3y(s%nl),l3z(s%nl)
    double precision l4x(s%nl),l4y(s%nl),l4z(s%nl)

    
    call s%calceachstep(m1,m2,m3,l1x,l1y,l1z)
    call s%calceachstep(m1+l1x*0.5d0,m2+l1y*0.5d0,m3+l1z*0.5d0,l2x,l2y,l2z)
    call s%calceachstep(m1+l2x*0.5d0,m2+l2y*0.5d0,m3+l2z*0.5d0,l3x,l3y,l3z)
    call s%calceachstep(m1+l3x,m2+l3y,m3+l3z,l4x,l4y,l4z)

    do i=1,s%nl
      dm1(i)=(l1x(i)+2.0d0*l2x(i)+2.0d0*l3x(i)+l4x(i))/6.0d0
      dm2(i)=(l1y(i)+2.0d0*l2y(i)+2.0d0*l3y(i)+l4y(i))/6.0d0
      dm3(i)=(l1z(i)+2.0d0*l2z(i)+2.0d0*l3z(i)+l4z(i))/6.0d0
    end do

  end subroutine

  subroutine mLLGS_calceachstep(s,m1,m2,m3,dm1,dm2,dm3)
    use func_LLGS
    implicit none
    class(mLLGS) s
    integer i,j
    double precision m1(s%nl),m2(s%nl),m3(s%nl),dm1(s%nl),dm2(s%nl),dm3(s%nl),H1(s%nl),H2(s%nl),H3(s%nl)
    double precision dir_uni(3),dir_D1(3),dir_D2(3)

    do i=1,s%nl

      if(i==1)then
        do j=1,3
          dir_D1(j)=s%dir_D(j,i)
          dir_D2(j)=s%dir_D(j,i)
          dir_uni(j)=s%dir_uni(j,i)
        end do
        call s%effectivefield(m1(i),m2(i),m3(i),H1(i),H2(i),H3(i),s%Hkeff(i),s%Huni(i),dir_uni,&
        0.0d0,0.0d0,0.0d0,dir_D1,m1(i),m2(i),m3(i),&
        s%Hex(i),s%Hex2(i),s%Hdmi(i),dir_D2,m1(i+1),m2(i+1),m3(i+1),s%hpx(i),s%hpy(i),s%hpz(i))
      elseif(i==s%nl)then
        do j=1,3
          dir_D1(j)=s%dir_D(j,i-1)
          dir_D2(j)=s%dir_D(j,i-1)
          dir_uni(j)=s%dir_uni(j,i)
        end do
        call s%effectivefield(m1(i),m2(i),m3(i),H1(i),H2(i),H3(i),s%Hkeff(i),s%Huni(i),dir_uni,&
        s%Hex(i-1),s%Hex2(i-1),s%Hdmi(i-1),dir_D1,m1(i-1),m2(i-1),m3(i-1),&
        0.0d0,0.0d0,0.0d0,dir_D2,m1(i),m2(i),m3(i),s%hpx(i),s%hpy(i),s%hpz(i))
      else
        do j=1,3
          dir_D1(j)=s%dir_D(j,i-1)
          dir_D2(j)=s%dir_D(j,i)
          dir_uni(j)=s%dir_uni(j,i)
        end do
        call s%effectivefield(m1(i),m2(i),m3(i),H1(i),H2(i),H3(i),s%Hkeff(i),s%Huni(i),dir_uni,&
        s%Hex(i-1),s%Hex2(i-1),s%Hdmi(i-1),dir_D1,m1(i-1),m2(i-1),m3(i-1),&
        s%Hex(i),s%Hex2(i),s%Hdmi(i),dir_D2,m1(i+1),m2(i+1),m3(i+1),s%hpx(i),s%hpy(i),s%hpz(i))
      end if

    end do
    
    do i=1,s%nl
      call LLGS_calc_tau(s%torque(i),s%sigma,m1(i),m2(i),m3(i),s%taux(i),s%tauy(i),s%tauz(i))
      
      dm1(i)=s%dt*LLGS_fx(m1(i),m2(i),m3(i),H1(i),H2(i),H3(i),s%alpha,s%gamma_p,s%taux(i),s%tauy(i),s%tauz(i))
      dm2(i)=s%dt*LLGS_fy(m1(i),m2(i),m3(i),H1(i),H2(i),H3(i),s%alpha,s%gamma_p,s%taux(i),s%tauy(i),s%tauz(i))
      dm3(i)=s%dt*LLGS_fz(m1(i),m2(i),m3(i),H1(i),H2(i),H3(i),s%alpha,s%gamma_p,s%taux(i),s%tauy(i),s%tauz(i))

    end do

  end subroutine

  subroutine mLLGS_effectivefield(s,mx,my,mz,Hx,Hy,Hz,Hkeff,Huni,dir_uni,Hex1,Hex21,Hdmi1,dir1,m1x,m1y,m1z,Hex2,Hex22,Hdmi2,dir2,&
    m2x,m2y,m2z,hpx,hpy,hpz)
    implicit none
    class(mLLGS) s
    double precision mx,my,mz,Hx,Hy,Hz,Hkeff,Huni,dir_uni(3),Hex1,Hex21,Hdmi1,dir1(3),m1x,m1y,m1z
    double precision Hex2,Hex22,Hdmi2,dir2(3),m2x,m2y,m2z,mdu,midmj1,midmj2,hpx,hpy,hpz

    mdu=mx*dir_uni(1)+my*dir_uni(2)+mz*dir_uni(3)
    midmj1=mx*m1x+my*m1y+mz*m1z
    midmj2=mx*m2x+my*m2y+mz*m2z

    Hx=s%Hext*sin(s%theta_h)*cos(s%phi_h)+Huni*mdu*dir_uni(1)+Hex1*m1x+Hex21*midmj1*m1x+Hex2*m2x+Hex22*midmj2*m2x&
    +Hdmi1*(dir1(2)*m1z-dir1(3)*m1y)-Hdmi2*(dir2(2)*m2z-dir2(3)*m2y)+hpx
    Hy=s%Hext*sin(s%theta_h)*sin(s%phi_h)+Huni*mdu*dir_uni(2)+Hex1*m1y+Hex21*midmj1*m1y+Hex2*m2y+Hex22*midmj2*m2y&
    +Hdmi1*(dir1(3)*m1x-dir1(1)*m1z)-Hdmi2*(dir2(3)*m2x-dir2(1)*m2z)+hpy
    Hz=s%Hext*cos(s%theta_h)+Hkeff*mz+Huni*mdu*dir_uni(3)+Hex1*m1z+Hex21*midmj1*m1z+Hex2*m2z+Hex22*midmj2*m2z&
    +Hdmi1*(dir1(1)*m1y-dir1(2)*m1x)-Hdmi2*(dir2(1)*m2y-dir2(2)*m2x)+hpz

  end subroutine

end module
