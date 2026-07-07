module mod_twoLLGS !module of LLG+spin-transfer torque for two macrospins coupled by Jex
  implicit none
  type twoLLGS
    integer iout,N,Neq
    double precision time,dt,Hext,theta_h,phi_h,hpx(2),hpy(2),hpz(2)
    double precision alpha(2),gamma(2),gamma_p(2),gfac(2),Hk(2),Hkeff(2),Ms(2),dF(2),Ku(2),Jex,Hex(2)
    double precision mx(2),my(2),mz(2),Hx(2),Hy(2),Hz(2),dmx(2),dmy(2),dmz(2),dm_norm(2),norm(2),mvec(3),nvec(3)
    double precision torque1(2),sigma1(3),taux1(2),tauy1(2),tauz1(2),torque2(2),sigma2(3),taux2(2),tauy2(2),tauz2(2)
    double precision ener_Z,ener_K,ener_J,ener_tot
    logical lout
  contains
    procedure::init=>twoLLGS_init
    procedure::calcparam=>twoLLGS_calcparam
    procedure::calcequiv=>twoLLGS_calcequiv
    procedure::calctimestep=>twoLLGS_calctimestep
    procedure::calcRK=>twoLLGS_calcRK
    procedure::calceachstep=>twoLLGS_calceachstep
    procedure::effectivefield=>twoLLGS_effectivefield
    procedure::calcenergy=>twoLLGS_calcenergy
  end type
  type(twoLLGS) type_twoLLGS

  contains

    subroutine twoLLGS_init(s)
    use constant
    implicit none
    class(twoLLGS) s

    s%dt=1.0d0
    s%N=50000
    s%Neq=50000

    s%iout=50
    s%gfac=2.1d0
    s%gamma=s%gfac*myuB/hbar*1.0d-12
    s%alpha=0.01
    s%gamma_p=s%gamma/(1.0d0+s%alpha**2.0d0)

    s%sigma1(1)=1.0d0
    s%sigma1(2)=0.0d0
    s%sigma1(3)=0.0d0
    s%torque1=0.0d0
    s%sigma2(1)=0.0d0
    s%sigma2(2)=0.0d0
    s%sigma2(3)=1.0d0
    s%torque2=0.0d0

    s%Hext=0.2d0 !T
    s%Ms=1.0d6
    s%Ku=0.0d0 !J/m2
    s%Hkeff=2.0d0*s%Ku/s%Ms-myu0*s%Ms

    s%dF=1.0d-9
    s%Jex=-0.9d-3

    s%theta_h=0.5d0*pi !rad
    s%phi_h=0.0d0 !rad

    s%mx(1)=0.0d0
    s%my(1)=0.0d0
    s%mz(1)=1.0d0
    s%mx(2)=0.0d0
    s%my(2)=0.0d0
    s%mz(2)=-1.0d0

    s%lout=.true.

    return
  end subroutine

  subroutine twoLLGS_calcparam(s)
    use constant
    implicit none
    class(twoLLGS) s

    s%gamma_p=s%gamma/(1.0d0+s%alpha**2.0d0)
    s%Hkeff=2.0d0*s%Ku/s%Ms-myu0*s%Ms
    s%Hex(1)=s%Jex/s%Ms(1)/s%dF(1)
    s%Hex(2)=s%Jex/s%Ms(2)/s%dF(2)
    
    return
  end subroutine

  subroutine twoLLGS_calcequiv(s)
    implicit none
    class(twoLLGS) s
    integer i,j

    if(s%lout .eqv. .true.)then
      write(*,'(A,4e12.4)')'initial magnetization 1 = ',s%mx(1),s%my(1),s%mz(1),s%norm(1)
      write(*,'(A,4e12.4)')'initial magnetization 2 = ',s%mx(2),s%my(2),s%mz(2),s%norm(2)
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
      write(*,'(A,7e12.4)')'Nan is found',s%time,(s%mx(j),s%my(j),s%mz(j),j=1,2)
      exit
      end if
               
      if(s%dm_norm(1)<1.0d-8 .and. s%dm_norm(2)<1.0d-8)then
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
      write(*,'(A,4e12.4)')'equilibrium magnetization 1 = ',s%mx(1),s%my(1),s%mz(1),s%norm(1)
      write(*,'(A,4e12.4)')'equilibrium magnetization 2 = ',s%mx(2),s%my(2),s%mz(2),s%norm(2)
    end if

    return
  end subroutine

  subroutine twoLLGS_calctimestep(s)
    implicit none
    class(twoLLGS) s
    integer i
    double precision dm1(2),dm2(2),dm3(2)
    
    call s%calcRK(s%mx,s%my,s%mz,dm1,dm2,dm3)
    s%mx=s%mx+dm1
    s%my=s%my+dm2
    s%mz=s%mz+dm3
    
    do i=1,2
      s%norm(i)=sqrt(s%mx(i)*s%mx(i)+s%my(i)*s%my(i)+s%mz(i)*s%mz(i))
    end do

    s%mx=s%mx/s%norm
    s%my=s%my/s%norm
    s%mz=s%mz/s%norm

    s%dm_norm(1)=sqrt(dm1(1)**2.0d0+dm2(1)**2.0d0+dm3(1)**2.0d0)
    s%dm_norm(2)=sqrt(dm1(2)**2.0d0+dm2(2)**2.0d0+dm3(2)**2.0d0)

    s%mvec(1)=0.5d0*(s%mx(1)+s%mx(2))
    s%mvec(2)=0.5d0*(s%my(1)+s%my(2))
    s%mvec(3)=0.5d0*(s%mz(1)+s%mz(2))

    s%nvec(1)=0.5d0*(s%mx(1)-s%mx(2))
    s%nvec(2)=0.5d0*(s%my(1)-s%my(2))
    s%nvec(3)=0.5d0*(s%mz(1)-s%mz(2))

    return
  end subroutine

  subroutine twoLLGS_calcRK(s,m1,m2,m3,dm1,dm2,dm3)
    implicit none
    class(twoLLGS) s
    integer i
    double precision m1(2),m2(2),m3(2),dm1(2),dm2(2),dm3(2),mp1(2),mp2(2),mp3(2)
    double precision l1x(2),l1y(2),l1z(2)
    double precision l2x(2),l2y(2),l2z(2)
    double precision l3x(2),l3y(2),l3z(2)
    double precision l4x(2),l4y(2),l4z(2)

    call s%calceachstep(m1,m2,m3,l1x,l1y,l1z)
    call s%calceachstep(m1+l1x*0.5d0,m2+l1y*0.5d0,m3+l1z*0.5d0,l2x,l2y,l2z)
    call s%calceachstep(m1+l2x*0.5d0,m2+l2y*0.5d0,m3+l2z*0.5d0,l3x,l3y,l3z)
    call s%calceachstep(m1+l3x,m2+l3y,m3+l3z,l4x,l4y,l4z)

    dm1(1)=(l1x(1)+2.0d0*l2x(1)+2.0d0*l3x(1)+l4x(1))/6.0d0
    dm2(1)=(l1y(1)+2.0d0*l2y(1)+2.0d0*l3y(1)+l4y(1))/6.0d0
    dm3(1)=(l1z(1)+2.0d0*l2z(1)+2.0d0*l3z(1)+l4z(1))/6.0d0

    dm1(2)=(l1x(2)+2.0d0*l2x(2)+2.0d0*l3x(2)+l4x(2))/6.0d0
    dm2(2)=(l1y(2)+2.0d0*l2y(2)+2.0d0*l3y(2)+l4y(2))/6.0d0
    dm3(2)=(l1z(2)+2.0d0*l2z(2)+2.0d0*l3z(2)+l4z(2))/6.0d0

    return
  end subroutine

  subroutine twoLLGS_calceachstep(s,m1,m2,m3,dm1,dm2,dm3)
    use func_LLGS
    implicit none
    class(twoLLGS) s
    double precision m1(2),m2(2),m3(2),dm1(2),dm2(2),dm3(2),H1(2),H2(2),H3(2)

    call s%effectivefield(m1,m2,m3,H1,H2,H3,s%Hext,s%theta_h,s%phi_h,s%Hkeff,s%Hex,s%hpx,s%hpy,s%hpz)

    call LLGS_calc_tau(s%torque1(1),s%sigma1,m1(1),m2(1),m3(1),s%taux1(1),s%tauy1(1),s%tauz1(1))
    call LLGS_calc_tau(s%torque2(1),s%sigma2,m1(1),m2(1),m3(1),s%taux2(1),s%tauy2(1),s%tauz2(1))
    dm1(1)=s%dt*LLGS_fx(m1(1),m2(1),m3(1),H1(1),H2(1),H3(1),s%alpha(1),s%gamma_p(1),&
    s%taux1(1)+s%taux2(1),s%tauy1(1)+s%tauy2(1),s%tauz1(1)+s%tauz2(1))
    dm2(1)=s%dt*LLGS_fy(m1(1),m2(1),m3(1),H1(1),H2(1),H3(1),s%alpha(1),s%gamma_p(1),&
    s%taux1(1)+s%taux2(1),s%tauy1(1)+s%tauy2(1),s%tauz1(1)+s%tauz2(1))
    dm3(1)=s%dt*LLGS_fz(m1(1),m2(1),m3(1),H1(1),H2(1),H3(1),s%alpha(1),s%gamma_p(1),&
    s%taux1(1)+s%taux2(1),s%tauy1(1)+s%tauy2(1),s%tauz1(1)+s%tauz2(1))

    call LLGS_calc_tau(s%torque1(2),s%sigma1,m1(2),m2(2),m3(2),s%taux1(2),s%tauy1(2),s%tauz1(2))
    call LLGS_calc_tau(s%torque2(2),s%sigma2,m1(2),m2(2),m3(2),s%taux2(2),s%tauy2(2),s%tauz2(2))
    dm1(2)=s%dt*LLGS_fx(m1(2),m2(2),m3(2),H1(2),H2(2),H3(2),s%alpha(2),s%gamma_p(2),&
    s%taux1(2)+s%taux2(2),s%tauy1(2)+s%tauy2(2),s%tauz1(2)+s%tauz2(2))
    dm2(2)=s%dt*LLGS_fy(m1(2),m2(2),m3(2),H1(2),H2(2),H3(2),s%alpha(2),s%gamma_p(2),&
    s%taux1(2)+s%taux2(2),s%tauy1(2)+s%tauy2(2),s%tauz1(2)+s%tauz2(2))
    dm3(2)=s%dt*LLGS_fz(m1(2),m2(2),m3(2),H1(2),H2(2),H3(2),s%alpha(2),s%gamma_p(2),&
    s%taux1(2)+s%taux2(2),s%tauy1(2)+s%tauy2(2),s%tauz1(2)+s%tauz2(2))

    return
  end subroutine

  subroutine twoLLGS_effectivefield(s,mx,my,mz,Hx,Hy,Hz,Hext,theta_h,phi_h,Hkeff,Hex,hpx,hpy,hpz)
    implicit none
    class(twoLLGS) s
    double precision mx(2),my(2),mz(2),Hx(2),Hy(2),Hz(2),Hext,theta_h,phi_h,Hkeff(2),Hex(2),hpx(2),hpy(2),hpz(2)

    Hx(1)=Hext*sin(theta_h)*cos(phi_h)+Hex(1)*mx(2)+hpx(1)
    Hy(1)=Hext*sin(theta_h)*sin(phi_h)+Hex(1)*my(2)+hpy(1)
    Hz(1)=Hext*cos(theta_h)+Hkeff(1)*mz(1)+Hex(1)*mz(2)+hpz(1)
    Hx(2)=Hext*sin(theta_h)*cos(phi_h)+Hex(2)*mx(1)+hpx(2)
    Hy(2)=Hext*sin(theta_h)*sin(phi_h)+Hex(2)*my(1)+hpy(2)
    Hz(2)=Hext*cos(theta_h)+Hkeff(2)*mz(2)+Hex(2)*mz(1)+hpz(2)

    return
  end subroutine

  subroutine twoLLGS_calcenergy(s,mx,my,mz)
    use constant
    implicit none
    class(twoLLGS) s
    double precision mx(2),my(2),mz(2),Hx,Hy,Hz

    Hx=s%Hext*sin(s%theta_h)*cos(s%phi_h)
    Hy=s%Hext*sin(s%theta_h)*sin(s%phi_h)
    Hz=s%Hext*cos(s%theta_h)

    s%ener_Z=-1.0d0*(s%Ms(1)*s%dF(1)*(mx(1)*Hx+my(1)*Hy+mz(1)*Hz)+s%Ms(2)*s%dF(2)*(mx(2)*Hx+my(2)*Hy+mz(2)*Hz))
    s%ener_K=-0.5d0*(s%Ms(1)*s%dF(1)*s%Hkeff(1)*mz(1)*mz(1)+s%Ms(2)*s%dF(2)*s%Hkeff(2)*mz(2)*mz(2))
    s%ener_J=-1.0d0*s%Jex*(mx(1)*mx(2)+my(1)*my(2)+mz(1)*mz(2))
    s%ener_tot=s%ener_Z+s%ener_K+s%ener_J

    return
  end subroutine
end module
