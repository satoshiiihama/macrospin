  module mod_rfptwoLLGS
    ! parameter list
    implicit none
    type rfptwoLLGS
      integer nparam,nlist
      character(200),allocatable::charlist(:),paramlist(:),paramvalues(:)
      character(200) filename
      integer N,Neq
      double precision dt,time_exc
      double precision Hext,theta_h,phi_h
      double precision alpha,Ms,Ku,dF,Jex
      double precision mx(2),my(2),mz(2),sigma(3),jc
      character(200) cfilename
      character(50) cn,cneq,cdt,cHext,ctheta_h,cphi_h,calpha,cMs,cKu,cdF,cJex,cmx1,cmx2,cmy1,cmy2,cmz1,cmz2
      character(50) csigma,cjc,ctime_exc
    contains
      procedure::init=>rfptwoLLGS_init
      procedure::help=>rfptwoLLGS_help
      procedure::set=>rfptwoLLGS_set
      procedure::read=>rfptwoLLGS_read
      procedure::ana=>rfptwoLLGS_ana
      procedure::print=>rfptwoLLGS_print
      procedure::close=>rfptwoLLGS_close
    end type
  contains

    subroutine rfptwoLLGS_init(s)
      use constant
      implicit none
      class(rfptwoLLGS) s
      character(200) cread

      if( iargc() == 1 )then
      else
        write(*,*)'***Error*** Invalid number of command line arguments, describe parameterfiles or "--help"'
        stop
      end if

      call getarg(1,cread)

      if(trim(cread)=='--help' .or. trim(cread)=='-h')then
        call s%help()
        stop
      else
        s%filename=trim(cread)
      end if

      s%nparam=17
      allocate(s%paramlist(s%nparam),s%paramvalues(s%nparam))
      s%paramlist(1)='dt'
      s%paramlist(2)='N'
      s%paramlist(3)='Neq'
      s%paramlist(4)='Hext'
      s%paramlist(5)='theta_h'
      s%paramlist(6)='phi_h'
      s%paramlist(7)='alpha'
      s%paramlist(8)='Ms'
      s%paramlist(9)='Ku'
      s%paramlist(10)='dF'
      s%paramlist(11)='Jex'
      s%paramlist(12)='m1_ini'
      s%paramlist(13)='m2_ini'
      s%paramlist(14)='sigma'
      s%paramlist(15)='jc'
      s%paramlist(16)='time_exc'
      s%paramlist(17)='filename'

      s%dt=1.0d0
      s%N=1000
      s%Neq=5000
      s%Hext=0.2d0
      s%theta_h=0.5d0*pi
      s%phi_h=0.0d0
      s%alpha=0.01d0
      s%Ms=1.0d6
      s%Ku=0.0d0
      s%dF=1.0d-9
      s%Jex=-1.0d-3
      s%mx(1)=0.0d0
      s%my(1)=1.0d0
      s%mz(1)=0.0d0
      s%mx(2)=0.0d0
      s%my(2)=-1.0d0
      s%mz(2)=0.0d0
      s%sigma(1)=0.0d0
      s%sigma(2)=0.0d0
      s%sigma(3)=1.0d0
      s%jc=1.0d-3/(pi*(20.0d-9)*(20.0d-9))
      s%time_exc=10.0d0

      s%cfilename='test.txt'
      
      return
    end subroutine

    subroutine rfptwoLLGS_help(s)
      implicit none
      class(rfptwoLLGS) s

        write(*,*)' This program is to calculate macrospin dynamics of two coupled magnetization of synthetic antiferromagnet '
        write(*,*)' Parameters for the calculation must be argued after .exe statement'
        write(*,*)' The meaning of each parameters are described below'
        write(*,*)' dt: time-step for Runge-Kutta simulation'
        write(*,*)' N: total time step number'
        write(*,*)' Neq: time step number for equilibration'
        write(*,*)' Hext: external magnetic field strength in unit of [T]'
        write(*,*)' theta_h: external magnetic field angle from film normal in unit of radian'
        write(*,*)' phi_h: external magnetic field in-plane angle from x-axis in unit of radian'
        write(*,*)' alpha: Gilbert damping parameter for both layer'
        write(*,*)' Ms: Saturation magnetization in unit of [A/m]'
        write(*,*)' Ku: Uniaxial magnetic anisotropy constant in unit of [J/m3]'
        write(*,*)' dF: ferromagnetic layer thickness in unit of [m]'
        write(*,*)' Jex: RKKY interlayer coupling strength in unit of [J/m2]'
        write(*,*)' m(1,2)ini: initial magnetization configuration for layer 1 or 2 (x, y, z)'
        write(*,*)' sigma: spin-polarization direction (x, y, z)'
        write(*,*)' jc: current density for spin-transfer torque excitation'
        write(*,*)' time_exc: time for spin-transfer torque excitation in unit of ps'
        write(*,*)' filename: filename for output'
       
      return
    end subroutine

    subroutine rfptwoLLGS_read(s)
      implicit none
      class(rfptwoLLGS) s
      integer i,j,ios,n
      character(80) dummy(20)

      call s%set()

      !write(*,*)'-----Reading a file-----'
      n=1
      open(20,file=trim(s%filename), status='old')
      do i=1,s%nlist
        read(20,'(a)',iostat=ios)s%charlist(n)
        !write(*,*)ios
        if(ios<0) exit
        n=n+1
      end do
      close(20)

      !do i=1,s%nlist
      !  write(*,'(a)')trim(s%charlist(i))
      !end do

      return
    end subroutine

    subroutine rfptwoLLGS_ana(s)
      implicit none
      class(rfptwoLLGS) s
      integer i,j,k,len_k,len_j,n
      character(200) file,file2
      character(50),allocatable::cparam(:),cnum(:)

      !write(*,*)'-----Analyzing file-----'
      n=0
      do i=1,s%nlist
        file=s%charlist(i)
        
        ! searching comment out
        j=index(file,"#")
        if(j==0)then
          file2=file
        else
          file2=file(1:j-1)
        end if

        ! serching equal statement
        k=index(file2,"=")
        len_k=len_trim(file2(1:k-1))

        if(len_k==0)then

        else
          n=n+1
        end if

      end do

      allocate(cparam(n),cnum(n))

      n=0
      do i=1,s%nlist
        file=s%charlist(i)
        
        ! searching comment out
        j=index(file,"#")
        if(j==0)then
          file2=file
        else
          file2=file(1:j-1)
        end if

        ! serching equal statement
        k=index(file2,"=")
        len_j=len_trim(file2)
        len_k=len_trim(file2(1:k-1))

        if(len_k==0)then

        else
          n=n+1
          cparam(n)=trim(file2(1:k-1))
          cnum(n)=trim(file2(k+1:len_j))
        end if

      end do

      !write(*,*)n
      !do i=1,n
      !  write(*,'(A)')cparam(i)
      !end do
      !do i=1,n
      !  write(*,'(A)')cnum(i)
      !end do

      do j=1,s%nparam

        do i=1,n

          if(cparam(i)==s%paramlist(j))then
            s%paramvalues(j)=trim(cnum(i))
          else
          end if
        end do

      end do

      deallocate(cparam,cnum)

      return
    end subroutine

    subroutine rfptwoLLGS_set(s)
      implicit none
      class(rfptwoLLGS) s
      integer i,n,nmax,ios
      character(200) dum

      nmax=1000
      write(*,*)'filename = ',trim(s%filename)
      open(20,file=trim(s%filename), status='old')
      n=0
      do i=1,nmax
        !read(20,*,end=999)dum
        read(20,'(a)',iostat=ios)dum
        !write(*,*)trim(dum)
        if(ios<0) exit
        n=n+1
      end do
999   close(20)
      s%nlist=n
      write(*,'(A,i8)')'Number of data = ',n

      !stop
      allocate(s%charlist(n))
      !write(*,*)'Allocating parameters in rfptwoLLGS'

      return
    end subroutine

    subroutine rfptwoLLGS_print(s)
      use constant
      implicit none
      class(rfptwoLLGS) s
      integer i
      character(20) cm1,cm2,csigma
      !integer ndis,ntra,ntes,nset,ntype,np,nv,ncomp,nccomp,ndim,ndelay,nzero,nloop,nnorm,nzsc,nprint,npnum
      !character(20) cdate,ci,cj 

        !do i=1,s%nparam
        !  write(*,*)s%paramvalues(i)
        !end do
      
        read(s%paramvalues(1),*)s%dt
        read(s%paramvalues(2),*)s%n
        read(s%paramvalues(3),*)s%neq
        read(s%paramvalues(4),*)s%Hext
        read(s%paramvalues(5),*)s%theta_h
        read(s%paramvalues(6),*)s%phi_h
        read(s%paramvalues(7),*)s%alpha
        read(s%paramvalues(8),*)s%Ms
        read(s%paramvalues(9),*)s%Ku
        read(s%paramvalues(10),*)s%dF
        read(s%paramvalues(11),*)s%Jex
        read(s%paramvalues(15),*)s%jc
        read(s%paramvalues(16),*)s%time_exc
        i=index(trim(s%paramvalues(17))," ",back=.true.)
        s%cfilename=trim(s%paramvalues(17)(1+i:len_trim(s%paramvalues(17))))
        
        s%theta_h = s%theta_h * pi / 180.0d0
        s%phi_h = s%phi_h * pi / 180.0d0

        i=index(trim(s%paramvalues(12))," ",back=.true.)
        cm1=trim(s%paramvalues(12)(1+i:len_trim(s%paramvalues(12))))
        i=index(trim(s%paramvalues(13))," ",back=.true.)
        cm2=trim(s%paramvalues(13)(1+i:len_trim(s%paramvalues(13))))
        i=index(trim(s%paramvalues(14))," ",back=.true.)
        csigma=trim(s%paramvalues(14)(1+i:len_trim(s%paramvalues(14))))
        if(cm1=='x')then
          s%mx(1)=1.0d0
          s%my(1)=0.0d0
          s%mz(1)=0.0d0
        elseif(cm1=='-x')then
          s%mx(1)=-1.0d0
          s%my(1)=0.0d0
          s%mz(1)=0.0d0
        elseif(cm1=='y')then
          s%mx(1)=0.0d0
          s%my(1)=1.0d0
          s%mz(1)=0.0d0
        elseif(cm1=='-y')then
          s%mx(1)=0.0d0
          s%my(1)=-1.0d0
          s%mz(1)=0.0d0
        elseif(cm1=='z')then
          s%mx(1)=0.0d0
          s%my(1)=0.0d0
          s%mz(1)=1.0d0
        elseif(cm1=='-z')then
          s%mx(1)=0.0d0
          s%my(1)=0.0d0
          s%mz(1)=-1.0d0
        else
        end if

        if(cm2=='x')then
          s%mx(2)=1.0d0
          s%my(2)=0.0d0
          s%mz(2)=0.0d0
        elseif(cm2=='-x')then
          s%mx(2)=-1.0d0
          s%my(2)=0.0d0
          s%mz(2)=0.0d0
        elseif(cm2=='y')then
          s%mx(2)=0.0d0
          s%my(2)=1.0d0
          s%mz(2)=0.0d0
        elseif(cm2=='-y')then
          s%mx(2)=0.0d0
          s%my(2)=-1.0d0
          s%mz(2)=0.0d0
        elseif(cm2=='z')then
          s%mx(2)=0.0d0
          s%my(2)=0.0d0
          s%mz(2)=1.0d0
        elseif(cm2=='-z')then
          s%mx(2)=0.0d0
          s%my(2)=0.0d0
          s%mz(2)=-1.0d0
        else
        end if

        if(csigma=='x')then
          s%sigma(1)=1.0d0
          s%sigma(2)=0.0d0
          s%sigma(3)=0.0d0
        elseif(csigma=='-x')then
          s%sigma(1)=-1.0d0
          s%sigma(2)=0.0d0
          s%sigma(3)=0.0d0
        elseif(csigma=='y')then
          s%sigma(1)=0.0d0
          s%sigma(2)=1.0d0
          s%sigma(3)=0.0d0
        elseif(csigma=='-y')then
          s%sigma(1)=0.0d0
          s%sigma(2)=-1.0d0
          s%sigma(3)=0.0d0
        elseif(csigma=='z')then
          s%sigma(1)=0.0d0
          s%sigma(2)=0.0d0
          s%sigma(3)=1.0d0
        elseif(csigma=='-z')then
          s%sigma(1)=0.0d0
          s%sigma(2)=0.0d0
          s%sigma(3)=-1.0d0
        else
        end if

        do i=1,s%nparam
          write(*,'(A,2x,A)')trim(s%paramlist(i)),trim(s%paramvalues(i))
        end do

      return
    end subroutine

    subroutine rfptwoLLGS_close(s)
      implicit none
      class(rfptwoLLGS) s

      deallocate(s%charlist,s%paramlist,s%paramvalues)
      !write(*,*)'Deallocating parameters in rfptwoLLGS'

      return
    end subroutine

  end module
