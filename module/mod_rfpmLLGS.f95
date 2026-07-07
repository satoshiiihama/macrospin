  module mod_rfpmLLGS
    ! parameter list
    implicit none
    type rfpmLLGS
      integer nparam,nlist
      character(200),allocatable::charlist(:),paramlist(:),paramvalues(:)
      character(200) filename
      integer N,Neq,nl,nint_mh
      double precision dt,time_exc,alpha_mh
      double precision Hext,theta_h,phi_h
      double precision alpha,Ms,Ku,dF,Jex,Jex2,DMI,Kuni
      double precision sigma(3),jc
      double precision,allocatable::mx(:),my(:),mz(:),dir_D(:,:),dir_uni(:,:)
      double precision Hmin,Hmax,dH
      character(200) cfilename,cfilename_mh
      character(50) cn,cneq,cdt,cHext,ctheta_h,cphi_h,calpha,cMs,cKu,cdF,cJex,cJex2,cDMI,cdir_D,cmx,cmy,cmz
      character(50) csigma,cjc,ctime_exc,cinit_m,cprint_mh_all,cdir_uni
      logical print_mh_all
    contains
      procedure::init=>rfpmLLGS_init
      procedure::help=>rfpmLLGS_help
      procedure::set=>rfpmLLGS_set
      procedure::read=>rfpmLLGS_read
      procedure::ana=>rfpmLLGS_ana
      procedure::print=>rfpmLLGS_print
      procedure::close=>rfpmLLGS_close
    end type
  contains

    subroutine rfpmLLGS_init(s)
      use constant
      implicit none
      class(rfpmLLGS) s
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

      s%nparam=26
      allocate(s%paramlist(s%nparam),s%paramvalues(s%nparam))
      s%paramvalues='0'
      s%paramlist(1)='nl'
      s%paramlist(2)='dt'
      s%paramlist(3)='N'
      s%paramlist(4)='Neq'
      s%paramlist(5)='Hext'
      s%paramlist(6)='theta_h'
      s%paramlist(7)='phi_h'
      s%paramlist(8)='alpha'
      s%paramlist(9)='Ms'
      s%paramlist(10)='Ku'
      s%paramlist(11)='dF'
      s%paramlist(12)='Jex'
      s%paramlist(13)='Jex2'
      s%paramlist(14)='DMI'
      s%paramlist(15)='dir_D'
      s%paramlist(16)='m_ini'
      s%paramlist(17)='filename'
      s%paramlist(18)='Hmin'
      s%paramlist(19)='Hmax'
      s%paramlist(20)='dH'
      s%paramlist(21)='alpha_mh'
      s%paramlist(22)='filename_mh'
      s%paramlist(23)='init_m'
      s%paramlist(24)='print_mh_all'
      s%paramlist(25)='Kuni'
      s%paramlist(26)='dir_uni'

      s%nl=2
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
      s%Jex2=0.0d-3
      s%DMI=1.0d-3
      s%alpha_mh=1.0d0

      s%Hmin=-2.0d0
      s%Hmax=2.0d0
      s%dH=0.1d0

      s%jc=1.0d-3/(pi*(20.0d-9)*(20.0d-9))
      s%time_exc=10.0d0

      s%cfilename='test.txt'
      
      return
    end subroutine

    subroutine rfpmLLGS_print(s)
      use constant
      implicit none
      class(rfpmLLGS) s
      integer i
      character(20) cm1,cm2,csigma
      !integer ndis,ntra,ntes,nset,ntype,np,nv,ncomp,nccomp,ndim,ndelay,nzero,nloop,nnorm,nzsc,nprint,npnum
      !character(20) cdate,ci,cj 

        !do i=1,s%nparam
        !  write(*,*)s%paramvalues(i)
        !end do
      
        read(s%paramvalues(1),*)s%nl
        read(s%paramvalues(2),*)s%dt
        read(s%paramvalues(3),*)s%n
        read(s%paramvalues(4),*)s%neq
        read(s%paramvalues(5),*)s%Hext
        read(s%paramvalues(6),*)s%theta_h
        read(s%paramvalues(7),*)s%phi_h
        read(s%paramvalues(8),*)s%alpha
        read(s%paramvalues(9),*)s%Ms
        read(s%paramvalues(10),*)s%Ku
        read(s%paramvalues(11),*)s%dF
        read(s%paramvalues(12),*)s%Jex
        read(s%paramvalues(13),*)s%Jex2
        read(s%paramvalues(14),*)s%DMI

        allocate(s%mx(s%nl),s%my(s%nl),s%mz(s%nl),s%dir_D(3,s%nl),s%dir_uni(3,s%nl))

        i=index(trim(s%paramvalues(17))," ",back=.true.)
        s%cfilename=trim(s%paramvalues(17)(1+i:len_trim(s%paramvalues(17))))
        
        s%theta_h = s%theta_h * pi / 180.0d0
        s%phi_h = s%phi_h * pi / 180.0d0

        i=index(trim(s%paramvalues(15))," ",back=.true.)
        cm1=trim(s%paramvalues(15)(1+i:len_trim(s%paramvalues(15))))
        i=index(trim(s%paramvalues(16))," ",back=.true.)
        cm2=trim(s%paramvalues(16)(1+i:len_trim(s%paramvalues(16))))

        read(s%paramvalues(18),*)s%Hmin
        read(s%paramvalues(19),*)s%Hmax
        read(s%paramvalues(20),*)s%dH
        read(s%paramvalues(21),*)s%alpha_mh

        i=index(trim(s%paramvalues(22))," ",back=.true.)
        s%cfilename_mh=trim(s%paramvalues(22)(1+i:len_trim(s%paramvalues(22))))

        i=index(trim(s%paramvalues(23))," ",back=.true.)
        s%cinit_m=trim(s%paramvalues(23)(1+i:len_trim(s%paramvalues(23))))
        if(s%cinit_m=='m_ini')then
          s%nint_mh=1
        else if(s%cinit_m=='cont')then
          s%nint_mh=0
        end if

        i=index(trim(s%paramvalues(24))," ",back=.true.)
        s%cprint_mh_all=trim(s%paramvalues(24)(1+i:len_trim(s%paramvalues(24))))
        if(s%cprint_mh_all=='true')then
          s%print_mh_all=.true.
        else
          s%print_mh_all=.false.
        end if

        read(s%paramvalues(25),*)s%Kuni

        i=index(trim(s%paramvalues(26))," ",back=.true.)
        s%cdir_uni=trim(s%paramvalues(26)(1+i:len_trim(s%paramvalues(26))))
        
        if(cm1=='x')then
          do i=1,s%nl
            s%dir_D(1,i)=1.0d0
            s%dir_D(2,i)=0.0d0
            s%dir_D(3,i)=0.0d0
          end do
        elseif(cm1=='-x')then
          do i=1,s%nl
            s%dir_D(1,i)=-1.0d0
            s%dir_D(2,i)=0.0d0
            s%dir_D(3,i)=0.0d0
          end do
        elseif(cm1=='x-alt')then
          do i=1,s%nl
            if(mod(i,2)==0)then
              s%dir_D(1,i)=1.0d0
            else
              s%dir_D(1,i)=-1.0d0
            end if
            s%dir_D(2,i)=0.0d0
            s%dir_D(3,i)=0.0d0
          end do
        elseif(cm1=='y')then
          do i=1,s%nl
            s%dir_D(1,i)=0.0d0
            s%dir_D(2,i)=1.0d0
            s%dir_D(3,i)=0.0d0
          end do
        elseif(cm1=='-y')then
          do i=1,s%nl
            s%dir_D(1,i)=0.0d0
            s%dir_D(2,i)=-1.0d0
            s%dir_D(3,i)=0.0d0
          end do
        elseif(cm1=='y-alt')then
          do i=1,s%nl
            s%dir_D(1,i)=0.0d0
            if(mod(i,2)==0)then
              s%dir_D(2,i)=1.0d0
            else
              s%dir_D(2,i)=-1.0d0
            end if
            s%dir_D(3,i)=0.0d0
          end do
        elseif(cm1=='z')then
          do i=1,s%nl
            s%dir_D(1,i)=0.0d0
            s%dir_D(2,i)=0.0d0
            s%dir_D(3,i)=1.0d0
          end do
        elseif(cm1=='-z')then
          do i=1,s%nl
            s%dir_D(1,i)=0.0d0
            s%dir_D(2,i)=0.0d0
            s%dir_D(3,i)=-1.0d0
          end do
        elseif(cm1=='z-alt')then
          do i=1,s%nl
            s%dir_D(1,i)=0.0d0
            s%dir_D(2,i)=0.0d0
            if(mod(i,2)==0)then
              s%dir_D(3,i)=1.0d0
            else
              s%dir_D(3,i)=-1.0d0
            end if
          end do
        else
        end if

        if(s%cdir_uni=='x')then
          do i=1,s%nl
            s%dir_uni(1,i)=1.0d0
            s%dir_uni(2,i)=0.0d0
            s%dir_uni(3,i)=0.0d0
          end do
        elseif(s%cdir_uni=='y')then
          do i=1,s%nl
            s%dir_uni(1,i)=0.0d0
            s%dir_uni(2,i)=1.0d0
            s%dir_uni(3,i)=0.0d0
          end do
        elseif(s%cdir_uni=='z')then
          do i=1,s%nl
            s%dir_uni(1,i)=0.0d0
            s%dir_uni(2,i)=0.0d0
            s%dir_uni(3,i)=1.0d0
          end do
        elseif(s%cdir_uni=='xy')then
          do i=1,s%nl
            s%dir_uni(3,i)=0.0d0
            if(mod(i,2)==0)then
              s%dir_uni(1,i)=1.0d0
              s%dir_uni(2,i)=0.0d0
            else
              s%dir_uni(1,i)=0.0d0
              s%dir_uni(2,i)=1.0d0
            end if
          end do
        else
        end if

        if(cm2=='x')then
          s%mx=1.0d0
          s%my=0.0d0
          s%mz=0.0d0
        elseif(cm2=='-x')then
          s%mx=-1.0d0
          s%my=0.0d0
          s%mz=0.0d0
        elseif(cm2=='x-alt')then
          s%my=0.0d0
          s%mz=0.0d0
          do i=1,s%nl
            if(mod(i,2)==0)then
              s%mx(i)=1.0d0
            else
              s%mx=(i)-1.0d0
            end if
          end do
        elseif(cm2=='y')then
          s%mx=0.0d0
          s%my=1.0d0
          s%mz=0.0d0
        elseif(cm2=='-y')then
          s%mx=0.0d0
          s%my=-1.0d0
          s%mz=0.0d0
        elseif(cm2=='y-alt')then
          s%mx=0.0d0
          s%mz=0.0d0
          do i=1,s%nl
            if(mod(i,2)==0)then
              s%my(i)=1.0d0
            else
              s%my(i)=-1.0d0
            end if
          end do
        elseif(cm2=='z')then
          s%mx=0.0d0
          s%my=0.0d0
          s%mz=1.0d0
        elseif(cm2=='-z')then
          s%mx=0.0d0
          s%my=0.0d0
          s%mz=-1.0d0
        elseif(cm2=='z-alt')then
          s%mx=0.0d0
          s%my=0.0d0
          do i=1,s%nl
            if(mod(i,2)==0)then
              s%mz(i)=1.0d0
            else
              s%mz(i)=-1.0d0
            end if
          end do
        else
        end if

        do i=1,s%nparam
          write(*,'(A,2x,A)')trim(s%paramlist(i)),trim(s%paramvalues(i))
        end do

      return
    end subroutine

    subroutine rfpmLLGS_help(s)
      implicit none
      class(rfpmLLGS) s

        write(*,*)' This program is to calculate macrospin dynamics of multilayer coupled magnetization '
        write(*,*)' Parameters for the calculation must be argued after .exe statement'
        write(*,*)' The meaning of each parameters are described below'
        write(*,*)' dt: time-step for Runge-Kutta simulation'
        write(*,*)' N: total time step number'
        write(*,*)' Neq: time step number for equilibration'
        write(*,*)' Hext: external magnetic field strength in unit of [T]'
        write(*,*)' theta_h: external magnetic field angle from film normal in unit of degree'
        write(*,*)' phi_h: external magnetic field in-plane angle from x-axis in unit of degree'
        write(*,*)' alpha: Gilbert damping parameter for both layer'
        write(*,*)' Ms: Saturation magnetization in unit of [A/m]'
        write(*,*)' Ku: Uniaxial magnetic anisotropy constant in unit of [J/m3]'
        write(*,*)' dF: ferromagnetic layer thickness in unit of [m]'
        write(*,*)' Jex: RKKY interlayer coupling strength in unit of [J/m2]'
        write(*,*)' DMI: interlayer DMI coupling strength in unit of [J/m2]'
        write(*,*)' dir_D: direction of the DMI vector (+- x, y, z)'
        write(*,*)' m_ini: initial magnetization configuration (+- x, y, z)'
        write(*,*)' Hmin: minimum Hext value for m-h calculation'
        write(*,*)' Hmax: maximum Hext value for m-h calculation'
        write(*,*)' dH: step Hext value for m-h calculation'
        write(*,*)' filename: filename for output'
       
      return
    end subroutine

    subroutine rfpmLLGS_read(s)
      implicit none
      class(rfpmLLGS) s
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

    subroutine rfpmLLGS_ana(s)
      implicit none
      class(rfpmLLGS) s
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

    subroutine rfpmLLGS_set(s)
      implicit none
      class(rfpmLLGS) s
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
      !write(*,*)'Allocating parameters in rfpmLLGS'

      return
    end subroutine

    subroutine rfpmLLGS_close(s)
      implicit none
      class(rfpmLLGS) s

      deallocate(s%charlist,s%paramlist,s%paramvalues)
      !write(*,*)'Deallocating parameters in rfpmLLGS'

      return
    end subroutine

  end module
