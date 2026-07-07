  module func_LLGS
    implicit none
  contains
    double precision function LLGS_fx(mx,my,mz,Hx,Hy,Hz,alpha,gammap,taux,tauy,tauz)
      implicit none
      double precision mx,my,mz,Hx,Hy,Hz,alpha,gammap,taux,tauy,tauz
      LLGS_fx=-1.0d0*gammap*((my*Hz-mz*Hy)+alpha*((mx*Hx+my*Hy+mz*Hz)*mx-Hx)-taux-alpha*(my*tauz-mz*tauy))
      return
    end function

    double precision function LLGS_fy(mx,my,mz,Hx,Hy,Hz,alpha,gammap,taux,tauy,tauz)
      implicit none
      double precision mx,my,mz,Hx,Hy,Hz,alpha,gammap,taux,tauy,tauz
      LLGS_fy=-1.0d0*gammap*((mz*Hx-mx*Hz)+alpha*((mx*Hx+my*Hy+mz*Hz)*my-Hy)-tauy-alpha*(mz*taux-mx*tauz))
      return
    end function

    double precision function LLGS_fz(mx,my,mz,Hx,Hy,Hz,alpha,gammap,taux,tauy,tauz)
      implicit none
      double precision mx,my,mz,Hx,Hy,Hz,alpha,gammap,taux,tauy,tauz
      LLGS_fz=-1.0d0*gammap*((mx*Hy-my*Hx)+alpha*((mx*Hx+my*Hy+mz*Hz)*mz-Hz)-tauz-alpha*(mx*tauy-my*taux))
      return
    end function

    subroutine LLGS_calc_tau(torque,sigma,mx,my,mz,tau_x,tau_y,tau_z)
      implicit none
      double precision torque,sigma(3),mx,my,mz,tau_x,tau_y,tau_z
      double precision mdsigma

      mdsigma=mx*sigma(1)+my*sigma(2)+mz*sigma(3)
      tau_x=torque*(sigma(1)-mdsigma*mx)
      tau_y=torque*(sigma(2)-mdsigma*my)
      tau_z=torque*(sigma(3)-mdsigma*mz)

      return
    end subroutine

  end module
