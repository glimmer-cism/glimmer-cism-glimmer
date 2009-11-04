module plume_functions

use plume_global

contains

real(kind=kdp) function freezing_temp_func(s,dep)

  implicit none

! calculate freezing temperature according to linear equation

! millero, f.j.
! annex 6: freezing point of seawater.  eighth report of the joint
! panel of oceanographic tables and standards.
! unesco technical papers in marine science; 28; 1978, pp29-31.

      real(kind=kdp) s,dep
      real(kind=kdp) a,b,c

! preliminary definitions
      a = -5.73d-2
      b = 8.32d-2
      c = -7.61d-4

! apply function definition      
      freezing_temp_func = a*s + b + c*dep
!    
      return
      end function

!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************
real(kind=kdp) function rho_func_linear(t,s) 
  
  implicit none

! calculate density according to linear equation of state
! (no thermobaricity by definition)

! jenkins, a. and bombosch, a.
! modeling the effects of frazil ice crystals on the dynamics and 
! thermodynamics of ice shelf water plumes
! journal of geophysical research; 100; c4; 1995; p.6967-6981

      real(kind=kdp) t,s
      real(kind=kdp) rho0,t0,s0,bett,bets

! preliminary definitions
      rho0 = 1030.d0
      t0 = -2.d0
      s0 = 34.5d0
      bett = 3.87d-5
      bets = 7.86d-4     

! apply function definition      
      rho_func_linear = rho0*(1.d0 + bets*(s - s0) - bett*(t - t0))
!    
      return
      end function

!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************
real(kind=kdp) function rho_func_nonlinear(t,s,p) 

  implicit none

! calculate density acoording to nonlinear equation of state

! fofonoff, n.p. and millard, r.c. 
! algorithms for computation of fundamental properties of seawater 
! unesco technical papers in marine science; 44; 1983; 53 p.

      real(kind=kdp) t,s,p
      real(kind=kdp) t2,t3,t4,t5,s2,s32,p2,rt,rst,kt,kst,kstp

! get powers of variables
      t2 = t*t
      t3 = t2*t
      t4 = t3*t
      t5 = t4*t
      s2 = s*s
      s32 = dsqrt(s2*s)
      p2 = p*p

! apply function definition      
      rt = 999.842594d0 + 6.793952d-2*t - 9.095290d-3*t2 &
                     + 1.001685d-4*t3 - 1.120083d-6*t4 + 6.536332d-9*t5

      rst = rt + s*(8.24493d-1 - 4.0899d-3*t + 7.6438d-5*t2 &
                                         - 8.2467d-7*t3 + 5.3875d-9*t4) &
              + s32*(-5.72466d-3 + 1.0227d-4*t - 1.6546d-6*t2) &
              + 4.8314d-4*s2 

      kt = 19652.21d0 + 148.4206d0*t - 2.327105d0*t2 + 1.360477d-2*t3  &
                                                       - 5.155288d-5*t4

      kst = kt + s*(54.6746d0 - 6.03459d-1*t + 1.09987d-2*t2 &
                                                        - 6.1670d-5*t3) &
              + s32*(7.944d-2 + 1.6483d-2*t - 5.3009d-4*t2)

      kstp = kst + p*(3.239908d0 + 1.43713d-3*t + 1.16092d-4*t2  &
                                                      - 5.77905d-7*t3) &
                + p*s*(2.2838d-3 - 1.0981d-5*t - 1.6078d-6*t2) &
                + 1.91075d-4*p*s32  &
                + p2*(8.50935d-5 - 6.12293d-6*t + 5.2787d-8*t2) &
                + p2*s*(-9.9348d-7 + 2.0816d-8*t + 9.1697d-10*t2) 

      rho_func_nonlinear = rst/(1.d0 - p/kstp)
!    
      return
      end function

!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************
real(kind=kdp) function tinsitu_func(tp,s,pin)  
  implicit none

! calculate in-situ temperature from potential temperature

! fofonoff, n.p. and millard, r.c. 
! algorithms for computation of fundamental properties of seawater 
! unesco technical papers in marine science; 44; publ: 1983; 53 p.

! uses potential temperature in celsius, pressure in kilobars 
! and salinity in psu) 
      integer it
      real(kind=kdp) pin,s,tp
      real(kind=kdp) p,t,z,tt,pit,sit,tit

      p = pin*1.0d-3
      t = tp + 2.d0                                                           
      z = s - 35.d0                                                           

      do it = 1,10                                                  
        tt = t*t                                                          
        pit = p*(8.9309d-1 - t*(3.1628d-2 - 2.1987d-4*t - 5.0484d-3*p) &
     &                                                    - 1.6056d-1*p)        
        sit = z*(1.7439d-2 - 2.9778d-4*t - 4.1057d-3*p)                            
        tit = t*(8.3198d-2 - 5.4065d-4*t + 4.0274d-6*tt)                        
        t = tp + p*(3.6504d-1 + pit + sit + tit)                                       
      end do

      tinsitu_func=t    

      return                                                            
end function

end module                                                                
