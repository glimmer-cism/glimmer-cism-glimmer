module plume_ambient

  use plume_global
  use plume_functions

contains

  subroutine ambient()

    ! set ambient profiles of temperture, salinity and density
    ! density is unaffected by frazil as the ambient is ice-free



    implicit none

    ! local variables

    integer i

    real(kind=kdp) depth,pressure,ttt,rhopot

    depth = 0.d0

    open(41,file='ambout')

    do i = 1,namb

       tamb(i) = temptop + depth*tgrad
       samb(i) = salttop + depth*sgrad

       if (rholinear) then
          rhovf(i) = rho_func_linear(tamb(i),samb(i))              
          rhopot = rhovf(i)
       else
          if (thermobar) then
             pressure = depth*1.0d-1
             ttt = tinsitu_func(tamb(i),samb(i),pressure)
          else
             pressure = 0.d0    
             ttt = tamb(i)
          end if
          rhovf(i) = rho_func_nonlinear(ttt,samb(i),pressure)              
          rhopot = rho_func_nonlinear(tamb(i),samb(i),0.d0)
       end if

       write(41,*) depth,tamb(i),ttt,samb(i),rhovf(i),rhopot

       depth = depth + dzincr

    end do

  end subroutine ambient

end module
