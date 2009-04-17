! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_ground.f90 - part of the GLIMMER ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

#include <glide_mask.inc>
module glide_ground
  use glide_types
  use glimmer_global
  implicit none
contains
!-------------------------------------------------------------------------
  subroutine glide_initialise_backstress(thck,backstress)
     implicit none
     
     real(dp), dimension(:,:), intent(in) :: thck !*FD Ice thickness
     real(sp), dimension(:,:), intent(inout) :: backstress !*FD Backstress
     backstress = 0.0
     where(thck > 0.0) 
         backstress = 1.0
     end where
  end subroutine glide_initialise_backstress
!-------------------------------------------------------------------------

  subroutine glide_marinlim(which,thck,relx,topg,flwa,levels,mask,mlimit,calving_fraction,eus,ablation_field,backstress, & 
                 tempanmly,dew,dns)


    !*FD Removes non-grounded ice, according to one of two altenative
    !*FD criteria, and sets upper surface of non-ice-covered points 
    !*FD equal to the topographic height, or sea-level, whichever is higher.

    use glimmer_global, only : dp, sp
    use glimmer_physcon, only : rhoi, rhoo, grav, gn
    use glide_vertint, only : vertint_output2d
    use glimmer_paramets, only: thk0
    use glide_mask
    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer,                intent(in)    :: which   !*FD Option to choose ice-removal method
                                                     !*FD \begin{description}
                                                     !*FD \item[0] Set thickness to zero if 
                                                     !*FD relaxed bedrock is below a given level.
                                                     !*FD \item[1] Set thickness to zero if
                                                     !*FD ice is floating.
                                                     !*FD \end{description}
    real(dp),dimension(:,:),intent(inout) :: thck    !*FD Ice thickness (scaled)
    real(dp),dimension(:,:),intent(in)    :: relx    !*FD Relaxed topography (scaled)
    real(dp),dimension(:,:),intent(in)    :: topg    !*FD Present bedrock topography (scaled)
    real(dp),dimension(:,:,:),intent(in)  :: flwa    !*FD Vertically averaged ice hardness
    real(dp),dimension(:)    ,intent(in)  :: levels    !*FD Vertically averaged ice hardness
    integer, dimension(:,:),pointer       :: mask    !*FD grid type mask
    real(dp)                              :: mlimit  !*FD Lower limit on topography elevation for
                                                     !*FD ice to be present (scaled). Used with 
                                                     !*FD $\mathtt{which}=0$.
    real(dp), intent(in) :: calving_fraction         !*FD fraction of ice lost when calving Used with 
                                                     !*FD $\mathtt{which}=3$.
    real, intent(inout) :: eus                       !*FD eustatic sea level
    real(sp),dimension(:,:),intent(inout) :: ablation_field !*FD this is passed as climate%calving

    real(dp), parameter :: con = - rhoi / rhoo
    real(dp), parameter :: sigmaxx = 0.5 * rhoi * grav * (1.0 - rhoi / rhoo)
    real(dp), parameter :: theta = 0.5
    real(dp), dimension(2,2) :: A
    
    real(sp), dimension(:,:), intent(inout) :: backstress
    real(sp) :: tempanmly
    real(dp), intent(in) ::  dew,dns

    integer ew,ns
    real(dp) :: sigmab
    !---------------------------------------------------------------------

    sigmab = 0.95

    ablation_field=0.0

    select case (which)
        
    case(1) ! Set thickness to zero if ice is floating
      where (GLIDE_IS_FLOAT(mask))
        ablation_field=thck
        thck = 0.0d0
      end where

    case(2) ! Set thickness to zero if relaxed bedrock is below a 
       ! given level
       where (relx <= mlimit+eus)
          ablation_field=thck
          thck = 0.0d0
       end where

    case(3) ! remove fraction of ice when floating
       do ns = 2,size(thck,2)-1
          do ew = 2,size(thck,1)-1
             if (GLIDE_IS_CALVING(mask(ew,ns))) then
                ablation_field(ew,ns)=(1.0-calving_fraction)*thck(ew,ns)
                thck(ew,ns) =  calving_fraction*thck(ew,ns)
                !mask(ew,ns) = glide_mask_ocean
             end if
          end do
       end do

    case(4) ! Set thickness to zero at marine edge if present
            ! bedrock is below a given level
       where (GLIDE_IS_MARINE_ICE_EDGE(mask).and.topg<mlimit+eus)
          ablation_field=thck
          thck = 0.0d0
       end where
    case(5) ! Relation based on computing the horizontal stretching
            ! of the unconfined ice shelf (\dot \epsilon_{xx}) and multiplying by H.
            ! 
    do ns = 2, size(backstress,2)-1
       do ew = 2, size(backstress,1)-1
            if(backstress(ew,ns) < 1.0) then
               if (tempanmly > 0.0) then
                  backstress(ew,ns) = 0.0
               else
                  backstress(ew,ns) = 1-0.85*exp(tempanmly) 
               end if
            end if
       end do
    end do 
       do ns = 2,size(thck,2)-1
          do ew = 2,size(thck,1)-1
             if (GLIDE_IS_GROUNDING_LINE(mask(ew,ns))) then
                call vertint_output2d(flwa(:,ew-1:ew,ns-1:ns),A, levels * thck(ew,ns))
                ablation_field(ew,ns)= ((dew*dns)/(50.0d3)**2)* theta * A(2,2) * (sigmaxx * &
                thck(ew,ns)  * (1 - backstress(ew,ns))) ** gn
                if ((thck(ew,ns) - ablation_field(ew,ns)) >= 0.0) then
                  thck(ew,ns) = thck(ew,ns) - ablation_field(ew,ns) 
                else 
                  thck(ew,ns) = 0.0d0
                end if
            end if
          end do
       end do
         
       where (GLIDE_IS_MARINE_ICE_EDGE(mask).and.topg<mlimit+eus)
          ablation_field=thck
          thck = 0.0d0
       end where

    end select

  end subroutine glide_marinlim
  
!-------------------------------------------------------------------------
  !This function returns the correct grounding line using the data given 
  ! the mask reference point.  dir is specifying 'ew' or 'ns', but can be 
  ! left null if there's only one option.
  real function get_ground_line(model,ew1,ns1,ew2,ns2)
     use glide_types
     implicit none
     type(glide_global_type) :: model        !*FD model instance
     integer ns1,ew1,ns2,ew2,slot_ns,slot_ew !grounding line in ns/ew direction
     real appr_ground !grounding line
     !this is assuming a greedy grounding line progression.  Always setting the
     ! mask by rounding up to the next grid point
     if (ns1 .eq. ns2) then
         slot_ew = min(ew1,ew2)
         appr_ground = model%ground%gl_ns(slot_ew,ns1)
     else if (ew1 .eq. ew2) then
         slot_ns = min(ns1,ns2)
         appr_ground = model%ground%gl_ew(ew1,slot_ns)
     end if
     return
  end function get_ground_line
    
!-------------------------------------------------------------------------
  subroutine set_ground_line(model,ew1,ns1,ew2,ns2,value)
     use glide_types
     implicit none
     type(glide_global_type) :: model        !*FD model instance
     integer, intent(in) :: ns1 !grounding line in ns direction
     integer, intent(in) :: ew1 !grounding line in ew direction
     integer, intent(in) :: ns2 !grounding line in ns direction
     integer, intent(in) :: ew2 !grounding line in ew direction
     real(sp), intent(in) :: value !grounding line in ew direction
     integer slot_ew, slot_ns !integers to compute the min
     !this is assuming a greedy grounding line progression.  Always setting the
     ! mask by rounding up to the next grid point
     if (ns1 .eq. ns2) then
         slot_ew = min(ew1,ew2)
         model%ground%gl_ew(slot_ew,ns1) = value
     else if (ew1 .eq. ew2) then
         slot_ns = min(ns1,ns2)
         model%ground%gl_ns(ew1,slot_ns) = value
     end if
  end subroutine set_ground_line
!-------------------------------------------------------------------------
  real function lin_reg_xg(model,ew,ns,j1ew,j1ns)
     use glide_types
     use glimmer_physcon, only : rhoi, rhoo
     type(glide_global_type) :: model        !*FD model instance
     integer, intent(in) :: ns !grounding line in ns direction
     integer, intent(in) :: ew !grounding line in ew direction
     integer, intent(in) :: j1ns !ice shelf in ns direction
     integer, intent(in) :: j1ew !ice shelf line in ew direction
     real(sp) ::  xg                     !grounding line
     real(sp) ::  dx                      !distance between gridpts
     real(sp) ::  xj                        !grounding line
     real(sp) :: fj                        !f at grid pnt j
     real(sp) :: fj_1                      !f evaluated at j (+/-) 1
     real(sp) :: df                        !delta f of fj,jf_1
     if (ew .eq. j1ew) then
        dx = model%numerics%dns
        xj = ns*dx
     else
        dx = model%numerics%dew
        xj = ew*dx
     end if
     !set the pattyn f function - assuming ocean water
     fj = (model%climate%eus - model%geometry%topg(ew,ns))*rhoo/(rhoi*model%geometry%thck(ew,ns))
     fj_1 = (model%climate%eus - model%geometry%topg(j1ew,j1ns))*rhoo/(rhoi*model%geometry%thck(j1ew,j1ns))
     df = (fj_1 - fj)/dx                
     xg = (1 - fj + df*xj)/df
     return 
  end function lin_reg_xg
!-------------------------------------------------------------------------
  !also need an update mask
  subroutine update_ground_line(model,mask)
     use glide_mask
     implicit none
     type(glide_global_type) :: model        !*FD model instance
     integer, dimension(:,:),pointer :: mask    !*FD grid type mask
     integer ew,ns,jns,jew,j1ns,j1ew
     real(sp) :: xg                        !grounding line
     !this is assuming the grounding line is the last grounded pt on the mask
     !reset grounding line data to zero
     model%ground%gl_ew = 0.0
     model%ground%gl_ns = 0.0
     do ns = 1,model%general%nsn
        do ew = 1,model%general%ewn
            if (GLIDE_IS_GROUNDING_LINE(mask(ew,ns))) then
                !the grounding line always rounds down so it is grounded.
                
                !northern grounding line
                if (GLIDE_IS_OCEAN(mask(ew,ns - 1)) &
                        .or. (GLIDE_IS_FLOAT(mask(ew,ns - 1)))) then
                    xg = lin_reg_xg(model,ew,ns,ew,ns-1)
                    call set_ground_line(model,ew,ns,ew,ns-1,xg)
                !southern grounding line
                else if (GLIDE_IS_OCEAN(mask(ew,ns + 1)) &
                        .or. (GLIDE_IS_FLOAT(mask(ew,ns + 1)))) then
                    xg = lin_reg_xg(model,ew,ns,ew,ns+1)
                    call set_ground_line(model,ew,ns,ew,ns+1,xg) 
                end if 
                
                !western grounding line
                if (GLIDE_IS_OCEAN(mask(ew - 1,ns)) &
                        .or. GLIDE_IS_FLOAT(mask(ew - 1,ns))) then
                    xg = lin_reg_xg(model,ew,ns,ew - 1,ns)
                    call set_ground_line(model,ew,ns,ew-1,ns,xg)
                !eastern grounding line
                else if (GLIDE_IS_OCEAN(mask(ew + 1,ns)) &
                        .or. GLIDE_IS_FLOAT(mask(ew + 1,ns))) then
                    xg = lin_reg_xg(model,ew,ns,ew + 1,ns)
                    call set_ground_line(model,ew,ns,ew + 1,ns,xg)
                !need to update mask
                end if
            end if 
        end do
     end do
  end subroutine update_ground_line
!---------------------------------------------------------------------------
end module glide_ground
