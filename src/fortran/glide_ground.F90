! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_ground.f90 - part of the GLIMMER ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module glide_ground
  use glide_types
  use glimmer_log  !do i need this?
  implicit none

!-------------------------------------------------------------------------

  subroutine glide_marinlim(which,thck,relx,topg,flwa,levels,mask,mlimit,calving_fraction,eus,ablation_field)

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
       do ns = 2,size(thck,2)-1
          do ew = 2,size(thck,1)-1
             if (GLIDE_IS_GROUNDING_LINE(mask(ew,ns))) then
                call vertint_output2d(flwa(:,ew-1:ew,ns-1:ns),A, levels * thck(ew,ns))
                ablation_field(ew,ns)=  theta * A(2,2) * (sigmaxx * &
                thck(ew,ns)  * (1 - sigmab)) ** gn
                if ((thck(ew,ns) - ablation_field(ew,ns)) >= 0.0) then
                thck(ew,ns) = thck(ew,ns) - ablation_field(ew,ns) 
                else 
                thck(ew,ns) = 0.0d0
                end if
            end if
          end do
       end do

    end select

  end subroutine glide_marinlim

end module glide_ground
