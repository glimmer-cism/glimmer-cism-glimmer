
!***********************************************************************
module remap_glamutils      
!***********************************************************************

    ! *sp* contains various subroutines needed when using LANL incremental remapping code
    ! for thickness evolution in glam/glimmer codes

    ! *sp* just stubs for now

    implicit none
    private

    public :: horizontal_remap_init, horizontal_remap_final, &
              horizontal_remap_in, horizontal_remap_out

    contains


    ! *sp* initialize variables for use in inc. remapping code   
    subroutine horizontal_remap_init( )
    end subroutine horizontal_remap_init


    ! *sp* get GLAM variables in order for use in inc. remapping code   
    subroutine horizontal_remap_in( ) 
    end subroutine horizontal_remap_in


    ! *sp* take output from inc. remapping and put back in GLAM format
    subroutine horizontal_remap_out( )
    end subroutine horizontal_remap_out


    ! *sp* deallocate variables for use in inc. remapping code   
    subroutine horizontal_remap_final( )
    end subroutine horizontal_remap_final



!***********************************************************************
end module remap_glamutils
!***********************************************************************

