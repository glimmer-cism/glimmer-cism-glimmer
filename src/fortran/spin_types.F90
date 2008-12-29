! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  spin_types.f90 - part of the GLIMMER ice model           + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module spin_types

  use glimmer_global, only: sp
  use glimmer_pdd
  use spin_mb
  use spin_temp
  use spin_slc
  implicit none
  type spin_climate_type
     type(spin_mb_type)  :: mb  !*FD mass balance forcing 
     type(spin_temp_type) :: temp !*FD temperature forcing
     type(spin_slc_type)  :: slc  !*FD sea level forcing
  end type spin_climate_type
end module spin_types
