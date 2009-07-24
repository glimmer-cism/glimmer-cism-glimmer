! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  spin_types.f90 - part of the GLIMMER ice model           + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! The Spin drivers are intended for model spin-up for Antarctica and Greenland.  The
! model_type variable differentiates between the two models, (0) for Antarctica,
! and (1) for Greenland.  These can be set inside the configuration file.  Also,
! included in the drivers are a earlier versions of the mass balance and
! temperature and are meant to be the same as the ones used in the EISMINT model
! intercomparison experiment.  More information can be found at
! http://homepages.vub.ac.be/~phuybrec/eismint.html.  These can be specified by
! setting the variable use_simple = 1 in the configuration file under temperature
! and massbalance separately.  These drivers also vary for Antarctica and
! Greenland.  The other drivers available in these drivers are based on the
! current papers of Philippe Huybrechts:  

! P. Huybrechts and J. de Wolde. The dynamic response of the greenland and
! antarctic ice sheets to multiple-century climatic warming. Journal of Climate,
! 12(8):2169-2188, 1999.

! I. Janssens and P. Huybrects. The treatment of meltwater retention in
! massbalance parameterisations of the greenland ice sheet. Annals of Glaciology,
! 31:133-140, 2000.

! P. Huybrechts. Sea-level changes at the LGM from ice-dynamic reconstructions
! of the Greenland and Antarctic ice sheets during the glacial cycles. Quater-
! nary Science Reviews, 21:203-231, 2002.

! P. Huybrechts, O. Rybak, F. Pattyn, U. Ruth, and D. Steinhage. Ice thinning,
! upstream advection, and non-climatic biases for the upper ice sheet. Climate
! of the Past, 3:577-589, 2007.

! The Spin series of drivers were developed by Brian Hand, at the University of
! Montana, 2008.
#ifdef HAVE_CONFIG_H
#include "config.inc"
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
