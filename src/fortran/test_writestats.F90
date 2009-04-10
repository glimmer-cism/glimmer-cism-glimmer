!> test glimmer_writestats module
!!
!! \author Magnus Hagdorn
!! \date April 2009

program test_writestats
  use glimmer_writestats_module
  implicit none

  call glimmer_writestats("results","model.conf")
end program test_writestats
