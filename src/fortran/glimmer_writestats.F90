module glimmer_writestats_module
  !> F90 wrapper to gc_writestats
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
contains
  subroutine glimmer_writestats(resname, cfgname)
    use glimmer_global, only : dp
    implicit none
    character(len=*), intent(in) :: resname !< name of the output result file
    character(len=*), intent(in) :: cfgname !< name of configuration file

    call gf_writestats(resname,cfgname)
  end subroutine glimmer_writestats

end module glimmer_writestats_module
