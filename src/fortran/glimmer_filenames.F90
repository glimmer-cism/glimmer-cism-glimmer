
module glimmer_filenames

  !*FD Module to handle setting a working directory for glimmer

  use glimmer_global,only: dirsep

  implicit none

  character(200) :: workingdir = ''
  !*FD Working directory for all file operations. Absolute paths are unaffected

contains

  subroutine glimmer_set_path(path)

    use glimmer_log

    character(*),intent(in) :: path

    workingdir=path
    call write_log('Set GLIMMER working dir to :'//trim(workingdir))

  end subroutine glimmer_set_path

  character(200) function process_path(path)

    character(*),intent(in) :: path

    character(200) :: alpath

    alpath=adjustl(path)

    if (alpath(1:1)/=dirsep.and.trim(workingdir)/='') then
       process_path=trim(workingdir)//'/'//alpath
    else
       process_path=alpath
    end if

  end function process_path

end module glimmer_filenames
