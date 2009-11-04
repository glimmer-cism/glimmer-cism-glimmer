program iswplume

  use plume

  implicit none

  integer::iarg_count
  integer::command_argument_count
  

  ! *******************************************************************
  ! *** process command-line arguments ********************************
  ! *******************************************************************

  iarg_count = command_argument_count()
	
  if (iarg_count > 0) then
     call get_command_argument(1,jobid)
     if (jobid == "-h") then
        write(*,*)'Usage: plume <jobid> <namelist_file> <output_dir>'
        stop
     end if
  else
     !get job id
     write(*,*) 'please enter 3-character job id:'
     read(*,*) jobid
  end if

  if (iarg_count > 1) then
     call get_command_argument(2,nl_filename)
  else
     nl_filename = 'plume.nl'
  end if

  if (iarg_count > 2) then
     call get_command_argument(3,output_dir)
  else
     output_dir = './'
  end if


  call plume_initialize()

  ! main loop 

  do istep = 1,nsteps

     call plume_runstep()

  end do

  call plume_finalize()


end program iswplume
