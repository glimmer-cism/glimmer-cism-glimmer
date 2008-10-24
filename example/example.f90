<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN"><HTML>
<HEAD>
<!-- hennerik CVSweb $Revision: 1.112 $ -->
<TITLE>glimmer-example/example.f90 - view - 1.5</TITLE></HEAD>
<BODY BGCOLOR="#eeeeee">
<table width="100%" border=0 cellspacing=0 cellpadding=1 bgcolor="#9999ee"><tr valign=bottom><td><a href="example.f90?cvsroot=glimmer#rev1.5"><IMG SRC="/icons/cvsweb/back.gif" ALT="[BACK]" BORDER="0" WIDTH="20" HEIGHT="22"></a> <b>Return to <A HREF="example.f90?cvsroot=glimmer#rev1.5">example.f90</A>
 CVS log</b> <IMG SRC="/icons/cvsweb/text.gif" ALT="[TXT]" BORDER="0" WIDTH="20" HEIGHT="22"></td><td align=right><IMG SRC="/icons/cvsweb/dir.gif" ALT="[DIR]" BORDER="0" WIDTH="20" HEIGHT="22"> <b>Up to  <a href="/cgi-bin/cvsweb.cgi/?cvsroot=glimmer#dirlist">[glimmer]</a> / <a href="/cgi-bin/cvsweb.cgi/glimmer-example/?cvsroot=glimmer#dirlist">glimmer-example</a></b></td></tr></table><HR noshade><table width="100%"><tr><td bgcolor="#ffffff">File:  <a href="/cgi-bin/cvsweb.cgi/?cvsroot=glimmer#dirlist">[glimmer]</a> / <a href="/cgi-bin/cvsweb.cgi/glimmer-example/?cvsroot=glimmer#dirlist">glimmer-example</a> / <a href="/cgi-bin/cvsweb.cgi/glimmer-example/example.f90?cvsroot=glimmer">example.f90</a>&nbsp;(<A HREF="/cgi-bin/cvsweb.cgi/~checkout~/glimmer-example/example.f90?rev=1.5&amp;cvsroot=glimmer" target="cvs_checkout" onClick="window.open('/cgi-bin/cvsweb.cgi/~checkout~/glimmer-example/example.f90?rev=1.5','cvs_checkout','resizeable,scrollbars');"><b>download</b></A>)<BR>
Revision <B>1.5</B>, <i>Mon Oct  3 13:48:31 2005 UTC</i> (3 years ago) by <i>magi</i>
<BR>CVS Tags: <b>RELEASE_0_4, HEAD</b><BR>Changes since <b>1.4: +3 -2
 lines</b><PRE>
use new glide API
</PRE>
</td></tr></table><HR noshade><PRE>! ******************************************************************************
! example.f90
! Magnus Hagdorn
!
! simple example climate driver demonstrating how to use the library
! ******************************************************************************
!
! ChangeLog
! 2004-11-10 Magnus Hagdorn
! &nbsp;* initial version

program example

 &nbsp;! load various modules
 &nbsp;use glimmer_global, only:rk ! precision of the model
 &nbsp;use glide &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; ! main glide module
 &nbsp;use glimmer_log &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; ! module for logging messages
 &nbsp;use glimmer_config &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;! module for handling configuration files
 &nbsp;implicit none

 &nbsp;! some variables
 &nbsp;type(glide_global_type) :: model &nbsp; &nbsp;
 &nbsp;! this variable holds all the data associated with this particular model
 &nbsp;! instance.
 &nbsp;type(ConfigSection), pointer :: config
 &nbsp;! this pointer points to the first element of a linked list which contains
 &nbsp;! all the configuration variables and their values
 &nbsp;character(len=50) :: fname &nbsp; 
 &nbsp;! name of paramter file
 &nbsp;real(kind=rk) time
 &nbsp;! current time
 &nbsp;
 &nbsp;integer ew,ns,ewct,nsct &nbsp;! loop variable and grid centre
 &nbsp;real grid, dist &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;! node spacing and distance

 &nbsp;! local ice sheet variables
 &nbsp;real, allocatable, dimension(:,:) :: ice_thickness
 &nbsp;real, allocatable, dimension(:,:) :: mass_balance
 &nbsp;real, allocatable, dimension(:,:) :: surface_temperature

 &nbsp;! Ask for configuration file
 &nbsp;write(*,*) 'Enter name of GLIDE configuration file to be read'
 &nbsp;read(*,*) fname
 &nbsp;
 &nbsp;! start logging
 &nbsp;call open_log(unit=50)
 &nbsp;
 &nbsp;! read configuration
 &nbsp;call ConfigRead(fname,config)
 &nbsp;call CheckSections(config)
 &nbsp;call glide_config(model,config)

 &nbsp;! initialise GLIDE
 &nbsp;call glide_initialise(model)
 &nbsp;! fill dimension variables
 &nbsp;call glide_nc_fillall(model)
 &nbsp;! get current time from start time
 &nbsp;time = get_tstart(model)
 &nbsp;
 &nbsp;! allocate variables
 &nbsp;allocate(ice_thickness(get_ewn(model),get_nsn(model)))
 &nbsp;allocate(mass_balance(get_ewn(model),get_nsn(model)))
 &nbsp;allocate(surface_temperature(get_ewn(model),get_nsn(model)))

 &nbsp;! setup some variables for BC
 &nbsp;ewct = real(get_ewn(model)+1) / 2.0 ! the grid centre (x)
 &nbsp;nsct = real(get_nsn(model)+1) / 2.0 ! and (y)
 &nbsp;grid = get_dew(model) &nbsp; &nbsp; &nbsp; ! node-spacing

 &nbsp;! loop over times
 &nbsp;do while(time.le.model%numerics%tend)
 &nbsp; &nbsp; ! setup boundary conditions
 &nbsp; &nbsp; ! this example implements the EISMINT-1 moving margin BC
 &nbsp; &nbsp; 
 &nbsp; &nbsp; ! mass balance
 &nbsp; &nbsp; ! loop over grid
 &nbsp; &nbsp; do ns = 1,model%general%nsn
 &nbsp; &nbsp; &nbsp; &nbsp;do ew = 1,model%general%ewn
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; ! calculate distance from centre
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; dist = grid * sqrt((real(ew) - ewct)**2 + (real(ns) - nsct)**2)
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; ! calculate mass balance
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; mass_balance(ew,ns) = min(0.5, 1.05e-5 * (450.0e3 - dist))
 &nbsp; &nbsp; &nbsp; &nbsp;end do
 &nbsp; &nbsp; end do
 &nbsp; &nbsp; ! and set it
 &nbsp; &nbsp; call glide_set_acab(model,mass_balance)
 &nbsp; &nbsp; 
 &nbsp; &nbsp; ! get ice thickness for temperature calculations
 &nbsp; &nbsp; call glide_get_thk(model,ice_thickness)
 &nbsp; &nbsp; ! calculate surface temperature
 &nbsp; &nbsp; surface_temperature = -3.150 &nbsp;-1.e-2 * ice_thickness
 &nbsp; &nbsp; ! and set it
 &nbsp; &nbsp; call glide_set_artm(model,surface_temperature)

 &nbsp; &nbsp; ! calculate temperature and velocity distribution
 &nbsp; &nbsp; call glide_tstep_p1(model,time)
 &nbsp; &nbsp; ! write to netCDF file, move ice
 &nbsp; &nbsp; call glide_tstep_p2(model)
 &nbsp; &nbsp; ! calculate isostatic adjustment
 &nbsp; &nbsp; call glide_tstep_p3(model)
 &nbsp; &nbsp; ! increment time counter
 &nbsp; &nbsp; time = time + get_tinc(model)
 &nbsp;end do

 &nbsp;! finalise GLIDE
 &nbsp;call glide_finalise(model)
end program example
</PRE>