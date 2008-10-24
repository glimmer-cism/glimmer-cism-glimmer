<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN"><HTML>
<HEAD>
<!-- hennerik CVSweb $Revision: 1.112 $ -->
<TITLE>glimmer-example/Makefile - view - 1.6</TITLE></HEAD>
<BODY BGCOLOR="#eeeeee">
<table width="100%" border=0 cellspacing=0 cellpadding=1 bgcolor="#9999ee"><tr valign=bottom><td><a href="Makefile?cvsroot=glimmer#rev1.6"><IMG SRC="/icons/cvsweb/back.gif" ALT="[BACK]" BORDER="0" WIDTH="20" HEIGHT="22"></a> <b>Return to <A HREF="Makefile?cvsroot=glimmer#rev1.6">Makefile</A>
 CVS log</b> <IMG SRC="/icons/cvsweb/text.gif" ALT="[TXT]" BORDER="0" WIDTH="20" HEIGHT="22"></td><td align=right><IMG SRC="/icons/cvsweb/dir.gif" ALT="[DIR]" BORDER="0" WIDTH="20" HEIGHT="22"> <b>Up to  <a href="/cgi-bin/cvsweb.cgi/?cvsroot=glimmer#dirlist">[glimmer]</a> / <a href="/cgi-bin/cvsweb.cgi/glimmer-example/?cvsroot=glimmer#dirlist">glimmer-example</a></b></td></tr></table><HR noshade><table width="100%"><tr><td bgcolor="#ffffff">File:  <a href="/cgi-bin/cvsweb.cgi/?cvsroot=glimmer#dirlist">[glimmer]</a> / <a href="/cgi-bin/cvsweb.cgi/glimmer-example/?cvsroot=glimmer#dirlist">glimmer-example</a> / <a href="/cgi-bin/cvsweb.cgi/glimmer-example/Makefile?cvsroot=glimmer">Makefile</a>&nbsp;(<A HREF="/cgi-bin/cvsweb.cgi/~checkout~/glimmer-example/Makefile?rev=1.6&amp;cvsroot=glimmer" target="cvs_checkout" onClick="window.open('/cgi-bin/cvsweb.cgi/~checkout~/glimmer-example/Makefile?rev=1.6','cvs_checkout','resizeable,scrollbars');"><b>download</b></A>)<BR>
Revision <B>1.6</B>, <i>Fri Jul 20 14:24:50 2007 UTC</i> (15 months ago) by <i>icrutt</i>
<BR>CVS Tags: <b>HEAD</b><BR>Changes since <b>1.5: +1 -1
 lines</b><PRE>
Prepared for 0.6 release
</PRE>
</td></tr></table><HR noshade><PRE># Makefile for building the example climate driver
# Magnus Hagdorn, November 2004

# setup compilation flags: tell the compiler where to look for f90 module files
GLIMMER_FFLAGS = `glimmer-config --fcflags`
# setup linker flags: tell the linker where to look for libraries and which libraries to link
GLIMMER_LIBRARIES = `glimmer-config --libs` -lglide -lglimmer -lslap -lnetcdf

# set variables FCFLAGS and LDFLAGS for extra libraries, e.g. netCDF

# how to build the example executable
example: &nbsp; &nbsp; &nbsp; &nbsp;example.o
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;$(FC) -o $@ $&lt; $(LDFLAGS) $(GLIMMER_LIBRARIES)

# implicit rule, telling the make system how to build an object file from a .f90 file
.f90.o:
 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;$(FC) -fpp $(FCFLAGS) $(GLIMMER_FFLAGS) -c $&lt; -o $@


.SUFFIXES: &nbsp; &nbsp; &nbsp;.f90 .o


# stuff below here is not too interesting. You might as well ignore it...
# build a distribution
DIST_FILES = ChangeLog example.config example.f90 Makefile README
VERSION = 0.6
DIST_NAME = glimmer-example-$(VERSION)
dist::
 &nbsp; &nbsp; &nbsp; &nbsp;mkdir $(DIST_NAME)
 &nbsp; &nbsp; &nbsp; &nbsp;cp $(DIST_FILES) $(DIST_NAME)
 &nbsp; &nbsp; &nbsp; &nbsp;tar cvzf $(DIST_NAME).tar.gz $(DIST_NAME)
 &nbsp; &nbsp; &nbsp; &nbsp;rm -r $(DIST_NAME)
</PRE>