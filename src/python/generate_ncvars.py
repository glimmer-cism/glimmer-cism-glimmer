#! /usr/bin/env python

# Copyright 2004, Magnus Hagdorn
# 
# This file is part of glimmer.
# 
# PyGMT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# PyGMT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyGMT; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# python script used to generate source code files given a variable definition file

import ConfigParser, sys, time, string

class Variables(dict):
    """Dictionary containing variable definitions."""

    def __init__(self,filename):
        """Initialise Variable class.

        filename: name of file containing variable definitions."""

        dict.__init__(self)

        # reading variable configuration file
        vars = ConfigParser.ConfigParser()
        vars.readfp(open(filename))

        for v in vars.sections():
            vardef = {}
            vardef['name'] = v
            for (name, value) in vars.items(v):
                vardef[name] = value
            self.__setitem__(v,vardef)

    def keys(self):
        """Reorder standard keys alphabetically."""
        dk = []
        vk = []
        for v in dict.keys(self):
            if v == self.__getitem__(v)['dimensions']:
                dk.append(v)
            else:
                vk.append(v)
        dk.sort()
        vk.sort()
        return dk+vk

class PrintVars:
    """Base class for printing variables."""
    canhandle = None
    comment = '!'

    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""
        if filename != '%s.in'%self.canhandle:
            raise NotImplementedError, 'Can only handle %s'%self.canhandle

        self.infile = open("%s.in"%self.canhandle,'r')
        self.stream = open(self.canhandle,'w')

        self.handletoken = {'!GENVARS!' : self.print_var}

    def print_warning(self):
        """Write a warning message to stream"""

        self.stream.write("%s\n"%(80*self.comment))
        self.stream.write("%s WARNING: this file was automatically generated on\n! %s\n! from %s.in\n"%(self.comment,time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()), self.canhandle))
        self.stream.write("%s\n\n"%(80*self.comment))
        
    def print_var(self, var):
        """Template for writing single variable"""

        raise NotImplementedError, 'You should use one of the derived classes'

    def write(self,vars):
        """Merge file with definitions"""

        self.print_warning()
        for l in self.infile.readlines():
            for token in self.handletoken:
                if string.find(l,token) is not -1:
                    break
            if string.find(l,token) is not -1:
                for v in vars.keys():
                    self.handletoken[token](vars[v])
            else:
                self.stream.write("%s"%l)
        self.infile.close()
        self.stream.close()

class PrintNCDF(PrintVars):
    """Process ncdf.f90"""
    canhandle = 'ncdf.f90'

    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        PrintVars.__init__(self,filename)

        self.handletoken['!GENVAR_TYPES!'] = self.print_var_types

    def print_var_types(self,var):
        """Write variable id to stream."""

        if var['dimensions'] != var['name']:
            self.stream.write("  integer, parameter :: %s = %d ! %s\n"%(var_type(var),self.thisvar,var['long_name']))
            self.thisvar = self.thisvar + 1

    def print_var(self, var):
        """Write single variable block to stream for ncdf."""

        # skip variables associated with dimension 
        if var['dimensions'] != var['name']:
            self.stream.write("    if (nc%%do_var(%s)) then\n"%var_type(var))
            self.stream.write("       write(unit,*) '%s'\n"%var['name'])
            self.stream.write("    end if\n")

    def write(self,vars):
        """Merge ncdf.F90.in with definitions."""

        numvars = 0
        for v in vars:
            if vars[v]['dimensions'] != v:
                numvars = numvars + 1

        self.thisvar = 1

        self.print_warning()
        for l in self.infile.readlines():
            for token in self.handletoken:
                if string.find(l,token) is not -1:
                    break
            if string.find(l,token) is not -1:
                for v in vars.keys():
                    self.handletoken[token](vars[v])
            elif string.find(l,'!GENVARS_NUMVARS!'):
                self.stream.write("%s"%l.replace('!GENVARS_NUMVARS!','%d'%numvars))
            else:
                self.stream.write("%s"%l)
        self.infile.close()
        self.stream.close()


class PrintNCDF_FILE(PrintVars):
    """Process ncdf_file.f90"""
    canhandle = 'ncdf_file.f90'

    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        PrintVars.__init__(self,filename)

        self.handletoken['!GENVAR_WRITE!'] = self.print_var_write

    def print_var(self, var):
        """Write single variable block to stream for ncdf_file."""

        dims = string.split(var['dimensions'],',')
        dims.reverse()
        dimstring = 'NC%%%sdim'%dims[0].strip()
        for d in dims[1:]:
            dimstring = '%s, NC%%%sdim'%(dimstring,d.strip())
        
        self.stream.write("    !     %s -- %s\n"%(var['name'],var['long_name']))
        spaces = 0
        if var['dimensions'] != var['name']:
            spaces=3
            self.stream.write("    if (NC%%do_var(%s)) then\n"%(var_type(var)))
            id = 'NC%%varids(%s)'%var_type(var)
        else:
            id = 'NC%%%svar'%var['name']
        self.stream.write("%s    write(*,*) 'Creating variable %s'\n"%(spaces*' ',var['name']))
        self.stream.write("%s    status = nf90_def_var(NC%%id,'%s',NF90_FLOAT,(/%s/),%s)\n"%(spaces*' ',
                                                                                             var['name'],
                                                                                             dimstring,
                                                                                             id
                                                                                             ))
        self.stream.write("%s    call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces*' '))
        for attrib in var:
            if attrib not in ['name','dimensions','data','factor']:
                self.stream.write("%s    status = nf90_put_att(NC%%id, %s, '%s', &\n%s           '%s')\n"%(spaces*' ',
                                                                                                           id,
                                                                                                           attrib,
                                                                                                           spaces*' ',
                                                                                                           var[attrib]))
        if var['dimensions'] != var['name']:
            self.stream.write("    end if\n")
        self.stream.write("\n")

    def print_var_write(self,var):
        """Write single variable block to stream for ncdf_file."""

        # skip variables associated with dimension 
        if var['dimensions'] != var['name']:
            dims = string.split(var['dimensions'],',')
            dims.reverse()
            for i in range(0,len(dims)):
                dims[i] = dims[i].strip()
            self.stream.write("    if (NC%%do_var(%s)) then\n"%(var_type(var)))

            dimstring = ''
            spaces = ''
            for i in range(0,len(dims)):
                if i > 0:
                    dimstring = dimstring + ','
                if dims[i] == 'time':
                    dimstring = dimstring + 'outfile%timecounter'
                elif dims[i] == 'level':
                    dimstring = dimstring + 'up'
                else:
                    dimstring = dimstring + '1'
            
            if  'level' in dims:
                # handle 3D fields
                spaces = ' '*3
                self.stream.write("       do up=1,model%general%upn\n")

                        
            if 'factor' in var:
                data = '(%s)*(%s)'%(var['factor'], var['data'])
            else:
                data = var['data']
            self.stream.write("%s       status = nf90_put_var(NC%%id, %s, &\n%s            %s, (/%s/))\n"%(spaces,'NC%%varids(%s)'%var_type(var),
                                                                                                         spaces,data, dimstring))
            self.stream.write("%s       call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces))

            if  'level' in dims:
                self.stream.write("       end do\n")
            self.stream.write("    end if\n\n")

            
def var_type(var):
    """Map variable to type parameter."""
    return 'NC_B_%s'%string.upper(var['name'])

def usage():
    """Short help message."""

    print 'Usage generate_ncvars.py vardef [outfile0.in [,... [,outfileN.in]]]'
    print 'generate source code files given a variable definition file'
    print ''
    print 'vardef: file containing variable definition'
    print 'outfile.in: output template to be processed'
    print 'print variables if no templates are given'

HandleFile = {'ncdf.f90.in' : PrintNCDF, 'ncdf_file.f90.in' : PrintNCDF_FILE}

if __name__ == '__main__':

    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    vars = Variables(sys.argv[1])

    if len(sys.argv) == 2:
        for v in vars.keys():
            print v
            for o in vars[v]:
                print '%s: %s'%(o, vars[v][o])
            print ''
    else:
        for f in sys.argv[2:]:
            if f in HandleFile:
                handle = HandleFile[f](f)
                handle.write(vars)
                     
