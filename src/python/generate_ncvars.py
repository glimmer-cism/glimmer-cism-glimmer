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

import ConfigParser, sys, time, string,re, os.path

NOATTRIB = ['name','dimensions','data','factor','load','f90file','hot','type']
VAR1D = ['eus']
hotvars = []

def spotname(name):
    """Return name of spotvariable."""
    return "%s_spot"%name

def isspot(var):
    """Return True if variable is a spot variable."""
    return '_spot' in var['name']

def is_dimvar(var):
    """Return True if variable is associated with a dimension.

    this is assumed to be the case if no time dim is present
    """

    if len(string.split(var['dimensions'],',')) == 1 and var['name'] not in VAR1D:
        return True
    else:
        return False

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
            if 'hot' in vardef:
                if vardef['hot'].lower() in ['1','true','t']:
                    hotvars.append(v)
                    vardef['load'] = '1'
            if 'type' not in vardef:
                vardef['type'] = 'float'
            self.__setitem__(v,vardef)
            self.__add_spot(vardef)

    def __add_spot(self,vdef):
        """Add spot variable.

        vname: name of variable
        vdef:  variable definition"""

        if 'time' not in vdef['dimensions']:
            return

        spdef = {}
        for k in vdef:
            if k=='name':
                spdef[k] = spotname(vdef[k])
            elif k=='dimensions':
                search = re.search('y[0-1]\s*,\s*x[0-1]', vdef[k])
                if search!=None:
                    spdef[k] = vdef[k][:search.start()] + 'spot' + vdef[k][search.end():]
                else:
                    return
            else:
                spdef[k] = vdef[k]
        if 'x0' in vdef['dimensions']:
            spdef['coordinates'] = 'y0_spot x0_spot'
        else:
            spdef['coordinates'] = 'y1_spot x1_spot'
        self.__setitem__(spdef['name'],spdef)

    def keys(self):
        """Reorder standard keys alphabetically."""
        dk = []
        vk = []
        for v in dict.keys(self):
            if is_dimvar(self.__getitem__(v)):
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
        self.stream.write("%s WARNING: this file was automatically generated on\n%s %s\n%s from %s.in\n"%(self.comment,
                                                                                                          self.comment,time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()),
                                                                                                          self.comment, self.canhandle))
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

class PrintNCDF_BASEIO(PrintVars):
    """Class for defining custom I/O blocks."""

    def check_write(self,var):
        """Return True if we should write the variable."""
        return False
    def check_read(self,var):
        """Return True if we should read the variable."""
        return False
    def print_var_write(self,var):
        """Write single variable block to stream for ncdf_file."""

        # skip variables associated with dimension 
        if not is_dimvar(var) and self.check_write(var):
            dims = string.split(var['dimensions'],',')
            dims.reverse()
            for i in range(0,len(dims)):
                dims[i] = dims[i].strip()
            self.stream.write("    if (NCO%%do_var(%s)) then\n"%(var_type(var)))
            
            if isspot(var):
                if len(dims)==2:
                    if var['name']=='btemp_spot':
                        dimstring_in='model%general%upn,outfile%spotx(spot),outfile%spoty(spot)'
                    else:
                        dimstring_in='outfile%spotx(spot),outfile%spoty(spot)'
                    dimstring_out='(/spot,outfile%timecounter/)'
                elif len(dims)==3:
                    dimstring_in=':,outfile%spotx(spot),outfile%spoty(spot)'
                    dimstring_out='(/spot,1,outfile%timecounter/),(/1,model%general%upn,1/)'

                self.stream.write("       do spot=1,size(outfile%spotx)\n")
                pos = var['data'].find('(')
                if pos>-1:
                    data=var['data'][:pos]
                else:
                    data=var['data']
                data = '%s(%s)'%(data,dimstring_in)
                if 'factor' in var:
                    data = '(%s)*(%s)'%(var['factor'], data)
                self.stream.write("          status = nf90_put_var(NCO%%id, NCO%%varids(%s), &\n"%var_type(var))
                self.stream.write("               %s, &\n"%data)
                self.stream.write("               %s)\n"%dimstring_out)
                self.stream.write("          call nc_errorhandle(__FILE__,__LINE__,status)\n")
                self.stream.write("       end do\n")
                
            else:
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
                self.stream.write("%s       status = nf90_put_var(NCO%%id, %s, &\n%s            %s, (/%s/))\n"%(spaces,'NCO%%varids(%s)'%var_type(var),
                                                                                                             spaces,data, dimstring))
                self.stream.write("%s       call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces))

                if  'level' in dims:
                    self.stream.write("       end do\n")
                # remove self since it's not time dependent
                if 'time' not in dims:
                    self.stream.write("       NCO%%do_var(%s) = .False.\n"%(var_type(var)))
            self.stream.write("    end if\n\n")    

    def print_var_read(self, var):
        """Write single variable block to stream for ncdf_infile."""

        if 'load' in var and not isspot(var) and not is_dimvar(var) and self.check_read(var):
            if var['load'].lower() in ['1','true','t']:
                dims = string.split(var['dimensions'],',')
                dims.reverse()
                for i in range(0,len(dims)):
                    dims[i] = dims[i].strip()
                self.stream.write("    if (NCI%%do_var(%s)) then\n"%(var_type(var)))
                self.stream.write("       call write_log('  Loading %s')\n"%var['name'])
                dimstring = ''
                spaces = ''
                for i in range(0,len(dims)):
                    if i > 0:
                        dimstring = dimstring + ','
                    if dims[i] == 'time':
                        dimstring = dimstring + 'infile%current_time'
                    elif dims[i] == 'level':
                        dimstring = dimstring + 'up'
                    else:
                        dimstring = dimstring + '1'

                if  'level' in dims:
                    # handle 3D fields
                    spaces = ' '*3
                    self.stream.write("       do up=1,model%general%upn\n")

                self.stream.write("%s       status = nf90_get_var(NCI%%id, %s, &\n%s            %s, (/%s/))\n"%(spaces,'NCI%%varids(%s)'%var_type(var),
                                                                                                               spaces,var['data'], dimstring))
                self.stream.write("%s       call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces))
                if 'factor' in var:
                    self.stream.write("%s       if (scale) %s = %s/(%s)\n"%(spaces,var['data'],var['data'],var['factor']))

                if  'level' in dims:
                    self.stream.write("       end do\n")
                # remove self since it's not time dependent
                if 'time' not in dims:
                    self.stream.write("       NCI%%do_var(%s) = .False.\n"%(var_type(var)))
                self.stream.write("    end if\n\n")
                

class PrintDoc(PrintVars):
    """Process varlist.tex"""
    canhandle = 'varlist.tex'
    comment = '%'

    def print_var(self, var):
        """Write single variable block to stream for ncdf_params."""

        # skip variables associated with dimension 
        if not is_dimvar(var) and not isspot(var):
            load = ''
            if 'load' in var:
                if var['load'].lower() in ['1','true','t']:
                    load = '$^\\ast$'

            self.stream.write("\\texttt{%s}%s & %s & %s\\\\\n"%(var['name'],load,var['long_name'],
                                                    var['units'].replace('_','\_')))
            if 'standard_name' in var:
                self.stream.write("&CF name: \\texttt{%s}&\\\\\n"%(var['standard_name'].replace('_','\_')))
            self.stream.write("\\hline\n")

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

        if not is_dimvar(var):
            self.stream.write("  integer, parameter :: %s = %d ! %s\n"%(var_type(var),self.thisvar,var['long_name']))
            self.thisvar = self.thisvar + 1

    def print_var(self, var):
        """Write single variable block to stream for ncdf."""

        # skip variables associated with dimension 
        if not is_dimvar(var):
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

class PrintNCDF_PARAMS(PrintVars):
    """Process ncdf_params.f90"""
    canhandle = 'ncdf_params.f90'

    def write(self,vars):
        """Merge ncdf.F90.in with definitions."""

        self.print_warning()
        for l in self.infile.readlines():
            for token in self.handletoken:
                if string.find(l,token) is not -1:
                    break
            if string.find(l,token) is not -1:
                for v in vars.keys():
                    self.handletoken[token](vars[v])
            elif string.find(l,'!GENVARS_HOT!') is not -1:
                self.stream.write("    if (index(vars,' hot ').ne.0) then\n")
                for v in hotvars:
                    self.stream.write("       handle_output%%nc%%do_var(%s) = .true.\n"%var_type(vars[v]))
                self.stream.write("    end if\n\n")
            else:
                self.stream.write("%s"%l)
        self.infile.close()
        self.stream.close()

    def print_var(self, var):
        """Write single variable block to stream for ncdf_params."""

        # skip variables associated with dimension 
        if not is_dimvar(var):
            if isspot(var):
                self.stream.write("    if (index(vars,' %s ').ne.0 .and. handle_output%%nc%%do_spot) then\n"%(var['name']))
            else:
                self.stream.write("    if (index(vars,' %s ').ne.0) then\n"%(var['name']))
            self.stream.write("       handle_output%%nc%%do_var(%s) = .true.\n"%var_type(var))
            self.stream.write("    end if\n\n")
        

class PrintNCDF_FILE(PrintNCDF_BASEIO):
    """Process ncdf_file.f90"""
    canhandle = 'ncdf_file.f90'

    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        PrintVars.__init__(self,filename)

        self.handletoken['!GENVAR_WRITE!'] = self.print_var_write

    def check_write(self,var):
        """Return True if we should write the variable."""

        dowrite = True
        if 'f90file' in var:
            dowrite = var['f90file'] == self.canhandle
        return dowrite

    def print_var(self, var):
        """Write single variable block to stream for ncdf_file."""

        dims = string.split(var['dimensions'],',')
        dims.reverse()
        dimstring = 'NCO%%%sdim'%dims[0].strip()
        for d in dims[1:]:
            dimstring = '%s, NCO%%%sdim'%(dimstring,d.strip())
        
        self.stream.write("    !     %s -- %s\n"%(var['name'],var['long_name']))
        spaces = 0
        if not is_dimvar(var):
            spaces=3
            self.stream.write("    if (NCO%%do_var(%s)) then\n"%(var_type(var)))
            id = 'NCO%%varids(%s)'%var_type(var)
        else:
            if 'spot' in var['dimensions']:
                spaces=3
                self.stream.write("    if (NCO%do_spot) then\n")
            id = 'NCO%%%svar'%var['name']
        self.stream.write("%s    call write_log('Creating variable %s')\n"%(spaces*' ',var['name']))
        self.stream.write("%s    status = nf90_def_var(NCO%%id,'%s',NF90_%s,(/%s/),%s)\n"%(spaces*' ',
                                                                                           var['name'],
                                                                                           var['type'].upper(),
                                                                                           dimstring,
                                                                                           id
                                                                                           ))
        self.stream.write("%s    call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces*' '))
        for attrib in var:
            if attrib not in NOATTRIB:
                self.stream.write("%s    status = nf90_put_att(NCO%%id, %s, '%s', &\n%s           '%s')\n"%(spaces*' ',
                                                                                                           id,
                                                                                                           attrib,
                                                                                                           spaces*' ',
                                                                                                           var[attrib]))
        if not is_dimvar(var) and 'spot' not in var['dimensions']:
            self.stream.write("%s    if (CFproj_allocated(model%%projection)) then\n"%(spaces*' '))
            self.stream.write("%s       status = nf90_put_att(NCO%%id, %s, 'grid_mapping',mapvarname)\n"%(spaces*' ',id))
            self.stream.write("%s    end if\n"%(spaces*' '))
        if not is_dimvar(var) or 'spot' in var['dimensions']:
            self.stream.write("    end if\n")
        self.stream.write("\n")

class PrintNCDF_INFILE(PrintNCDF_BASEIO):
    """Process ncdf_infile.f90"""
    canhandle = 'ncdf_infile.f90'

    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        PrintVars.__init__(self,filename)

        self.handletoken['!GENVAR_READ!'] = self.print_var_read

    def check_read(self,var):
        """Return True if we should read the variable."""

        doread = True
        if 'f90file' in var:
            doread = var['f90file'] == self.canhandle
        return doread

    def print_var(self, var):
        """Write single variable block to stream for ncdf_infile."""

        if 'load' in var and not isspot(var) and not is_dimvar(var):
            if var['load'].lower() in ['1','true','t']:
                self.stream.write("       case('%s')\n"%var['name'])
                self.stream.write("          NCI%%do_var(%s) = .true.\n"%var_type(var))
                self.stream.write("          NCI%%varids(%s) = i\n"%var_type(var))

class PrintNCDF_IO(PrintNCDF_BASEIO):
    """Process customised I/O files."""

    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        self.fname = os.path.splitext(filename)[0]
        self.infile = open("%s.in"%self.fname,'r')
        self.stream = open(self.fname,'w')

        self.handletoken = {}
        self.handletoken['!GENVAR_WRITE!'] = self.print_var_write
        self.handletoken['!GENVAR_READ!'] = self.print_var_read

    def check_write(self,var):
        """Return True if we should write the variable."""

        dowrite = False
        if 'f90file' in var:
            dowrite = var['f90file'] == self.fname
        return dowrite    

    def check_read(self,var):
        """Return True if we should read the variable."""

        doread = False
        if 'f90file' in var:
            doread = var['f90file'] == self.fname
        return doread    
            
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

HandleFile={}
HandleFile['ncdf.f90.in'] = PrintNCDF
HandleFile['ncdf_file.f90.in'] = PrintNCDF_FILE
HandleFile['ncdf_infile.f90.in'] = PrintNCDF_INFILE
HandleFile['ncdf_params.f90.in'] = PrintNCDF_PARAMS
HandleFile['varlist.tex.in'] = PrintDoc

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
            else:
                handle = PrintNCDF_IO(f)
                handle.write(vars)
                     
