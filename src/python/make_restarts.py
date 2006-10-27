
import getopt,sys,time
from fort_interp import *
from fort_parser import *

marker='make_restart'
restart_marker='!'+marker

usefiles='      use glimmer_restart_common\n'
usefiles=usefiles+'      use glimmer_restart_statscal\n'
usefiles=usefiles+'      use glimmer_restart_statarr\n'
usefiles=usefiles+'      use glimmer_restart_pointarr\n'

# Top-level code for generating restarts

def check_restart(fname):
    f=open(fname)
    lines=f.readlines()
    for l in lines:
        if l.lower().lstrip().find(restart_marker)!=-1:
            return True
    return False

def get_re_fname(mod):
    for l in mod.comments:
        pos=l.lower().lstrip().find(marker)
        if pos!=-1:
            return l[pos+len(marker):].lstrip()
    print 'ERROR: no restart output name found'
    sys.exit(2)

def write_type_name(var):
    if var.pointer:
        return 'rsw_'+var.type+'_point'
    else:
        return 'rsw_'+var.type

def read_type_name(var):
    if var.pointer:
        return 'rsr_'+var.type+'_point'
    else:
        return 'rsr_'+var.type

def get_no_restart(a):
    out=[]
    for c in a.comments:
        if c.lstrip().find('no_restart')!=-1:
            out=out+c.split()[1:]
    return out

def var_write(var,file,prefix,dtype=''):
    varname=dtype+''+var.name
    if var.parameter: return
    if var.type in types:
        if var.array:
            if var.pointer:
                file.write("      call write_pointarr(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
            elif var.allocatable:
                print 'ERROR: Allocatable arrays not supported yet'
                sys.exit(2)
                #file.write("      call write_allocarr(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
            else:
                file.write("      call write_statarr(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
        else:
            if var.pointer:
                print 'ERROR: Pointer scalars not supported yet'
                sys.exit(2)
                #file.write("      call write_pointscalvar(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
            else:
                file.write("      call write_statscalvar(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
    else:
        if var.array:
            if var.pointer:
                file.write("      call rsw_%s_pointarr(file,'%s','%s',%s)\n" % (var.type,prefix,var.name,varname))
            else:
                print 'ERROR: Derived-type arrays not supported yet',var.Print()
                sys.exit(2)
        else:
            if var.pointer:
                file.write("      call write_pointer(file,'%s','%s',associated(%s))\n" % (prefix,var.name,varname))
                file.write("      if (associated(%s)) call %s(file,%s)\n" % (varname,write_type_name(var),varname))
            else:
                file.write("      call %s(file,%s)\n" % (write_type_name(var),varname))

def var_read(var,file,prefix,dtype=''):
    varname=dtype+''+var.name
    if var.parameter: return
    if var.type in types:
        if var.array:
            if var.pointer:
                file.write("      call read_pointarr(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
            elif var.allocatable:
                print 'ERROR: Allocatable arrays not supported yet'
                sys.exit(2)
                #file.write("      call read_allocarr(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
            else:
                file.write("      call read_statarr(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
        else:
            if var.pointer:
                print 'ERROR: Pointer scalars not supported yet'
                sys.exit(2)
                #file.write("      call read_pointscalvar(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
            else:
                file.write("      call read_statscalvar(file,'%s','%s',%s)\n" % (prefix,var.name,varname))
    else:
        if var.array:
            if var.pointer:
                file.write("      call rsr_%s_pointarr(file,'%s','%s',%s)\n" % (var.type,prefix,var.name,varname))
            else:
                print 'ERROR: Derived-type arrays not supported yet',var.Print()
                sys.exit(2)
        else:
            if var.pointer:
                file.write("      call read_pointer(file,'%s','%s',assoc)\n" % (prefix,var.name))
                file.write("      if (assoc) call %s(file,%s)\n" % (read_type_name(var),varname))
            else:
                file.write("      call %s(file,%s)\n" % (read_type_name(var),varname))

def write_common_restart(headfile,bodyfile):
    headfile.write('!++++++++++++++++++++++++++++++++++++++++\n')
    headfile.write('! WARNING: this file was automatically generated on\n! %s\n'
                   %(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())))
    headfile.write('!++++++++++++++++++++++++++++++++++++++++\n')
    bodyfile.write('!++++++++++++++++++++++++++++++++++++++++\n')
    bodyfile.write('! WARNING: this file was automatically generated on\n! %s\n'
                   %(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())))
    bodyfile.write('!++++++++++++++++++++++++++++++++++++++++\n')

def write_mod_restart(mod,headfile,bodyfile):
    writesubname=mod.name+'_modrsw'
    readsubname=mod.name+'_modrsr'
    no_restart=get_no_restart(mod)
    # Make public
    headfile.write('   public :: %s,%s\n' % (writesubname,readsubname))
    # Write code
    bodyfile.write('   subroutine %s(file)\n' % (writesubname))
    bodyfile.write('\n')
    bodyfile.write(usefiles)
    bodyfile.write('      implicit none\n')
    bodyfile.write('\n')
    bodyfile.write('      type(restart_file),intent(inout) :: file\n')
    bodyfile.write('\n')
    bodyfile.write('      call check_write(file)\n')    
    bodyfile.write('\n')
    for var in mod.modvars:
        if not (var.name in no_restart):
            var_write(var,bodyfile,mod.name+'mod')
    bodyfile.write('\n')
    bodyfile.write('   end subroutine %s\n' % (writesubname))
    bodyfile.write('\n')
    # Read code
    bodyfile.write('   subroutine %s(file)\n' % (readsubname))
    bodyfile.write('\n')
    bodyfile.write(usefiles)
    bodyfile.write('      implicit none\n')
    bodyfile.write('\n')
    bodyfile.write('      type(restart_file),intent(inout) :: file\n')
    bodyfile.write('\n')
    bodyfile.write('      call check_read(file)\n')    
    bodyfile.write('\n')
    for var in mod.modvars:
        if not (var.name in no_restart):
            var_read(var,bodyfile,mod.name+'mod')
    bodyfile.write('\n')
    bodyfile.write('   end subroutine %s\n' % (readsubname))
    bodyfile.write('\n')

def write_type_restart(t,headfile,bodyfile):
    writesubname = 'rsw_'+t.name
    readsubname  = 'rsr_'+t.name
    wpointn=writesubname+'_point'
    rpointn=readsubname+ '_point'
    wparrn_1d=writesubname+'_pointarr_1d'
    rparrn_1d=readsubname +'_pointarr_1d'
    wparrn_2d=writesubname+'_pointarr_2d'
    rparrn_2d=readsubname +'_pointarr_2d'
    wparrn_3d=writesubname+'_pointarr_3d'
    rparrn_3d=readsubname +'_pointarr_3d'
    no_restart=get_no_restart(t)
    # Make public
    headfile.write('   public :: %s,%s,%s,%s,%s,%s\n' % (writesubname,readsubname,wpointn
                                                      ,rpointn,writesubname+'_pointarr',readsubname+'_pointarr'))
    headfile.write('\n')
    headfile.write('interface %s\n'%(writesubname+'_pointarr'))
    headfile.write('    module procedure %s\n'%(wparrn_1d))
    headfile.write('end interface\n')
    headfile.write('\n')
    headfile.write('interface %s\n'%(readsubname+'_pointarr'))
    headfile.write('    module procedure %s\n'%(rparrn_1d))
    headfile.write('end interface\n')
    headfile.write('\n')
    # Write code (scalar)
    bodyfile.write('   recursive subroutine %s(file,dat)\n' % (writesubname))
    bodyfile.write('\n')
    bodyfile.write(usefiles)
    bodyfile.write('      implicit none\n')
    bodyfile.write('\n')
    bodyfile.write('      type(restart_file),intent(inout) :: file\n')
    bodyfile.write('      type(%s)    :: dat\n'%(t.name))
    bodyfile.write('      logical :: assoc\n')
    bodyfile.write('\n')
    bodyfile.write('      call check_write(file)\n')    
    bodyfile.write('\n')
    for var in t.components:
        if not (var.name in no_restart):
            var_write(var,bodyfile,t.name,'dat%')
    bodyfile.write('\n')
    bodyfile.write('   end subroutine %s\n' % (writesubname))
    bodyfile.write('\n')
    # Write code (scalar pointer)
    bodyfile.write('   recursive subroutine %s(file,dat)\n' % (wpointn))
    bodyfile.write('\n')
    bodyfile.write(usefiles)
    bodyfile.write('      implicit none\n')
    bodyfile.write('\n')
    bodyfile.write('      type(restart_file),intent(inout) :: file\n')
    bodyfile.write('      type(%s),pointer    :: dat\n'%(t.name))
    bodyfile.write('      logical :: assoc\n')
    bodyfile.write('\n')
    bodyfile.write('      call check_write(file)\n')    
    bodyfile.write('\n')
    for var in t.components:
        if not (var.name in no_restart):
            var_write(var,bodyfile,t.name,'dat%')
    bodyfile.write('\n')
    bodyfile.write('   end subroutine %s\n' % (wpointn))
    bodyfile.write('\n')
    # Write code (1D pointer array)
    bodyfile.write('recursive subroutine %s(file,prefix,name,values)\n'%(wparrn_1d))
    bodyfile.write('\n')
    bodyfile.write(usefiles)
    bodyfile.write('      implicit none\n')
    bodyfile.write('\n')
    bodyfile.write('    type(restart_file),intent(inout) :: file\n')
    bodyfile.write('    character(*),      intent(in)    :: prefix\n')
    bodyfile.write('    character(*),      intent(in)    :: name\n')
    bodyfile.write('    type(%s),dimension(:),pointer    :: values\n'%(t.name))
    bodyfile.write('\n')
    bodyfile.write('    integer :: status,varid\n')
    bodyfile.write('    character(varnamelen) :: varname\n')
    bodyfile.write('    integer,dimension(1) :: sh\n')
    bodyfile.write('    integer :: i\n')
    bodyfile.write('\n')
    bodyfile.write('    if (associated(values)) then\n')
    bodyfile.write('       sh=shape(values)\n')
    bodyfile.write('       call write_array_pointer(file,prefix,name,sh)\n')
    bodyfile.write('       do i=1,sh(1)\n')
    bodyfile.write('          call %s(file,values(i))\n'%(writesubname))
    bodyfile.write('       end do\n')
    bodyfile.write('    else\n')
    bodyfile.write('       call write_null_array_pointer(file,prefix,name)\n')
    bodyfile.write('    end if\n')
    bodyfile.write('\n')
    bodyfile.write('end subroutine %s\n'%(wparrn_1d))
    # Read code (scalar)
    bodyfile.write('   recursive subroutine %s(file,dat)\n' % (readsubname))
    bodyfile.write('\n')
    bodyfile.write(usefiles)
    bodyfile.write('      implicit none\n')
    bodyfile.write('\n')
    bodyfile.write('      type(restart_file),intent(inout) :: file\n')
    bodyfile.write('      type(%s)   :: dat\n'%(t.name))
    bodyfile.write('      logical :: assoc\n')
    bodyfile.write('\n')
    bodyfile.write('      call check_read(file)\n')    
    bodyfile.write('\n')
    for var in t.components:
        if not (var.name in no_restart):
            var_read(var,bodyfile,t.name,'dat%')
    bodyfile.write('\n')
    bodyfile.write('   end subroutine %s\n' % (readsubname))
    bodyfile.write('\n')
    # Read code (scalar pointer)
    bodyfile.write('   recursive subroutine %s(file,dat)\n' % (rpointn))
    bodyfile.write('\n')
    bodyfile.write(usefiles)
    bodyfile.write('      implicit none\n')
    bodyfile.write('\n')
    bodyfile.write('      type(restart_file),intent(inout) :: file\n')
    bodyfile.write('      type(%s),pointer   :: dat\n'%(t.name))
    bodyfile.write('      logical :: assoc\n')
    bodyfile.write('\n')
    bodyfile.write('      call check_read(file)\n')    
    bodyfile.write('\n')
    bodyfile.write('      allocate(dat)\n')
    for var in t.components:
        if not (var.name in no_restart):
            var_read(var,bodyfile,t.name,'dat%')
    bodyfile.write('\n')
    bodyfile.write('   end subroutine %s\n' % (rpointn))
    bodyfile.write('\n')
    # Read code (1D pointer array)
    bodyfile.write('recursive subroutine %s(file,prefix,name,values)\n'%(rparrn_1d))
    bodyfile.write('\n')
    bodyfile.write(usefiles)
    bodyfile.write('      implicit none\n')
    bodyfile.write('\n')
    bodyfile.write('      type(restart_file),intent(inout) :: file\n')
    bodyfile.write('      character(*),      intent(in)    :: prefix\n')
    bodyfile.write('      character(*),      intent(in)    :: name\n')
    bodyfile.write('      type(%s),dimension(:),pointer    :: values\n'%(t.name))
    bodyfile.write('\n')
    bodyfile.write('    integer :: status,varid\n')
    bodyfile.write('    character(varnamelen) :: varname\n')
    bodyfile.write('    integer,dimension(1) :: sh\n')
    bodyfile.write('    logical :: assoc\n')
    bodyfile.write('    integer :: i\n')
    bodyfile.write('\n')
    bodyfile.write('    if (associated(values)) then\n')
    bodyfile.write('       deallocate(values)\n')
    bodyfile.write('       nullify(values)\n')
    bodyfile.write('    end if\n')
    bodyfile.write('\n')
    bodyfile.write('    call read_array_pointer(file,prefix,name,sh,assoc)\n')
    bodyfile.write('    if (assoc) then\n')
    bodyfile.write('           allocate(values(sh(1)))\n')
    bodyfile.write('           do i=1,sh(1)\n')
    bodyfile.write('              call %s(file,values(i))\n'%(readsubname))
    bodyfile.write('           end do\n')
    bodyfile.write('        end if\n')
    bodyfile.write('\n')
    bodyfile.write('end subroutine %s\n'%(rparrn_1d))
    bodyfile.write('\n')

def generate_restart(fort_fname,headfile,bodyfile):
    print 'Generating restart code for',fort_fname
    parsed_lines=ParseFortran(fort_fname)
    filestruct=InterpFortFile(parsed_lines)
    # Go through each module in turn
    for mod in filestruct.modules:
        defname='RST_'+mod.name.upper()
        headfile.write('#ifdef %s\n' %(defname))
        bodyfile.write('#ifdef %s\n' %(defname))
        write_mod_restart(mod,headfile,bodyfile)
        for t in mod.typedefs:
            write_type_restart(t,headfile,bodyfile)
        headfile.write('#endif\n')
        bodyfile.write('#endif\n')

def usage():
    "short help message"
    print 'Usage: make_restarts.py [OPTIONS] f90files'
    print 'generate restart include files for a set of f90/95 files'
    print ''
    print '  -h, --help\n\tthis message'
    print '  -o <filename>\n\tstem of output file'

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],'ho:',['help','outstem'])
    except getopt.GetoptError:
        # print usage and exit
        usage()
        sys.exit(2)
   
    if len(args) < 1:
        # print usage and exit
        usage()
        sys.exit(2)

    outstem=''

    for o,a in opts:
        if o in ('-h', '--help'):
            usage()
            sys.exit(0)
        if o in ('-o', '--outstem'):
            outstem=a

    if outstem=='':
        print 'ERROR: No output stem supplied'
        usage()
        sys.exit(0)

    headfile=open(outstem+'_head.inc','w')
    bodyfile=open(outstem+'_body.inc','w')
    write_common_restart(headfile,bodyfile)
    
    for arg in args:
        if check_restart(arg): generate_restart(arg,headfile,bodyfile)
        
    headfile.close()
    bodyfile.close()

