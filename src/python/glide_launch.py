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

# a handy python script which will launch glide and record the execution time of the model

import os,sys,fcntl,time, ConfigParser, getopt,os.path

def usage():
    """Print short usage message."""

    print 'Useage: %s [options] config_file'%sys.argv[0]
    print 'launch glide and record execution time of the model.'
    print 'the model binary can be either set using the environment variable'
    print 'GLIDE_MODEL or is automatically determined from the config file.'
    print ''
    print 'config_file is the name of the model configuration to be used.'
    print ''
    print 'options'
    print '--help\tthis message'
    print '-r NAME --results=NAME\tname of file where timing info is stored'
    print '--prefix PFX\tGLIMMER prefix'
    print '--src\tassume prefix is source directory'

def get_runtype(config):
    """Determine which model to run given configuration.

    First we check for the presence of the environment variable GLIDE_MODEL, if
    that fails we try to determine which binary to run given the config file.

    config: ConfigParser object"""

    try:
        model = os.environ['GLIDE_MODEL']
    except KeyError:
        # simple_glide
        if ('EISMINT-1 fixed margin' in config.sections() 
            or 'EISMINT-1 moving margin' in config.sections()
            or 'EISMINT-2' in config.sections()):
            model = 'simple_glide'
            # eis_glide
        elif ('EIS ELA' in config.sections()):
            model = 'eis_glide'
            # no idea
        else:
            raise KeyError, 'no idea what model I should start'
    return model

def get_lock(f):
    """Wait until we get lock for file f."""

    got_lock = False
    while not got_lock:
        try:
            fcntl.lockf(f,(fcntl.LOCK_EX|fcntl.LOCK_NB))
        except IOError:
            got_lock=False
            time.sleep(1)
        else:
            got_lock=True

def release_lock(f):
    """Release lock for file f."""

    fcntl.lockf(f,fcntl.LOCK_UN)

def get_gmtdate():
    """Get current date."""

    d = time.gmtime()
    return '%d-%d-%d_%d:%d'%(d[0],d[1],d[2],d[3],d[4])

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],'hr:',['help','prefix=','results=','src'])
    except getopt.GetoptError,error:
        # print usage and exit
        print error
        usage()
        sys.exit(1)

    if len(args) == 1:
        configname = args[0]
    else:
        usage()
        sys.exit(1)        

    prefix=None
    is_src = False
    results_name = 'results'
    for o,a in opts:
        if o in ('-h', '--help'):
            usage()
            sys.exit(0)
        if o in ('-r', '--results'):
            results_name = a
        if o == '--prefix':
            prefix=a
        if o == '--src':
            is_src = True


    config = ConfigParser.ConfigParser()
    config.readfp(open(configname))

    model = get_runtype(config)
    if prefix!=None:
        if is_src:
            model = os.path.join(prefix,'src','fortran',model)
        else:
            model = os.path.join(prefix,'bin',model)

    # open results file
    if os.path.exists(results_name):
        status = open(results_name,'a')
    else:
        status = open(results_name,'a')
        get_lock(status)
        status.write('#cfg_file\tusr_time\tsys_time\tdate\t\tversion\tfcflags\n')
        release_lock(status)
        
    prog = os.popen(model,'w')
    prog.write(configname)
    prog.close()

    t = os.times()

    # extract some info
    if prefix!=None:
        if is_src:
            glimmer_cfg = 'sh %s'%os.path.join(prefix,'glimmer-config')
        else:
            glimmer_cfg = os.path.join(prefix,'bin','glimmer-config')
    else:
        glimmer_cfg = 'glimmer-config'
    prog = os.popen('%s --version'%glimmer_cfg)
    version = prog.readline().strip()
    prog.close()
    prog = os.popen('%s --fcflags'%glimmer_cfg)
    data = prog.readline()
    prog.close()
    fcflags=''
    for f in data.split():
        interesting = True
        for i in ('-I','-fpp'):
            if f.startswith(i):
                interesting = False
        if interesting:
            fcflags='%s %s'%(fcflags,f)
    fcflags='\"%s\"'%fcflags.strip()
    
    # get lock
    get_lock(status)
    status.write('%s\t%f\t%f\t%s\t%s\t%s\n'%(configname,t[2],t[3],get_gmtdate(),version,fcflags))

    # release lock
    release_lock(status)
    status.close()
