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
    print '--prefix PFX\tGLIMMER prefix'

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

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],'h',['help','prefix='])
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
    for o,a in opts:
        if o in ('-h', '--help'):
            usage()
            sys.exit(0)
        if o == '--prefix':
            prefix=os.path.join(a,'bin')


    config = ConfigParser.ConfigParser()
    config.readfp(open(configname))

    model = get_runtype(config)
    if prefix!=None:
        model = os.path.join(prefix,model)

    prog = os.popen(model,'w')
    prog.write(configname)
    prog.close()

    t = os.times()

    # get lock
    status = open('results','a')
    got_lock = False
    while not got_lock:
        try:
            fcntl.lockf(status,(fcntl.LOCK_EX|fcntl.LOCK_NB))
        except IOError:
            got_lock=False
            time.sleep(1)
        else:
            got_lock=True

    status.write('%s\t%f\t%f\n'%(configname,t[2],t[3]))

    fcntl.lockf(status,fcntl.LOCK_UN)
    status.close()



