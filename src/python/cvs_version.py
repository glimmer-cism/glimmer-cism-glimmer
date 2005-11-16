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

# python script that parses the CVS/Entries file to get date and branch of current version

import time,sys

if len(sys.argv)<2:
    print 'Error, specify Entries file to be read'
    sys.exit(1)
inname = sys.argv[1]


# check if file is present
try:
    infile = open(inname,'r')
except:
    infile = None
    sys.exit(1)

# get array of times and tags
times=[]
cvstag=''
for line in infile.readlines():
    l=line.split('/')
    if len(l)>5:
        try:
            t = time.strptime(l[3].strip())
        except:
            t = None
            pass
        if t!=None:
            times.append(t)
            tag=l[5].strip()
            if len(tag)>0:
                tag=tag[1:]
                if cvstag=='':
                    cvstag=tag

infile.close()

# get largest time
times.sort()
cvstime = time.strftime('%Y-%m-%d %H:%M:%S',times[-1])

if len(cvstag)>0:
    cvstag=' '+cvstag

# reading input file and writing to output
if len(sys.argv) == 4:
    infile = open(sys.argv[2],'r')
    outfile = open(sys.argv[3],'w')

    for l in infile.readlines():
        l = l.replace('cvs_vers_string','CVS: %s%s'%(cvstime,cvstag))
        outfile.write(l)
else:
    print cvstime, cvstag
