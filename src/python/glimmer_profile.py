#! /usr/bin/env python

# Copyright 2004, Magnus Hagdorn
# 
# This file is part of GLIMMER.
# 
# GLIMMER is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GLIMMER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GLIMMER; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""plot profiling data produced by profile module"""

import PyGMT,sys,string

Colours = ['255/0/0','0/255/0','0/0/255','0/255/255','255/0/255','255/255/0','127/0/0','0/127/0','0/0/127','0/127/127','127/0/127','127/127/0']

def read_profile(inname):
    """read a profile and return data"""

    # read data file
    profile_data = {}
    infile = open(inname)
    for l in infile.readlines():
        l = l.strip()
        if l[0] != '#':
            l = l.split()
            if len(l)>4:
                id = l[2]
                # insert new profile type into dir
                if id not in profile_data:
                    profile_data[id] = {}
                    profile_data[id]['name'] = string.join(l[4:])
                    profile_data[id]['time'] = []
                    profile_data[id]['duration'] = []
                    profile_data[id]['model_t'] = []
                profile_data[id]['time'].append(float(l[0]))
                profile_data[id]['duration'].append(float(l[1]))
                profile_data[id]['model_t'].append(float(l[3]))
    return profile_data


data = read_profile('glide.profile')

plot = PyGMT.Canvas('profile.ps',size='A4')
plot.defaults['LABEL_FONT_SIZE']='12p'
plot.defaults['ANOT_FONT_SIZE']='10p'

area = PyGMT.AutoXY(plot,size=[18.,10.],pos=[0.,4.])
keyarea = PyGMT.KeyArea(plot,size=[18.,2.])
keyarea.num=[4,4]
i = 0
for id in data.keys():
    area.line('-W1/%s'%Colours[i],data[id]['time'], data[id]['duration'])
    keyarea.plot_line(data[id]['name'],'1/%s'%Colours[i])
    i = i+1
area.finalise()
area.xlabel='total elapsed time [CPU time]'
area.ylabel='duration of operation [CPU time]'
area.axis='WeSn'
area.coordsystem()
plot.close()
