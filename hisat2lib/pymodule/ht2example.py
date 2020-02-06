#!/usr/bin/python

#
# Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.

#
# ./setup.py build
# ./setup.py install
#
import ht2py 

# Path to index
ht2_index = '../../evaluation/indexes/HISAT2_22/22_rep'

# Get default options
ht2_options = ht2py.get_options()

print ht2_options
ht2_options['gVerbose'] = 1
ht2_options['startVerbose'] = 1
# or
ht2_options = {}

handle = ht2py.init(ht2_index, ht2_options)

print ht2py.index_getrefnamebyid(handle, 0)

#print ht2py.index_getrefnamebyid(handle, 0, 1, 3, 5, 7, 9)
# outofindex
#print ht2py.index_getrefnamebyid(handle, 1)
#print ht2py.index_getrefnamebyid(handle, -1)

refnames = ht2py.index_getrefnames(handle)

#for name in refnames:
#    print name

# repeat expansion
positions = ht2py.repeat_expand(handle, 'rep100-300', 8308, 100) 

for pos in positions:
    chr_id = pos[0]
    direction = pos[1]
    chr_pos = pos[2]

    chr_dir = '+'
    if direction == 1:
        chr_dir = '-'

    print refnames[chr_id].split()[0] + ":" + str(chr_pos) + ':' + chr_dir

# close handle
ht2py.close(handle)
