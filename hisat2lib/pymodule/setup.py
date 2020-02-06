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


from distutils.core import setup, Extension

module1 = Extension('ht2py',
#                    define_macros = [('DEBUG', '1')],
                    include_dirs=['../'],
                    libraries=['stdc++'],
                    extra_objects = ['../../libhisat2lib.a'],
                    sources = ['ht2module.c'])

setup(name = 'ht2py',
        version = '1.0',
        ext_modules = [module1])
