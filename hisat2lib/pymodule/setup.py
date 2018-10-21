#!/usr/bin/python
from distutils.core import setup, Extension

module1 = Extension('ht2py',
                    define_macros = [('DEBUG', '1')],
                    include_dirs=['../'],
                    libraries=['hisat2lib'],
                    library_dirs = ['/home/parkch/work/04-hisat2bench/hisat2.api/cmake-build-debug'],
                    sources = ['ht2module.c'])

setup(name = 'ht2py',
        version = '1.0',
        ext_modules = [module1])
