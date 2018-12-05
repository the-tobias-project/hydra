# distutils: language = c++
from distutils.core import setup, Extension
from Cython.Build import cythonize

#setup(ext_modules=cythonize("corr.pyx"))
setup(ext_modules = cythonize(Extension('corr',
  sources=['corr.pyx'],
  language = "c++", extra_compile_args = ['-std=gnu++0x']
  )))
#setup(ext_modules = cythonize(Extension('corr', 
#  sources=['corr.pyx',  'plink-ng/2.0/plink2_cmdline.cc', 'plink-ng/2.0/plink2_string.cc','plink-ng/2.0/plink2_base.cc',
#    ],
#  language = "c++",
#  )))
