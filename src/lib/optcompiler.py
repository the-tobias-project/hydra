# distutils: language = c++
from distutils.core import setup, Extension
from Cython.Build import cythonize

# setup(ext_modules=cythonize("corr.pyx"))
setup(ext_modules=cythonize(Extension('opt',
      sources=['opt.pyx', "../../../L-BFGS-B-C/src/lbfgsb.c", "../../../L-BFGS-B-C/src/subalgorithms.c",
               "../../../L-BFGS-B-C/src/linpack.c", "../../../L-BFGS-B-C/src/timer.c",
               "../../../L-BFGS-B-C/src/linesearch.c",
               "../../../L-BFGS-B-C/src/miniCBLAS.c", "../../../L-BFGS-B-C/src/print.c"],
  language="c", extra_compile_args=['-std=gnu++0x']
  )))
# setup(ext_modules = cythonize(Extension('corr',
#  sources=['corr.pyx',  'plink-ng/2.0/plink2_cmdline.cc', 'plink-ng/2.0/plink2_string.cc',
#           'plink-ng/2.0/plink2_base.cc',
#    ],
#  language = "c++",
#  )))
