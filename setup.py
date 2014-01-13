import sys
if 'build_sphinx' in sys.argv or 'develop' in sys.argv:
    from setuptools import setup,Extension
else:
    from distutils.core import setup,Extension
import distutils.util
#from scipy_distutils.core import Extension as scipyExtension
#from scipy_distutils.core import setup as scipysetup
from numpy.distutils.core import Extension as numpyExtension
from numpy.distutils.core import setup as numpysetup
#from numpy.distutils.core import build_ext
from numpy.distutils.command import build_src
import Cython
import Cython.Compiler.Main
build_src.Pyrex = Cython
build_src.have_pyrex = True
from Cython.Distutils import build_ext
import Cython
import numpy
import os
import shutil

# removed following lines as per http://www.mail-archive.com/numpy-discussion@scipy.org/msg19932.html
# OLD from numpy.distutils.core import setup
# OLD from numpy.distutils.core import Extension

# copied from pyspeckit's...
with open('README.rst') as file:
    long_description = file.read()

with open('CHANGES') as file:
    long_description += file.read()

with open('REQUIREMENTS') as file:
    requirements = file.readlines()


# print "To create cplfit.so (for importing), call command: "
# print "python setup.py build_ext --inplace"
# print "If this fails, make sure c_numpy.pxd is in the path somewhere (e.g. this directory)"

try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
    numpy_include_dirs = get_numpy_include_dirs()
except AttributeError:
    numpy_include_dirs = numpy.get_include()


dirs = list(numpy_include_dirs)
dirs.extend(Cython.__path__)
dirs.append('.')

ext_cplfit = Extension("plfit/cplfit",
                       ["plfit/cplfit.pyx"],
                       include_dirs=dirs,
                       extra_compile_args=['-O3'])

#ext_fplfit = numpyExtension(name="fplfit",
#                    sources=["fplfit.f"])

if __name__=="__main__":

    # can't specify fcompiler if numpysetup is included
    # therefore, run this command separately
    # gfortran = OK.  g77, g95 NOT ok
    # also, this is kind of a ridiculous hack...
    if any([x in sys.argv for x in ['build','install']]):
        fortran_compile_command = "cd plfit && f2py -c fplfit.f -m fplfit --fcompiler=gfortran && cd .."
        os.system(fortran_compile_command)
    # do this first so it gets copied (in principle...)
    # in practice, see hack cont'd
    if os.path.exists('plfit/fplfit.so'):
        build_dir = 'build/lib.{0}-{1}.{2}/plfit/'.format(distutils.util.get_platform(),
                                                          sys.version_info[0],
                                                          sys.version_info[1])
        if not os.path.exists(build_dir):
            os.makedirs(build_dir)
        shutil.copy2('plfit/fplfit.so',build_dir+"/fplfit.so")

    S = setup(
        name="plfit",
        version="1.0.2",
        description="Python implementation of Aaron Clauset's power-law distribution fitter",
        long_description=long_description,
        author="Adam Ginsburg",
        author_email="adam.g.ginsburg@gmail.com",
        url="https://github.com/keflavich/plfit",
        download_url="https://github.com/keflavich/plfit/archive/master.zip",
        license="MIT",
        platforms=["Linux","MacOS X"],
        packages=['plfit','plfit.tests'],
        # obsolete package_dir={'plfit':'.'},
        install_requires=["numpy","cython"],
        ext_modules=[ext_cplfit],
        cmdclass={'build_ext': build_ext}
    )

    #numpysetup(name = 'fplfit',
    #      ext_modules = [ext_fplfit]
    #      )




#print "I can't get numpy.distutils to compile the fortran.  To do it yourself, run some variant of:"
#print 'f2py -c fplfit.f -m fplfit'
# keep an eye on this: http://stackoverflow.com/questions/7932028/setup-py-for-packages-that-depend-on-both-cython-and-f2py

# try:
#     os.system('f2py -c fplfit.f -m fplfit')
# except:
#     print "Could not build fplfit"

