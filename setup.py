from setuptools import setup, find_packages, Extension

import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), 'src/'))
import py_qcode as pq

#qcode_dist = Extension('qcode_dist', 
#						sources = ['src/c/qcode_dist.c'])

setup(
    name='py_qcode',
    version='{0}.{1}.{2}'.format(*pq.__version__),
    url='http://bcriger.github.com/py_qcode/',
    author='Ben Criger',
    author_email='bcriger@gmail.com',
    package_dir={'': 'src'},
    packages=['py_qcode'],
    include_package_data=True,
#    ext_modules=[qcode_dist]
)
