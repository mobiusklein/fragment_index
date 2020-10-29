import os
import sys
import traceback

from setuptools import find_packages, Extension, setup


def make_extensions():
    try:
        import numpy
    except ImportError:
        print("Installation requires `numpy`")
        raise

    try:
        from Cython.Build import cythonize
        cython_directives = {
            'embedsignature': True,
            "profile": False
        }
        extensions = cythonize([
            Extension(name='fragment_index.fragment_index', sources=["fragment_index/fragment_index.pyx"],
                    include_dirs=[numpy.get_include()]),
        ], compiler_directives=cython_directives)
    except ImportError:
        extensions = [
            Extension(name='fragment_index.fragment_index', sources=["fragment_index/fragment_index.c"],
                      include_dirs=[numpy.get_include()]),
        ]
    return extensions


setup(name='fragment_index',
      packages=find_packages(),
      ext_modules=make_extensions(),
      include_package_data=True)