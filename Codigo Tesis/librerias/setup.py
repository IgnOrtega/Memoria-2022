# Archivo: setup.py

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension("compute_c_ab_sample_7_8_cython", ["compute_c_ab_sample_7_8_cython.pyx"], include_dirs=[np.get_include()])
]


setup(
    ext_modules = cythonize(extensions)
)
