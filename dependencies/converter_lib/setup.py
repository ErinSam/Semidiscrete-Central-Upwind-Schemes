from setuptools import Extension, setup
from Cython.Build import cythonize


extensions = [
    Extension("converter", ["converter.pyx"])
]


setup(
    ext_modules = cythonize(extensions)
)
