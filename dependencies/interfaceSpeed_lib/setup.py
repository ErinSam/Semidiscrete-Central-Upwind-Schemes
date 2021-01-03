from setuptools import Extension, setup
from Cython.Build import cythonize


extensions = [
    Extension("interfaceSpeed", ["interfaceSpeed.pyx"])
]


setup(
    name = "Function definitions to calculating cell interface speed", 
    ext_modules = cythonize(extensions)
)
