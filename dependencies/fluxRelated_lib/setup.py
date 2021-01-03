from setuptools import setup, Extension
from Cython.Build import cythonize


extns = [
    Extension("fluxRelated", ["fluxRelated.pyx"])
]


setup(
    name = "Function definitions that are related to flux terms",
    ext_modules = cythonize(extns)
)
