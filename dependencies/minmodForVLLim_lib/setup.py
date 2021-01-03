from setuptools import Extension, setup
from Cython.Build import cythonize


extns = [
    Extension("minmodForVLLim", ["minmodForVLLim.pyx"])
]


setup(
    name = "Definition of minmod function",
    ext_modules = cythonize(extns)
)
