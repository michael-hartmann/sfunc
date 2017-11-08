from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

extension = Extension("bessel", sources=["bessel.pyx", "bessel_c.c"],
        extra_compile_args = ["-std=c99"],
        include_dirs = [".."],
        libraries = ["m"],
)

setup(
    ext_modules = cythonize(extension)
)
