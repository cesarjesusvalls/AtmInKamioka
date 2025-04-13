from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension

extensions = [
    Extension(
        "EarthDensity",
        sources=["EarthDensity.pyx"],
        language="c++",
        include_dirs=["/Users/cjesus/Documents/Prob3plusplus/"],
        extra_objects=["/Users/cjesus/Documents/Prob3plusplus/libThreeProb_3.20.a"],
        extra_compile_args=["-std=c++14"]
    ),
    Extension(
        "BargerPropagator",
        sources=["BargerPropagator.pyx"],
        language="c++",
        include_dirs=["/Users/cjesus/Documents/Prob3plusplus/"],
        extra_objects=["/Users/cjesus/Documents/Prob3plusplus/libThreeProb_3.20.a"],
        extra_compile_args=["-std=c++14"]
        ),
]

setup(
    ext_modules = cythonize(extensions, language_level="3")
)