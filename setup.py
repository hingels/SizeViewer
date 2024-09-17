from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext): # https://stackoverflow.com/a/21621689
    def finalize_options(self):
        _build_ext.finalize_options(self)
        __builtins__['__NUMPY_SETUP__'] = False
        import numpy
        self.include_dirs.append(numpy.get_include())

setup(
    cmdclass={'build_ext':build_ext},
    setup_requires=['numpy'],
    ext_modules=[
        Extension(
            name = "nanotracking.data_handler",
            sources = ["src/nanotracking/data_handler.c"]
        ),
    ],
)
# setup(
#     ext_modules=[
#         Extension(
#             name = "nanotracking.data_handler",
#             sources = ["src/nanotracking/data_handler.c"]
#         ),
#     ],
# )