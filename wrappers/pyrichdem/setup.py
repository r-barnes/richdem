import os
import re
import subprocess
import sys
from typing import Optional

import setuptools
from pybind11.setup_helpers import Pybind11Extension
from setuptools.command.build_ext import build_ext as _build_ext

richdem_compile_time: Optional[str] = None
richdem_git_hash: Optional[str] = None

# Compiler specific arguments
BUILD_ARGS = {
    "msvc": ["-std=c++17", "-g", "-fvisibility=hidden", "-O3"],
    "gcc": ["-std=c++17", "-g", "-fvisibility=hidden", "-O3", "-Wno-unknown-pragmas"],
    "unix": ["-std=c++17", "-g", "-fvisibility=hidden", "-O3", "-Wno-unknown-pragmas"],
}

library_dirs = []
if sys.platform.startswith("win"):
    library_dirs.extend(
        [
            os.path.join(sys.prefix, "Library", "lib"),
            os.path.join(sys.prefix, "Library", "bin"),
        ]
    )

# Magic that hooks compiler specific arguments up with the compiler
class build_ext_compiler_check(_build_ext):
    def build_extensions(self):
        compiler = self.compiler.compiler_type
        print(f"COMPILER {compiler}")
        args = BUILD_ARGS[compiler]
        for ext in self.extensions:
            ext.extra_compile_args = args
            print(f"COMPILER ARGUMENTS: {ext.extra_compile_args}")
        _build_ext.build_extensions(self)


if richdem_git_hash is None:
    try:
        shash = (
            subprocess.Popen(
                ["git log --pretty=format:'%h' -n 1"],
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
            )
            .stdout.readlines()[0]
            .decode("utf8")
            .strip()
        )
        sdate = (
            subprocess.Popen(
                ["git log -1 --pretty='%ci'"],
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE,
            )
            .stdout.readlines()[0]
            .decode("utf8")
            .strip()
        )
        if re.match(r"^[0-9a-z]+$", shash) and re.match(
            r"^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}.*$", sdate
        ):
            richdem_compile_time = sdate
            richdem_git_hash = shash
    except:
        print(
            "Warning! Could not find RichDEM version. Software will still work, but reproducibility will be compromised."
        )

if richdem_git_hash is None:
    richdem_compile_time = "Unknown"
    richdem_git_hash = "Unknown"

print("Using RichDEM hash={0}, time={1}".format(richdem_git_hash, richdem_compile_time))

ext_modules = [
    Pybind11Extension(
        "_richdem",
        ["src/pywrapper.cpp"],
        include_dirs=["lib/"],
        library_dirs=library_dirs,
        libraries=["richdem"],
        define_macros=[
            ("DOCTEST_CONFIG_DISABLE", None),
            ("RICHDEM_COMPILE_TIME", f'"\\"{richdem_compile_time}\\""'),
            ("RICHDEM_GIT_HASH", f'"\\"{richdem_git_hash}\\""'),
            # ("RICHDEM_LOGGING", None),
            (
                "_USE_MATH_DEFINES",
                None,
            ),  # To ensure that `#include <cmath>` imports `M_PI` in MSVC
        ],
    ),
]

long_description = """RichDEM is a set of digital elevation model (DEM) hydrologic analysis tools.

RichDEM uses parallel processing and state of the art algorithms to quickly process even very large DEMs.

RichDEM offers a variety of flow metrics, such as D8 and D-infinity.

It can flood or breach depressions, as well as calculate flow accumulation, slopes, curvatures, &c."""


# TODO: https://packaging.python.org/tutorials/distributing-packages/#configuring-your-project
setuptools.setup(
  name              = 'richdem',
  version           = '0.3.5',
  description       = 'High-Performance Terrain Analysis',
  long_description  = long_description,
  url               = 'https://github.com/r-barnes/richdem',
  author            = 'Richard Barnes',
  author_email      = 'rbarnes@umn.edu',
  license           = 'GPLv3',
  packages          = setuptools.find_packages(),
  #scripts           = glob.glob('bin/*'),
  entry_points = {'console_scripts': [
    'rd_depression_filling=richdem.cli:DepressionFilling',
    'rd_breach_depressions=richdem.cli:BreachDepressions',
    'rd_flow_accumulation=richdem.cli:FlowAccumulation',
    'rd_terrain_attribute=richdem.cli:TerrainAttribute',
    'rd_info=richdem.cli:RdInfo',
    'rd_compare=richdem.cli:RdCompare'
  ]},
  ext_modules      = ext_modules,
  cmdclass         = {'build_ext': build_ext_compiler_check},
  keywords         = 'GIS terrain hydrology geomorphology raster',
  #packages        = find_packages(exclude=['contrib', 'docs', 'tests*']),
  install_requires = [
    "numpy>=1.7,<2; python_version > '3.4' or python_version < '3.0'",
    "numpy>=1.7,<1.12; python_version < '3.4' and python_version > '3.0'"
  ],
  # extras_require    = {
  #   ':python_version > "3.4"': [
  #     'numpy>=1.7,<2'
  #   ],
  #   ':python_version < "3.4"': [
  #     'numpy>=1.7,<1.12'
  #   ]
  # },
  #python_requires  = ' >= 2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',

  #TODO: https://pypi.python.org/pypi?%3Aaction=list_classifiers
  classifiers      = [
      'Development Status :: 4 - Beta',

      'Environment :: Console',

      'Intended Audience :: Developers',
      'Intended Audience :: End Users/Desktop',
      'Intended Audience :: Education',
      'Intended Audience :: Science/Research',
      'Intended Audience :: Other Audience',

      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

      'Natural Language :: English',

      'Topic :: Scientific/Engineering :: GIS',
      'Topic :: Scientific/Engineering :: Information Analysis',
      'Topic :: Scientific/Engineering :: Visualization',
      'Topic :: Software Development :: Libraries',

      # Specify the Python versions you support here. In particular, ensure
      # that you indicate whether you support Python 2, Python 3 or both.
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 2.6',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3',
      'Programming Language :: Python :: 3.2',
      'Programming Language :: Python :: 3.3',
      'Programming Language :: Python :: 3.4',
  ]
)
