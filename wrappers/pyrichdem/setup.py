import setuptools
import glob
import datetime
import subprocess

RICHDEM_COMPILE_TIME = '"'+datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')+'"'

try:
  RICHDEM_GIT_HASH = subprocess.check_output(["git", "describe"]).strip()
  RICHDEM_GIT_HASH = '"'+RICHDEM_GIT_HASH.decode("utf-8")+'"'
except:
  RICHDEM_GIT_HASH = ''

ext_modules = [
  setuptools.Extension(
    "_richdem",
    glob.glob('src/*.cpp') + ['lib/richdem/common/random.cpp', 'lib/richdem/richdem.cpp'],
    include_dirs       = ['lib/'],
    language           = 'c++',
    extra_compile_args = ['-std=c++11','-g','-fvisibility=hidden','-O3','-Wno-unknown-pragmas'], #Figure out if we can do '-flto'
    define_macros      = [
      ('DOCTEST_CONFIG_DISABLE', None                ),
      ('RICHDEM_COMPILE_TIME',   RICHDEM_COMPILE_TIME),
      ('RICHDEM_GIT_HASH',       RICHDEM_GIT_HASH    ),
      ('RICHDEM_LOGGING',        None                )
    ]
  )
]

long_description = """RichDEM is a set of digital elevation model (DEM) hydrologic analysis tools.

RichDEM uses parallel processing and state of the art algorithms to quickly process even very large DEMs.

RichDEM offers a variety of flow metrics, such as D8 and D∞.

It can flood or breach depressions, as well as calculate flow accumulation, slopes, curvatures, &c."""


#TODO: https://packaging.python.org/tutorials/distributing-packages/#configuring-your-project
setuptools.setup(
  name              = 'richdem',
  version           = '0.0.3',
  description       = 'High-Performance Terrain Analysis',
  long_description  = long_description,
  url               = 'https://github.com/r-barnes/richdem',
  author            = 'Richard Barnes',
  author_email      = 'rbarnes@umn.edu',
  license           = 'GPLv3',
  packages          = setuptools.find_packages(),
  #scripts           = glob.glob('bin/*'),
  entry_points = {'console_scripts': [
    'rd_depression_remove=richdem.cli:DepressionFilling'
  ]},
  ext_modules       = ext_modules,
  keywords          = 'GIS terrain hydrology geomorphology raster',
  #packages         = find_packages(exclude=['contrib', 'docs', 'tests*']),
  install_requires  = ['numpy>=1.10,<2'],
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
