import setuptools
import glob
import datetime

RICHDEM_COMPILE_TIME = '"'+datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')+'"'
print(RICHDEM_COMPILE_TIME)

ext_modules = [
  setuptools.Extension(
    "_richdem",
    glob.glob('src/*.cpp'),
    include_dirs       = ['lib/'],
    language           = 'c++',
    extra_compile_args = ['-std=c++11','-O3','-fvisibility=hidden', '-flto'],
    define_macros      = [
      ('DOCTEST_CONFIG_DISABLE', None                ),
      ('RICHDEM_COMPILE_TIME',   RICHDEM_COMPILE_TIME)
    ]
  )
]


#TODO: https://packaging.python.org/tutorials/distributing-packages/#configuring-your-project
setuptools.setup(
  name              = 'richdem',
  version           = '0.0.1',
  description       = 'High-Performance Terrain Analysis',
  long_description  = 'TODO',
  url               = 'https://github.com/r-barnes/richdem',
  author            = 'Richard Barnes',
  author_email      = 'rbarnes@umn.edu',
  license           = 'GPLv3', #TODO
  packages          = setuptools.find_packages(),
  #scripts          = ['bin/mander'],
  ext_modules       = ext_modules,
  keywords          = 'TODO TODO',
  #packages         = find_packages(exclude=['contrib', 'docs', 'tests*']),
  #install_requires = ['peppercorn'],
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
