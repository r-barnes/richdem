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
    extra_compile_args = ['-std=c++11','-O3'],
    define_macros      = [
      ('DOCTEST_CONFIG_DISABLE', None                ),
      ('RICHDEM_COMPILE_TIME',   RICHDEM_COMPILE_TIME)
    ]
  )
]

setuptools.setup(name='richdem',
  version      = '0.0.1',
  description  = 'High-Performance Terrain Analysis',
  url          = 'http://github.com/r-barnes/richdem',
  author       = 'Richard Barnes',
  author_email = 'rbarnes@umn.edu',
  #license      = 'MIT', #No licensing yet
  packages     = setuptools.find_packages(),
  #scripts      = ['bin/mander'],
  ext_modules  = ext_modules
)
