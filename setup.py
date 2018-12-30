from setuptools import setup, find_packages

DESCRIPTION = "WebSTR: a population-wide database of short tandem repeat variation in humans"
LONG_DESCRIPTION = DESCRIPTION
NAME = "WebSTR"
AUTHOR = "Richard Yanicky"
AUTHOR_EMAIL = "richyanicky@gmail.com"
MAINTAINER = "Richrad Yanicky"
MAINTAINER_EMAIL = "richyanicky@gmail.com"
DOWNLOAD_URL = 'http://github.com/gymreklab/dbstr'
LICENSE = 'MIT'

VERSION = '1.0.0'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=find_packages(),
      package_data={
          'WebSTR': ['static/css/*.css',
                     'static/*.js',
                     'templates/*.html']
      },
      entry_points={
        'console_scripts': [
          'WebSTR = WebSTR.WebSTR:main',
        ],
      },
      install_requires=['argparse', 'flask','dash','pandas','plotly',\
                        'numpy', 'pyfaidx',
                        'dash_core_components', 'dash_html_components',\
                        'dash_table_experiments'],
      classifiers=['Development Status :: 4 - Beta',\
                       'Programming Language :: Python :: 3.4',\
                       'License :: OSI Approved :: MIT License',\
                       'Operating System :: OS Independent',\
                       'Intended Audience :: Science/Research',\
                       'Topic :: Scientific/Engineering :: Bio-Informatics',\
                       'Topic :: Scientific/Engineering :: Visualization']
     )
