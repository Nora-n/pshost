#! /usr/bin/env python3
#

DESCRIPTION = "pshost: Basic tools for astrophysics and cosmology"
LONG_DESCRIPTION = """\
This module gathers the basic tools one usually needs for astrophysics and
cosmology usage.
"""

DISTNAME = 'pshost'
AUTHOR = 'pshost Developers'
MAINTAINER = 'Nora Nicolas'
MAINTAINER_EMAIL = 'nora.nicolas@ens-lyon.org'
URL = 'https://github.com/Nora-n/pshost/'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/Nora-n/pshost/'
VERSION = '0.1'

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup

# def check_dependencies():
#     install_requires = []
#
#     # Just make sure dependencies exist, I haven't rigorously
#     # tested what the minimal versions that will work are
#     # (help on that would be awesome)
#     try:
#         import propobject
#     except ImportError:
#         install_requires.append('propobject')
#     try:
#         import astropy
#     except ImportError:
#         install_requires.append('astropy')
#
#     return install_requires

if __name__ == "__main__":

    #    install_requires = check_dependencies()

    if _has_setuptools:
        packages = find_packages()
        print(packages)
    else:
        # This should be updated if new submodules are added
        packages = [
            'pshost']

    setup(name=DISTNAME,
          author=AUTHOR,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          #          install_requires=install_requires,
          packages=packages,
          classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 3.7',
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'],
          )
