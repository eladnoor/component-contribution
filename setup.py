#!/usr/bin/env python


import versioneer
from setuptools import setup


# Most arguments are set in the `setup.cfg`.
setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
