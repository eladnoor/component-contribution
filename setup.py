#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import absolute_import

import versioneer
from setuptools import setup


# Most arguments are set in the `setup.cfg`.
setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass()
)
