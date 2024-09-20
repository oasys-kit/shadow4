# coding: utf-8
# /*##########################################################################
#
# Copyright (c) 2016-2019 European Synchrotron Radiation Facility
# Copyright (c) 2019-2020 Lawrence Berkeley National Laboratory
# Copyright (c) 2020-now European Synchrotron Radiation Facility
# Copyright (c) 2023-now Argonne National Laboratory
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/

#
# Memorandum: 
#
# Install from sources: 
#     git clone https://github.com/oasys-kit/shadow4
#     cd shadow4
#     python -m pip install -e . --no-deps --no-binary :all:
#
# Upload to pypi (when uploading, increment the version number):
#     python setup.py register (only once, not longer needed)
#     python setup.py sdist
#     python -m twine upload dist/...
#          
# Install from pypi:
#     pip install shadow4
#

__authors__ = ["M Sanchez del Rio, Luca Rebuffi"]
__license__ = "MIT"
__date__ = "2016-now"

import os
from setuptools import setup

PACKAGES = [
    "shadow4",
    "shadow4.beam",
    "shadow4.sources",
    "shadow4.sources.source_geometrical",
    "shadow4.sources.bending_magnet",
    "shadow4.sources.undulator",
    "shadow4.sources.wiggler",
    "shadow4.optical_surfaces",
    "shadow4.physical_models.prerefl",
    "shadow4.physical_models.mlayer",
    "shadow4.beamline",
    "shadow4.beamline.optical_elements",
    "shadow4.beamline.optical_elements.ideal_elements",
    "shadow4.beamline.optical_elements.absorbers",
    "shadow4.beamline.optical_elements.mirrors",
    "shadow4.beamline.optical_elements.multilayers",
    "shadow4.beamline.optical_elements.crystals",
    "shadow4.beamline.optical_elements.gratings",
    "shadow4.beamline.optical_elements.refractors",
    "shadow4.tools",
]

INSTALL_REQUIRES = (
    'setuptools',
    'numpy',
    'scipy',
    'syned>=1.0.40',
    'srxraylib>=1.0.54',
    'crystalpy>=0.0.23',
)

README_FILE = os.path.join(os.path.dirname(__file__), 'README.rst')
LONG_DESCRIPTION = open(README_FILE).read()

setup(name='shadow4',
      version='0.1.44',
      description='shadow implementation in python',
      long_description=LONG_DESCRIPTION,
      long_description_content_type='text/x-rst',
      author='Manuel Sanchez del Rio, Luca Rebuffi',
      author_email='srio@esrf.eu',
      url='https://github.com/oasys-kit/shadow4/',
      download_url='https://github.com/oasys-kit/shadow4/',
      packages=PACKAGES,
      install_requires=INSTALL_REQUIRES,
     )

