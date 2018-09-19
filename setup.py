# coding: utf-8
# /*##########################################################################
#
# Copyright (c) 2016 European Synchrotron Radiation Facility
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
#     git clone https://github.com/SophieTh/und_Sophie_2016.git pySRU
#     cd pySRU
#     python setup.py install
#
# Upload to pypi (when uploading, increment the version number):
#     python setup.py register (only once, not longer needed)
#     python setup.py sdist
#     python setup.py upload
#          
# Install from pypi:
#     pip install pySRU
#

__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "25/11/2016"

from setuptools import setup

setup(name='minishadow',
      version='0.0.1',
      description='implementation in python of some basic ray-tracing that mimic shadow',
      author='Manuel Sanchez del Rio',
      author_email='srio@esrf.eu',
      url='https://github.com/srio/minishadow/',
      packages=['minishadow'],
      install_requires=[
                        'numpy',
                        'scipy'
                       ],
      # test_suite='tests',
     )

