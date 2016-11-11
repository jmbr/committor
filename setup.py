#!/usr/bin/env python

from setuptools import setup, find_packages

import committor.version


classifiers = '''\
Development Status :: 3 - Alpha
Environment :: Console
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
Natural Language :: English
Operating System :: OS Independent
Programming Language :: Python :: 3.5
Topic :: Scientific/Engineering :: Mathematics
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Chemistry
'''

setup(name='committor',
      version=committor.version.version.v_short,
      description='Compute committor function',
      long_description=('Compute committor probability using the gradient '
                        'of the potential of mean force.'),
      license='MIT License',
      classifiers=classifiers.splitlines(),
      author='Juan M. Bello-Rivas',
      author_email='jmbr@superadditive.com',
      url='https://github.com/jmbr/committor/',
      packages=find_packages(),
      package_data={'': ['LICENSE']},
      install_requires=['scipy', 'numpy', 'fenics'],
      entry_points={
          'console_scripts': 'committor = committor.cli:main'
      })
