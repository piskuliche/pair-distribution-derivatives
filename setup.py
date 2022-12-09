#!/usr/bin/env python

from setuptools import find_packages, setup

setup(name="pair_distribution_derivatives",
      version="0.1",
      description="A code that calculates pair distribution derivatives",
      author="Zeke Piskulich",
      author_email='piskulichz@gmail.com',
      platforms=["any"],  # or more specific, e.g. "win32", "cygwin", "osx"
      license="BSD",
      url="https://github.com/piskuliche/pair-distribution-derivatives",
      packages=find_packages('pair_distribution_derivatives'),
      )
