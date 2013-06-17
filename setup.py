#!/usr/bin/env python
from distutils.core import setup,Extension
module1 = Extension('ndfit',
                    include_dirs=['/usr/local/include'],
                    libraries=['pthread'],
                    sources=['ndfitmodule.c','ndfitstruct.c'])


setup(name="ndfit",
      version="0.3",
      description="Non-linear curve fitting module",
      author="Tanhauser",
      url="Place",
      ext_modules=[module1])
