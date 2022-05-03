#!/usr/bin/env python3
# -*- coding:utf8 -*-

import setuptools
from kmerator import info

setuptools.setup(
    name = 'kmerator2',
    version = info.VERSION,
    author = info.AUTHOR,
    author_email = info.AUTHOR_EMAIL,
    description = info.SHORTDESC,
    long_description = open('README.md').read(),
    long_description_content_type = "text/markdown",
    # ~ url="https://github.com/bio2m/kmerator",
    packages = setuptools.find_packages(),
    classifiers = [
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
    ],
    entry_points = {
        'console_scripts': [
            'kmerator2 = kmerator.kmerator:main',
            'ktools = kmerator.ktools:main',
        ],
    },
    include_package_data = True,
    install_requires=['PyYAML', 'bs4'],
    python_requires = ">=3.6",
    licence = "GPLv3"
)
