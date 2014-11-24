# -*- coding: utf-8 -*-
__revision__ = "$Id$"
import sys
import os
from setuptools import setup, find_packages
import glob


_MAJOR               = 0
_MINOR               = 0
_MICRO               = 5
version              = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
release              = '%d.%d' % (_MAJOR, _MINOR)

metainfo = {
    'authors': {
        'Cokelaer':('Thomas Cokelaer','cokelaer@ebi.ac.uk'),
        },
    'version': version,
    'license' : 'BSD',
    'download_url' : ['http://pypi.python.org/pypi/cno'],
    'url' : ['http://pypi.python.org/pypi/cno'],
    'description': "CNO (Cell Net Optimiser): Manipulate, Visualise and Optimise Biological Networks to Perturbation Data." ,
    'platforms' : ['Linux', 'Unix', 'MacOsX', 'Windows'],
    'keywords' : ['CellNOpt', 'CellNOptR', 'CNO', 'Logical model', 'SBML', 'SIF', 
        'Boolean', 'ODE', 'Fuzzy', 'Optimisation', 'Protein Network'],
    'classifiers' : [
          'Development Status :: 1 - Planning',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics']
    }



with open('README.rst') as f:
    readme = f.read()


setup(
    name             = 'cno',
    version          = version,
    maintainer       = metainfo['authors']['Cokelaer'][0],
    maintainer_email = metainfo['authors']['Cokelaer'][1],
    author           = metainfo['authors']['Cokelaer'][0],
    author_email     = metainfo['authors']['Cokelaer'][1],
    long_description = readme,
    keywords         = metainfo['keywords'],
    description = metainfo['description'],
    license          = metainfo['license'],
    platforms        = metainfo['platforms'],
    url              = metainfo['url'],      
    download_url     = metainfo['download_url'],
    classifiers      = metainfo['classifiers'],

    # package installation
    #package_dir = {'':''},
    #packages = ['cno'],
    #package_dir  = package_dir,
    install_requires = ['numpy', 'matplotlib', 'pandas', 'networkx', 'bioservices',
    'colormap>=0.9.3',  'pygraphviz'],
    zip_safe=False,


    packages = find_packages(),
    
    # is this required ?
    package_data={
        'cno.datasets.ToyMMB': ['*'],
        'cno.datasets.ToyPB': ['*'],
    },

    entry_points = {
        'console_scripts': [
            'cno_plotmodel=cno.apps.plotmodel:plotmodel',
            ]
        },

    )


