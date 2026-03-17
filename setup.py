from skbuild import setup  # This line replaces 'from setuptools import setup'
import argparse

import io,os,sys
this_directory = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="edep2supera",
    version="2.0.0",
    cmake_source_dir='src/',
    include_package_data=True,
    cmake_args=[
        #'-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON',
        '-DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.9',
    ],
    author=['Kazuhiro Terao, Zach Hulcher'],
    author_email='kterao@slac.stanford.edu, zhulcher@slac.stanford.edu',
    description='Supera interface for EDepSim input data files',
    license='MIT',
    keywords='supera edep-sim larcv larcv3 neutrinos deep learning lartpc_mlreco3d',
    project_urls={
        'Source Code': 'https://github.com/DeepLearnPhysics/edep2supera'
    },
    url='https://github.com/DeepLearnPhysics/edep2supera',
    scripts=['bin/run_edep2supera.py'],
    packages=['edep2supera'],
    package_dir={'': 'python'},
    package_data={'edep2supera': ['config_data/*.yaml']},
    install_requires=[
        'numpy',
        'scikit-build',
        'supera',
    ],
    long_description=long_description,
    long_description_content_type='text/markdown',
)
