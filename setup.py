import sys
import os
from setuptools import setup, Extension, find_packages

if sys.version_info < (3, 11):
	sys.stdout.write("At least Python 3.11 is required.\n")
	sys.exit(1)

# Get the long description from the README file
with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'README.md')) as f:
	long_description = f.read()

setup(
	name='MSInspector',
	version='2.1.0',
	author='Yin Lu',
	author_email='yin.lu@icf.com',
	description='Quality control of the five experiments in CPTAC assay portal',
	long_description=long_description,
	license='MIT',
	package_dir={'': 'src'},
	packages=find_packages('src', exclude=['.txt']),
	package_data = {'MSInspector':['rScripts/*.R', 'skyrTemps/*.skyr', 'htmlTemps/*.html']},
	install_requires=['pandas>=2.0.0', 'Jinja2>=3.1.2',
	],
	entry_points={'console_scripts': ['MSInspector = MSInspector.__main__:main']},
	classifiers=[
		"Development Status :: 1 - Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 3.11",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
