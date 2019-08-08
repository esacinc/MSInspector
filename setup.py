import sys
import os
from setuptools import setup, Extension, find_packages

if sys.version_info < (2, 7):
	sys.stdout.write("At least Python 2.7 is required.\n")
	sys.exit(1)

# Get the long description from the README file
with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'README.md')) as f:
	long_description = f.read()

setup(
	name='MSInspector',
	version='1.1',
	author='Yin Lu',
	author_email='yin.lu@esacinc.com',
	description='Quality control of the five experiments in CPTAC assay portal',
	long_description=long_description,
	license='MIT',
	package_dir={'': 'src'},
	packages=find_packages('src', exclude=['.txt']),
	package_data = {'MSInspector':['rScripts/*.R', 'skyrTemps/*.skyr', 'htmlTemps/*.html']},
	install_requires=['pandas>=0.21.0', 'Jinja2>=2.9.6',
	],
	entry_points={'console_scripts': ['MSInspector = MSInspector.__main__:main']},
	classifiers=[
		"Development Status :: 1 - Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 2.7",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
