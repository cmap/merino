# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))


setup(name="merino", packages=find_packages(exclude=['plate_tracking', 'functional_tests', 'requirements_artifacts', 'test_output']))

