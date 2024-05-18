from setuptools import setup, find_packages
import sys
import os

# os.chdir(os.path.abspath(os.path.dirname(__file__)))

setup(
    name="orbdot",
    version="0.1.0",
    description="A set of tools for studying secular evolution of exoplanet orbits",
    url="https://github.com/simonehagey/orbdot",
    author="Simone R. Hagey",
    author_email="shagey@phas.ubc.ca",
    packages=["orbdot", "orbdot.tools", "orbdot.models", "orbdot.defaults"]

)