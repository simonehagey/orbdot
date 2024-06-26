from setuptools import setup

setup(
    name="orbdot",
    version="0.1.0",
    description="A set of tools for studying secular evolution of exoplanet orbits.",
    url="https://github.com/simonehagey/orbdot",
    author="Simone R. Hagey",
    author_email="shagey@phas.ubc.ca",
    packages=["orbdot", "orbdot.tools", "orbdot.models", "orbdot.defaults"],
    package_data={"orbdot": ['orbdot/defaults/*']}
)