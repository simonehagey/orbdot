from setuptools import setup, find_packages

with open("requirements.txt", "r") as fh:
    requires = fh.readlines()

setup(
    name="orbdot",
    version="0.1.0",
    packages=find_packages(),
    install_requires=requires,
    description="A set of tools for studying secular evolution of exoplanet orbits.",
    url="https://github.com/simonehagey/orbdot",
    author="Simone R. Hagey",
    author_email="shagey@phas.ubc.ca",
    package_data={"orbdot": ["defaults/*"]},
    include_package_data=True
)