"""orbdot.

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    defaults
    models
    tools
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("orbdot")
except PackageNotFoundError:
    # package is not installed
    pass

del version, PackageNotFoundError
