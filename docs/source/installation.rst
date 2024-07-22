.. _installation:

************
Installation
************

To install the OrbDot package, you have a few options depending on your environment and preferences.

Using pip
---------
The easiest way to install OrbDot is via pip. You can install the package directly from PyPI (Python Package Index) by running:

.. code-block:: bash

    pip install orbdot

This will download and install the latest release of OrbDot and its dependencies.

Using a Local Copy
------------------
If you have a local copy of the OrbDot source code or you want to contribute to the development, you can install it directly from the source. First, clone the repository:

.. code-block:: bash

    git clone https://github.com/simonehagey/orbdot.git

Navigate into the project directory:

.. code-block:: bash

    cd orbdot

Then, install the package using:

.. code-block:: bash

    pip ???

This will allow you to make changes to the code and see them reflected immediately without needing to reinstall.

Dependencies
------------
OrbDot requires Python 3.9 or higher. The package depends on the following libraries:

- `corner <https://github.com/dfm/corner.py>`_
- `matplotlib <https://github.com/matplotlib/matplotlib>`_
- `numpy <https://github.com/numpy/numpy>`_
- `scipy <https://github.com/scipy/scipy>`_
- `astropy <https://github.com/astropy/astropy>`_

and one of:

- `Nestle <https://github.com/kbarbary/nestle>`_ by Kyle Barbary.
- `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_ by Johannes Buchner, a Python interface
  for `MultiNest <https://github.com/JohannesBuchner/MultiNest>`_.

.. note::
    To perform the nested sampling methods the user may choose between two packages: ``Nestle`` and ``PyMultiNest``. ``PyMultiNest`` is generally faster and more robust, but it can be tricky to install, thus it is not a requirement to use this code. ``Nestle`` is included in the ``requirements.txt`` file.