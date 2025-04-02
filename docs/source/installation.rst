.. _installation:

************
Installation
************

To install the OrbDot package, you have different options depending on your environment and preferences.

Using `pip <http://www.pip-installer.org/>`_
--------------------------------------------
The recommended way to install the stable version of OrbDot is using ``pip``:

.. code-block:: bash

    pip install orbdot

This will download and install the latest release and its dependencies.

Installing a Local Copy
-----------------------
You can also install OrbDot directly from a local copy of the source code.

First, clone the OrbDot repository:

.. code-block:: bash

    git clone https://github.com/simonehagey/orbdot.git

Next, navigate into the project directory and install the package:

.. code-block:: bash

    cd orbdot
    python -m pip install .

Test the Installation
---------------------
To ensure that the installation worked, download or navigate (if local) to the ``examples/`` directory and run the ``example_wasp12.py`` script:

.. code-block:: bash

    cd orbdot/examples
    python example_wasp12.py

If the model fitting process starts successfully, the installation was successful!

Dependencies
------------
OrbDot requires Python 3.10 or higher and depends on the following libraries:

- `numpy <https://github.com/numpy/numpy>`_
- `scipy <https://github.com/scipy/scipy>`_
- `matplotlib <https://github.com/matplotlib/matplotlib>`_
- `astropy <https://github.com/astropy/astropy>`_
- `corner <https://github.com/dfm/corner.py>`_
- `nestle <https://github.com/kbarbary/nestle>`_

and optionally:

- `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_ by Johannes Buchner, a Python interface
  for `MultiNest <https://github.com/JohannesBuchner/MultiNest>`_.

.. note::
    When using the nested sampling methods the user may choose between two packages: Nestle and PyMultiNest. The latter is generally faster and more robust, but it can be tricky to install, and thus it is not a requirement to use this code.

    Nestle is included in the ``requirements.txt`` file and will be installed automatically. In order to use PyMultiNest, you will have to follow their installation instructions `here <https://johannesbuchner.github.io/PyMultiNest/install.html>`_ after setting up OrbDot.