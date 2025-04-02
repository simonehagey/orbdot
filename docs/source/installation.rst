.. _installation:

************
Installation
************

To install the OrbDot package, you have different options depending on your environment and preferences. Before installing OrbDot, it is recommended to create and activate a virtual environment using ``venv`` to prevent dependency conflicts.

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
To ensure that the installation was successful, run an example script.

If you installed OrbDot from a local copy, navigate to the ``examples/`` directory and execute the ``example_wasp12.py`` script:

.. code-block:: bash

    cd orbdot/examples
    python example_wasp12.py

If you installed via ``pip``, you will need to download the ``examples/`` directory from the repository. The example will not work if you only download the ``example_wasp12.py`` file, as the data and config files are also required to run the script.

If the model fitting process starts successfully, the installation was successful!

Dependencies
------------
OrbDot requires Python 3.10 or higher and depends on the following libraries:

- `numpy <https://github.com/numpy/numpy>`_ (>=1.24)
- `scipy <https://github.com/scipy/scipy>`_ (>=1.13)
- `astropy <https://github.com/astropy/astropy>`_ (>=5.1.1)
- `corner <https://github.com/dfm/corner.py>`_ (>=2.2.1)
- `nestle <https://github.com/kbarbary/nestle>`_ (>=0.2.0)
- `matplotlib <https://github.com/matplotlib/matplotlib>`_ (>=3.6.2)

Additional optional dependencies:

- `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_ by Johannes Buchner, a Python interface
  for `MultiNest <https://github.com/JohannesBuchner/MultiNest>`_.

.. note::
    When using the nested sampling methods the users can choose between two packages: Nestle and PyMultiNest. PyMultiNest is generally faster and more robust, but it can be difficult to install. Therefore, it is not required for using OrbDot.

    Nestle is included as a dependency and will be installed automatically. In order to use PyMultiNest, you will have to follow their installation instructions `here <https://johannesbuchner.github.io/PyMultiNest/install.html>`_ after setting up OrbDot.