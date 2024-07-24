.. _installation:

************
Installation
************

To install the OrbDot package, you have a few options depending on your environment and preferences.

Using `pip <http://www.pip-installer.org/>`_
--------------------------------------------

.. important::
    OrbDot is not yet available to pip! Coming soon...

The recommended way to install the stable version of ``orbdot`` is using ``pip``:

.. code-block:: bash

    pip install orbdot

This will download and install the latest release of OrbDot and its dependencies.

Using a Local Copy
------------------
You can also install it directly from a local copy of the source code.

First, clone the repository:

.. code-block:: bash

    git clone https://github.com/simonehagey/orbdot.git

Next, navigate into the project directory and install the package with the ``setup.py`` file:

.. code-block:: bash

    cd orbdot
    python setup.py install

Test the Installation
---------------------
To ensure that the installation worked, navigate to the ``examples/`` directory and run the ``example_wasp12.py`` script:

.. code-block:: bash

    cd orbdot/examples
    python example_wasp12.py

If the first model fit begins to run, you're good to go.

Dependencies
------------
OrbDot requires Python 3.9 or higher and depends on the following libraries:

- `numpy <https://github.com/numpy/numpy>`_
- `scipy <https://github.com/scipy/scipy>`_
- `matplotlib <https://github.com/matplotlib/matplotlib>`_
- `astropy <https://github.com/astropy/astropy>`_
- `corner <https://github.com/dfm/corner.py>`_

and one of:

- `Nestle <https://github.com/kbarbary/nestle>`_ by Kyle Barbary.
- `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_ by Johannes Buchner, a Python interface
  for `MultiNest <https://github.com/JohannesBuchner/MultiNest>`_.

.. note::
    When using the nested sampling methods the user may choose between two packages: Nestle and PyMultiNest. The latter is generally faster and more robust, but it can be tricky to install, and thus it is not a requirement to use this code.

    Nestle is included in the ``requirements.txt`` file and will be installed automatically. In order to use PyMultiNest, you will have to follow their installation instructions `here <https://johannesbuchner.github.io/PyMultiNest/install.html>`_.