.. _installation:

************
Installation
************

You can pip install it:

...

Or you can clone the repository:


Dependencies
------------

OrbDot is dependent on the following packages:

- `corner <https://github.com/dfm/corner.py>`_
- `matplotlib <https://github.com/matplotlib/matplotlib>`_
- `numpy <https://github.com/numpy/numpy>`_
- `scipy <https://github.com/scipy/scipy>`_
- `astropy <https://github.com/astropy/astropy>`_

and one of:

- `Nestle <https://github.com/kbarbary/nestle>`_ by Kyle Barbary.
- `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_ by Johannes Buchner, a Python interface
  for `MultiNest <https://github.com/JohannesBuchner/MultiNest>`_.

Nestle or PyMultiNest?
----------------------
To perform the nested sampling methods the user may choose between two packages: Nestle [1]_
and PyMultiNest [2]_. PyMultiNest is generally faster and more robust, but it can be tricky to
install, thus it is not a requirement to use this code. The desired sampler is specified in the
settings file as 'nestle' or 'multinest'.

The Nestle package is imported within this function so that it does not need to be installed if the user already uses PyMultiNest.

.. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
.. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest/
