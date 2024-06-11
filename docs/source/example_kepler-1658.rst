.. _example-kepler1568:

**********************************
The Orbital Decay of Kepler-1658 b
**********************************

This example is a reproduction of the results from XXX using OrbDot. In this study, the authors found that...

Setup
=====
To set this up, we need to take the data from the paper, replicate their priors, and populate the system info file.

Data
----

.. code-block:: text

    Epoch BJD Err_BJD Source
    -12 2454959.7314 0.0015 "Kepler LC"
    -6 2454982.82835 0.00061 "Kepler LC"
    11 2455048.26751 0.00022 "Kepler SC"
    35 2455140.65189 0.00042 "Kepler LC"
    59 2455233.03736 0.00035 "Kepler LC"
    83 2455325.42133 0.00036 "Kepler LC"
    131 2455510.19192 0.00023 "Kepler SC"
    155 2455602.57708 0.00027 "Kepler SC"
    178 2455691.11211 0.00031 "Kepler LC"
    228 2455883.58121 0.00033 "Kepler LC"
    252 2455975.96583 0.00036 "Kepler LC"
    275 2456064.50087 0.00036 "Kepler LC"
    325 2456256.97026 0.00037 "Kepler LC"
    349 2456349.35438 0.00039 "Kepler LC"
    365 2456410.94385 0.00064 "Kepler LC"
    1063 2459097.8002 0.0015 "Palomar/WIRC"
    1151 2459436.5407 0.0023 "TESS"
    1243 2459790.6819 0.0013 "Palomar/WIRC"
    1242 2459786.8359 0.0030 "TESS"
    1249 2459813.7791 0.0027 "TESS"

Priors
------

System Info File
----------------


Running Model Fits
==================

Constant-Period TTV Model
-------------------------
.. code-block:: python

 ttv_fit_c = sp.run_ttv_fit(['t0', 'P0'], model='constant')

.. code-block:: text

    Stats
    -----
    Sampler: nestle
    Free parameters: ['t0' 'P0']
    log(Z) = -43.38 ± 0.09
    Run time (s): 2.16
    Num live points: 1000
    Evidence tolerance: 0.1
    Eff. samples per second: 1810

    Results
    -------
    t0 = 2455005.924893011 + 0.00012903939932584763 - 0.00012579653412103653
    P0 = 3.8493666077754116 + 5.376794542932828e-07 - 5.947995269650619e-07


Orbital Decay TTV Model
-----------------------

.. code-block:: python

 ttv_fit_d = sp.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')

.. code-block:: text

    Stats
    -----
    Sampler: nestle
    Free parameters: ['t0' 'P0' 'PdE']
    log(Z) = -21.15 ± 0.11
    Run time (s): 3.49
    Num live points: 1000
    Evidence tolerance: 0.1
    Eff. samples per second: 1345

    Results
    -------
    t0 = 2455005.9242031113 + 0.0001594102941453457 - 0.00016830721870064735
    P0 = 3.8493733034167033 + 1.1344747679054024e-06 - 1.1537980562081884e-06
    PdE = -1.60906722485075e-08 + 2.283352281141653e-09 - 2.2936184207406122e-09
    dPdt (ms/yr) = -131.9131605393512 + 18.719181608964227 - 18.803344588616152

Apsidal Precession TTV Model
----------------------------

.. code-block:: python

 ttv_fit_a = sp.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

.. code-block:: text

    Stats
    -----
    Sampler: nestle
    Free parameters: ['t0' 'P0' 'e0' 'w0' 'wdE']
    log(Z) = -23.39 ± 0.12
    Run time (s): 208.68
    Num live points: 1000
    Evidence tolerance: 0.1
    Eff. samples per second: 23

    Results
    -------
    t0 = 2455005.9246068085 + 0.0028900164179503918 - 0.0064444104209542274
    P0 = 3.8493588745373857 + 5.480888950248897e-06 - 2.3021433830372473e-05
    e0 = 0.007656003736720481 + 0.021499437949721097 - 0.005828913924036578
    w0 = 1.5433150534526858 + 0.5575896087334574 - 0.6313063693238424
    wdE = 0.0016103164922047814 + 0.001830668502334377 - 0.0006607950175296083


Interpretation
==============
We can see that


