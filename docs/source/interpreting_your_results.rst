.. _interpreting-results:

************************
Interpreting the Results
************************

The Analysis Class
==================
The :class:`~orbdot.analysis.Analysis` class is designed to facilitate and interpret various analyses related to the model fits. It combines the results, star-planet system info, and data together to compute and summarize effects such as proper motion, orbital decay, and apsidal precession.

To use the :class:`~orbdot.analysis.Analysis`  class, you need an instance of a StarPlanet class and a dictionary containing the results of the model fit. the dictionary can either be passed in directly from the model fit in the script, or it can be read from a preexisting file. Either way, however, you still need to hvae a planet instance.

In the script right after a model fit:

.. code-block:: python

    analysis = Analysis(planet_instance, results_dic)

From a pre-existing results file:

.. code-block:: python

    analysis = Analysis(planet_instance, results_dic)

As soon as you make an analysis object a file is made to summarize what you do with it. This file is named after the model and whatever suffix you chose. For example...

Also an analysis directory is made.

The following methods will add to the file and print to the console if the argument ``printout=True``.

Key Methods
------------

1. **Proper Motion Analysis**:

 Computes and summarizes predicted transit timing variations (TTVs) and transit duration variations (TDVs) due to systemic proper motion, use the `proper_motion` method:

Running this:

 .. code-block:: python

    ttv_c = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')
    a = Analysis(wasp12, ttv_c)
    proper_motion()

yields the following output in the file

.. code-block:: text

    Systemic Proper Motion Effects
    -----------------------------------------------------------------
     * Apparent apsidal precession rate due to proper motion:
          maximum: dw/dt = 3.46E-08 rad/yr
          minimum: dw/dt = -2.12E-24 rad/yr
     * Apparent rate of change of the inclination due to proper motion:
          maximum: di/dt = -3.48E-08 rad/yr
          minimum: di/dt = -4.27E-24 rad/yr
     * Transit duration variation due to proper motion:
          maximum: dT/dt = -4.71E-01 ms/yr
          minimum: dT/dt = -5.77E-17 ms/yr
     * Apparent orbital period drift due to proper motion:
          maximum: dP/dt = -5.37E-11 ms/yr
     * Apparent orbital period drift due to the Shklovskii effect:
          maximum: dP/dt = 1.58E-04 ms/yr


2. **Orbital Decay Prediction**:

 Computes and summarizes predicted orbital decay parameters based on an empirical law for the stellar tidal quality factor, use the `orbital_decay_predicted` method:

 .. code-block:: python

    analysis.orbital_decay_predicted()

3. **Orbital Decay Model Fit Interpretation**:

 Provides a summary of various interpretations of the results of an orbital decay model fit.

4. **Apsidal Precession Model Fit Interpretation**:

 Provides a summary of various interpretations of the results of an apsidal precession model fit, use the `apsidal_precession_fit` method:

 .. code-block:: python

    analysis.apsidal_precession_fit(printout=True)

4. **Model Comparison**:

 To determine the preferred model, the Baye's factor is compared to the thresholds established by Kass and Raftery (1995), tabulated below.

 .. math::

    \log{B_{12}} = \log{\mathrm{Z}}_{1} - \log{\mathrm{Z}}_{2}

 To compare two models, this method calculate the Bayes factor, denoted as:

uses Bayesian evidence, denoted as $\log{\mathrm{Z}}$, as a fundamental metric for comparing the outcomes of various model fits. A lower $\log{\mathrm{Z}}$ value signifies a superior fit to the observed data.

 .. table::
   :name: tab:bayesian_evidence
   :width: 80%
   :align: center

   +----------------------------------+---------------------------------------------------+
   | Condition                        | Evidence for Model 1 (Model 1)                    |
   +==================================+===================================================+
   | :math:`B_{12} \leq 1`            | Model 1 is not supported over Model 2             |
   +----------------------------------+---------------------------------------------------+
   | :math:`1 < B_{12} \leq 3`        | Evidence for Model 1 barely worth mentioning      |
   +----------------------------------+---------------------------------------------------+
   | :math:`3 < B_{12} \leq 20`       | Positive evidence for Model 1                     |
   +----------------------------------+---------------------------------------------------+
   | :math:`20 < B_{12} \leq 150`     | Strong evidence for Model 1                       |
   +----------------------------------+---------------------------------------------------+
   | :math:`150 < B_{12}`             | Very strong evidence for Model 1                  |
   +----------------------------------+---------------------------------------------------+

.. code-block:: text

    """

    The only star-planet parameters that are necessary to fully
    utilize the :class::class:`~orbdot.analysis.Analysis` class are:

        --> The star mass 'M_s' in solar masses.
        --> The star radius 'R_s' in solar radii.
        --> The planet mass 'M_p' in earth masses.
        --> The planet radius 'R_p' in earth radii.

    System Parameters
    ------------------
        - 'RA': Right ascension of the system.
        - 'DEC': Declination of the system.
        - 'mu': Proper motion of the system (mas/yr).
        - 'mu_RA': Proper motion in right ascension (mas/yr).
        - 'mu_DEC': Proper motion in declination (mas/yr).
        - 'parallax': Parallax of the system (mas).
        - 'distance': Distance to the system (pc).
        - 'v_r': Radial velocity of the system (km/s).
        - 'gaia_dr3_ID': Gaia DR3 identifier of the system.
        - 'discovery_year': Year of discovery of the system.

    Stellar Parameters
    ------------------
        - 'spectral_type': spectral type of the star.
        - 'age': age of the star (Gyr).
        - 'm_v': visual magnitude of the star.
        - 'Teff': effective temperature of the star (K).
        - 'M_s': mass of the star (Solar masses).
        - 'R_s': radius of the star (Solar radii).
        - 'metallicity': metallicity of the star ([Fe/H]).
        - 'log_g': suresultsace gravity of the star (log10(cm/s^2)).
        - 'rho_s': density of the star (g/cm^3).
        - 'k2_s': dimensionless tidal Love number k2 of the star.
        - 'vsini': projected rotational velocity of the star (km/s).
        - 'P_rot_s': rotation period of the star (days).
        - 'luminosity': luminosity of the star (W).

    Planet Parameters
    -----------------
        - 'sm_axis': semi-major axis of the planet's orbit (AU).
        - 'M_p': mass of the planet (Earth masses).
        - 'R_p': radius of the planet (Earth radii).
        - 'rho_p': density of the planet (g/cm^3).
        - 'P_rot_p': rotation period of the planet (days).
        - 'k2_p': dimensionless tidal Love number k2 of the planet.
        - 'T_eq': equilibrium temperature of the planet (K).
        - 'lambda': obliquity of the planet (degrees).
        - 'Psi': longitude of the ascending node of the planet (degrees).

    General Orbit Parameters
    ------------------------
        - 't0': time of periastron passage (BJD_TDB).
        - 'P0': orbital period (days).
        - 'e0': orbital eccentricity.
        - 'w0': argument of periastron (radians).
        - 'i0': orbital inclination (degrees).
        - 'O0': longitude of the ascending node (radians).

    Time-Dependent Parameters
    -------------------------
        - 'PdE': derivative of orbital period with respect to eccentric anomaly (days/E).
        - 'wdE': derivative of argument of periastron with respect to eccentric anomaly (radians/E).
        - 'edE': derivative of eccentricity with respect to eccentric anomaly (/E).
        - 'idE': derivative of inclination with respect to eccentric anomaly (degrees/E).
        - 'OdE': derivative of longitude of the ascending node with respect to eccentric anomaly (radians/E).

    Radial Velocity Parameters
    --------------------------
        - 'K': radial velocity semi-amplitude (m/s).
        - 'v0': zero-point offset in radial velocity (m/s).
        - 'jit': radial velocity jitter (m/s).
        - 'dvdt': linear trend in radial velocity (m/s/day).
        - 'ddvdt': quadratic trend in radial velocity (m/s^2/day).

    """

------------

Theoretical Background
======================
The :py:mod:`~orbdot.models.theory` module provides several analytical models that can be used to investigate long-term variations in planetary orbits and their causes, such as :ref:`orbital decay <orbital_decay_theory>`, :ref:`apsidal precession <_apsidal_precession_theory>`, nonresonant :ref:`companion planets <planet_companion_theory>`, bound :ref:`stellar companions <binary_star_theory>`, and systemic :ref:`proper motion <proper_motion_theory>`.

The following sections briefly introduce these things...

------------

.. _orbital_decay_theory:
Orbital Decay
-------------
Orbital decay refers to a transfer of angular momentum from the planet to the host star that results in a shrinking of the orbital period, which we denote :math:`\dot{P}_{\mathrm{decay}}`, eventually leading to planetary engulfment.

------------

**Equilibrium Tides**

 .. autofunction:: orbdot.models.theory.decay_quality_factor_from_pdot
 .. autofunction:: orbdot.models.theory.decay_pdot_from_quality_factor

**Empirical Quality Factor**

 .. autofunction:: orbdot.models.theory.decay_empirical_quality_factor

------------

**Decay Timescale**

 .. autofunction:: orbdot.models.theory.decay_timescale

------------

**Energy and Angular Momentum Loss**

 .. autofunction:: orbdot.models.theory.decay_energy_loss
 .. autofunction:: orbdot.models.theory.decay_angular_momentum_loss

------------

.. _apsidal_precession_theory:
Apsidal Precession
------------------
Apsidal precession is the gradual increase of the argument of pericentre :math:`\omega` of a planet's orbit over time, meaning the line connecting the pericentre and apocentre of the orbit rotates through :math:`2\pi` in one precession period.

This can result from several factors, including components due to general relativistic effects :cite:p:`pal_periastron_2008,jordan_observability_2008`, perturbations from other planets :cite:p:`heyl_using_2007`, and gravitational moments arising from both the host star's rotation and planetary tidal bulges :cite:p:`greenberg_apsidal_1981`. The following sections describe the equations and OrbDot methods that are relevant to these effects.

------------

**Transit Variations**
 .. autofunction:: orbdot.models.theory.get_tdot_from_wdot
 .. autofunction:: orbdot.models.theory.get_pdot_from_wdot

------------

**General Relativity**

 .. autofunction:: orbdot.models.theory.precession_gr

------------

**Rotational Flattening**

 .. autofunction:: orbdot.models.theory.precession_rotational_planet
 .. autofunction:: orbdot.models.theory.precession_rotational_planet_k2
 .. autofunction:: orbdot.models.theory.precession_rotational_star
 .. autofunction:: orbdot.models.theory.precession_rotational_star_k2

------------

**Tidal Bulges**

 .. autofunction:: orbdot.models.theory.precession_tidal_planet
 .. autofunction:: orbdot.models.theory.precession_tidal_planet_k2
 .. autofunction:: orbdot.models.theory.precession_tidal_star
 .. autofunction:: orbdot.models.theory.precession_tidal_star_k2

.. _proper_motion_theory:
Proper Motion
-------------
The apparent secular evolution of exoplanet transit signatures that are induced by the systemic proper motion, which is the movement of the star-planet system with respect to reference frame of the Solar System. This motion in 3D space is partially constrained with measurements of the proper motion on the sky-plane :math:`\mu` and radial velocity :math:`v_r`.

.. autofunction:: orbdot.models.theory.proper_motion_wdot
.. autofunction:: orbdot.models.theory.proper_motion_idot
.. autofunction:: orbdot.models.theory.proper_motion_pdot
.. autofunction:: orbdot.models.theory.proper_motion_tdot
.. autofunction:: orbdot.models.theory.proper_motion_shklovskii

.. _planet_companion_theory:
Planetary Companion
-------------------

.. _binary_star_theory:
Resolved Stellar Companion
--------------------------
