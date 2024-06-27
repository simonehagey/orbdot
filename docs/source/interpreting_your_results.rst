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

Theoretical Background
======================
The ``theory`` class provides several analytical models that are needed to investigate the existence and cause of long-term variations in the orbit of a planet. these include equations for quantifying the effect of energy dissipation, tidal and rotational deformation, massive
outer companions, general relativity, and more. See here for more information.

Orbital Decay
-------------
Orbital decay refers to a transfer of angular momentum from the planet to the host star that results in a shrinking of the orbital period, which we denote :math:`\dot{P}_{\rm decay}`, eventually leading to planetary engulfment. :cite:p:`levrard2009, matsumura2010`

If equilibrium tides dominate the evolution of the HJ system, the rate of orbital decay depends on the efficiency of tidal energy dissipation within the star :cite:p:`levrard2009, matsumura2010, tejada_arevalo_further_2021`, typically parameterized by the dimensionless modified stellar tidal quality factor :math:`Q_\star^{'}`.

the star's "modified tidal quality factor," defined as the quality factor  divided by 2/3 of the Love number k2.

In the "constant phase lag" model of :cite:author:`Goldreich1966`, assuming that the planet's mass stays constant, the decay rate is

:math:`\dot{P}_{\rm decay} = -\frac{27\pi}{2Q_\star^{'}}\left(\frac{M_p}{M_*}\right)\left(\frac{R_*}{a}\right)^5\rm`

.. autofunction:: orbdot.models.theory.quality_factor_from_decay
.. autofunction:: orbdot.models.theory.decay_from_quality_factor

------------

We estimate a value for :math:`Q_\star^{'}`, and thus :math:`\dot{P}_{\rm decay}` using an empirical model derived by \citet{penev_empirical_2018}, which was derived by studying a population with: :math:`$`M_p>0.1M_{\text{Jup}}`, :math:`P<3.5` days, and :math:`T_{eff,*}<6100\,K`, the parameter space of HJ systems. This derived model is given by:

:math:`Q^{'}_* = \max\left[{\frac{10^6}{(P_{\rm tide}/\text{day}})^{3.1}}\,,\,{10^5}\right]`

where :math:`P_{\rm tide}` is the tidal forcing period of the star-planet system, given by:

:math:`P_{\rm tide} = \frac{1}{2\left({P_{\rm {orb}}}^{-1} -{P_{\rm rot}}^{-1} \right)}`.

\citet{penev_empirical_2018} derived an empirical model for :math:`Q^{'}_\star` given the tidal forcing period of the star-planet system :math:`P_{\mathrm{tide}}` from an analysis of all known exoplanet systems with :math:`M_p>0.1` M:math:`_{\text{Jup}}`, :math:`P<3.5` days, and :math:`T_{eff,\star}<6100\,K`, the parameter space in which the TrES-1 system resides. They found that

:math:`Q^{'}_\star = \max\left[{\frac{10^6}{(P_{\text{tide}}/\text{day})^{3.1}}},{10^5}\right]`

where,

:math:`P_{\text{tide}} = \frac{1}{{2\left(P_{\text{orb}}^{-1} - P_{\text{spin}}^{-1}\right)}}`.

.. autofunction:: orbdot.models.theory.empirical_quality_factor

------------

The timescale over which the orbit is shrinking is:

:math:`\tau=\frac{P}{|\dot{P}|}`

.. autofunction:: orbdot.models.theory.remaining_lifetime

------------

The orbital energy and angular momentum are both decreasing, at rates of

:math:`\frac{d E}{d t}=\frac{(2 \pi)^{2 / 3} M_{\mathrm{p}}}{3}\left(\frac{G M_{\star}}{P}\right)^{2 / 3} \frac{1}{P} \frac{d P}{d t}`

.. autofunction:: orbdot.models.theory.tidal_energy_loss


:math:`\frac{d L}{d t}=\frac{M_{\mathrm{p}}}{3(2 \pi)^{1 / 3}}\left(\frac{G M_{\star}}{P}\right)^{2 / 3} \frac{d P}{d t}`

.. autofunction:: orbdot.models.theory.tidal_angular_momentum_loss

------------

Apsidal Precession
------------------
The apsidal precession rate of a planetary orbit can result from several factors, including components due to general relativistic effects, perturbations from other planets, and gravitational moments arising from both the host star's rotation and planetary tidal bulges.

General Relativity
^^^^^^^^^^^^^^^^^^

Due to their close proximity to their host star, general relativity can contribute to HJ apsidal precession rates. The lowest order of the relativistic contribution is given by:

.. math::

    \dot{\omega}_{\mathrm GR} = \frac{3 \eta G M_s}{ac^2(1-e^2)}

where :math:`G` is the gravitational constant, :math:`M_s` is the host star mass, :math:`a` is the planet's semi-major axis, :math:`c` is the speed of light in a vacuum, :math:`e` is the eccentricity, and :math:`\eta = 2\pi/P` is the mean motion.

.. autofunction:: orbdot.models.theory.precession_gr

Rotational Flattening
^^^^^^^^^^^^^^^^^^^^^

Another source of apsidal precession is the rotational flattening of host stars and their planets, as the resulting oblate distortion perturbs the gravitational potential.

Use the :math:`k_2` formulation from [Ragozzine2009]_ equation X, the rotation-induced precession rate is [Ragozzine2009]_:

.. math::

     \dot{\omega}_{\mathrm rot,p} = \frac{\eta {R_p}^5\,k_{2p}\,{\dot{\theta}_p}^2}{2 a^2\,G M_p} \,g_2(e)

.. math::

     \dot{\omega}_{\mathrm rot,s} = \frac{\eta {R_s}^5\,k_{2s}\,{\dot{\theta}_s}^2}{2 a^2\,G M_s} \,g_2(e)

where,

.. math::

    g_2(e) = (1-e^2)^{-2}

:math:`\dot{\theta}_p` and :math:`\dot{\theta}_s` represent the rotation speed of the planet and star, respectively.

.. note::

 The Love number represents how centrally condensed the body is, and is a fixed property of the body. The lower the :math:`k_2`, the more centrally condensed the planetary interior structure, which in turn leads to a slower precession rate. The theoretical upper limit of :math:`k_2` is :math:`3/2`, which corresponds to a uniform density sphere [lissauer2019]_. Note that :math:`k_2` is generally much lower for main-sequence stars [claret_love_num]_ (:math:`\sim 0.03`) than planets :cite:p:`Ragozzine2009` (0.1 -- 0.3).

.. autofunction:: orbdot.models.theory.precession_rotational_planet
.. autofunction:: orbdot.models.theory.precession_rotational_star

Tidal Bulges
^^^^^^^^^^^^

Due to the close proximity of HJs to their host stars, significant tidal bulges -- an ellipsoidal distortion -- are raised in both the planet and star. Both pairs of tidal bulges induce apsidal precession, but for Hot Jupiters, the planet's bulge is again expected to dominate :cite:p:`Ragozzine2009`. The precession rate itself depends on the internal density distribution of the HJ, which affects the extent to which the planet is elongated. This is again parameterized by the planetary Love number :math:`k_{2,p}`. For completeness, we also consider the effect of the star's tidal bulge. :cite:t:`Ragozzine2009` formulate the tides-induced precession as:

.. math::

    \begin{aligned}
    \dot{\omega}_{\mathrm tide} = \dot{\omega}_{\mathrm tide,p} +  \dot{\omega}_{\mathrm tide,s} \,\,
    &= \frac{15}{2}k_{2,p}n\left(\frac{R_p}{a}\right)^5\left(\frac{M_s}{M_p}\right)f_2(e) \\
    &+ \frac{15}{2}k_{2,s}n\left(\frac{R_s}{a}\right)^5\left(\frac{M_p}{M_s}\right)f_2(e),
    \end{aligned}

where,

.. math::

    \begin{aligned}
    f_2(e)= (1-e^2)^{-5} \left(1 + \frac{3}{2}e^2 + \frac{1}{8}e^4 \right).
    \end{aligned}

.. autofunction:: orbdot.models.theory.precession_tidal_planet
.. autofunction:: orbdot.models.theory.precession_tidal_star

Proper Motion
-------------

.. autofunction:: orbdot.models.theory.get_wdot_pm
.. autofunction:: orbdot.models.theory.get_idot_pm
.. autofunction:: orbdot.models.theory.get_pdot_pm
.. autofunction:: orbdot.models.theory.get_tdot_pm
.. autofunction:: orbdot.models.theory.shklovskii_effect

Planetary Companion
-------------------


Resolved Stellar Companion
--------------------------
