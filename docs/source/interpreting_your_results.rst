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

 If equilibrium tides dominate the evolution of the HJ system, the rate of orbital decay depends on the efficiency of tidal energy dissipation within the star :cite:p:`Goldreich1966, barker_tidal_2020`, which is typically parameterized by the star's "modified" tidal quality factor :math:`Q_\star^{'}`. From the "constant phase lag" model of :cite:t:`Goldreich1966`, the decay rate is:

 .. math::

    \dot{P}_{\mathrm{decay}} = -\frac{27\pi}{2Q_\star^{'}}\left(\frac{M_p}{M_\star}\right)\left(\frac{R_\star}{a}\right)^5

 The method :meth:`~orbdot.models.theory.quality_factor_from_decay` calculates :math:`Q_\star^{'}` from a given :math:`\dot{P}_{\mathrm{decay}}`, and the method :meth:`~orbdot.models.theory.decay_from_quality_factor` calculates :math:`\dot{P}_{\mathrm{decay}}` from a given :math:`Q_\star^{'}`.

**Empirical Quality Factors**

 The method :meth:`~orbdot.models.theory.empirical_quality_factor` estimates a value for :math:`Q_\star^{'}` with an empirical law derived by :cite:t:`Penev2018`, given by:

 .. math::

    Q^{'}_* = \max\left[{\frac{10^6}{(P_{\mathrm{tide}}/\mathrm{day}})^{3.1}}\,,\,{10^5}\right]

 where :math:`P_{\mathrm{tide}}` is the tidal forcing period of the star-planet system, given by:

 .. math::

    P_{\mathrm{tide}} = \frac{1}{2\left({P_{\mathrm{orb}}}^{-1} -{P_{\mathrm{rot}}}^{-1} \right)}

 where :math:`P_{\mathrm{orb}}` is the orbital period and :math:`P_{\mathrm{rot}}` is the rotational period of the host star.

------------

**Decay Timescale**

 The method :meth:`~orbdot.models.theory.remaining_lifetime` computes the timescale over which the orbit is shrinking for any given decay rate :math:`\dot{P}_{\mathrm{decay}}` and initial orbital period :math:`P_0` using the equation:

 .. math::

    \tau=\frac{P_0}{|\dot{P}_{\mathrm{decay}}|}

 .. autofunction:: orbdot.models.theory.remaining_lifetime

------------

**Energy and Angular Momentum Loss**

 As a planet experiences orbital decay, both the orbital energy and angular momentum will decrease over time. The methods :meth:`~orbdot.models.theory.tidal_energy_loss` and :meth:`~orbdot.models.theory.tidal_angular_momentum_loss` calculate these loss rates for any given orbital decay rate, using the following equations from :cite:t:`yee2020`:

 .. math::

    \frac{d E}{d t}=\frac{(2 \pi)^{2 / 3} M_{\mathrm{p}}}{3}\left(\frac{G M_{\star}}{P}\right)^{2 / 3} \frac{1}{P} \frac{d P}{d t}

    \frac{d L}{d t}=\frac{M_{\mathrm{p}}}{3(2 \pi)^{1 / 3}}\left(\frac{G M_{\star}}{P}\right)^{2 / 3} \frac{d P}{d t}

 where :math:`\frac{d E}{d t}` and :math:`\frac{d L}{d t}` are the loss rates of orbital energy and angular momentum, respectively.

 .. autofunction:: orbdot.models.theory.tidal_energy_loss

 .. autofunction:: orbdot.models.theory.tidal_angular_momentum_loss

------------

.. _apsidal_precession_theory:
Apsidal Precession
------------------
Apsidal precession is the gradual increase of the argument of pericentre :math:`\omega` of a planet's orbit over time, meaning the line connecting the pericentre and apocentre of the orbit rotates through :math:`2\pi` in one precession period.

This can result from several factors, including components due to general relativistic effects :cite:p:`pal_periastron_2008,jordan_observability_2008`, perturbations from other planets :cite:p:`heyl_using_2007`, and gravitational moments arising from both the host star's rotation and planetary tidal bulges :cite:p:`greenberg_apsidal_1981`. The following sections describe the equations and OrbDot methods that are relevant to these effects.

------------

**General Relativity**

 The lowest order of the relativistic contribution to apsidal precession is given by:

 .. math::

    \dot{\omega}_{\mathrm{GR}} = \frac{3 \eta G M_s}{ac^2(1-e^2)}

The method :meth:`~orbdot.models.theory.precession_gr` calculates the expected precession rate :math:`\dot{\omega}_{\mathrm{GR}}` for any given system:

 .. autofunction:: orbdot.models.theory.precession_gr

------------

**Rotational Flattening**
 Another source of apsidal precession is the rotational flattening of host stars and their planets, which is an oblate distortion that perturbs the gravitational potential.



The Love number represents how centrally condensed the body is, and is a fixed property of the body. The lower the :math:`k_2`, the more centrally condensed the planetary interior structure, which in turn leads to a slower precession rate. The theoretical upper limit of :math:`k_2` is :math:`3/2`, which corresponds to a uniform density sphere [lissauer2019]_. Note that :math:`k_2` is generally much lower for main-sequence stars [claret_love_num]_ (:math:`\sim 0.03`) than planets :cite:p:`Ragozzine2009` (0.1 -- 0.3).

OrbDott the :math:`k_2` formulation from [Ragozzine2009]_ equation X, the rotation-induced precession rate is [Ragozzine2009]_:

 .. math::

     \dot{\omega}_{\mathrm{rot,p}} = \frac{\eta {R_p}^5\,k_{2p}\,{\dot{\theta}_p}^2}{2 a^2\,G M_p} \,g_2(e)

 .. math::

     \dot{\omega}_{\mathrm{rot,s}} = \frac{\eta {R_s}^5\,k_{2s}\,{\dot{\theta}_s}^2}{2 a^2\,G M_s} \,g_2(e)

 where,

 .. math::

    g_2(e) = (1-e^2)^{-2}

 :math:`\dot{\theta}_p` and :math:`\dot{\theta}_s` represent the rotation speed of the planet and star, respectively.

.. autofunction:: orbdot.models.theory.precession_rotational_planet
  :noindex:

.. autofunction:: orbdot.models.theory.k2p_from_wdot_rot_p
  :noindex:

.. autofunction:: orbdot.models.theory.precession_rotational_star
  :noindex:

.. autofunction:: orbdot.models.theory.k2s_from_wdot_rot_s
  :noindex:

**Tidal Bulges**
 Due to the close proximity of HJs to their host stars, significant tidal bulges -- an ellipsoidal distortion -- are raised in both the planet and star. Both pairs of tidal bulges induce apsidal precession, but for Hot Jupiters, the planet's bulge is again expected to dominate :cite:p:`Ragozzine2009`. The precession rate itself depends on the internal density distribution of the HJ, which affects the extent to which the planet is elongated. This is again parameterized by the planetary Love number :math:`k_{2,p}`. For completeness, we also consider the effect of the star's tidal bulge. :cite:t:`Ragozzine2009` formulate the tides-induced precession as:

 .. math::

    \begin{aligned}
    \dot{\omega}_{\mathrm{tide}} = \dot{\omega}_{\mathrm{tide,p}} +  \dot{\omega}_{\mathrm{tide,s}} \,\,
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

.. _proper_motion_theory:
Proper Motion
-------------

.. autofunction:: orbdot.models.theory.get_wdot_pm
.. autofunction:: orbdot.models.theory.get_idot_pm
.. autofunction:: orbdot.models.theory.get_pdot_pm
.. autofunction:: orbdot.models.theory.get_tdot_pm
.. autofunction:: orbdot.models.theory.shklovskii_effect

.. _planet_companion_theory:
Planetary Companion
-------------------

.. _binary_star_theory:
Resolved Stellar Companion
--------------------------
