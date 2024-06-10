.. _interpreting-results:

*************************
Interpreting Your Results
*************************

The Analysis Class
==================
The :class:`orbdot.analysis.Analysis` class is designed to facilitate and interpret various analyses related to the transit timing and orbital
dynamics of exoplanets. It processes data from a given planetary system and provides functionalities to compute and
summarize effects such as proper motion, orbital decay, and apsidal precession.

The `Analysis` class provides a comprehensive toolset for analyzing and interpreting the orbital dynamics and
transit timing variations of exoplanets. By leveraging empirical laws and model fits, it helps researchers understand
the intricate gravitational interactions within exoplanetary systems and predict future behavior. The class's methods
are designed to be user-friendly, offering detailed summaries of results and the underlying physics.

For further details and examples, refer to the source code and the referenced paper by Penev et al. (2018).

To use the `Analysis` class, you need an instance of a planet class containing system parameters and observational
data. Additionally, you need a dictionary containing the results of the model fit. Here is an example of how to
initialize the `Analysis` class:

.. code-block:: python

    planet_instance = ...  # Your planet instance
    results_dic = ...      # Your results dictionary

    analysis = Analysis(planet_instance, results_dic)



Key Methods
------------

1. **Proper Motion Analysis**:

 Computes and summarizes predicted transit timing variations (TTVs) and transit duration variations (TDVs) due to systemic proper motion, use the `proper_motion` method:

 .. code-block:: python

    analysis.proper_motion(printout=True)


2. **Orbital Decay Prediction**:

 Computes and summarizes predicted orbital decay parameters based on an empirical law for the stellar tidal quality factor, use the `orbital_decay_predicted` method:

 .. code-block:: python

    analysis.orbital_decay_predicted(printout=True)

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
   | Condition                        | Evidence for Model 1 (M:math:`_1`)                |
   +==================================+===================================================+
   | :math:`B_{12} \leq 1`            | M:math:`_1` is not supported over M:math:`_2`     |
   +----------------------------------+---------------------------------------------------+
   | :math:`1 < B_{12} \leq 3`        | Evidence for M:math:`_1` barely worth mentioning  |
   +----------------------------------+---------------------------------------------------+
   | :math:`3 < B_{12} \leq 20`       | Positive evidence for M:math:`_1`                 |
   +----------------------------------+---------------------------------------------------+
   | :math:`20 < B_{12} \leq 150`     | Strong evidence for M:math:`_1`                   |
   +----------------------------------+---------------------------------------------------+
   | :math:`150 < B_{12}`             | Very strong evidence for M:math:`_1`              |
   +----------------------------------+---------------------------------------------------+


Theoretical Background
======================
The ``theory`` class provides several analytical models that are needed to investigate the existence and cause of long-term variations in the orbit of a planet. these include equations for quantifying the effect of energy dissipation, tidal and rotational deformation, massive
outer companions, general relativity, and more. See here for more information.

Orbital Decay
-------------


Apsidal Precession
------------------
The apsidal precession rate of a planetary orbit can result from several factors, including components due to general relativistic effects, perturbations from other planets, and gravitational moments arising from both the host star's rotation and planetary tidal bulges.

**General Relativity**:

 Due to their close proximity to their host star, general relativity can contribute to HJ apsidal precession rates. The lowest order of the relativistic contribution is given by:

 .. math::

    \dot{\omega}_{\mathrm GR} = \frac{3 \eta G M_s}{ac^2(1-e^2)}

 where :math:`G` is the gravitational constant, :math:`M_s` is the host star mass, :math:`a` is the planet's semi-major axis, :math:`c` is the speed of light in a vacuum, :math:`e` is the eccentricity, and :math:`\eta = 2\pi/P` is the mean motion.

**Rotational Flattening**:

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

 The Love number represents how centrally condensed the body is, and is a fixed property of the body. The lower the :math:`k_2`, the more centrally condensed the planetary interior structure, which in turn leads to a slower precession rate. The theoretical upper limit of :math:`k_2` is :math:`3/2`, which corresponds to a uniform density sphere [lissauer2019]_. Note that :math:`k_2` is generally much lower for main-sequence stars [claret_love_num]_ (:math:`\sim 0.03`) than planets [Ragozzine2009]_ (0.1 -- 0.3).

**Tidal Bulges**:

 Due to the close proximity of HJs to their host stars, significant tidal bulges -- an ellipsoidal distortion -- are raised in both the planet and star. Both pairs of tidal bulges induce apsidal precession, but for Hot Jupiters, the planet's bulge is again expected to dominate [Ragozzine2009]_. The precession rate itself depends on the internal density distribution of the HJ, which affects the extent to which the planet is elongated. This is again parameterized by the planetary Love number :math:`k_{2,p}`. For completeness, we also consider the effect of the star's tidal bulge. [Ragozzine2009]_ formulate the tides-induced precession as:

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


Proper Motion
-------------


Planetary Companion
-------------------


Resolved Stellar Companion
--------------------------

    """

    The only star-planet parameters that are necessary to fully
    utilize the :class:`Analysis` class are:

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