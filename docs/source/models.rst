.. _models:

******
Theory
******
This section describes the models that can be fit, and also all the equations for theoretical calculations of various effects.

Section XX goes over the coordinate system, section XX reviews the transit and eclipse timing models, section XX discusses radial velocity models, and section XX details the theoretical calculations for long-term orbital variations and their causes.
.. _coordinate_system:

Coordinate System
=================
It is critical to be consistent in the definition of the argument of pericenter :math:`\omega` when simultaneously fitting transit and eclipse mid-times and radial velocities. In the OrbDot coordinate system, the argument of pericenter is determined from the positive x-axis, such that a transit occurs when the true anomaly :math:`\phi` is equal to:

.. math::
 \phi_{\mathrm I}\,=\,\frac{\pi}{2} - \omega_{\mathrm p}

and an eclipse occurs when:

.. math::
 \phi_{\mathrm II} = \frac{3\pi}{2} - \omega_{\mathrm p}

The OrbDot models are written in a coordinate system in which the sky plane lies on the x-z plane and the y-axis points toward the observer along the line of sight.

.. image:: _static/coordinate_system.png

------------

.. _model_parameters:

Model Parameters
================
To keep the use of this class clean and concise, there is a main set of "allowed" model parameters. In every model fit, the list of free parameters is compared to this set and the given order is recorded. This means that any number of free parameters can be provided in any order as long as they are part of the physical model(s) defined in the log-likelihoods.

It is important for the user to familiarize themselves with the parameter symbols and definitions, which are described in the following table. Some of the parameters listed below are not currently implemented in the OrbDot models, but are planned for future integration. These parameters are indicated with a footnote [*].

**Orbital Elements**
 .. list-table::
   :header-rows: 1
   :widths: 30 25 25 100

   * - Parameter
     - Symbol
     - Unit
     - Description
   * - :math:`t_0`
     - ``t0``
     - BJD
     - The reference transit mid-time.
   * - :math:`P_0`
     - ``P0``
     - days
     - The observed (sidereal) orbital period.
   * - :math:`e`
     - ``e0``
     -
     - The eccentricity of the orbit.
   * - :math:`\omega_p`
     - ``w0``
     - radians
     - The argument of pericenter of the planet's orbit.
   * - :math:`i`
     - ``i0``
     - degrees
     - The line-of-sight inclination of the orbit.
   * - :math:`\Omega`
     - ``O0``
     - radians
     - The longitude of the ascending node. [*]_

**Coupled Parameters**
 .. list-table::
   :header-rows: 1
   :widths: 30 25 20 100

   * - Parameter
     - Symbol
     - Unit
     - Description
   * - :math:`e\,\cos{\,\omega_p}`
     - ``ecosw``
     -
     - The eccentricity :math:`e` multiplied by the cosine of :math:`\omega_p`.
   * - :math:`e\,\sin{\,\omega_p}`
     - ``esinw``
     -
     - the eccentricity :math:`e` multiplied by the sine of :math:`\omega_p`.
   * - :math:`\sqrt{e}\,\cos{\,\omega_p}`
     - ``sq_ecosw``
     -
     - The square root of :math:`e` multiplied by the cosine of :math:`\omega_p`.
   * - :math:`\sqrt{e}\,\sin{\,\omega_p}`
     - ``sq_esinw``
     -
     - The square root of :math:`e` multiplied by the sine of :math:`\omega_p`.

**Time-Dependent Parameters**
 .. list-table::
    :header-rows: 1
    :widths: 25 15 25 100

    * - Parameter
      - Symbol
      - Unit
      - Description
    * - :math:`\frac{dP}{dE}`
      - ``PdE``
      - days :math:`E^{-1}`
      - A constant change of the orbital period.
    * - :math:`\frac{d \omega}{dE}`
      - ``wdE``
      - rad :math:`E^{-1}`
      - A constant change of the argument of pericenter.
    * - :math:`\frac{de}{dE}`
      - ``edE``
      - :math:`E^{-1}`
      - A constant change of the orbital eccentricity. [*]_
    * - :math:`\frac{di}{dE}`
      - ``idE``
      - deg :math:`E^{-1}`
      - A constant change of the line-of-sight inclination. [*]_
    * - :math:`\frac{d \Omega}{dE}`
      - ``OdE``
      - rad :math:`E^{-1}`
      - A constant change of the long. of the ascending node. [*]_

**Radial Velocity Parameters**
 .. list-table::
   :header-rows: 2
   :widths: 20 15 25 70

   * - Parameter
     - Symbol
     - Unit
     - Description
   * - :math:`K`
     - ``K``
     - m :math:`{\mathrm s}^{-1}`
     - The radial velocity semi-amplitude.
   * - :math:`\gamma_j`
     - ``v0``
     - m :math:`{\mathrm s}^{-1}`
     - An instrument specific systemic radial velocity.
   * - :math:`\sigma_j`
     - ``jit``
     - m :math:`{\mathrm s}^{-1}`
     - An instrument-specific radial velocity "jitter" term.
   * - :math:`\dot{\gamma}`
     - ``dvdt``
     - m :math:`{\mathrm s}^{-1}` :math:`{\mathrm day}^{-1}`
     - A linear radial velocity trend.
   * - :math:`\ddot{\gamma}`
     - ``ddvdt``
     - m :math:`{\mathrm s}^{-1}` :math:`{\mathrm day}^{-2}`
     - A second order radial velocity trend.
   * - :math:`K_{\mathrm{tide}}`
     - ``K_tide``
     - m :math:`{\mathrm s}^{-1}`
     - The amplitude of a tidal radial velocity signal from the star. [*]_

.. [*] Not currently implemented in the OrbDot models.

------------

Transit and Eclipse Timing Models
=================================
Secular variations in the orbits of transiting exoplanets may be detected by measuring deviations in the observed transit and eclipse timing from what is expected of a circular, unchanging orbit. These deviations are commonly referred to as Transit Timing Variations (TTVs).

OrbDot currently supports model fitting for three evolutionary cases:

    1. An unchanging orbit that is circular or eccentric.
    2. A constant evolution of the orbital period, :math:`\dot{P}` (orbital decay).
    3. A constant evolution of the argument of pericenter, :math:`\dot{\omega}` (apsidal precession).

Constant-Period
---------------
For a planet on a circular orbit with a constant orbital period, we expect a linear increase in the center times of both transits :math:`t_{\mathrm{I}}` and eclipses :math:`t_{\mathrm{II}}`:

.. math:: t_{\mathrm{I}} = t_0 + PE

and,

.. math:: t_{\mathrm{II}} = t_0 + PE + \frac{P}{2}

where :math:`t_0` is the reference transit time, :math:`P` is the orbital period, and  :math:`E` is the epoch, which represents the number of orbits that have passed since time :math:`t_0`.

If the orbit is eccentric, we add must add an offset to the eclipse times to account for the variable speed of the planet:

.. math:: t_{\mathrm{I}} = t_0 + PE

.. math:: t_{\mathrm{II}} = t_0 + PE + \frac{P}{2} + \frac{P_a\,e}{\pi}\,\cos{\,\omega_p}

The :meth:`~orbdot.models.ttv_models.ttv_constant` method implements this model:

.. autofunction:: orbdot.models.ttv_models.ttv_constant

Orbital Decay
-------------
For a planet on a decaying circular orbit with a constant orbital period, we expect the mid-times of the transits (:math:`t_{\mathrm{I}}`) and eclipses (:math:`t_{\mathrm{II}}`) to be:

.. math:: t_{\mathrm{I}} = t_0 + PE + \frac{1}{2}\,\frac{dP}{dE}\,E^2

and,

.. math:: t_{\mathrm{II}} = t_0 + PE + \frac{P}{2} + \frac{1}{2}\,\frac{dP}{dE}\,E^2

where :math:`t_0` is the reference transit time, :math:`P` is the orbital period, :math:`dP/dE` is the rate of change of the period in units of days per epoch, and :math:`E` is the epoch, which represents the number of orbits that have passed since time :math:`t_0`.

If the orbit is eccentric, an offset of :math:`\frac{P_a\,e}{\pi}\,\cos{\,\omega_p}` is added to the eclipse times.

The :meth:`~orbdot.models.ttv_models.ttv_decay` method implements this model:

.. autofunction:: orbdot.models.ttv_models.ttv_decay

Apsidal Precession
------------------
For a planet on an elliptical orbit undergoing apsidal precession, we expect the
mid-times of the transits (:math:`t_{\mathrm{I}}`) and eclipses (:math:`t_{\mathrm{II}}`)
to be :cite:p:`Gimenez1995, Patra2017`:

.. math:: t_{\mathrm{I}} = t_0 + P_s E - \frac{e P_a}{\pi}\cos{\omega_p}

and,

.. math:: t_{\mathrm{II}} = t_0 + P_s E+ \frac{P_a}{2} + \frac{eP_a}{\pi}\cos{\omega_p}

where :math:`t_0` is the reference transit time, :math:`e` is the orbit eccentricity,
:math:`\omega_p` is the argument of pericentre, and :math:`E` is the epoch, which represents
the number of orbits that have passed since time :math:`t_0`. We assume that :math:`\omega_p`
evolves at a constant rate, denoted as :math:`d\omega/dE`, such that any given epoch
:math:`\omega_p` is given by

.. math:: \omega_{p}\left(E\right) = \omega_0 + \frac{d\omega}{dE}\,E

where :math:`\omega_0` is the value of :math:`\omega_p` at time :math:`t_0`.

In the equations above :math:`P_a` represents the anomalistic orbital period -- ie. the elapsed
time between subsequent pericentre passages, which characterizes the osculating orbit -- and
:math:`P_s` is the sidereal period. The latter represents the observed orbital period of the
system, and is related to the anomalistic period by:

.. math:: P_s = P_a\left(1-\frac{d\omega/{dE}}{2\pi}\right)

The :meth:`~orbdot.models.ttv_models.ttv_precession` method implements this model:

.. autofunction:: orbdot.models.ttv_models.ttv_precession

------------

Radial Velocity Models
======================
This module provides functions for modelling exoplanet radial velocity observations. OrbDot
currently supports model fitting for three evolutionary cases:

    1. An unchanging orbit that is circular or eccentric ("constant-period").
    2. A constant evolution of the orbital period, :math:`\dot{P}` ("orbital decay").
    3. A constant evolution of the argument of pericenter, :math:`\dot{\omega}` ("apsidal precession").

The following section provides model-independent background information, with details of the
implementation provided in the docstrings of the individual methods.

Background
----------
Many exoplanet host stars exhibit noticeable periodic radial velocity (RV) variations as they
"wobble" around the center of mass of the star-planet system. The amplitude of this effect
can be expressed as [1]_:

.. math::
    K=\left(\frac{2 \pi G}{P}\right)^{1/3} \frac{M_p \sin i}{\left(M_{
    \star}+M_p\right)^{2/3}} \frac{1}{\left(1-e^2\right)^{1/2}}

where :math:`M_p` is the planet mass, :math:`M_{\star}` is the star mass, :math:`i` is
the line-of-sight inclination of the orbit, and :math:`G` is the universal gravitational
constant.

At any given time, the periodic signal depends on the planet's position in its orbit,
defined by the true anomaly :math:`\phi`, and the systemic velocity along the line-of-sight,
denoted as :math:`\gamma`. When RV observations from multiple instruments are combined,
it is standard practice to fit :math:`\gamma` individually for each source to account for
instrumental variations. Additionally, we consider first and second-order acceleration terms,
:math:`\dot{\gamma}` and :math:`\ddot{\gamma}` respectively, to accommodate potential
perturbations from an outer, non-resonant companion planet with an orbital period longer than
the observational baseline. The total observed radial velocity signal is:

.. math::
    v_r = K[\cos{(\phi\left(t\right)+\omega_p)}+e\cos{\omega_p}] + \gamma_j + \dot{
    \gamma} \left(t-t_0\right) + \ddot{\gamma}\left(t-t_0\right)^2

where :math:`\omega_\star` is the argument of pericenter of the star's orbit, defined as
:math:`\omega_\star = \omega_p + \pi`, and :math:`t_0` is a transit mid-time. OrbDot
requires the reference time to that of a transit so that the phase at any time :math:`t` can be
determined by the knowledge that when :math:`t=t_0` the true anomaly is :math:`\phi_0 =
\frac{\pi}{2} - \omega_p`.

To avoid underestimating the uncertainties of the parameters in the RV model fit,
OrbDot implements an instrument-dependent 'jitter' parameter instrument into the log
likelihood (when ``jit`` is given as a free parameter), to account for systematic noise.
These terms are added in quadrature with the individual measurement errors:

.. math:: \sigma = \sqrt{\sigma_i^2 + \sigma_j^2}

where :math:`\sigma_j` is the individual measurement error and :math:`\sigma_j` is the
instrument-specific jitter term.

Constant-Period
---------------
.. autofunction:: orbdot.models.rv_models.rv_constant

Orbital Decay
-------------
.. autofunction:: orbdot.models.rv_models.rv_decay

Apsidal Precession
------------------
.. autofunction:: orbdot.models.rv_models.rv_precession

------------

Transit Duration Models
=======================

.. attention::

   Under active development and testing.

------------

.. _theory_module:

The ``theory`` Module
=====================
The :py:mod:`~orbdot.models.theory` module provides several analytical models that can be used to investigate long-term variations in planetary orbits and their causes, some of which are used in the :class:`~orbdot.analysis.Analyzer` class (see :ref:`interpreting-results`). The key methods of this class are listed below, organized by subject.

.. _orbital_decay_theory:

Equilibrium Tides
-----------------
Due to the close proximity of HJs to their host stars, significant tidal bulges -- an ellipsoidal distortion -- are raised in both the planet and star. In the case of orbital decay, the planet orbital rate is faster than the star's rotational rate. As a result, the star's tidal bulge lags behind the HJ, creating a net torque that spins up the star at the expense of the planet's orbital angular momentum :cite:p:`levrard2009, Matsumura2010`. The tidal forces raised by the misaligned tidal bulges are known as `equilibrium tides' and are believed to be the most significant process governing the future evolution of HJ orbits :cite:p:`Ma2021, Barker2020`. The following are OrbDot methods that are relevant to these effects.

.. autofunction:: orbdot.models.theory.decay_pdot_from_quality_factor
.. autofunction:: orbdot.models.theory.decay_quality_factor_from_pdot
.. autofunction:: orbdot.models.theory.decay_empirical_quality_factor
.. autofunction:: orbdot.models.theory.decay_timescale
.. autofunction:: orbdot.models.theory.decay_energy_loss
.. autofunction:: orbdot.models.theory.decay_angular_momentum_loss

------------

.. _apsidal_precession_theory:

Apsidal Precession
------------------
Apsidal precession is the gradual increase of the argument of pericenter of a planet's orbit :math:`\omega_p` over time, meaning the line connecting the pericenter and apocenter of the orbit rotates through :math:`2\pi` in one precession period. This can result from several factors, including components due to general relativistic effects, perturbations from other planets, and gravitational moments arising from both the host star's rotation and planetary tidal bulges. The following sections describe the equations and OrbDot methods that are relevant to these effects.

General Relativity
^^^^^^^^^^^^^^^^^^
.. autofunction:: orbdot.models.theory.precession_gr

Rotation
^^^^^^^^
.. autofunction:: orbdot.models.theory.precession_rotational_planet
.. autofunction:: orbdot.models.theory.precession_rotational_planet_k2
.. autofunction:: orbdot.models.theory.precession_rotational_star
.. autofunction:: orbdot.models.theory.precession_rotational_star_k2

Tides
^^^^^
.. autofunction:: orbdot.models.theory.precession_tidal_planet
.. autofunction:: orbdot.models.theory.precession_tidal_planet_k2
.. autofunction:: orbdot.models.theory.precession_tidal_star
.. autofunction:: orbdot.models.theory.precession_tidal_star_k2

Transit Variations
^^^^^^^^^^^^^^^^^^
.. autofunction:: orbdot.models.theory.get_pdot_from_wdot
.. autofunction:: orbdot.models.theory.get_tdot_from_wdot

------------

.. _proper_motion_theory:

Proper Motion
-------------
Proper motion is the movement of a star through space with respect to reference frame of the Solar System, and is partially constrained with measurements of the proper motion on the sky-plane :math:`\mu` and radial velocity :math:`v_r`. The systemic proper motion of an exoplanet host star changes the orientation of the orbit in the sky-plane, resulting in an apparent variation of the argument of pericenter :math:`\omega_p` and the line-of-sight inclination :math:`\i` of the orbit that can raise apparent secular transit variations. The OrbDot methods that are described below estimate the upper limit of these proper motion effects by applying equations derived in :cite:t:`Rafikov2009`.

.. autofunction:: orbdot.models.theory.proper_motion_idot
.. autofunction:: orbdot.models.theory.proper_motion_wdot
.. autofunction:: orbdot.models.theory.proper_motion_pdot
.. autofunction:: orbdot.models.theory.proper_motion_tdot
.. autofunction:: orbdot.models.theory.proper_motion_shklovskii

.. _planet_companion_theory:

Companion Planets
-----------------
Blurb...

Long-Term Radial Velocity Trends
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: orbdot.models.theory.companion_rv_trend_from_mass
.. autofunction:: orbdot.models.theory.companion_mass_from_rv_trend
.. autofunction:: orbdot.models.theory.companion_from_quadratic_rv
.. autofunction:: orbdot.models.theory.companion_doppler_pdot_from_rv_trend
.. autofunction:: orbdot.models.theory.companion_doppler_rv_trend_from_pdot

Companion-Induced Precession
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: orbdot.models.theory.companion_precession
.. autofunction:: orbdot.models.theory.companion_mass_from_precession

.. _binary_star_theory:

Resolved Stellar Binary
-----------------------
Blurb...

.. autofunction:: orbdot.models.theory.resolved_binary_mass_from_rv_trend
.. autofunction:: orbdot.models.theory.resolved_binary_rv_trend_from_mass

