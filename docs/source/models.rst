.. _models:

**********
The Models
**********
This section describes the models...

.. _model_parameters:

Model Parameters
================
To keep the use of this class clean and concise, there is a main set of "allowed" model parameters. For every model fit, the list of free parameters is compared to this set and the given order is recorded. This means that any number of free parameters can be provided in any order as long as they are part of the physical model(s), as defined by the log-likelihood. It is important for the user to familiarize themselves with the parameter symbols and definitions.

The OrbDot parameter set also includes time derivatives of other orbital elements, such as :math:`\dot{e}` and :math:`\dot{\Omega}`, that are planned for integration in future iterations but are not currently implemented in any OrbDot models.

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
     - --
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
     - --
     - The eccentricity :math:`e` multiplied by the cosine of :math:`\omega_p`.
   * - :math:`e\,\sin{\,\omega_p}`
     - ``esinw``
     - --
     - the eccentricity :math:`e` multiplied by the sine of :math:`\omega_p`.
   * - :math:`\sqrt{e}\,\cos{\,\omega_p}`
     - ``sq_ecosw``
     - --
     - The square root of :math:`e` multiplied by the cosine of :math:`\omega_p`.
   * - :math:`\sqrt{e}\,\sin{\,\omega_p}`
     - ``sq_esinw``
     - --c
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

.. _coordinate_system:

Coordinate System
=================
In this work, we establish the sky plane to lie within
the x-z plane, with the y-axis pointing towards the observer along the line of sight (refer to Figure X).

It is critical to be consistent in the definition of the argument of pericenter :math:`\omega` when simultaneously
fitting transit timing and RV data to an eccentric orbit model. The results can become inaccurate if the pericenter
angles for the planet and star are measured from different axes.

.. image:: _static/coordinate_system.png

The argument of pericenter is determined from the positive x-axis, such that a transit occurs when the true anomaly
:math:`\phi` is equal to:

.. math::
 \phi_{\mathrm I}\,=\,\frac{\pi}{2} - \omega_{\mathrm p}

and an eclipse occurs when:

.. math::
 \phi_{\mathrm II} = \frac{3\pi}{2} - \omega_{\mathrm p}

------------

Transit and Eclipse Timing Models
=================================
Secular variations in the orbits of transiting exoplanets may be detected by measuring deviations
in the observed transit and eclipse timing from what is expected of a circular, unchanging orbit.
These deviations are commonly referred to as Transit Timing Variations (TTVs).

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
For a planet on a decaying circular orbit with a constant orbital period, we expect the mid-times of the transits (:math:`t_{\\mathrm{I}}`) and eclipses (:math:`t_{\\mathrm{II}}`) to be:

.. math:: t_{\\mathrm{I}} = t_0 + PE + \\frac{1}{2}\\,\\frac{dP}{dE}\\,E^2

and,

.. math:: t_{\\mathrm{II}} = t_0 + PE + \\frac{P}{2} + \\frac{1}{2}\\,\\frac{dP}{dE}\\,E^2

where :math:`t_0` is the reference transit time, :math:`P` is the orbital period, :math:`dP/dE` is the rate of change of the period in units of days per epoch, and :math:`E` is the epoch, which represents the number of orbits that have passed since time :math:`t_0`.

If the orbit is eccentric, an offset of :math:`\\frac{P_a\\,e}{\\pi}\\,\\cos{\\,\\omega_p}` is added to the eclipse times.

The :meth:`~orbdot.models.ttv_models.ttv_decay` method implements this model:

.. autofunction:: orbdot.models.ttv_models.ttv_decay

Apsidal Precession
------------------
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

.. automodule:: orbdot.models.tdv_models

------------

.. _theory_module:

The ``theory`` Module
=====================
The :py:mod:`~orbdot.models.theory` module provides several analytical models that can be used to investigate long-term variations in planetary orbits and their causes, which are used in the :class:`~orbdot.analysis.Analyzer`. The key methods of this class are listed below, organized by subjects: :ref:`orbital decay <orbital_decay_theory>`, :ref:`apsidal precession <_apsidal_precession_theory>`, nonresonant :ref:`companion planets <planet_companion_theory>`, bound :ref:`stellar companions <binary_star_theory>`, and systemic :ref:`proper motion <proper_motion_theory>`.

.. _orbital_decay_theory:

Orbital Decay
-------------
Orbital decay refers to a transfer of angular momentum from the planet to the host star that results in a shrinking of the orbital period, eventually leading to planetary engulfment.

Due to the close proximity of HJs to their host stars, significant tidal bulges -- an ellipsoidal distortion -- are raised in both the planet and star. In the case of orbital decay, the planet orbital rate is faster than the star's rotational rate. As a result, the star's tidal bulge lags behind the HJ, creating a net torque that spins up the star at the expense of the planet's orbital angular momentum \citep{levrard_falling_2009, penev_empirical_2018, ma_orbital_2021}. The tidal forces raised by the misaligned tidal bulges are known as `equilibrium tides' and are believed to be the most significant process governing the future evolution of HJ orbits \citep{ma_orbital_2021, barker_tidal_2020}.

.. autosummary::
   :nosignatures:

   orbdot.models.theory.decay_pdot_from_quality_factor
   orbdot.models.theory.decay_quality_factor_from_pdot
   orbdot.models.theory.decay_empirical_quality_factor
   orbdot.models.theory.decay_timescale
   orbdot.models.theory.decay_energy_loss
   orbdot.models.theory.decay_angular_momentum_loss

------------

.. _apsidal_precession_theory:

Apsidal Precession
------------------
Apsidal precession is the gradual increase of the argument of pericentre :math:`\omega` of a planet's orbit over time, meaning the line connecting the pericentre and apocentre of the orbit rotates through :math:`2\pi` in one precession period.

This can result from several factors, including components due to general relativistic effects :cite:p:`pal_periastron_2008,jordan_observability_2008`, perturbations from other planets :cite:p:`heyl_using_2007`, and gravitational moments arising from both the host star's rotation and planetary tidal bulges :cite:p:`greenberg_apsidal_1981`. The following sections describe the equations and OrbDot methods that are relevant to these effects.

.. autosummary::

   orbdot.models.theory.precession_gr
   orbdot.models.theory.precession_rotational_planet
   orbdot.models.theory.precession_rotational_planet_k2
   orbdot.models.theory.precession_rotational_star
   orbdot.models.theory.precession_rotational_star_k2
   orbdot.models.theory.precession_tidal_planet
   orbdot.models.theory.precession_tidal_planet_k2
   orbdot.models.theory.precession_tidal_star
   orbdot.models.theory.precession_tidal_star_k2
   orbdot.models.theory.get_tdot_from_wdot
   orbdot.models.theory.get_pdot_from_wdot

------------

.. _proper_motion_theory:

Proper Motion
-------------
The apparent secular evolution of exoplanet transit signatures that are induced by the systemic proper motion, which is the movement of the star-planet system with respect to reference frame of the Solar System. This motion in 3D space is partially constrained with measurements of the proper motion on the sky-plane :math:`\mu` and radial velocity :math:`v_r`.

.. autosummary::

        orbdot.models.theory.proper_motion_wdot
        orbdot.models.theory.proper_motion_idot
        orbdot.models.theory.proper_motion_pdot
        orbdot.models.theory.proper_motion_tdot
        orbdot.models.theory.proper_motion_shklovskii

.. _planet_companion_theory:

Planetary Companion
-------------------

.. autosummary::

        orbdot.models.theory.companion_doppler_pdot_from_rv_trend
        orbdot.models.theory.companion_doppler_rv_trend_from_pdot
        orbdot.models.theory.companion_from_quadratic_rv
        orbdot.models.theory.companion_mass_from_rv_trend
        orbdot.models.theory.companion_rv_trend_from_mass
        orbdot.models.theory.companion_precession
        orbdot.models.theory.companion_mass_from_precession

.. _binary_star_theory:

Resolved Stellar Companion
--------------------------

.. autosummary::

        orbdot.models.theory.resolved_binary_mass_from_rv_trend
        orbdot.models.theory.resolved_binary_rv_trend_from_mass
