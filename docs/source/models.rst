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

.. automodule:: orbdot.models.ttv_models

------------

.. automodule:: orbdot.models.rv_models

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
   :nosignatures:

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
   :nosignatures:

        orbdot.models.theory.proper_motion_wdot
        orbdot.models.theory.proper_motion_idot
        orbdot.models.theory.proper_motion_pdot
        orbdot.models.theory.proper_motion_tdot
        orbdot.models.theory.proper_motion_shklovskii

.. _planet_companion_theory:

Planetary Companion
-------------------

.. autosummary::
   :nosignatures:

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
   :nosignatures:

        orbdot.models.theory.resolved_binary_mass_from_rv_trend
        orbdot.models.theory.resolved_binary_rv_trend_from_mass
