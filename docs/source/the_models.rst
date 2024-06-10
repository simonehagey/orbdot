.. _models:

******
Models
******

The Model Parameters
====================

.. list-table:: **Orbital Elements**
   :header-rows: 2
   :widths: 30 25 25 100

   * -
     -
     -
     -
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

.. [*] Not currently implemented.

|

.. list-table:: **Coupled Parameters**
   :header-rows: 2
   :widths: 30 25 20 100

   * -
     -
     -
     -
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

|

.. list-table:: **Time-Dependent Parameters**
   :header-rows: 2
   :widths: 25 15 25 100

   * -
     -
     -
     -
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
     - A constant change of the orbital eccentricity.
   * - :math:`\frac{di}{dE}`
     - ``idE``
     - deg :math:`E^{-1}`
     - A constant change of the line-of-sight inclination [*]_.
   * - :math:`\frac{d \Omega}{dE}`
     - ``OdE``
     - rad :math:`E^{-1}`
     - A constant change of the long. of the ascending node.

.. [*] Not currently implemented.

|

.. list-table:: **Radial Velocity Parameters**
   :header-rows: 2
   :widths: 20 15 25 70

   * -
     -
     -
     -
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
     - m :math:`{\mathrm s}^{-2}` :math:`{\mathrm day}^{-1}`
     - A second order radial velocity trend.


Coordinate System
=================

It is critical to be consistent in the definition of the argument of pericenter :math:`\omega` when simultaneously
fitting transit timing and RV data to an eccentric orbit model. The results can become inaccurate if the pericenter
angles for the planet and star are measured from different axes. In this work, we establish the sky plane to lie within
the x-z plane, with the y-axis pointing towards the observer along the line of sight (refer to Figure X).

.. image:: coordinate_system.png
   :width: 40pt

The argument of pericenter is determined from the positive x-axis, such that a transit occurs when the true anomaly
:math:`\phi` is equal to:

.. math::
    \phi_{\mathrm I}\,=\,\frac{\pi}{2} - \omega_{\mathrm p}

and an eclipse occurs when:

.. math::
    \phi_{\mathrm II} = \frac{3\pi}{2} - \omega_{\mathrm p}


Transit and Eclipse Timing Models
=================================
Secular variations in the orbits of transiting exoplanets may be detected by measuring deviations in the observed
transit and eclipse timing from what is expected of a circular, unchanging orbit. These deviations are commonly
referred to as Transit Timing Variations (TTVs). We examine two such deviations:

1. A quadratic TTV trend, which suggests orbital decay.
2. A sinusoidal TTV trend, indicative of apsidal precession.

Constant-Period
---------------

For a planet on a **circular** orbit with a constant orbital period, we expect a linear increase in the center times of both transits :math:`t_{\mathrm I}` and eclipses :math:`t_{\mathrm II}`:

.. math::

     &t_{\mathrm I} = t_0 + PE \\
     &t_{\mathrm II} = t_0 + PE + \frac{P}{2}

where :math:`t_0` is the reference transit time, :math:`P` is the orbital period, and  :math:`E` is the epoch, which represents the number of orbits that have passed since time :math:`t_0`.

If the orbit is **eccentric**, we add must add an offset to the eclipse times to account for the variable speed of the planet:

 .. math::

     &t_{\mathrm I} = t_0 + PE \\
     &t_{\mathrm II} = t_0 + PE + \frac{P}{2} + \frac{P_a\,e}{\pi}\,\cos{\,\omega_p}

.. autofunction:: orbdot.models.ttv_models.ttv_constant

Orbital Decay
-------------

The first case of secular evolution that we consider is a planet on a circular orbit experiencing **orbital decay**. This gradual shrinking of the orbital period results in a progressive shift in the measured transit and eclipse mid-times. The equations for this scenario are:

.. math::

     &t_{\mathrm I} = t_0 + PE + \frac{1}{2}\,\frac{dP}{dE}\,E^2 \\
     &t_{\mathrm II} = t_0 + PE + \frac{P}{2} + \frac{1}{2}\,\frac{dP}{dE}\,E^2

where :math:`dP/dE` is the rate of change of the period in units of days per epoch, and the value of which will be negative in the case of orbital decay.

.. autofunction:: orbdot.models.ttv_models.ttv_decay

Apsidal Precession
------------------

In the case of **apsidal precession**, which requires an eccentric orbit, we assume that the argument of pericenter of the planet's orbit evolves at a constant rate, denoted as :math:`\frac{d\omega}{dE}`.

 So, at any given :math:`E` we have: :math:`\, \omega_{p}\,(E)\, = \omega_0 + \frac{d\omega}{dE}\,E`

This secular evolution induces long-term sinusoidal trends in the transit and eclipse timing data. The equations for this scenario are:

.. math::
  &t_{\mathrm I} = t_0 + P_s\,E - \frac{e P_a}{\pi}\,\cos\,{\omega_{p}} \\
  &t_{\mathrm II} = t_0 + P_s E + \frac{P_a}{2} + \frac{eP_a}{\pi}\,\cos{\,\omega_{p}}

where :math:`P_a` is the 'anomalistic' orbital period and :math:`P_s` is the sidereal period, related to the anomalistic period by:

.. math::
  P_s = P_a\left(1-\frac{d\omega/{dE}}{2\pi}\right)

.. autofunction:: orbdot.models.ttv_models.ttv_precession

------------

Radial Velocity Models
======================

Many host stars of transiting Hot Jupiters exhibit noticeable periodic radial velocity (RV) variations as they
'wobble' around the center of mass of the star-planet system. The amplitude of this effect can be expressed as:

.. math::

    K=\left(\frac{2 \pi G}{P}\right)^{1/3} \frac{M_p \sin i}{\left(M_{\star}+M_p\right)^{2/3}} \frac{1}{\left(1-e^2\right)^{1/2}}

where :math:`M_p` is the planet mass, :math:`M_{\star}` is the stellar mass, :math:`i` is the inclination of the orbit with respect to the line-of-sight, and :math:`G` is the universal gravitational constant [1]_.

At any given time, the periodic signal depends on the planet's position in its orbit, defined by the true anomaly :math:`\phi`, and the systemic velocity along the line-of-sight, denoted as :math:`\gamma`. When RV observations from multiple instruments are combined, as is done in this study, it is standard practice to fit $\gamma$ individually for each source to account for instrumental variations. Additionally, we consider first and second-order acceleration terms, :math:`\dot{\gamma}` and :math:`\ddot{\gamma}` respectively, to accommodate potential perturbations from an outer, non-resonant companion planet with an orbital period longer than the observational baseline. The total observed RV signal is:

.. math::

    v_r = -K[\cos{(\phi\left(t\right)+\omega_\star)}+e\cos{\omega_\star}] + \gamma_j + \dot{\gamma} \left(t-t_{_{\rm ref}}\right)

where :math:`\omega_\star` is the argument of pericentre of the star's orbit, defined as :math:`\omega_\star = \omega_{p} + \pi`. The sign in Equation X is needed in this coordinate system so that blue-shifted spectral lines correspond to a negative radial velocity. Additionally, to avoid underestimating the uncertainties of the parameters in the RV model fit, we incorporate an instrument-dependent 'jitter' parameter instrument into the log likelihood, :math:`\sigma_j`. These terms are added in quadrature with the individual measurement errors to accommodate additional noise that may arise from stellar activity and instrument systematics. SHOW EQUATION

.. [1] Heyl and Gladman (2007).

**Orbital Decay**
To perform a joint fit, the analytical models described in Sections :ref:`sec:timing_models`
and :ref:`sec:rv_model` must share as many free parameters as possible.

**Apsidal Precession**
To perform a joint fit, the analytical models described in Sections :ref:`sec:timing_models`
and :ref:`sec:rv_model` must share as many free parameters as possible.

The orbital period $P$
------------------------
First, one must be careful when defining the orbital period $P$. For a precessing orbit, there is a distinction
between the sidereal, or observed, period $P_s$ and the anomalistic, or true, period $P_a$. These are related by:

.. math::

    P_s = P_a\left(1-\frac{d\omega/{dE}}{2\pi}\right)

To simultaneously fit the apsidal precession and radial velocity models, the radial velocity model must use calculate
$P_a$. If the orbit is not precessing, $d\omega/dE=0$ and $P_a = P_s = P$.

The reference time :math:`t{_{\mathrm ref}}`
--------------------------------------------
The reference time :math:`t{_{\mathrm ref}}` in the radial velocity model can be directly substituted with the reference transit
time, :math:`t_0`. At any time :math:`t`, the epoch number :math:`E` can be calculated by:

.. math::
    E = \frac{{(t - t_0)}}{{P_a}}

It is notable here that the variable :math:`E` is not restricted to being an integer number, and it only is for transits
because of its definition.

Argument of Pericenter :math:`\omega_\star`
-------------------------------------------
At any position in the planet's orbit, the argument of pericenter of the star $\omega_\star$ can be expressed in
terms of :math:`\omega_p`:

.. math::
    \omega_\star = \omega_{p} + \pi

For a precessing orbit, :math:`\omega_p` is a function of time. In Equation \ref{eq:omega_p}, this time-dependence expressed
in terms of the epoch:

.. math::
    \omega_{p}\left(E\right) = \omega_0 + \frac{d\omega}{dE}E.

In this way, :math:`\omega_\star` is written as a function of $\omega_0$ and :math:`d\omega/dE`. Apsidal precession is
continuous, and thus the variable $E$ is not restricted to being an integer number and this substitution can be
made directly in the RV equations, regardless of the position in the orbit. If the orbit is not precessing, $d\omega/dE=0$.

True anomaly
-------------
The most involved part is figuring out the true anomaly :math:`\phi` at any given time :math:`t`. This involves solving Kepler's
equation numerically:

.. math::
    \textrm{M} = \textrm{E} - e \sin \textrm{E}

Where :math:`\textrm{E}` is the eccentric anomaly (not to be confused with the epoch number :math:`E`) and :math:`\textrm{M}` is the
mean anomaly. The true anomaly :math:`\phi` is connected to the eccentric anomaly through:

.. math::
    \tan \left(\frac{\phi}{2}\right) = \sqrt{\frac{1+e}{1-e}} \tan \left(\frac{\textrm{E}}{2}\right)

and the mean anomaly is related to time through:

.. math::
    \textrm{M} = \eta \left(t - t_p\right)

where :math:`\eta = 2 \pi/P_a` is the mean motion of the planet's orbit and :math:`t_p` is the time marking the passage of pericentre.

Time of Pericenter Passage
---------------------------
But first we must calculate the time elapsed since pericenter passage. This is straightforward, as we know that at
time :math:`t_0` the true anomaly is :math:`\phi_0 = \frac{\pi}{2} - \omega_{p}`. So we can calculate the true anomaly at
transit $\phi_0$, and then work backwards through Kepler's equation to calculate the mean anomaly at transit
:math:`\textrm{M}_0`, and then solve for the time of pericentre passage :math:`t_p`. After all that, we can calculate the true
anomaly at any time, :math:`\phi\left(t\right)`.

Now, by re-parameterizing the RV model, we have successfully established shared parameters (:math:`t_0`, :math:`P`, :math:`e`,
:math:`\omega_0`, and :math:`d\omega/dE`) for a more robust analysis. To integrate the models in Sections
:ref:`sec:timing_models` and :ref:`sec:rv_model`, the log-likelihoods are summed and this joint log-likelihood is
given to the nested sampling algorithm. This effectively captures the interplay between the timing variations and
radial velocity signals in the models, better constraining the long-term behaviour of the data.


Transit Duration Models
=======================

