"""
Theory
======
This module defines analytical models for investigating the source of long-term variations of
exoplanet orbits. These include equations for assessing the effects of tidal energy dissipation,
apsidal precession, systemic proper motion, and more.
"""

import numpy as np
import scipy.special as sci

# define constants
c = 2.99792458e8    # m/s
G = 6.6743e-11      # m^3 / kg / s^2

# define unit conversions
AU = 1.495978707e11   # 1 AU = 1.496 x 10^11 m
parsec = 3.0857e16    # 1 pc = 3.0857 x 10^16 m
R_earth = 6.371e6     # earth radius = 6,371.000 km
M_earth = 5.9722e24   # earth mass = 5.9722 x 10^24 kg
R_jup = 6.9911e7      # jupiter radius = 69911 km
M_jup = 1.89813e27    # jupiter mass = 1,898.13 x 10^24 kg
R_sun = 6.957e8       # solar radius = 695,700 km
M_sun = 1.9885e30     # solar mass = 1,988,500 x 10^30 kg


def companion_doppler_pdot_from_rv_trend(P, dvdt):
    """Calculate the apparent variation of an orbital period due to a line-of-sight acceleration.

    This method returns the time derivative of the observed orbital period due to the Doppler
    effect induced by an acceleration along the line-of-sight (i.e. a linear RV trend). It uses
    equation (6) from Bouma et al. (2020) [1]_, which is derived in convenient units.

    Parameters
    ----------
    P : float
        The orbital period in days.
    dvdt : float
        The linear radial velocity trend in m/s/day.

    Returns
    -------
    float
        The apparent time derivative of the orbital period in milliseconds per year.

    Notes
    -----
    When radial velocity data of a transiting exoplanet's host star show a linear trend (i.e.,
    acceleration), and assuming this trend is a real effect, the Doppler effect can cause the
    observed period between transits to vary. This change in the transit period is produced by
    the line-of-sight acceleration of the system. Equation (6) from Bouma et al. (2020) [1]_
    expresses the expected period derivative in very convenient units:

    .. math::
        \\dot{P}_{\\mathrm{RV}} = 105.3 \\mathrm{\\,ms\\,yr^{-1}}\\left(\\frac{P}{\\mathrm{
        day}}\\right)\\left(\\frac{\\dot{\\gamma}}{\\mathrm{m\\,s^{-1}\\,day^{-1}}\\right)

    which is derived from,

    .. math::
        \\dot{P}_{\\mathrm{RV}} = \\frac{\\dot{v}_r P}{c}

    where :math:`\\dot{P}_{\\mathrm{RV}}` is the time derivative of the observed orbital period,
    :math:`\\dot{v}_r` is the linear radial velocity trend, :math:`P` is the orbital period of
    the transiting planet, and :math:`c` is the speed of light in a vacuum.

    References
    ----------
    .. [1] :cite:t:`Bouma2020`. https://doi.org/10.3847/2041-8213/ab8563.

    """
    return 105.3 * P * dvdt


def companion_doppler_rv_trend_from_pdot(P, dPdt):
    """Calculate the line-of-sight acceleration given the apparent variation of an orbital period.

    This method returns the linear radial velocity trend (acceleration) that corresponds to a given
    time derivative of the observed orbital period due to the Doppler effect. It uses equation (6)
    from Bouma et al. (2020) [1]_, which is derived in convenient units.

    Parameters
    ----------
    P : float
        The orbital period in days.
    dPdt : float
        The apparent time derivative of the orbital period in milliseconds per year (ms/yr).

    Returns
    -------
    float
        The linear radial velocity trend (acceleration) in m/s/day.

    Notes
    -----
    For a transiting exoplanet, if a variation in the observed period between transits is measured,
    it may be due to the Doppler effect induced by a line-of-sight acceleration. Assuming this
    trend is a real effect, Equation (6) from Bouma et al. (2020) [1]_ can be used to determine
    the corresponding radial velocity trend from the period derivative:

    .. math::
        \\dot{\\gamma} = \\frac{\\dot{P}_{\\mathrm{RV}}}{105.3 \\mathrm{\\,ms\\,yr^{-1}}}
        \\left(\\frac{\\mathrm{day}}{P}\\right)

    which is derived from:

    .. math:: \\dot{\\gamma} = \\frac{\\dot{P}_{\\mathrm{RV}} c}{P}

    where math:`\\dot{\\gamma}` is the linear radial velocity trend, :math:`\\dot{P}_{\\mathrm{
    RV}}` is the time derivative of the observed orbital period, :math:`P` is the orbital period of
    the transiting planet, and :math:`c` is the speed of light in a vacuum.

    References
    ----------
    .. [1] :cite:t:`Bouma2020`. https://doi.org/10.3847/2041-8213/ab8563.

    """
    return dPdt / (105.3 * P)


def companion_mass_from_rv_trend(tau, dvdt, M_s):
    """Calculate the minimum mass of an outer companion given the slope of a linear RV trend.

    This method computes the minimum mass for an unseen outer companion planet given a linear
    trend in radial velocity data (i.e., an acceleration). While this approach to determining
    the minimum companion mass was originally described in Feng et al. (2015) [1]_ (see eq. 1),
    this method implements the version given by equation (8) in Bouma et al. (2020) [2]_.

    Parameters
    ----------
    tau : float
        The time span of the observations in years.
    dvdt : float
        The linear radial velocity trend in m/s/day.
    M_s : float
        The mass of the host star in solar masses.

    Returns
    -------
    float
        The companion planet mass in Earth masses.

    Notes
    -----
    Given a measured acceleration in radial velocity observations and the time span over which
    there are RV measurements, this method calculates the minimum mass of an outer companion that
    could induce such an acceleration. While this approach to determining the minimum companion
    mass was originally described in Feng et al. (2015) [1]_ (see eq. 1), this method implements
    the version given by equation (8) in Bouma et al. (2020) [2]_:

    .. math::
        M_{c} \\approx 5.99 M_{\\mathrm{Jup}}\\left(\\frac{\\tau}{\\mathrm{yr}}\\right)^{
        4/3}\\left|\\frac{\\dot{\\gamma}}{\\mathrm{m\\,s^{-1}\\,day^{-1}}}\\right|\\left(
        \\frac{M_{\\star}}{M_{\\odot}}\\right)^{2 / 3}

    where :math:`\\tau` is the time span of the observations in years, :math:`\\dot{\\gamma}` is
    the linear radial velocity trend, and :math:`M_\\star` is the mass of the host star.

    As a means of constraining the absolute minimum possible companion mass, the above
    equation was derived with the assumption that the companion's orbit has an eccentricity of 0.5,
    an argument of pericenter that is exactly 90 degrees, and a period that is 1.25 times the time
    span of the observations. In this case, the observed linear trend in the radial velocity
    model is effectively a section of the sawtooth-like curve of the companion planet (see
    Figure 1 of [1]_), for which the semi-amplitude can be approximated as half of
    the baseline multiplied by the acceleration (:math:`0.5 * \\tau * \\frac{dv}{dt}`).

    References
    ----------
    .. [1] :cite:t:`Feng2015`. https://doi.org/10.1088/0004-637X/800/1/22
    .. [2] :cite:t:`Bouma2020`. https://doi.org/10.3847/2041-8213/ab8563

    """
    # calculate the companion planet mass (Jupiter masses)
    M_c = 5.99 * tau ** (4/3) * np.abs(dvdt) * M_s ** (2/3)

    # return the mass in Earth masses
    return M_c * M_jup / M_earth


def companion_rv_trend_from_mass(tau, M_c, M_s):
    """Calculate the slope of a radial velocity trend given a lower limit on the companion mass.

    This method computes the expected linear radial velocity trend (i.e., acceleration) induced
    by an unseen outer companion planet given its minimum mass. This is calculated using equation
    (8) from Bouma et al. (2020) [1]_, which was adapted from Feng et al. (2015) [2]_.

    Parameters
    ----------
    tau : float
        The time span of the observations in years.
    M_c : float
        The companion planet mass in earth masses.
    M_s : float
        The mass of the host star in solar masses.

    Returns
    -------
    float
        The linear radial velocity trend in m/s/day.

    Notes
    -----
    Given a companion planet's mass, the time span over which radial velocity measurements have
    been taken, and the mass of the host star, this method calculates the expected linear trend
    in the radial velocity data using equation (8) from Bouma et al. (2020) [1]_:

    .. math::
        \\dot{\\gamma} \\approx \\frac{M_c}{5.99 M_{\\mathrm{Jup}}} \\left(\\frac{\\mathrm{yr}}{
        \\tau}\\right)^{4/3} \\left(\\frac{M_{\\odot}}{M_{\\star}}\\right)^{2 / 3} \\mathrm{m\\,
        s^{-1}\\,day^{-1}}

    where :math:`\\tau` is the time span of the observations in years, :math:`M_c` is the mass of
    the companion planet, and :math:`M_\\star` is the mass of the host star.

    The above equation assumes that the companion's orbit has an eccentricity of 0.5, an argument
    of pericenter that is exactly 90 degrees, and a period that is :math:`1.25\\times` times the
    timespan of the observations. In this scenario, the observed linear trend in the radial
    velocity model is a segment of the sawtooth-like curve of the companion planet (see Figure 1
    of [1]_), for which the semi-amplitude can be approximated as half of the baseline multiplied
    by the acceleration, ie. :math:`0.5 \\tau \\dot{\\gamma}`.

    References
    ----------
    .. [1] :cite:t:`Feng2015`. https://doi.org/10.1088/0004-637X/800/1/22
    .. [2] :cite:t:`Bouma2020`. https://doi.org/10.3847/2041-8213/ab8563

    """
    # convert earth masses to jupiter masses
    M_c *= M_earth / M_jup

    # calculate and return the RV slope in m/s/day
    return M_c / (5.99 * tau ** (4 / 3) * M_s ** (2 / 3))


def companion_from_quadratic_rv(P_min, t_pivot, dvdt, ddvdt, M_s):
    """Constrain properties of the orbit of an outer companion given a quadratic RV trend.

    Given first and second order acceleration terms that describe a quadratic radial velocity
    trend, this method follows the formulation from Kipping et al. (2011) [1]_ (see Equations 1,
    3, and 4) to constrain properties of a possible outer companion. It requires a minimum
    possible orbital period of the companion, which should be informed by the timespan of the
    radial velocity measurements (see the :meth:`~orbdot.analysis.Analyzer.rv_trend_quadratic`
    method from the :class:`~orbdot.analysis.Analyzer` class).

    Parameters
    ----------
    P_min : float
        The minimum possible orbital period of the companion in days.
    t_pivot : float
        The 'pivot' point in days. This is often fixed as the mean time of the RV observations,
        but in the case of OrbDot joint fitting, it is the reference mid-time of the transiting
        planet.
    dvdt : float
        Linear radial velocity trend in m/s/day.
    ddvdt : float
        Quadratic radial velocity trend in m/s/day^2.
    M_s : float
        The mass of the host star in solar masses.

    Returns
    -------
    tuple
        The minimum possible orbital period in days, the associated lower limit on the RV
        semi-amplitude in m/s, the companion mass in Earth masses, and the time when the outer
        companion RV signal is at a minimum.

    References
    ----------
    .. [1] :cite:t:`Kipping2011`. https://doi.org/10.1088/0004-6256/142/3/95

    """
    # unit conversions
    M_s *= M_sun    # solar masses to kg

    # calculate the lower limit of the RV semi-amplitude from the quadratic RV term in m/s
    K_min = np.abs(ddvdt) * P_min ** 2 / (4 * np.pi ** 2)

    # calculate the lower limit of the mass from K_min
    M_min = K_min * M_s ** (2/3) * ((P_min * 86400) / (2 * np.pi * G)) ** (1/3)

    # solve for the time when the outer companion RV signal is at a minimum
    tau_c = (-dvdt + ddvdt * t_pivot) / ddvdt - P_min/4

    return P_min, K_min, M_min/M_earth, tau_c


def companion_precession(P, M2, P2, M_s):
    """
    Calculates the rate of apsidal precession driven by a nonresonant planetary companion.

    This method returns the rate of apsidal precession expected to be induced by a nonresonant
    companion planet in the system, calculated using equation (8) from Heyl and Gladman (2007) [
    1]_. This method is valid for a companion planet that is either interior or exterior to the
    observed (i.e., transiting) planet, but it assumes that the companion mass is much less than
    that of the host star.

    Parameters
    ----------
    P : float
        The orbital period of the transiting planet in days.
    M2 : float
        The mass of the outer planet in Earth masses.
    P2 : float
        The orbital period of the outer planet in days.
    M_s : float
        The mass of the host star in solar masses.

    Returns
    -------
    float
        The precession rate in radians per orbit.

    Notes
    -----
    This method enables the exploration of the possibility that a nonresonant companion in the
    system is driving apsidal precession of the observed planet. It assumes that the companion
    object is on a circular, coplanar, and nonresonant orbit, and that its mass is far less than
    the host star. In the formulation of Heyl and Gladman (2007) [1]_, the induced precession
    rate in arcseconds per year is:

    .. math::
        \\delta \\varpi = \\frac{m_2}{M_\\star} \\frac{\\alpha}{(\\alpha+1)(\\alpha-1)^2}\\left[
        \\left(\\alpha^2+1\\right) E\\left(\\frac{2 \\alpha^{1 / 2}}{\\alpha+1}\\right)
        -(\\alpha-1)^2 K\\left(\\frac{2 \\alpha^{1/2}}{\\alpha+1}\\right)\\right]

    where :math:`m_2` is the mass of the perturbing planet, :math:`\\alpha = a_1/a_2` is the
    semi-major axis ratio, and :math:`\\mathcal{K}` and :math:`\\mathcal{E}` are the complete
    elliptic integrals of the first and second kind, respectively.

    .. important::
        The ``scipy`` implementation of the complete elliptic integral of the first kind
        :math:`\\mathcal{K}` expects an integrand of the form :math:`(1 - m \\sin^2 t)^{-1/2}`,
        but we expect to have :math:`(1 - m^2 \\sin^2 t)^{-1/2}`, which can be seen via the
        expansion of :math:`\\mathcal{K}`. To keep the code implementation consistent with the
        equation, the argument for :math:`\\mathcal{K}` must be squared when using the ``scipy``
        version, such that :math:`\\mathcal{K}( 4 \\alpha/( \\alpha+1)^2)`.

    References
    ----------
    .. [1] :cite:t:`Heyl2007`. https://doi.org/10.1111/j.1365-2966.2007.11697.x

    """
    # derive parameters
    a1 = get_semi_major_axis_from_period(P, M_s)
    a2 = get_semi_major_axis_from_period(P2, M_s)

    # unit conversions
    M_s *= M_sun    # solar masses to kg
    M2 *= M_earth   # earth masses to kg

    # define the ratio of the semi major axes
    alpha = a1 / a2

    # calculate the rate of apsidal precession
    t1 = alpha / ((alpha + 1) * (alpha - 1) ** 2)

    # exterior orbit
    if alpha < 1.0:
        t2 = (alpha ** 2 + 1) * sci.ellipe(2 * alpha ** (1/2) / (alpha + 1))
        t3 = (alpha - 1) ** 2 * sci.ellipk(2 * alpha ** (1/2) / (alpha + 1))

    # interior orbit
    else:
        t2 = (alpha ** 2 + 1) * sci.ellipe(4 * alpha / (alpha + 1) ** 2)
        t3 = (alpha - 1) ** 2 * sci.ellipk(4 * alpha / (alpha + 1) ** 2)

    dwdE = (M2 / M_s) * t1 * (t2 - t3)

    # return the apsidal precession rate in radians per orbit (rad/E)
    return dwdE


def companion_mass_from_precession(P, P2, dwdE, M_s):
    """Calculate the mass of a nonresonant companion given a precession rate and orbital period.

    Given the orbital period of a nonresonant companion planet in the system, this method returns
    the necessary mass for it to be responsible for a measured apsidal precession rate of the
    studied orbit. It uses equation (8) from Heyl and Gladman (2007) [1]_, and is appropriate for
    a companion planet that is either interior or exterior to the observed (i.e., transiting)
    planet. However, it assumes that the companion mass is much less than that of the host star.

    Parameters
    ----------
    P : float
        The orbital period of the transiting planet in days.
    P2 : float
        The orbital period of the outer planet days.
    dwdE : float
        The rate of apsidal precession in radians per orbit.
    M_s : float
        The mass of the host star in solar masses.

    Returns
    -------
    float
        The mass of the outer companion in Earth masses.

    Notes
    -----
    This method enables the exploration of the possibility that a nonresonant companion in the
    system is driving apsidal precession of the observed planet. It assumes that the companion
    object is on a circular, coplanar, and nonresonant orbit, and that its mass is far less than
    the host star. In the formulation of Heyl and Gladman (2007) [1]_, the mass of the companion
    planet can be calculated from the induced precession rate in arcseconds per year using:

    .. math::
        m_2 = \\delta \\varpi \\frac{M_\\star (\\alpha+1)(\\alpha-1)^2}{\\alpha \\left[
        \\left(\\alpha^2+1\\right) E\\left(\\frac{2 \\alpha^{1 / 2}}{\\alpha+1}\\right)
        - (\\alpha-1)^2 K\\left(\\frac{2 \\alpha^{1/2}}{\\alpha+1}\\right)\\right]}

    where :math:`\\delta \\varpi` is the observed precession rate, :math:`\\alpha = a_1/a_2` is
    the semimajor axis ratio, and :math:`\\mathcal{K}` and :math:`\\mathcal{E}` are the complete
    elliptic integrals of the first and second kind, respectively.

    .. important::
        The ``scipy`` implementation of the complete elliptic integral of the first kind
        :math:`\\mathcal{K}` expects an integrand of the form :math:`(1 - m \\sin^2 t)^{-1/2}`,
        but we expect to have :math:`(1 - m^2 \\sin^2 t)^{-1/2}`, which can be seen via the
        expansion of :math:`\\mathcal{K}`. To keep the code implementation consistent with the
        equation, the argument for :math:`\\mathcal{K}` must be squared when using the ``scipy``
        version, such that :math:`\\mathcal{K}(4 \\alpha/(\\alpha+1)^2)`.

    References
    ----------
    .. [1] :cite:t:`Heyl2007`. https://doi.org/10.1111/j.1365-2966.2007.11697.x

    """
    # derive parameters
    a1 = get_semi_major_axis_from_period(P, M_s)
    a2 = get_semi_major_axis_from_period(P2, M_s)

    # unit conversions
    M_s *= M_sun    # solar masses to kg

    # define the ratio of the semi major axes
    alpha = a1 / a2

    # calculate the mass of the companion planet
    t1 = alpha / ((alpha + 1) * (alpha - 1) ** 2)
    t2 = (alpha ** 2 + 1) * sci.ellipe(4 * alpha / (alpha + 1) ** 2)
    t3 = (alpha - 1) ** 2 * sci.ellipk(4 * alpha / (alpha + 1) ** 2)

    M2 = dwdE * M_s / (t1 * (t2 - t3))

    # return the companion mass in Earth masses
    return M2 / M_earth


def decay_quality_factor_from_pdot(P, dPdE, M_s, M_p, R_s):
    """Calculate the modified stellar quality factor given the orbital decay rate.

    This method returns a predicted modified stellar quality factor for any given orbital decay
    rate, assuming that the star-planet system is coplanar and that equilibrium tides dominate
    the dynamical evolution.

    Parameters
    ----------
    P : float
        Orbital period in days.
    dPdE : float
        Orbital decay rate in days per orbit.
    M_s : float
        Host star mass in solar masses.
    M_p : float
        Planet mass in earth masses.
    R_s : float
        Host star radius in solar radii.

    Returns
    -------
    float
        The modified stellar quality tidal factor.

    Notes
    -----
    Assuming that equilibrium tides dominate the evolution, the rate of orbital decay depends on
    the efficiency of tidal energy dissipation within the star, which is parameterized by the
    tidal quality factor. In the constant-phase lag model for equilibrium tides, the modified
    stellar tidal quality factor is [1]_:

    .. math::
        Q_\\star^{'} = -\\frac{27\\pi}{2\\dot{P}_{\\mathrm{decay}}}\\left(\\frac{M_p}{
        M_\\star}\\right)\\left(\\frac{R_\\star}{a}\\right)^5

    where :math:`M_\\star` is the host star mass, :math:`M_p` is the planet mass,
    :math:`R_\\star` is the host star radius, :math:`a` is the orbital semi-major axis,
    and :math:`\\dot{P}_{\\mathrm{decay}}` is the orbital decay rate.

    The modified tidal quality factor is defined as:

    .. math:: Q_\\star^{'} = \\frac{3}{2} \\frac{Q_\\star}{k_{2,\\star}}

    where :math:`Q_\\star` is the tidal quality factor and :math:`k_{2,\\star}` is the second-order
    potential Love number [2]_. The theoretical upper limit is :math:`k_{2,\\star}` for a uniform
    density sphere, in which case :math:`Q_\\star^{'} = Q_\\star`.

    References
    ----------
    .. [1] :cite:t:`Goldreich1966`. https://doi.org/10.1016/0019-1035(66)90051-0
    .. [2] :cite:t:`Ogilvie2007`. https://doi.org/10.1017/9781108304061

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)

    # unit conversions
    M_p *= M_earth    # earth masses to kg
    M_s *= M_sun      # solar masses to kg
    R_s *= R_sun      # solar radii to m
    dPdt = dPdE / P   # days/E to days/day

    # compute and return the stellar tidal quality factor
    Q_star = -27 * np.pi / (2 * dPdt) * (M_p / M_s) * (R_s / a) ** 5

    return Q_star


def decay_pdot_from_quality_factor(P, M_s, M_p, R_s, Q_star):
    """Calculate the orbital decay rate given the host star's modified tidal quality factor.

    This method returns a predicted orbital decay rate for any given "modified" stellar tidal
    quality factor, assuming that the star-planet system is coplanar and that equilibrium
    tides dominate the dynamical evolution.

    Parameters
    ----------
    P : float
        Orbital period in days.
    M_s : float
        Host star mass in solar masses.
    M_p : float
        Planet mass in Earth masses.
    R_s : float
        Host star radius in solar radii.
    Q_star : float
        The modified stellar tidal quality factor.

    Returns
    -------
    float
        Rate of change of the planet's orbital period in days per orbit.

    Notes
    -----
    Assuming that equilibrium tides dominate the evolution, the rate of orbital decay depends on
    the efficiency of tidal energy dissipation within the star, which is parameterized by the
    tidal quality factor. In the constant-phase lag model for equilibrium tides, the orbital
    decay rate is [1]_:

    .. math::
        \\dot{P}_{\\mathrm{decay}} = -\\frac{27\\pi}{2 Q_\\star^{'}}\\left(\\frac{M_p}{
        M_\\star}\\right)\\left(\\frac{R_\\star}{a}\\right)^5

    where :math:`M_\\star` is the host star mass, :math:`M_p` is the planet mass,
    :math:`R_\\star` is the host star radius, :math:`a` is the orbital semi-major axis,
    and :math:`Q_\\star^{'} ` is the star's "modified" tidal quality factor.

    The modified tidal quality factor is defined as:

    .. math:: Q_\\star^{'} = \\frac{3}{2} \\frac{Q_\\star}{k_{2,\\star}}

    where :math:`Q_\\star` is the tidal quality factor and :math:`k_{2,\\star}` is the second-order
    potential Love number [2]_. The theoretical upper limit is :math:`k_{2,\\star}` for a uniform
    density sphere, in which case :math:`Q_\\star^{'} = Q_\\star`.

    References
    ----------
    .. [1] :cite:t:`Goldreich1966`. https://doi.org/10.1016/0019-1035(66)90051-0
    .. [2] :cite:t:`Ogilvie2007`. https://doi.org/10.1017/9781108304061

    """

    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)

    # unit conversions
    M_p *= M_earth  # earth masses to kg
    M_s *= M_sun    # solar masses to kg
    R_s *= R_sun    # solar radii to m

    # calculate the predicted orbital decay rate (days/day)
    pdot = -27 * np.pi / (2 * Q_star) * (M_p / M_s) * (R_s / a) ** 5

    return pdot * P   # days/E


def decay_empirical_quality_factor(P_orb, P_rot_s):
    """Calculate the modified stellar tidal quality factor from an empirical law.

    This method calculates the tidal forcing period of a star-planet system and uses it to
    estimate the host star's modified tidal quality factor with an empirical law derived
    by Penev et al. (2018) [1]_.

    Parameters
    ----------
    P_orb : float
        Orbital period in days.
    P_rot_s : float
        Period of the host star's rotation in days.

    Returns
    -------
    tuple
        The star's modified tidal quality factor and the tidal forcing period of the system in days.

    Notes
    -----
    The "modified" tidal quality factor is a typical parameterization of the efficiency of tidal
    energy dissipation within a body. It is defined as:

    .. math:: Q_\\star^{'} = \\frac{3}{2} \\frac{Q_\\star}{k_{2,\\star}}

    where :math:`Q_\\star` is the tidal quality factor and :math:`k_{2,\\star}` is the second-order
    potential Love number [2]_. The theoretical upper limit is :math:`k_{2,\\star}` for a uniform
    density sphere, in which case :math:`Q_\\star^{'} = Q_\\star`.

    This method estimates a value of :math:`Q_\\star^{'} ` with an empirical law derived by
    Penev et al. (2018) [1]_, given as:

    .. math::
        Q^{'}_\\star = \\max\\left(\\frac{10^6}{(P_{\\mathrm{tide}}/\\mathrm{day})^{3.1}},
        10^5\\right)

    where :math:`P_{\\mathrm{tide}}` is the tidal forcing period of the star-planet system,
    defined as:

    .. math::
        P_{\\mathrm{tide}} = \\frac{1}{2\\left(P_{\\mathrm{orb}}^{-1} - P_{\\mathrm{rot}}^{
        -1}\\right)}

    where :math:`P_{\\mathrm{orb}}` is the planet's orbital period and :math:`P_{\\mathrm{rot}}`
    is the rotational period of the host star.

    References
    ----------
    .. [1] :cite:t:`Penev2018`. https://doi.org/10.3847/1538-3881/aaaf71
    .. [2] :cite:t:`Ogilvie2007`. https://doi.org/10.1017/9781108304061

    """
    # determine the tidal forcing period
    t1 = 1 / P_orb - 1 / P_rot_s
    P_tide = 1 / (2 * t1)

    # use the empirical law for Q'
    Q_star = 10 ** 6.0 / P_tide ** 3.1

    # return the tidal quality factor Q'
    return Q_star, P_tide


def decay_timescale(P, dPdE):
    """Calculate the remaining lifetime of a planet on a decaying orbit.

    Parameters
    ----------
    P : float
        Orbital period in days.
    dPdE : float
        Orbital decay rate in days per orbit.

    Returns
    -------
    float
        The remaining lifetime in millions of years (Myr).

    Notes
    -----
    This method calculates the timescale over which a planetary orbit is shrinking using the
    equation:

    .. math:: \\tau = \\frac{P_0}{|\\dot{P}_{\\mathrm{decay}}|}

    where :math:`\\dot{P}_{\\mathrm{decay}}` is the orbital decay rate and :math:`P_0` is the
    initial orbital period.

    """

    # convert days/orbit to days/yr
    dPdt = dPdE / P * 365.25

    # calculate the remaining lifetime
    tau = P / np.abs(dPdt)

    # return the remaining lifetime in Myr
    return tau / 1e6


def decay_energy_loss(P, dPdE, M_s, M_p):
    """Calculate the rate of orbital energy loss due to tidal forces causing orbital decay.

    This function uses Equation (9) from Yee et al. (2020) [1]_ to compute the rate at which the
    orbital energy of a planetary orbit decreases as the orbit decays due to tidal interactions
    between the planet and its host star.

    Parameters
    ----------
    P : float
        Orbital period in days.
    dPdE : float
        Orbital decay rate in days per orbit.
    M_s : float
        Mass of the host star in solar masses.
    M_p : float
        Mass of the planet in Earth masses.

    Returns
    -------
    float
        The rate of orbital energy loss in watts.

    Notes
    -----
    As a planet's orbit decays, the orbital energy and angular momentum will decrease over time.
    Equation (9) of Yee et al. (2020) [1]_ defines the rate of orbital energy loss as:

    .. math::
        \\frac{dE}{dt} = \\frac{(2\\pi)^{2/3} M_{\\mathrm{p}}}{3} \\left(\\frac{G M_{\\star}}{
        P}\\right)^{2/3} \\frac{1}{P} \\frac{dP}{dt}

    where :math:`G` is the gravitational constant, :math:`M_{\\star}` is the host star mass,
    :math:`M_{\\mathrm{p}}` is the planet mass, :math:`P` is the orbital period,
    and :math:`\\frac{dP}{dt}` is the orbital decay rate.

    References
    ----------
    .. [1] :cite:t:`Yee2020`. https://doi.org/10.3847/2041-8213/ab5c16

    """
    # unit conversions
    M_p *= M_earth     # earth masses to kg
    M_s *= M_sun       # solar masses to kg
    dPdt = dPdE / P    # days/E to days/day
    P *= 86400         # days to seconds

    # calculate the time derivative of the orbital energy
    t1 = (2 * np.pi) ** (2 / 3) * M_p / 3
    t2 = (G * M_s / P) ** (2 / 3)
    t3 = (1 / P) * dPdt
    dEdt = t1 * t2 * t3

    # return the rate of energy loss in Watts
    return dEdt


def decay_angular_momentum_loss(P, dPdE, M_s, M_p):
    """Calculate the rate of angular momentum loss due to tidal forces causing orbital decay.

    This function uses Equation (10) from Yee et al. (2020) [1]_ to compute the rate at which
    the angular momentum of a planetary orbit decreases as the orbit decays due to tidal
    interactions between the planet and its host star.

    Parameters
    ----------
    P : float
        Orbital period in days.
    dPdE : float
        Orbital decay rate in days per orbit.
    M_s : float
        Host star mass in solar masses.
    M_p : float
        Planet mass in Earth masses.

    Returns
    -------
    float
        The rate of orbital angular momentum loss in kg m^2 / s^2.

    Notes
    -----
    As a planet's orbit decays, the orbital energy and angular momentum will decrease over time.
    Equation (10) of Yee et al. (2020) [1]_ defines the rate of orbital angular momentum loss as:

    .. math::
        \\frac{dL}{dt} = \\frac{M_{\\mathrm{p}}}{3(2\\pi)^{1/3}} \\left(\\frac{G M_{\\star}}{
        P}\\right)^{2/3} \\frac{dP}{dt}

    where :math:`G` is the gravitational constant, :math:`M_{\\star}` is the host star mass,
    :math:`M_{\\mathrm{p}}` is the planet mass, :math:`P` is the orbital period,
    and :math:`\\frac{dP}{dt}` is the orbital decay rate.

    References
    ----------
    .. [1] :cite:t:`Yee2020`. https://doi.org/10.3847/2041-8213/ab5c16

    """
    # unit conversions
    M_p *= M_earth     # earth masses to kg
    M_s *= M_sun       # solar masses to kg
    dPdt = dPdE / P    # days/E to days/day
    P *= 86400         # days to seconds

    # calculate the time derivative of the orbit angular momentum
    t1 = M_p / (3 * (2 * np.pi) ** (1 / 3))
    t2 = (G * M_s / P) ** (2 / 3)
    dLdt = t1 * t2 * dPdt

    # return the rate of angular momentum loss in kg m^2 / s^2
    return dLdt


def get_tdot_from_wdot(P, e, w, i, T, dwdE, M_s, R_s):
    """Calculate the expected transit duration variation (TDV) due to apsidal precession.

    This method determines the expected TDV signal for apsidal precession in general, independent
    of the physical mechanism. It uses Equation (9) from Rafikov (2009) [1]_, assuming there is
    no change in the line-of-sight inclination over time.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    w : float
        The argument of pericenter in radians.
    i : float
        The line-of-sight inclination in degrees.
    T : float
        The transit duration in minutes.
    dwdE : float
        The precession rate in radians per orbit.
    M_s : float
        The host star mass in solar masses.
    R_s : float
        The host star radius in solar radii.

    Returns
    -------
    float
        The time derivative of the transit duration in milliseconds per year (ms/yr).

    Notes
    -----
    For transiting exoplanets, apsidal precession causes the duration of the transits :math:`T`
    to vary over time. Assuming that the line-of-sight inclination :math:`i` of the orbit is
    constant, the expected transit duration variation (TDV) is given by [1]_:

    .. math::
        \\dot{T} = -\\frac{T}{1 + e\\sin{\\omega}} \\left[e\\dot{\\omega}\\cos{\\omega} -
        g \\dot{\\omega}\\cos{i}\\frac{e\\cos{\\omega}}{1+e\\sin{\\omega}}\\right]

    where,

    .. math:: g = \\frac{a}{R_\\star}\\frac{b}{1 - b^2}

    and where :math:`e` is the orbit eccentricity, :math:`\\omega` is the argument of pericenter,
    :math:`\\dot{\\omega}` is the apsidal precession rate, :math:`b` is the impact parameter,
    :math:`a` is the semi-major axis, and :math:`R_\\star` is the radius of the host star.

    References
    ----------
    .. [1] :cite:t:`Rafikov2009`. https://doi.org/10.1088/0004-637X/700/2/965

    """
    # unit conversions
    i *= np.pi / 180     # degrees to radians
    R_s *= R_sun         # solar radii to m
    T *= 60 * 1000       # minutes to ms
    wdot = dwdE * (1 / P) * 365.25   # rad/E to rad/yr

    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    b = (a / R_s) * np.cos(i) * (1 - e ** 2) / (1 + e * np.sin(w))
    g = (a / R_s) * b / (1 - b ** 2)

    # calculate the transit duration variation
    t1 = - T / (1 + e * np.sin(w))
    t2 = e * wdot * np.cos(w)
    t3 = - g * (wdot * np.cos(i) * e * np.cos(w)/(1 + e * np.sin(w)))

    # return the TDV rate in ms/yr
    return t1 * (t2 + t3)


def get_pdot_from_wdot(P, e, w, dwdE):
    """Calculate the apparent time derivative of the orbital period due to apsidal precession.

    This method determines the apparent time derivative of the orbital period due to apsidal
    precession in general, independent of the physical mechanism. It uses Equation (17) from
    Rafikov (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    w : float
        The argument of pericenter in radians.
    dwdE : float
        The precession rate in radians per orbit.

    Returns
    -------
    float
        The time derivative of the orbital period in milliseconds per year (ms/yr).

    Notes
    -----
    For transiting exoplanets, apsidal precession causes the sidereal period, i.e., the length of
    time between transits, to vary over time. Because apsidal precession periods are typically
    far greater than observational baselines available at this time, this effect yields an
    *apparent* constant time derivative :math:`\\dot{P}`. Equation (17) of Rafikov (2009) [1]_
    expresses this as:

    .. math::
        \\dot{P} = \\frac{4\\pi\\left(\\dot{\\omega}\\right)^2 e\\cos\\omega}{
        \\eta^2}\\frac{\\left(1-e^2\\right)^{3/2}}{(1+e\\sin\\omega)^3}

    where :math:`e` is the orbit eccentricity, :math:`\\omega` is the argument of pericenter,
    :math:`\\dot{\\omega}` is the apsidal precession rate, and :math:`\\eta = 2\\pi/P` is the
    orbital mean motion.

    References
    ----------
    .. [1] :cite:t:`Rafikov2009`. https://doi.org/10.1088/0004-637X/700/2/965

    """
    # derive parameters
    nu = 2 * np.pi / P

    # unit conversions
    dwdt = dwdE * (1 / P) * 365.25   # rad/E to rad/yr

    # calculate the apparent period derivative
    t1 = e * np.cos(w) * (1 - e ** 2) ** (3 / 2) / (1 + e * np.sin(w)) ** 3
    t2 = 4 * np.pi * (dwdt / nu) ** 2
    pdot = t1 * t2   # days^2/yr^2

    # return the apparent period derivative in ms/yr
    return pdot * (1 / 365.25) * 8.64e+7


def precession_gr(P, e, M_s):
    """Calculate the rate of apsidal precession predicted by general relativity.

    This method returns the expected apsidal precession rate of the planet's orbit due to general
    relativistic (GR) effects, using Equation (12) from Ragozzine and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.

    Returns
    -------
    float
        The precession rate in radians per orbit.

    Notes
    -----
    This method applies the lowest order of the relativistic contribution to apsidal precession,
    given in Equation (12) of Ragozzine and Wolf (2009) [1]_ as:

    .. math:: \\dot{\\omega}_{\\mathrm{GR}} = \\frac{3 \\eta G M_{\\star}}{a c^2(1 - e^2)}

    where :math:`G` is the gravitational constant, :math:`c` is the speed of light in a vacuum,
    :math:`M_{\\star}` is the host star mass, :math:`a` is the planet's semi-major axis,
    :math:`e` is the eccentricity, and :math:`\\eta = 2\\pi/P` is the mean motion.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # unit conversions
    M_s *= M_sun    # solar masses to kg
    nu *= 1 / 86400      # 1/days to 1/s

    # calculate the precession rate
    wdot = 3 * G * M_s * nu / (a * c ** 2 * (1 - e ** 2))

    # convert from rad/s to rad/E
    wdot *= 86400 * P

    # return the apsidal precession rate in radians per orbit
    return wdot


def precession_rotational_star(P, e, M_s, R_s, k2_s, P_rot_s):
    """Calculate the rate of apsidal precession driven by stellar rotation.

    This method returns the expected apsidal precession rate of a planet's orbit due to the
    rotational bulge of the host star, calculated using Equations (10) and (11) from Ragozzine
    and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.
    R_s : float
        The host star radius in solar radii.
    k2_s : float
        The host star's Love number.
    P_rot_s : float
        The period of the host star's rotation in days.

    Returns
    -------
    float
        The precession rate in radians per orbit.

    Notes
    -----
    The rotation of a fluid body raises an oblate distortion, i.e., flattening, that perturbs its
    gravitational potential. For close-in planets and their host stars, this effect is expected
    to drive apsidal precession of the planetary orbit. The rate of precession induced by the
    star's rotational bulge is given in Equation (10) of Ragozzine and Wolf (2009) [1]_ as:

    .. math::
        \\dot{\\omega}_{\\mathrm{rot,\\star}} = \\frac{\\eta {R_\\star}^5 k_{2,\\star} {\\dot{
        \\theta}_\\star}^2}{2 a^2 G M_\\star} g_2(e)

    where :math:`\\dot{\\theta}_\\star` is the rotation speed of the star, :math:`k_{2,\\star}` is
    its second-order potential Love number, :math:`G` is the gravitational constant, :math:`a` is
    the semi-major axis of the orbit, :math:`\\eta = 2\\pi/P` is the mean motion,
    and :math:`M_\\star` and :math:`R_\\star` are the stellar mass and radius, respectively.

    The parameter :math:`g_2(e)` represents an expansion in eccentricity :math:`e`, defined in
    Equation (11) of Ragozzine and Wolf (2009) [1]_ as:

    .. math:: g_2(e) = (1-e^2)^{-2}

    The Love number represents how centrally condensed the star is and is a fixed property of
    the body [2]_. The lower the :math:`k_{2,\\star}`, the more centrally condensed the star's
    interior structure, which in turn leads to a slower precession rate.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778
    .. [2] :cite:t:`Lissauer2019`. https://doi.org/10.1017/9781108304061

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # calculate the eccentricity expansion and rotational velocity
    g = (1 - e ** 2) ** (-2)
    v_s = 2 * np.pi / P_rot_s  # rad/day

    # unit conversions
    M_s *= M_sun        # solar masses to kg
    R_s *= R_sun        # solar radii to m
    nu *= 1 / 86400     # 1/days to 1/s
    v_s *= 1 / 86400    # rad/day to rad/s

    # calculate precession rate
    wdot = (k2_s / 2) * (R_s / a) ** 5 * (v_s ** 2 * a ** 3 / (G * M_s)) * g * nu
    wdot *= (86400 * P)  # convert rad/s to rad/E

    return wdot  # rad/E


def precession_rotational_planet(P, e, M_s, M_p, R_p, k2_p, P_rot_p):
    """
    Calculates the rate of apsidal precession driven by planetary rotation.

    This method returns the expected apsidal precession rate of a planet's orbit due to its
    rotational bulge, calculated using Equations (10) and (11) from Ragozzine and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.
    M_p : float
        The planet mass in earth masses.
    R_p : float
        The planet radius in earth radii.
    k2_p : float
        The planet's Love number.
    P_rot_p : float
        The period of the planet's rotation in days.

    Returns
    -------
    float
        The precession rate in radians per orbit.

    Notes
    -----
    The rotation of a fluid body raises an oblate distortion (i.e., flattening) that perturbs its
    gravitational potential. For close-in planets and their host stars, this effect is predicted
    to drive apsidal precession of the planetary orbit. The rate of precession induced by the
    planet's rotational bulge is given in Equation (10) of Ragozzine and Wolf (2009) [1]_ as:

    .. math::
        \\dot{\\omega}_{\\mathrm{rot,p}} = \\frac{\\eta {R_p}^5 k_{2,p} {\\dot{\\theta}_p}^2}{2
        a^2 G M_p} g_2(e)

    where :math:`\\dot{\\theta}_p` is the rotation speed of the planet, :math:`k_{2,p}` is its
    second-order potential Love number, :math:`G` is the gravitational constant, :math:`\\eta =
    2\\pi/P` is the orbital mean motion, :math:`a` is the semi-major axis, and :math:`M_p` and
    :math:`R_p` are the planet mass and radius, respectively.

    The parameter :math:`g_2(e)` represents an expansion in eccentricity :math:`e`, defined in
    Equation (11) of Ragozzine and Wolf (2009) [1]_ as:

    .. math:: g_2(e) = (1-e^2)^{-2}

    The Love number represents how centrally condensed the planet is and is a fixed property of
    the body [2]_. The lower the :math:`k_{2,p}`, the more centrally condensed the planetary
    interior structure, which in turn leads to a slower precession rate.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778
    .. [2] :cite:t:`Lissauer2019`. https://doi.org/10.1017/9781108304061

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # calculate the eccentricity expansion and rotational velocity
    g = (1 - e ** 2) ** (-2)
    v_p = 2 * np.pi / P_rot_p  # rad/day

    # unit conversions
    M_p *= M_earth  # earth masses to kg
    R_p *= R_earth  # earth radii to m
    nu *= 1 / 86400  # 1/days to 1/s
    v_p *= 1 / 86400  # rad/day to rad/s

    # calculate precession rate in rad/E
    wdot = (k2_p / 2) * (R_p / a) ** 5 * (v_p ** 2 * a ** 3 / (G * M_p)) * g * nu  # rad/s

    wdot *= (86400 * P)  # convert rad/s to rad/E

    return wdot  # rad/E


def precession_rotational_star_k2(P, e, M_s, R_s, P_rot_s, dwdE):
    """Calculate the Love number given a rate of apsidal precession driven by stellar rotation.

    Given an apsidal precession rate for a planetary orbit, this method returns the second-order
    potential Love number of the host star :math:`k_{2,\\star}` under the assumption that the
    precession is driven entirely by the rotational bulge of the host star. See the
    :meth:`~orbdot.models.theory.precession_rotational_star` docstring for a more detailed
    description of the relevant equations, which are adopted from Ragozzine and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.
    R_s : float
        The host star radius in solar radii.
    P_rot_s : float
        The period of the host star's rotation in days.
    dwdE : float
        The precession rate in radians per orbit.

    Returns
    -------
    float
        The host star's potential Love number of the second order.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # calculate the eccentricity expansion and rotational velocity
    g = (1 - e ** 2) ** (-2)
    v_s = 2 * np.pi / P_rot_s  # rad/day

    # unit conversions
    M_s *= M_sun        # solar masses to kg
    R_s *= R_sun        # solar radii to m
    nu *= 1 / 86400     # 1/days to 1/s
    v_s *= 1 / 86400    # rad/day to rad/s
    dwdt = dwdE / (86400 * P)  # rad/E to rad/s

    # compute and return the Love number
    t1 = (1 / 2) * (R_s / a) ** 5 * (v_s ** 2 * a ** 3 / (G * M_s)) * g * nu
    k2_s = dwdt * t1 ** (-1)

    return k2_s


def precession_rotational_planet_k2(P, e, M_s, M_p, R_p, P_rot_p, dwdE):
    """Calculate the Love number given a rate of apsidal precession driven by planetary rotation.

    Given an apsidal precession rate for a planetary orbit, this method returns the second-order
    potential Love number of the planet :math:`k_{2,p}` under the assumption that the precession
    is driven entirely by the rotational bulge of the planet. See the
    :meth:`~orbdot.models.theory.precession_rotational_planet` docstring for a more detailed
    description of the relevant equations, which are adopted from Ragozzine and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.
    M_p : float
        The planet mass in earth masses.
    R_p : float
        The planet radius in earth radii.
    P_rot_p : float
        The period of the planet's rotation in days.
    dwdE : float
        The precession rate in radians per orbit.

    Returns
    -------
    float
        The planet's potential Love number of the second order.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # calculate the eccentricity expansion and rotational velocity
    g = (1 - e ** 2) ** (-2)
    v_p = 2 * np.pi / P_rot_p  # rad/day

    # unit conversions
    M_p *= M_earth      # earth masses to kg
    R_p *= R_earth      # earth radii to m
    nu *= 1 / 86400     # 1/days to 1/s
    v_p *= 1 / 86400    # rad/day to rad/s
    dwdt = dwdE / (86400 * P)  # rad/E to rad/s

    # compute and return the Love number
    t1 = (1 / 2) * (R_p / a) ** 5 * (v_p ** 2 * a ** 3 / (G * M_p)) * g * nu
    k2_p = dwdt * t1 ** (-1)

    return k2_p


def precession_tidal_star(P, e, M_s, M_p, R_s, k2_s):
    """Calculate the rate of apsidal precession driven by the star's tidal bulge.

    This method returns the expected apsidal precession rate of the planet's orbit due to the
    star's tidal bulge, calculated using Equations (6) and (7) from Ragozzine and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.
    M_p : float
        The planet mass in earth masses.
    R_s : float
        The host star radius in solar radii.
    k2_s : float
        The host star's Love number.

    Returns
    -------
    float
        The precession rate in radians per orbit.

    Notes
    -----
    A massive, close-in planet raises a "tidal bulge" on its host star, which is an ellipsoidal
    distortion caused by the varying gravitational force experienced across the extent of
    the body. This physical distortion changes the gravitational potential, driving apsidal
    precession of the orbit. The rate of precession induced by the star's tidal bulge is given in
    Equation (6) of Ragozzine and Wolf (2009) [1]_ as:

    .. math::
        \\dot{\\omega}_{\\mathrm{tide,\\star}} = \\frac{15}{2}k_{2,\\star} \\eta \\left(\\frac{
        R_\\star}{a}\\right)^5\\left(\\frac{M_p}{M_\\star}\\right)f_2(e)

    where :math:`k_{2,\\star}` is the star's second-order potential Love number, :math:`\\eta =
    2\\pi/P` is the orbital mean motion, :math:`a` is the semi-major axis, :math:`M_p` is the
    planet mass, and :math:`M_\\star` and :math:`R_\\star` are the stellar mass and radius.

    The parameter :math:`f_2(e)` represents an expansion in eccentricity :math:`e`, defined in
    Equation (7) of Ragozzine and Wolf (2009) [1]_ as:

    .. math:: f_2(e)= (1-e^2)^{-5} \\left(1 + \\frac{3}{2}e^2 + \\frac{1}{8}e^4 \\right).

    The Love number represents how centrally condensed the star is and is a fixed property of
    the body [2]_. The lower the :math:`k_{2,\\star}`, the more centrally condensed the stellar
    interior structure, which in turn leads to a slower precession rate.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # calculate the eccentricity expansion
    f = (1 - e ** 2) ** (-5) * (1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4)

    # unit conversions
    M_p *= M_earth      # earth masses to kg
    M_s *= M_sun        # solar masses to kg
    R_s *= R_sun        # solar radii to m
    nu *= 1 / 86400     # 1/days to 1/s

    # calculate precession rate in rad/E
    wdot = (15 / 2) * k2_s * nu * f * (R_s / a) ** 5 * (M_p / M_s)  # rad/s

    wdot *= (86400 * P)  # convert rad/s to rad/E

    return wdot  # rad/E


def precession_tidal_planet(P, e, M_s, M_p, R_p, k2_p):
    """Calculate the rate of apsidal precession driven by the planet's tidal bulge.

    This method returns the expected apsidal precession rate of the planet's orbit due to its tidal
    bulge, calculated using Equations (6) and (7) from Ragozzine and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.
    M_p : float
        The planet mass in earth masses.
    R_p : float
        The planet radius in earth radii.
    k2_p : float
        The planet's Love number.

    Returns
    -------
    float
        The precession rate in radians per orbit.

    Notes
    -----
    The host star of a massive, close-in planet raises a "tidal bulge" on the planet, which is an
    ellipsoidal distortion caused by the varying gravitational force experienced across
    the extent of the body. This physical distortion changes the gravitational potential,
    driving apsidal precession of the orbit. The rate of precession induced by the planet's tidal
    bulge is given in Equation (6) of Ragozzine and Wolf (2009) [1]_ as:

    .. math::
        \\dot{\\omega}_{\\mathrm{tide,p}} = \\frac{15}{2}k_{2,p} \\eta \\left(\\frac{R_p}{
        a}\\right)^5\\left(\\frac{M_\\star}{M_p}\\right)f_2(e)

    where :math:`k_{2,p}` is the planet's second-order potential Love number, :math:`\\eta =
    2\\pi/P` is the orbital mean motion, :math:`a` is the semi-major axis, :math:`M_\\star` is
    the host star mass, and :math:`M_p` and :math:`R_p` are the planet mass and radius.

    The parameter :math:`f_2(e)` represents an expansion in eccentricity :math:`e`, defined in
    Equation (7) of Ragozzine and Wolf (2009) [1]_ as:

    .. math:: f_2(e)= (1-e^2)^{-5} \\left(1 + \\frac{3}{2}e^2 + \\frac{1}{8}e^4 \\right).

    The Love number represents how centrally condensed the planet is and is a fixed property of
    the body [2]_. The lower the :math:`k_{2,p}`, the more centrally condensed the planetary
    interior structure, which in turn leads to a slower precession rate.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # calculate the eccentricity expansion
    f = (1 - e ** 2) ** (-5) * (1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4)

    # unit conversions
    M_p *= M_earth      # earth masses to kg
    M_s *= M_sun        # solar masses to kg
    R_p *= R_earth        # earth radii to m
    nu *= 1 / 86400     # 1/days to 1/s

    # calculate precession rate in rad/E
    wdot = (15 / 2) * k2_p * nu * f * (R_p / a) ** 5 * (M_s / M_p)  # rad/s

    wdot *= (86400 * P)  # convert rad/s to rad/E

    return wdot  # rad/E


def precession_tidal_star_k2(P, e, M_s, M_p, R_s, dwdE):
    """Calculate a host star's Love number given a rate of apsidal precession driven by tides.

    Given an apsidal precession rate for a planetary orbit, this method returns the second-order
    potential Love number of the host star :math:`k_{2,\\star}` under the assumption that the
    precession is driven entirely by the tidal bulge of the star. See the
    :meth:`~orbdot.models.theory.precession_tidal_star` docstring for a more detailed description
    of the relevant equations, which are adopted from Ragozzine and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.
    M_p : float
        The planet mass in earth masses.
    R_s : float
        The host star radius in solar radii.
    dwdE : float
        The precession rate in radians per orbit.

    Returns
    -------
    float
        The host star's potential Love number of the second order.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # calculate the eccentricity expansion
    f = (1 - e ** 2) ** (-5) * (1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4)

    # unit conversions
    M_p *= M_earth      # earth masses to kg
    M_s *= M_sun        # solar masses to kg
    R_s *= R_sun        # solar radii to m
    nu *= 1 / 86400     # 1/days to 1/s
    dwdt = dwdE / (86400 * P)  # rad/E to rad/s

    # compute and return the Love number
    t1 = (15 / 2) * nu * f * (R_s / a) ** 5 * (M_p / M_s)
    k2_s = dwdt * t1 ** (-1)

    return k2_s


def precession_tidal_planet_k2(P, e, M_s, M_p, R_p, dwdE):
    """Calculate a planet's Love number given a rate of apsidal precession driven by tides.

    Given an apsidal precession rate for a planetary orbit, this method returns the second-order
    potential Love number of the planet :math:`k_{2,p}` under the assumption that the precession
    is driven entirely by the tidal bulge of the planet. See the
    :meth:`~orbdot.models.theory.precession_tidal_planet` docstring for a more detailed
    description of the relevant equations, which are adopted from Ragozzine and Wolf (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    M_s : float
        The host star mass in solar masses.
    M_p : float
        The planet mass in earth masses.
    R_p : float
        The planet radius in earth radii.
    dwdE : float
        The precession rate in radians per orbit.

    Returns
    -------
    float
        The planet's potential Love number of the second order.

    References
    ----------
    .. [1] :cite:t:`Ragozzine2009`. https://doi.org/10.1088/0004-637X/698/2/1778

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)
    nu = 2 * np.pi / P

    # calculate the eccentricity expansion
    f = (1 - e ** 2) ** (-5) * (1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4)

    # unit conversions
    M_p *= M_earth  # earth masses to kg
    M_s *= M_sun  # solar masses to kg
    R_p *= R_earth  # earth radii to m
    nu *= 1 / 86400  # 1/days to 1/s
    dwdt = dwdE / (86400 * P)  # convert rad/E to rad/s

    # compute and return the Love number
    t1 = (15 / 2) * nu * f * (R_p / a) ** 5 * (M_s / M_p)  # rad/s
    k2_p = dwdt * t1 ** (-1)

    return k2_p


def proper_motion_wdot(mu, i, beta):
    """Calculate the rate of the apparent apsidal precession due to systemic proper motion.

    This method returns the rate of the apparent apsidal precession induced by systemic proper
    motion using Equation (4) from Rafikov (2009) [1]_.

    Parameters
    ----------
    mu : float
        The proper motion of the system in mas/yr (milliarcseconds per year).
    i : float
        The line-of-sight inclination of the orbit in degrees.
    beta : float
        The angle in radians between the proper motion vector and the angular momentum vector
        projected onto the sky-plane.

    Returns
    -------
    float
        The apparent apsidal precession rate in radians per year.

    Notes
    -----
    The proper motion of a star-planet system alters its appearance on the sky-plane,
    which causes an apparent precession of the argument of pericenter :math:`\\dot{
    \\omega}_\\mu`. Rafikov (2009) [1]_ derives a simple expression for this effect in terms of
    the magnitude of the proper motion vector, :math:`\\mu = |\\vec{\\mu}|`:

    .. math:: \\dot{\\omega}_\\mu = -\\mu\\frac{\\sin{\\beta}}{\\sin{i}}

    where :math:`\\beta` is defined (in the coordinate system of [1]_) as the angle in the
    sky-plane between the proper motion vector :math:`\\vec{\\mu}` and the projection of the orbital
    angular momentum vector :math:`\\vec{L}`. The magnitude of the apparent precession rate
    :math:`|\\dot{\\omega}_\\mu|` is maximized when :math:`\\beta = 90^\\circ, 270^\\circ`.

    References
    ----------
    .. [1] :cite:t:`Rafikov2009`. https://doi.org/10.1088/0004-637X/700/2/965

    """
    # unit conversions
    beta *= np.pi / 180         # degrees to radians
    i *= np.pi / 180            # degrees to radians
    mu *= 1 / (1000 * 206265)   # mas/yr to rad/yr

    # return the apparent precession rate in rad/yr
    return - mu * np.sin(beta) / np.sin(i)


def proper_motion_idot(mu, beta):
    """Calculate the apparent rate of change of the inclination due to systemic proper motion.

    This method returns the rate of the variation of the line-of-sight inclination due to
    systemic proper motion using Equation (3) from Rafikov (2009) [1]_.

    Parameters
    ----------
    mu : float
        The proper motion of the system in mas/yr (milliarcseconds per year).
    beta : float
        The angle in radians between the proper motion vector and the angular momentum vector
        projected onto the sky-plane.

    Returns
    -------
    float
        The apparent time derivative of the line-of-sight inclination in radians per year.

    Notes
    -----
    The proper motion of a star-planet system alters its appearance on the sky-plane,
    which causes an evolution of the line-of-sight inclination :math:`\\dot{i}_\\mu`. Rafikov
    (2009) [1]_ derives a simple expression for this effect in terms of the magnitude of the
    proper motion vector, :math:`\\mu = |\\vec{\\mu}|`:

    .. math:: \\dot{i}_\\mu = -\\mu\\cos{\\beta}

    where :math:`\\beta` is defined (in the coordinate system of [1]_) as the angle in the
    sky-plane between the proper motion vector :math:`\\vec{\\mu}` and the projection of the
    orbital angular momentum vector :math:`\\vec{L}`. Thus, the magnitude of the inclination
    variation :math:`|\\dot{i}_\\mu|` is maximized when :math:`\\beta = 0^\\circ, 180^\\circ`.

    References
    ----------
    .. [1] :cite:t:`Rafikov2009`. https://doi.org/10.1088/0004-637X/700/2/965

    """
    # unit conversions
    beta *= np.pi / 180         # degrees to radians
    mu *= 1 / (1000 * 206265)   # mas/yr to rad/yr

    # return the time derivative of the inclination in rad/yr
    return - mu * np.cos(beta)


def proper_motion_pdot(P, e, w, mu):
    """Calculate the apparent time derivative of the orbital period due to systemic proper motion.

    This method returns the apparent time derivative of the orbital period that is expected as a
    result of the systemic proper motion-induced (apparent) apsidal precession. It uses Equation
    (15) from Rafikov (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    w : float
        The argument of pericenter in radians.
    mu : float
        The proper motion of the system in milliarcseconds per year (mas/yr).

    Returns
    -------
    float
        The apparent time derivative of the orbital period in milliseconds per year (ms/yr).

    Notes
    -----
    The systemic proper motion of a transiting star-planet system can have a significant
    influence on the observed period between transits due to the *apparent* apsidal precession
    :math:`\\dot{\\omega}_{\\mu}`. The following expression from Rafikov (2009) assumes a constant
    change in the observed period, denoted :math:`\\dot{P}_{\\omega,\\mu}`:

    .. math::
        \\dot{P}_{\\omega,\\mu} = -\\frac{2\\pi}{\\eta^2} \\frac{(1-e^2)^{3/2}}{(1+e\\sin{
        \\omega})^2} \\times \\left[\\ddot{\\omega}_{\\mu}-2(\\dot{\\omega}_{\\mu})^2\\frac{
        e\\cos{\\omega}}{1+e\\sin{\\omega}}\\right]

    where :math:`e` is the orbit eccentricity, :math:`\\omega` is the argument of pericenter,
    :math:`\\eta` is the orbital mean motion, :math:`\\dot{\\omega}_{\\mu}` is the apparent
    precession rate, and :math:`\\ddot{\\omega}_{\\mu}` is the time derivative of the latter.

    Rafikov [1]_ suggests that for transiting systems the following approximation is valid:

    .. math:: \\ddot{\\omega}_{\\mu} \\sim \\mu^2 \\sim (\\dot{\\omega}_{\\mu})^2

    which effectively maximizes the first equation.

    References
    ----------
    .. [1] :cite:t:`Rafikov2009`. https://doi.org/10.1088/0004-637X/700/2/965

    """
    # unit conversions
    mu *= 1 / (1000 * 206265)  # mas/yr to rad/yr

    # derive parameters
    nu = 2 * np.pi / P
    dwdt = mu
    ddwdt = mu ** 2

    # calculate apparent period derivative
    t1 = - 2 * np.pi * (1 / nu) ** 2
    t2 = (1 - e ** 2) ** (3/2) / (1 + e * np.sin(w)) ** 2
    t3 = ddwdt - 2 * dwdt ** 2 * (e * np.cos(w)) / (1 + e * np.sin(w))
    pdot_pm = (t1 * t2 * t3)  # days^2/yr^2

    # return the apparent period derivative in ms/yr
    return pdot_pm * (1 / 365.25) * 8.64e+7


def proper_motion_tdot(P, e, w, i, T, wdot_pm, idot_pm, M_s, R_s):
    """Calculate the time derivative of the transit duration due to systemic proper motion.

    This method returns the time derivative of the transit duration due to systemic proper motion
    using Equation (9) from Rafikov (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    e : float
        The eccentricity of the orbit.
    w : float
        The argument of pericenter in radians.
    i : float
        The line-of-sight inclination of the orbit in degrees.
    T : float
        The transit duration in minutes.
    wdot_pm : float
        The apparent apsidal precession rate in radians per year.
    idot_pm : float
        The apparent time derivative of the line-of-sight inclination in radians per year.
    M_s : float
        The host star mass in solar masses.
    R_s : float
        The host star radius in solar radii.

    Returns
    -------
    float
        A constant time derivative of the transit duration in milliseconds per year (ms/yr).

    Notes
    -----
    The systemic proper motion of a transiting star-planet system can have a significant
    influence on the observed transit duration, :math:`T`, due to a combination of the proper
    motion-induced variation of the line-of-sight inclination :math:`\\dot{i}_\\mu`, as well as
    the *apparent* apsidal precession :math:`\\dot{\\omega}_{\\mu}`. The transit duration
    variation (TDV) :math:`\\dot{T}_{tr,\\mu}` is governed by the followingn expression [1]_:

    .. math::
        \\dot{T}_{\\mu} = - \\frac{T}{1 + e\\sin{\\omega}} \\times \\left[e\\dot{\\omega}_{
        \\mu}\\cos{\\omega} - g \\left(\\dot{i}_{\\mu}\\sin{i} + \\dot{\\omega}_{\\mu}\\cos{
        i}\\frac{e\\cos{\\omega}}{1+e\\sin{\\omega}}\\right)\\right]

    where,

    .. math:: g = \\frac{a}{R_\\star}\\frac{b}{1 - b^2}

    and where :math:`e` is the orbit eccentricity, :math:`\\omega` is the argument of pericenter,
    :math:`b` is the impact parameter, :math:`a` is the semi-major axis, and :math:`R_\\star` is
    the radius of the host star.

    References
    ----------
    .. [1] :cite:t:`Rafikov2009`. https://doi.org/10.1088/0004-637X/700/2/965

    """
    # derive parameters
    a = get_semi_major_axis_from_period(P, M_s)

    # unit conversions
    i *= np.pi / 180      # degrees to radians
    R_s *= R_sun          # solar radii to m
    T *= 60 * 1000        # minutes to ms

    b = (a / R_s) * np.cos(i) * (1 - e ** 2) / (1 + e * np.sin(w))
    g = (a / R_s) * b / (1 - b ** 2)

    # calculate the transit duration variation
    t1 = - T / (1 + e * np.sin(w))
    t2 = e * wdot_pm * np.cos(w)
    t3 = - g * (idot_pm * np.sin(i) + wdot_pm * np.cos(i) * e * np.cos(w)/(1 + e * np.sin(w)))

    # return the transit duration variation in ms/yr
    return t1 * (t2 + t3)


def proper_motion_shklovskii(P, mu, D):
    """Calculate the apparent time derivative of observed period due to the Shklovskii effect.

    This method returns the apparent rate of change of the period between transits due to the
    Shklovskii effect, using Equation (21) from Rafikov (2009) [1]_.

    Parameters
    ----------
    P : float
        The orbital period in days.
    mu : float
        The proper motion of the system in mas/yr (milliarcseconds per year).
    D : float
        The distance to the system in parsecs.

    Returns
    -------
    float
        The apparent time derivative of the orbital period in milliseconds per year (ms/yr).

    Notes
    -----
    As a star-planet system moves through space, the radial component of its velocity :math:`v_r`
    varies in time as a result of the nonzero transverse velocity component :math:`v_t = \\mu D`,
    where :math:`D` is the distance to the system and :math:`\\mu = |\\vec{\\mu}|` is the
    magnitude of the proper motion vector.

    Consequently, the radial velocity component is time-dependent, which leads to a secular
    variation in the observed period between transits, denoted :math:`\\dot{P}_{ \\mathrm{Shk}}`.
    This is known as the Shklovskii Effect [1]_, expressed in the following convenient form
    by Rafikov (2009) [2]_:

    .. math::
        \\dot{P}_{\\mathrm{Shk}} = 20\\left(\\frac{\\mu}{100\\,{\\mathrm{mas\\,
        yr}^{-1}}}\\right)^2\\frac{D}{100}\\frac{P}{3\\,\\mathrm{day}}\\,\\mu {\\mathrm{s\\,
        yr}^{-1}}

    References
    ----------
    .. [1] :cite:t:`Rafikov2009`. https://doi.org/10.1088/0004-637X/700/2/965
    .. [2] :cite:t:`Shklovskii1970`. https://ui.adsabs.harvard.edu/abs/1970SvA....13..562S

    """
    # calculate period change in microseconds per year
    pdot = 20 * (mu / 100) ** 2 * (D / 100) * (P / 3)
    pdot *= 0.001  # convert microseconds to milliseconds

    # return the apparent period derivative in ms/yr
    return pdot


def resolved_binary_rv_trend_from_mass(theta, D, M_B):
    """Calculate a minimum RV trend (acceleration) given properties of a resolved stellar companion.

    This method returns the minimum possible acceleration, induced by a bound stellar companion,
    that may be observed in radial velocity observations of the primary. The secondary star (ie.
    the 'companion') must have been resolved through imaging or astrometric measurements such
    that the angular separation is known. This uses equation (6) from Torres (1999) [1]_,
    under the assumption that the radial velocity trend is due entirely to the gravitational
    influence of the secondary.

    Parameters
    ----------
    theta : float
        The angular separation of the binary in arcseconds.
    D : float
        The distance to the system in parsecs.
    dvdt : float
        The mass of the stellar companion in solar masses.

    Returns
    -------
    float
        The predicted linear RV trend (acceleration) in m/s/day.

    Notes
    -----
    Given the angular separation of a known companion object, the distance to the system, and the
    mass of the companion, this method uses Equation (6) from Torres et al. (1999) [1]_ to estimate
    the minimum possible acceleration in radial velocity observations:

    .. math::
        \\left|\\frac{d(R V)}{d t}\\right| = \\frac{M_B}{5.341 \\times 10^{-6}(D \\rho)^2 \\Phi}

    where :math:`\\rho` is the angular separation of the binary in arcseconds, :math:`D` is the
    distance to the system in parsecs, and :math:`M_B` is the mass of the stellar companion in
    solar masses. In the derivation of this equation, the author makes no assumptions about the mass
    or brightness of the secondary [1]_.

    The parameter :math:`\\Phi` is a function of the eccentricity, longitude of pericenter,
    and inclination of the companion's orbit. Assuming these parameters are unconstrained, this
    method uses the minimum value :math:`\\Phi = \\sqrt{3}\\frac{3}{2}` to determine the minimum
    companion mass.

    References
    ----------
    .. [1] :cite:t:`Torres1999`. https://doi.org/10.1086/316313.

    """
    # calculate the minimum value of the unknown component 'Phi'
    Phi = 3 * np.sqrt(3) / 2

    # calculate the expected linear RV trend (m/s/yr)
    dvdt = M_B / (5.341e-6 * (D * theta) ** 2 * Phi)

    # return the slope in m/s/day
    return dvdt / 365.25


def resolved_binary_mass_from_rv_trend(theta, D, dvdt):
    """Calculate the minimum mass of a resolved secondary star given a measured radial acceleration.

    This method applies to the case where an acceleration has been observed in radial velocity
    observations of a star, and there is a known secondary object for which an angular separation
    has been measured. It uses Equation (6) from Torres (1999) [1]_ to estimate a lower limit for
    the mass of the secondary object, under the assumption that the radial velocity trend is due
    entirely to its gravitational influence.

    Parameters
    ----------
    theta : float
        The angular separation of the binary in arcseconds.
    D : float
        The distance to the system in parsecs.
    dvdt : float
        The observed linear RV trend (acceleration) in m/s/day.

    Returns
    -------
    float
        The mass of the secondary star in solar masses.

    Notes
    -----
    Given the angular separation of a known companion object, the distance to the system, and a
    measured acceleration in radial velocity observations, this method uses Equation (6) from
    Torres (1999) [1]_ to estimate the minimum mass of the companion:

    .. math:: M_B = 5.341 \\times 10^{-6}(D \\rho)^2\\left|\\frac{d(R V)}{d t}\\right| \\Phi

    where :math:`\\rho` is the angular separation of the binary in arcseconds, :math:`D` is the
    distance to the system in parsecs, and :math:`\\frac{d(R V)}{d t}` is the measured radial
    acceleration. In the derivation of this equation, the author makes no assumptions about the
    mass or brightness of the secondary [1]_.

    The parameter :math:`\\Phi` is a function of the eccentricity, longitude of pericenter,
    and inclination of the companion's orbit. Assuming these parameters are unconstrained, this
    method uses the minimum value :math:`\\varphi = \\sqrt{3}\\frac{3}{2}` to determine the minimum
    companion mass.

    References
    ----------
    .. [1] :cite:t:`Torres1999`. https://doi.org/10.1086/316313.

    """
    # convert the slope from m/s/day to ms/s/yr
    dvdt *= 365.25

    # calculate the minimum value of the unknown component 'Phi'
    Phi = 3 * np.sqrt(3) / 2

    # return the mass of the resolved companion star in solar masses
    return 5.341e-6 * (D * theta) ** 2 * dvdt * Phi


def get_semi_major_axis_from_period(P, M_s):
    """Calculate the semi major axis in meters given the orbital period and host star mass.

    Parameters
    ----------
    P : float
        The orbital period in days.
    M_s : float
        The mass of the host star in solar masses.

    Returns
    -------
    float
        The semi major axis of the orbit in meters.

    """
    # unit conversions
    P *= 86400      # days to seconds
    M_s *= M_sun    # solar masses to kg

    # return the semi major axis in meters
    return (G * M_s * P ** 2 / (4 * np.pi ** 2)) ** (1/3)
