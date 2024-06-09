# orbdot
A set of tools for studying secular evolution of exoplanet orbits


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