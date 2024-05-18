# TODO: complete! Except for TDV plot maybe
import corner
import json
import numpy as np
from astropy.time import Time
from matplotlib import pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.gridspec as gridspec
import orbdot.models.ttv_models as ttv
import orbdot.tools.utilities as utl
from orbdot.models.rv_models import radial_velocity


def make_ttv_plot(plot_settings, outfile):
    """ Generate a TTV plot (observed minus calculated)

    This function generates a TTV plot using the provided plot settings and saves it to
    the specified output file.

    Notes
    -----
    This function relies on the :mod:'models.ttv_models' for TTV calculations.

    Parameters
    ----------
    plot_settings : dict
        A dictionary containing plot settings, including data file paths, results files,
        and plot configurations.

    outfile : str
        The path to save the generated TTV plot.

    Raises
    ------
    FileNotFoundError
        If the necessary files for a constant-period TTV fit are missing.

    Returns
    -------
    None

    """
    print('-' * 100)
    print('Generating TTV plot...')
    print('-' * 100)

    plt.rcParams.update(plot_settings["TTV_PLOT"]['rcParams'])

    data_colors = plot_settings['TTV_PLOT']['data_colors']
    dfmt = plot_settings['TTV_PLOT']['data_fmt']
    dms = plot_settings['TTV_PLOT']['data_markersize']
    delw = plot_settings['TTV_PLOT']['data_err_linewidth']
    decap = plot_settings['TTV_PLOT']['data_err_capsize']

    s_alpha = plot_settings['TTV_PLOT']['samples_alpha']
    s_lw = plot_settings['TTV_PLOT']['samples_linewidth']
    m_alpha = plot_settings['TTV_PLOT']['model_alpha']
    m_lw = plot_settings['TTV_PLOT']['model_linewidth']

    bbox = plot_settings['TTV_PLOT']['bbox_inches']
    dpi = plot_settings['TTV_PLOT']['dpi']
    pad_inches = plot_settings['TTV_PLOT']['pad_inches']

    # load full dataset
    data_file = plot_settings['TTV_PLOT']['data_file']
    data = utl.read_TTV_data(data_file)

    try:
        # load constant-period fit results
        with open(plot_settings['TTV_PLOT']['constant_results_file']) as jf:
            rf_c = json.load(jf)
            res_c = rf_c['params']

        # load constant-period samples
        with open(plot_settings['TTV_PLOT']['constant_samples_file']) as jf:
            s_c = json.load(jf)
            s_orb_c = np.array(s_c['orbital_elements'])  # t0, P0, e0, w0, i0, O0

    except FileNotFoundError:
        print('Error: Missing files for constant-period TTV fit. Please ensure that valid '
              'results and sample files are provided in the plot_settings configuration file '
              'to generate the O-C plot.')
        return

    # initialize figure
    if data['epoch_ecl'].size:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
    else:
        fig, ax1 = plt.subplots(1, 1)

    # define a full range of epochs
    preE = plot_settings['TTV_PLOT']['num_epochs_pre_data']
    postE = plot_settings['TTV_PLOT']['num_epochs_post_data']
    epochs_full = np.arange(min(data['epoch']) - preE, max(data['epoch']) + postE, 1)

    # generate best-fit constant period model (transits) over full epoch range
    cmod_full = ttv.constant_period(res_c['t0'][0], res_c['P0'][0], res_c['e0'][0],
                                    res_c['w0'][0], epochs_full, primary=True)

    # plot best-fit constant period model (transits) over full epoch range
    ax1.plot(epochs_full, np.array(cmod_full - cmod_full) * 1440,
             color='darkgrey', label='_', linewidth=m_lw, alpha=m_alpha)

    try:
        # generate best-fit constant-period model (eclipses) over full epoch range
        cmod_full_ecl = ttv.constant_period(res_c['t0'][0], res_c['P0'][0], res_c['e0'][0],
                                        res_c['w0'][0], epochs_full, primary=False)

        # plot best-fit constant period model (eclipses) over full epoch range
        ax2.plot(epochs_full, np.array(cmod_full_ecl - cmod_full_ecl) * 1440,
                 color='darkgrey', label='Constant Period', linewidth=m_lw, alpha=1)

    except UnboundLocalError:
        pass

    # plot 300 random samples from constant-period model fit
    for i in range(np.shape(s_orb_c)[0]):

        # generate constant-period model (transits)
        smod_c = ttv.constant_period(s_orb_c[i, 0], s_orb_c[i, 1], s_orb_c[i, 2],
                                   s_orb_c[i, 3], epochs_full, primary=True)

        # plot constant-period model (transits)
        ax1.plot(epochs_full, np.array(smod_c - cmod_full) * 1440,
                 color='darkgrey', label='_', linewidth=s_lw, alpha=s_alpha)
        try:
            # generate constant-period model (eclipses)
            smod_c_ecl = ttv.constant_period(s_orb_c[i, 0], s_orb_c[i, 1], s_orb_c[i, 2],
                                             s_orb_c[i, 3], epochs_full, primary=False)

            # plot constant-period model (eclipses)
            ax2.plot(epochs_full, np.array(smod_c_ecl - cmod_full_ecl) * 1440,
                     color='darkgrey', label='_', linewidth=s_lw, alpha=s_alpha)

        except UnboundLocalError:
            pass

    try:
        # load orbital decay fit results
        with open(plot_settings['TTV_PLOT']['decay_results_file']) as jf:
            rf_d = json.load(jf)
            res_d = rf_d['params']

        # load orbital decay samples
        with open(plot_settings['TTV_PLOT']['decay_samples_file']) as jf:
            s_d = json.load(jf)
            s_orb_d = np.array(s_d['orbital_elements'])  # t0, P0, e0, w0, i0, O0
            s_dt_d = np.array(s_d['time_dependent'])  # PdE, wdE, edE, idE, OdEc

        # generate best-fit orbital decay model (transits) over full epoch range
        dmod_full = ttv.orbital_decay(res_d['t0'][0], res_d['P0'][0], res_d['PdE'][0],
                                      res_d['e0'][0], res_d['w0'][0], epochs_full, primary=True)

        # plot best-fit orbital decay model (transits) over full epoch range
        ax1.plot(epochs_full, np.array(dmod_full - cmod_full) * 1440,
                 color='#9c3535', label='_', linewidth=m_lw, alpha=m_alpha)
        try:
            # generate best-fit orbital decay model (eclipses) over full epoch range
            dmod_full_ecl = ttv.orbital_decay(res_d['t0'][0], res_d['P0'][0], res_d['PdE'][0],
                                              res_d['e0'][0], res_d['w0'][0], epochs_full,
                                              primary=False)

            # plot best-fit orbital decay model (eclipses) over full epoch range
            ax2.plot(epochs_full, np.array(dmod_full_ecl - cmod_full_ecl) * 1440,
                     color='#9c3535', label='Orbital Decay', linewidth=m_lw, alpha=m_alpha)

        except UnboundLocalError:
            pass

        # plot 300 random samples orbital decay model fit
        for i in range(np.shape(s_orb_d)[0]):

            # generate orbital decay model (transits)
            smod_d = ttv.orbital_decay(s_orb_d[i, 0], s_orb_d[i, 1], s_dt_d[i, 0],
                                        s_orb_d[i, 2], s_orb_d[i, 3], epochs_full, primary=True)

            # generate orbital decay model (transits)
            ax1.plot(epochs_full, np.array(smod_d - cmod_full) * 1440,
                     color='#9c3535', label='_', linewidth=s_lw, alpha=s_alpha)
            try:
                # generate orbital decay model (eclipses)
                smod_d_ecl = ttv.orbital_decay(s_orb_d[i, 0], s_orb_d[i, 1], s_dt_d[i, 0],
                                        s_orb_d[i, 2], s_orb_d[i, 3], epochs_full, primary=False)

                # generate orbital decay model (eclipses)
                ax2.plot(epochs_full, np.array(smod_d_ecl - cmod_full_ecl) * 1440,
                         color='#9c3535', label='_', linewidth=s_lw, alpha=s_alpha)

            except UnboundLocalError:
                pass
    except:
        print(' --> No orbital decay fit results detected.')

    try:
        # load apsidal precession fit results
        with open(plot_settings['TTV_PLOT']['precession_results_file']) as jf:
            rf_p = json.load(jf)
            res_p = rf_p['params']

        # load apsidal precession samples
        with open(plot_settings['TTV_PLOT']['precession_samples_file']) as jf:
            s_p = json.load(jf)
            s_orb_p = np.array(s_p['orbital_elements'])  # t0, P0, e0, w0, i0, O0
            s_dt_p = np.array(s_p['time_dependent'])     # PdE, wdE, edE, idE, OdE

        # generate best-fit apsidal precession model (transits) over full epoch range
        pmod_full = ttv.apsidal_precession(res_p['t0'][0], res_p['P0'][0],
                                           res_p['e0'][0], res_p['w0'][0],
                                           res_p['wdE'][0], epochs_full, primary=True)

        # plot best-fit apsidal precession model (transits) over full epoch range
        ax1.plot(epochs_full, np.array(pmod_full - cmod_full) * 1440,
                 color='cadetblue', label='_', linewidth=m_lw, alpha=m_alpha)

        try:
            # generate best-fit apsidal precession model (eclipses) over full epoch range
            pmod_full_ecl = ttv.apsidal_precession(res_p['t0'][0], res_p['P0'][0],
                                                   res_p['e0'][0], res_p['w0'][0],
                                                   res_p['wdE'][0], epochs_full, primary=False)

            # plot best-fit apsidal precession model (eclipses) over full epoch range
            ax2.plot(epochs_full, np.array(pmod_full_ecl - cmod_full_ecl) * 1440,
                     color='cadetblue', label='Apsidal Precession',
                     linewidth=m_lw, alpha=m_alpha)

        except UnboundLocalError:
            pass

        # plot 300 random samples apsidal precession model fit
        for i in range(np.shape(s_orb_p)[0]):

            # generate apsidal precession model (transits)
            smod_p = ttv.apsidal_precession(s_orb_p[i, 0], s_orb_p[i, 1],
                                            s_orb_p[i, 2], s_orb_p[i, 3],
                                            s_dt_p[i, 1], epochs_full, primary=True)

            # plot apsidal precession model (transits)
            ax1.plot(epochs_full, np.array(smod_p - cmod_full) * 1440,
                     color='cadetblue', label='_', linewidth=s_lw, alpha=s_alpha)
            try:
                # generate apsidal precession model (eclipses)
                smod_p_ecl = ttv.apsidal_precession(s_orb_p[i, 0], s_orb_p[i, 1], s_orb_p[i, 2],
                                                    s_orb_p[i, 3], s_dt_p[i, 1], epochs_full,
                                                    primary=False)
                # plot apsidal precession model (eclipses)
                ax2.plot(epochs_full, np.array(smod_p_ecl - cmod_full_ecl) * 1440,
                         color='cadetblue', label='_', linewidth=s_lw, alpha=s_alpha)

            except UnboundLocalError:
                pass

    except:
        print(' --> No apsidal precession fit results detected.\n')

    # generate best-fit constant-period model (transits)
    cmod_obs = ttv.constant_period(res_c['t0'][0], res_c['P0'][0],
                                   res_c['e0'][0], res_c['w0'][0],
                                   data['epoch'], primary=True)

    # calculate O-C values for transit data
    oc = data['bjd'] - cmod_obs

    try:
        # generate best-fit constant-period model (eclipses)

        cmod_obs_ecl = ttv.constant_period(res_c['t0'][0], res_c['P0'][0],
                                           res_c['e0'][0], res_c['w0'][0],
                                           data['epoch_ecl'], primary=False)

        # calculate O-C values for eclipse data
        oc_ecl = data['bjd_ecl'] - cmod_obs_ecl

    except UnboundLocalError:
        pass

    # plot data, separating sources into different colours and labels
    sources_unique = np.unique(list(data['ss']) + list(data['ss_ecl']))
    num_sources = len(sources_unique)

    # iterate through transit data sources
    for i in range(num_sources):

        # plot transit data
        indx = np.where(data['ss'] == sources_unique[i])[0]
        ax1.errorbar(np.take(data['epoch'], indx), np.array(np.take(oc, indx)) * 1440,
                     yerr=np.take(data['err'], indx) * 1440,
                     label=sources_unique[i], color=data_colors[i], ecolor=data_colors[i],
                     fmt=dfmt, markersize=dms, elinewidth=delw, capsize=decap)

    try:
        # iterate through eclipse data sources
        for i in range(num_sources):

            # plot eclipse data
            indx = np.where(data['ss_ecl'] == sources_unique[i])[0]
            ax2.errorbar(np.take(data['epoch_ecl'], indx), np.array(np.take(oc_ecl, indx)) * 1440,
                         yerr=np.take(data['err_ecl'], indx) * 1440,
                         label=sources_unique[i], color=data_colors[i], ecolor=data_colors[i],
                         fmt=dfmt, markersize=dms, elinewidth=delw, capsize=decap)

    except UnboundLocalError:
        pass

    # plot data removed in sigma clipping process
    try:
        clipped_data = utl.read_TTV_data(plot_settings['TTV_PLOT']['clipped_data_file'])
    except FileNotFoundError:
        pass

    try:
        # generate best-fit constant-period model for clipped epochs (transits)
        cmod_obs_clipped = ttv.constant_period(res_c['t0'][0], res_c['P0'][0],
                                               res_c['e0'], res_c['w0'], clipped_data['epoch'])

        # calculate O-C values for clipped transit data
        clipped_OCs = clipped_data['bjd'] - cmod_obs_clipped

        # plot clipped data O-C
        ax1.errorbar(clipped_data['epoch'], np.array(clipped_OCs) * 1440,
                     yerr=clipped_data['err'] * 1440, label='Excluded',
                     fmt='x', markersize='5', ecolor='r', elinewidth=.8, color='r')

    except UnboundLocalError:
        pass

    # plot vertical lines for reference dates
    trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)

    date_refs = plot_settings['TTV_PLOT']['reference_dates']
    t_ref = Time(date_refs, format='iso', in_subfmt='date')
    dates_jd = t_ref.to_value('jd', subfmt='float')
    dates_e = utl.calculate_epochs(res_c['t0'][0], res_c['P0'][0], dates_jd, primary=True)

    t_t0 = Time([str(res_c['t0'][0])], format='jd', scale='tdb')  # add t0
    date_refs = list(date_refs)
    date_refs.append(t_t0.to_value('iso', subfmt='date')[0])
    dates_e.append(0)

    dates_year = [s.split('-')[0] for s in date_refs]  # remove month and day
    for i, date in enumerate(dates_e):
        ax1.axvline(x=date, linewidth=0.5, linestyle='-', color='grey')
        ax1.text(date, 0.05, dates_year[i], rotation=90, fontsize=12,
                 ha='right', va='bottom', transform=trans)
        try:
            ax2.axvline(x=date, linewidth=0.5, linestyle='-', color='grey', label='_')
        except UnboundLocalError:
            pass

    # finish plot
    ax1.set_title('{}'.format(plot_settings['TTV_PLOT']['title'] + ' Transits'))
    plt.xlabel('Epoch')
    ax1.set_ylabel('Timing Deviation (minutes)')

    plt.xlim(epochs_full[0], epochs_full[-1])
    plt.ylim(plot_settings['TTV_PLOT']['y_axis_limits'][0],
             plot_settings['TTV_PLOT']['y_axis_limits'][1])

    ax1.legend()

    try:
        ax2.set_title('{}'.format(plot_settings['TTV_PLOT']['title']) + ' Eclipses')
        ax2.set_ylabel('Timing Deviation (minutes)')
        legend2 = ax2.legend()
        for line in legend2.get_lines():
            line.set_linewidth(6)
            line.set_alpha(1)
    except UnboundLocalError:
        pass

    plt.savefig(outfile, bbox_inches=bbox, dpi=dpi, pad_inches=pad_inches)
    plt.close()

    print('Done!\n')


def make_rv_plots(plot_settings, outfile):
    """Radial velocity plot.

    Generates a 3-part plot of the RV models and data:

    1. The radial velocity measurements vs. time, with an option to plot the
       best-fit RV curve. The data are shifted by the zero velocity for each instrument.

    2. The residuals of the data after subtracting the signal of only the planet
       on the star. That is, subtracting the best-fit RV model without including the
       time-dependent radial velocity terms ('dvdt' and 'ddvdt').

    3. A phase-folded plot of the planet's orbit from 300 random samples from the model fit.

    Parameters
    ----------
    outfile : str
        The filename to save the plot.
    plot_settings : dict
        A dictionary containing plot settings.

    Raises
    ------
    FileNotFoundError
        If the specified results file or samples file does not exist.

    Returns
    -------
    None

    """
    data_colors = plot_settings['RV_PLOT']['data_colors']
    dfmt = plot_settings['RV_PLOT']['data_fmt']
    dms = plot_settings['RV_PLOT']['data_markersize']
    delw = plot_settings['RV_PLOT']['data_err_linewidth']
    decap = plot_settings['RV_PLOT']['data_err_capsize']

    s_alpha = plot_settings['RV_PLOT']['samples_alpha']
    s_lw = plot_settings['RV_PLOT']['samples_linewidth']
    m_alpha = plot_settings['RV_PLOT']['model_alpha']
    m_lw = plot_settings['RV_PLOT']['model_linewidth']

    bbox = plot_settings['RV_PLOT']['bbox_inches']
    dpi = plot_settings['RV_PLOT']['dpi']
    pad_inches = plot_settings['RV_PLOT']['pad_inches']

    results_file = plot_settings['RV_PLOT']['rv_results_file']
    samples_file = plot_settings['RV_PLOT']['rv_samples_file']
    data_file = plot_settings['RV_PLOT']['data_file']
    data = utl.read_RV_data(data_file)

    print('Generating RV plot...')
    plt.rcParams.update(plot_settings["RV_PLOT"]['rcParams'])

    # load fit results
    with open(results_file) as jf:
        rf = json.load(jf)

    res = rf['params']

    # load random samples
    with open(samples_file) as jf:
        s = json.load(jf)
        s_orb = np.array(s['orbital_elements'])  # t0, P0, e0, w0, i0, O0
        s_dt = np.array(s['time_dependent'])  # PdE, wdE, edE, idE, OdE
        s_rv = np.array(s['radial_velocity'])  # K, v0, jit, dvdt, ddvdt

    # start subplots
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1], width_ratios=[4, 5])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)
    ax3 = fig.add_subplot(gs[:, 1])

    # plot the best-fit radial velocity model
    if plot_settings["RV_PLOT"]['show_RV_curve'] == 'True':
        # define a range of time values spanning the observational baseline
        times_observed = data['trv_all']
        times_all = np.arange(min(times_observed),
                              max(times_observed), 0.01 * res['P0'][0])

        # compute the model with the best-fit parameters (without systemic velocity 'v0')
        mod_all = radial_velocity(t0=res['t0'][0], P0=res['P0'][0],
                                         e=res['e0'][0], w0=res['w0'][0],
                                         K=res['K'][0], v0=0.0, dvdt=res['dvdt'][0],
                                         ddvdt=res['ddvdt'][0], time=times_all,
                                         dPdE=res['PdE'][0], dwdE=res['wdE'][0])

        # plot the best-fit model
        ax1.plot(times_all - res['t0'][0], mod_all, label='Best-Fit Model',
                 color='darkgrey', alpha=m_alpha, linewidth=m_lw, zorder=0)

    # plot data, separating instruments into different colours and labels
    for i in data['ss_order']:
        # shift the RV measurement times by t0
        x = data['trv'][i] - res['t0'][0]

        # plot the RV data, adjusted by the best-fit systemic velocity 'v0'
        y = data['rvs'][i] - res['v0_' + data['ss_tags'][i]][0]
        ax1.errorbar(x, y, yerr=data['err'][i], label=data['ss_names'][i], color=data_colors[i],
                     fmt=dfmt, markersize=dms, elinewidth=delw, capsize=decap, zorder=1)

        # generate a best-fit model for each time step
        mod = radial_velocity(t0=res['t0'][0], P0=res['P0'][0],
                                     e=res['e0'][0], w0=res['w0'][0], K=res['K'][0],
                                     v0=res['v0_' + data['ss_tags'][i]][0],
                                     dvdt=res['dvdt'][0], ddvdt=res['ddvdt'][0],
                                     time=data['trv'][i], dPdE=res['PdE'][0],
                                     dwdE=res['wdE'][0])

        # calculate the residual (data - model)
        residuals = data['rvs'][i] - mod

        # plot the residuals
        ax2.errorbar(x, residuals, yerr=data['err'][i], label=data['ss_names'][i],
                     color=data_colors[i], fmt=dfmt, elinewidth=delw, markersize=dms,
                     capsize=decap, capthick=1)

        # phase-fold the data (subtracting long-term trends)
        x_fold = ((data['trv'][i] - res['t0'][0]) % res['P0'][0]) / res['P0'][0]
        y_fold = data['rvs'][i] - \
                 res['v0_' + data['ss_tags'][i]][0] - \
                 res['dvdt'][0] * (data['trv'][i] - res['t0'][0]) - \
                 0.5 * res['ddvdt'][0] * (data['trv'][i] - res['t0'][0]) ** 2

        # plot the phase-folded data
        ax3.errorbar(x_fold, y_fold, yerr=data['err'][i], label=data['ss_names'][i],
                     color=data_colors[i], fmt='.', markersize=dms, elinewidth=delw, capsize=decap)

    # phase-fold best-fit model (planet orbit only)
    times_fold = np.arange(res['t0'][0], res['t0'][0] + res['P0'][0], 0.01)
    mod_fold = radial_velocity(t0=res['t0'][0], P0=res['P0'][0], e=res['e0'][0],
                                      w0=res['w0'][0], K=res['K'][0], v0=0.0,
                                      dvdt=0.0, ddvdt=0.0, time=times_fold,
                                      dPdE=res['PdE'][0], dwdE=res['wdE'][0])

    # plot the phase-folded model
    ax3.plot((times_fold - res['t0'][0]) / res['P0'][0], mod_fold,
             label='RV Model', color='darkgrey', linewidth=m_lw)

    # repeat all of the above for the 300 random posterior samples
    for i in range(np.shape(s_orb)[0]):
        # define a time array
        time_s = np.arange(s_orb[i, 0], s_orb[i, 0] + s_orb[i, 1], 0.01)

        # calculate the best-fit model for the planet signal only
        smod = radial_velocity(t0=s_orb[i, 0], P0=s_orb[i, 1], e=s_orb[i, 2],
                                      w0=s_orb[i, 3], K=s_rv[i, 0], v0=0.0,
                                      dvdt=s_rv[i, 3], ddvdt=s_rv[i, 4], time=time_s,
                                      dPdE=s_dt[i, 0], dwdE=s_dt[i, 1])

        # plot the phase-folded model
        ax3.plot((time_s - s_orb[i, 0]) / s_orb[i, 1], smod, alpha=s_alpha,
                 label='_', color='darkgrey', linewidth=s_lw)

    # plot zero lines for reference
    ax1.axhline(y=0, linestyle='--', color='grey', linewidth=1)
    ax2.axhline(y=0, linestyle='--', color='grey', linewidth=1)
    ax1.axvline(x=0, linestyle='-', color='grey', linewidth=0.5)
    ax2.axvline(x=0, linestyle='-', color='grey', linewidth=0.5)

    # add text for t0 date
    from astropy.time import Time
    t = Time([str(res['t0'][0])], format='jd', scale='tdb')
    ax2.text(.99, .01, 't$_0$ = ' + t.to_value('iso', subfmt='date')[0],
             ha='right', va='bottom', transform=ax2.transAxes, fontsize=11)

    # finish plots
    plt.suptitle(plot_settings['RV_PLOT']['title'], fontsize=18)
    ax1.set_title('All Time')
    ax1.set_ylabel('Radial Velocity (m/s)')
    ax2.set_ylabel('Residuals (m/s)')
    ax2.set_xlabel('BJD $-$ t$_0$ (days)')
    ax3.set_ylabel('Radial Velocity (m/s)')
    ax3.set_xlabel('Phase (from t$_c$)'.format(res['t0'][0]))
    ax3.set_title('Phase Folded')
    ax3.legend(loc='lower right')

    plt.savefig(outfile, bbox_inches=bbox, dpi=dpi, pad_inches=pad_inches)
    plt.close()


def corner_plot(dic, samples, params, outfile):
    """Generate a corner plot from the weighted posterior samples.

    This method generates a corner plot from the weighted posterior samples using the
    corner.py package written by Daniel Foreman-Mackey [1].

    Parameters
    ----------
    dic : dict
        Dictionary containing the results of the sampler.
    samples : array_like
        Array containing the weighted posterior samples.
    params : list
        List of the parameter names.
    outfile : str
        File name of the saved plot.

    Returns
    -------
    None

    References
    ----------
    .. [1] Daniel Foreman-Mackey (2016). https://doi.org/10.21105/joss.00024

    """
    ndim = len(params)

    # generate corner plot
    figure = corner.corner(samples, labels=params, color='k', top_ticks=True,
                           how_titles=False, plot_contours=True, labelpad=0.1)
    axes = np.array(figure.axes).reshape((ndim, ndim))

    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(dic[params[i]][0], color="darkred")

    plt.savefig(outfile, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    return