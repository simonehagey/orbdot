"""
Plots
=====
This module contains the methods that generate the TTV, RV, and corner plots.
"""

import csv
import json
import corner
import numpy as np
import scipy.signal as sci
from astropy.time import Time
from astropy.timeseries import LombScargle
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as transforms
import orbdot.tools.utilities as utl
import orbdot.models.rv_models as rv
import orbdot.models.ttv_models as ttv

def make_ttv_plot(plot_settings, outfile, suffix=''):
    """Generate a TTV plot (observed minus calculated).

    This function generates a TTV plot using the provided plot settings and saves it to
    the specified output file name.

    Parameters
    ----------
    plot_settings : dict
        A dictionary containing plot settings, such as data file paths, results files,
        and plot configurations.
    outfile : str
        The path to save the generated TTV plot.
    suffix : str, optional
        Optional string that matches to particular model fit(s).

    Returns
    -------
    None

    """
    print('-' * 100)
    print('Generating TTV plot...')
    print('-' * 100)

    # load plot settings
    plt.rcParams.update(plot_settings['TTV_PLOT']['rcParams'])

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
    data_file = plot_settings['TTV_PLOT']['data_file'+suffix]
    data = utl.read_ttv_data(data_file)

    try:

        # load constant-period fit results
        with open(plot_settings['TTV_PLOT']['ttv_constant_results_file'+suffix]) as jf:
            rf_c = json.load(jf)
            res_c = rf_c['params']

        # load constant-period samples
        s_orb_c, s_tdp_c, s_rv_c = \
            read_random_samples(plot_settings['TTV_PLOT']['ttv_constant_samples_file'+suffix])

    except KeyError:
        print('ERROR: Missing \'*_results{}.json\' file for constant-period TTV fit. The O-C plot '
              'cannot be generated without first fitting the constant-period TTV model.\n\n'
              .format(suffix))
        return

    transits = False
    eclipses = False

    if data['epoch'].size > 0 and data['epoch_ecl'].size > 0:
        transits = True
        eclipses = True

    elif data['epoch'].size > 0:
        transits = True

    elif data['epoch_ecl'].size > 0:
        eclipses = True

    if transits and eclipses:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
    else:
        fig, ax1 = plt.subplots(1, 1)

    # define a full range of epochs
    preE = plot_settings['TTV_PLOT']['num_epochs_pre_data']
    postE = plot_settings['TTV_PLOT']['num_epochs_post_data']

    if transits:
        epochs_full = np.arange(min(data['epoch']) - preE, max(data['epoch']) + postE, 1)

    else:
        epochs_full = np.arange(min(data['epoch_ecl']) - preE, max(data['epoch_ecl']) + postE, 1)

    if transits and eclipses:
        # generate best-fit constant period model (transits) over full epoch range
        cmod_full = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0], res_c['e0'][0],
                                     res_c['w0'][0], epochs_full, primary=True)

        # plot best-fit constant period model (transits) over full epoch range
        ax1.plot(epochs_full, np.array(cmod_full - cmod_full) * 1440,
                 color='darkgrey', label='_', linewidth=m_lw, alpha=m_alpha)

        # generate best-fit constant-period model (eclipses) over full epoch range
        cmod_full_ecl = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0], res_c['e0'][0],
                                        res_c['w0'][0], epochs_full, primary=False)

        # plot best-fit constant period model (eclipses) over full epoch range
        ax2.plot(epochs_full, np.array(cmod_full_ecl - cmod_full_ecl) * 1440,
                 color='darkgrey', label='Constant Period', linewidth=m_lw, alpha=1)

    elif transits and not eclipses:
        # generate best-fit constant period model (transits) over full epoch range
        cmod_full = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0], res_c['e0'][0],
                                     res_c['w0'][0], epochs_full, primary=True)

        # plot best-fit constant period model (transits) over full epoch range
        ax1.plot(epochs_full, np.array(cmod_full - cmod_full) * 1440,
                 color='darkgrey', label='Constant Period', linewidth=m_lw, alpha=m_alpha)

    elif eclipses and not transits:
        # generate best-fit constant-period model (eclipses) over full epoch range
        cmod_full_ecl = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0], res_c['e0'][0],
                                        res_c['w0'][0], epochs_full, primary=False)

        # plot best-fit constant period model (eclipses) over full epoch range
        ax1.plot(epochs_full, np.array(cmod_full_ecl - cmod_full_ecl) * 1440,
                 color='darkgrey', label='Constant Period', linewidth=m_lw, alpha=1)

    # plot 300 random samples from constant-period model fit
    for i in range(np.shape(s_orb_c)[0]):

        if transits and eclipses:

            # generate constant-period model (transits)
            smod_c = ttv.ttv_constant(s_orb_c[i][0], s_orb_c[i][1], s_orb_c[i][2],
                                       s_orb_c[i][3], epochs_full, primary=True)

            # plot constant-period model (transits)
            ax1.plot(epochs_full, np.array(smod_c - cmod_full) * 1440,
                     color='darkgrey', label='_', linewidth=s_lw, alpha=s_alpha)

            # generate constant-period model (eclipses)
            smod_c_ecl = ttv.ttv_constant(s_orb_c[i][0], s_orb_c[i][1], s_orb_c[i][2],
                                             s_orb_c[i][3], epochs_full, primary=False)

            # plot constant-period model (eclipses)
            ax2.plot(epochs_full, np.array(smod_c_ecl - cmod_full_ecl) * 1440,
                     color='darkgrey', label='_', linewidth=s_lw, alpha=s_alpha)

        elif transits and not eclipses:
            # generate constant-period model (transits)
            smod_c = ttv.ttv_constant(s_orb_c[i][0], s_orb_c[i][1], s_orb_c[i][2],
                                      s_orb_c[i][3], epochs_full, primary=True)

            # plot constant-period model (transits)
            ax1.plot(epochs_full, np.array(smod_c - cmod_full) * 1440,
                     color='darkgrey', label='_', linewidth=s_lw, alpha=s_alpha)

        elif eclipses and not transits:
            # generate constant-period model (eclipses)
            smod_c_ecl = ttv.ttv_constant(s_orb_c[i][0], s_orb_c[i][1], s_orb_c[i][2],
                                             s_orb_c[i][3], epochs_full, primary=False)

            # plot constant-period model (eclipses)
            ax1.plot(epochs_full, np.array(smod_c_ecl - cmod_full_ecl) * 1440,
                     color='darkgrey', label='_', linewidth=s_lw, alpha=s_alpha)

    try:

        # load orbital decay fit results
        with open(plot_settings['TTV_PLOT']['ttv_decay_results_file'+suffix]) as jf:
            rf_d = json.load(jf)
            res_d = rf_d['params']

        # load orbital decay samples
        s_orb_d, s_tdp_d, s_rv_d \
            = read_random_samples(plot_settings['TTV_PLOT']['ttv_decay_samples_file'+suffix])

        if transits and eclipses:
            # generate best-fit orbital decay model (transits) over full epoch range
            dmod_full = ttv.ttv_decay(res_d['t0'][0], res_d['P0'][0], res_d['PdE'][0],
                                      res_d['e0'][0], res_d['w0'][0], epochs_full, primary=True)

            # plot best-fit orbital decay model (transits) over full epoch range
            ax1.plot(epochs_full, np.array(dmod_full - cmod_full) * 1440,
                     color='#9c3535', label='_', linewidth=m_lw, alpha=m_alpha)

            # generate best-fit orbital decay model (eclipses) over full epoch range
            dmod_full_ecl = ttv.ttv_decay(res_d['t0'][0], res_d['P0'][0], res_d['PdE'][0],
                                              res_d['e0'][0], res_d['w0'][0], epochs_full,
                                              primary=False)

            # plot best-fit orbital decay model (eclipses) over full epoch range
            ax2.plot(epochs_full, np.array(dmod_full_ecl - cmod_full_ecl) * 1440,
                     color='#9c3535', label='Orbital Decay', linewidth=m_lw, alpha=m_alpha)

        elif transits and not eclipses:

            # generate best-fit orbital decay model (transits) over full epoch range
            dmod_full = ttv.ttv_decay(res_d['t0'][0], res_d['P0'][0], res_d['PdE'][0],
                                      res_d['e0'][0], res_d['w0'][0], epochs_full, primary=True)

            # plot best-fit orbital decay model (transits) over full epoch range
            ax1.plot(epochs_full, np.array(dmod_full - cmod_full) * 1440,
                     color='#9c3535', label='Orbital Decay', linewidth=m_lw, alpha=m_alpha)

        elif eclipses and not transits:
            # generate best-fit orbital decay model (eclipses) over full epoch range
            dmod_full_ecl = ttv.ttv_decay(res_d['t0'][0], res_d['P0'][0], res_d['PdE'][0],
                                              res_d['e0'][0], res_d['w0'][0], epochs_full,
                                              primary=False)

            # plot best-fit orbital decay model (eclipses) over full epoch range
            ax1.plot(epochs_full, np.array(dmod_full_ecl - cmod_full_ecl) * 1440,
                     color='#9c3535', label='Orbital Decay', linewidth=m_lw, alpha=m_alpha)


        # plot 300 random samples orbital decay model fit
        for i in range(np.shape(s_orb_d)[0]):

            if transits and eclipses:

                # generate orbital decay model (transits)
                smod_d = ttv.ttv_decay(s_orb_d[i][0], s_orb_d[i][1], s_tdp_d[i][0],
                                            s_orb_d[i][2], s_orb_d[i][3], epochs_full, primary=True)

                # generate orbital decay model (transits)
                ax1.plot(epochs_full, np.array(smod_d - cmod_full) * 1440,
                         color='#9c3535', label='_', linewidth=s_lw, alpha=s_alpha)

                # generate orbital decay model (eclipses)
                smod_d_ecl = ttv.ttv_decay(s_orb_d[i][0], s_orb_d[i][1], s_tdp_d[i][0],
                                        s_orb_d[i][2], s_orb_d[i][3], epochs_full, primary=False)

                # generate orbital decay model (eclipses)
                ax2.plot(epochs_full, np.array(smod_d_ecl - cmod_full_ecl) * 1440,
                         color='#9c3535', label='_', linewidth=s_lw, alpha=s_alpha)

            elif transits and not eclipses:
                # generate orbital decay model (transits)
                smod_d = ttv.ttv_decay(s_orb_d[i][0], s_orb_d[i][1], s_tdp_d[i][0],
                                       s_orb_d[i][2], s_orb_d[i][3], epochs_full, primary=True)

                # generate orbital decay model (transits)
                ax1.plot(epochs_full, np.array(smod_d - cmod_full) * 1440,
                         color='#9c3535', label='_', linewidth=s_lw, alpha=s_alpha)

            elif eclipses and not transits:

                # generate orbital decay model (eclipses)
                smod_d_ecl = ttv.ttv_decay(s_orb_d[i][0], s_orb_d[i][1], s_tdp_d[i][0],
                                        s_orb_d[i][2], s_orb_d[i][3], epochs_full, primary=False)

                # generate orbital decay model (eclipses)
                ax1.plot(epochs_full, np.array(smod_d_ecl - cmod_full_ecl) * 1440,
                         color='#9c3535', label='_', linewidth=s_lw, alpha=s_alpha)

    except KeyError:
        print(' --> No orbital decay fit results detected.')

    try:

        # load apsidal precession fit results
        with open(plot_settings['TTV_PLOT']['ttv_precession_results_file'+suffix]) as jf:
            rf_p = json.load(jf)
            res_p = rf_p['params']

        # load apsidal precession samples
        s_orb_p, s_tdp_p, s_rv_p = \
            read_random_samples(plot_settings['TTV_PLOT']['ttv_precession_samples_file'+suffix])

        if transits and eclipses:

            # generate best-fit apsidal precession model (transits) over full epoch range
            pmod_full = ttv.ttv_precession(res_p['t0'][0], res_p['P0'][0],
                                           res_p['e0'][0], res_p['w0'][0],
                                           res_p['wdE'][0], epochs_full, primary=True)

            # plot best-fit apsidal precession model (transits) over full epoch range
            ax1.plot(epochs_full, np.array(pmod_full - cmod_full) * 1440,
                     color='cadetblue', label='_', linewidth=m_lw, alpha=m_alpha)

            # generate best-fit apsidal precession model (eclipses) over full epoch range
            pmod_full_ecl = ttv.ttv_precession(res_p['t0'][0], res_p['P0'][0],
                                                   res_p['e0'][0], res_p['w0'][0],
                                                   res_p['wdE'][0], epochs_full, primary=False)

            # plot best-fit apsidal precession model (eclipses) over full epoch range
            ax2.plot(epochs_full, np.array(pmod_full_ecl - cmod_full_ecl) * 1440,
                     color='cadetblue', label='Apsidal Precession',
                     linewidth=m_lw, alpha=m_alpha)

        elif transits and not eclipses:
            # generate best-fit apsidal precession model (transits) over full epoch range
            pmod_full = ttv.ttv_precession(res_p['t0'][0], res_p['P0'][0],
                                           res_p['e0'][0], res_p['w0'][0],
                                           res_p['wdE'][0], epochs_full, primary=True)

            # plot best-fit apsidal precession model (transits) over full epoch range
            ax1.plot(epochs_full, np.array(pmod_full - cmod_full) * 1440,
                     color='cadetblue', label='Apsidal Precession', linewidth=m_lw, alpha=m_alpha)

        elif eclipses and not transits:
            # generate best-fit apsidal precession model (eclipses) over full epoch range
            pmod_full_ecl = ttv.ttv_precession(res_p['t0'][0], res_p['P0'][0],
                                                   res_p['e0'][0], res_p['w0'][0],
                                                   res_p['wdE'][0], epochs_full, primary=False)

            # plot best-fit apsidal precession model (eclipses) over full epoch range
            ax1.plot(epochs_full, np.array(pmod_full_ecl - cmod_full_ecl) * 1440,
                     color='cadetblue', label='Apsidal Precession',
                     linewidth=m_lw, alpha=m_alpha)

        # plot 300 random samples apsidal precession model fit
        for i in range(np.shape(s_orb_p)[0]):

            if transits and eclipses:
                # generate apsidal precession model (transits)
                smod_p = ttv.ttv_precession(s_orb_p[i][0], s_orb_p[i][1],
                                            s_orb_p[i][2], s_orb_p[i][3],
                                            s_tdp_p[i][1], epochs_full, primary=True)

                # plot apsidal precession model (transits)
                ax1.plot(epochs_full, np.array(smod_p - cmod_full) * 1440,
                         color='cadetblue', label='_', linewidth=s_lw, alpha=s_alpha)

                # generate apsidal precession model (eclipses)
                smod_p_ecl = ttv.ttv_precession(s_orb_p[i][0], s_orb_p[i][1], s_orb_p[i][2],
                                                    s_orb_p[i][3], s_tdp_p[i][1], epochs_full,
                                                    primary=False)

                # plot apsidal precession model (eclipses)
                ax2.plot(epochs_full, np.array(smod_p_ecl - cmod_full_ecl) * 1440,
                         color='cadetblue', label='_', linewidth=s_lw, alpha=s_alpha)

            elif transits and not eclipses:
                # generate apsidal precession model (transits)
                smod_p = ttv.ttv_precession(s_orb_p[i][0], s_orb_p[i][1],
                                            s_orb_p[i][2], s_orb_p[i][3],
                                            s_tdp_p[i][1], epochs_full, primary=True)

                # plot apsidal precession model (transits)
                ax1.plot(epochs_full, np.array(smod_p - cmod_full) * 1440,
                         color='cadetblue', label='_', linewidth=s_lw, alpha=s_alpha)

            elif eclipses and not transits:
                # generate apsidal precession model (eclipses)
                smod_p_ecl = ttv.ttv_precession(s_orb_p[i][0], s_orb_p[i][1], s_orb_p[i][2],
                                                    s_orb_p[i][3], s_tdp_p[i][1], epochs_full,
                                                    primary=False)

                # plot apsidal precession model (eclipses)
                ax1.plot(epochs_full, np.array(smod_p_ecl - cmod_full_ecl) * 1440,
                         color='cadetblue', label='_', linewidth=s_lw, alpha=s_alpha)

    except KeyError:
        print(' --> No apsidal precession fit results detected.\n')

    if transits and eclipses:
        # generate best-fit constant-period model (transits)
        cmod_obs = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0],
                                    res_c['e0'][0], res_c['w0'][0],
                                    data['epoch'], primary=True)

        # calculate O-C values for transit data
        oc = data['bjd'] - cmod_obs

        # generate best-fit constant-period model (eclipses)
        cmod_obs_ecl = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0],
                                           res_c['e0'][0], res_c['w0'][0],
                                           data['epoch_ecl'], primary=False)

        # calculate O-C values for eclipse data
        oc_ecl = data['bjd_ecl'] - cmod_obs_ecl

    elif transits and not eclipses:
        # generate best-fit constant-period model (transits)
        cmod_obs = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0],
                                    res_c['e0'][0], res_c['w0'][0],
                                    data['epoch'], primary=True)

        # calculate O-C values for transit data
        oc = data['bjd'] - cmod_obs

    elif eclipses and not transits:
        # generate best-fit constant-period model (eclipses)
        cmod_obs_ecl = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0],
                                           res_c['e0'][0], res_c['w0'][0],
                                           data['epoch_ecl'], primary=False)

        # calculate O-C values for eclipse data
        oc_ecl = data['bjd_ecl'] - cmod_obs_ecl

    # plot data, separating sources into different colours and labels
    sources_unique = np.unique(list(data['src']) + list(data['src_ecl']))
    num_sources = len(sources_unique)

    if transits and eclipses:

        # iterate through transit data sources
        for i in range(num_sources):
            # plot transit data
            indx = np.where(data['src'] == sources_unique[i])[0]
            ax1.errorbar(np.take(data['epoch'], indx), np.array(np.take(oc, indx)) * 1440,
                         yerr=np.take(data['err'], indx) * 1440,
                         label=sources_unique[i], color=data_colors[i], ecolor=data_colors[i],
                         fmt=dfmt, markersize=dms, elinewidth=delw, capsize=decap)

        # iterate through eclipse data sources
        for i in range(num_sources):
            # plot eclipse data
            indx = np.where(data['src_ecl'] == sources_unique[i])[0]
            ax2.errorbar(np.take(data['epoch_ecl'], indx), np.array(np.take(oc_ecl, indx)) * 1440,
                         yerr=np.take(data['err_ecl'], indx) * 1440,
                         label='_', color=data_colors[i], ecolor=data_colors[i],
                         fmt=dfmt, markersize=dms, elinewidth=delw, capsize=decap)

    elif transits and not eclipses:
        # iterate through transit data sources
        for i in range(num_sources):
            # plot transit data
            indx = np.where(data['src'] == sources_unique[i])[0]
            ax1.errorbar(np.take(data['epoch'], indx), np.array(np.take(oc, indx)) * 1440,
                         yerr=np.take(data['err'], indx) * 1440,
                         label=sources_unique[i], color=data_colors[i], ecolor=data_colors[i],
                         fmt=dfmt, markersize=dms, elinewidth=delw, capsize=decap)

    elif eclipses and not transits:
        # iterate through eclipse data sources
        for i in range(num_sources):
            # plot eclipse data
            indx = np.where(data['src_ecl'] == sources_unique[i])[0]
            ax1.errorbar(np.take(data['epoch_ecl'], indx), np.array(np.take(oc_ecl, indx)) * 1440,
                         yerr=np.take(data['err_ecl'], indx) * 1440,
                         label='_', color=data_colors[i], ecolor=data_colors[i],
                         fmt=dfmt, markersize=dms, elinewidth=delw, capsize=decap)

    # plot data removed in sigma clipping process
    try:
        clipped_data = utl.read_ttv_data(plot_settings['TTV_PLOT']['clipped_data_file'])

    except FileNotFoundError:
        pass

    try:
        # generate best-fit constant-period model for clipped epochs (transits)
        cmod_obs_clipped = ttv.ttv_constant(res_c['t0'][0], res_c['P0'][0],
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

    # remove month and day
    dates_year = [s.split('-')[0] for s in date_refs]
    for i, date in enumerate(dates_e):

        ax1.axvline(x=date, linewidth=0.5, linestyle='-', color='grey')
        ax1.text(date, 0.05, dates_year[i], rotation=90, fontsize=12,
                 ha='right', va='bottom', transform=trans)

        try:
            ax2.axvline(x=date, linewidth=0.5, linestyle='-', color='grey', label='_')

        except UnboundLocalError:
            pass

    # finish plot
    plt.xlabel('Epoch')
    ax1.set_ylabel('Timing Deviation (minutes)')

    plt.xlim(epochs_full[0], epochs_full[-1])
    plt.ylim(plot_settings['TTV_PLOT']['y_axis_limits'][0],
             plot_settings['TTV_PLOT']['y_axis_limits'][1])

    if eclipses and not transits:
        ax1.set_title('{}'.format(plot_settings['TTV_PLOT']['title'] + ' Eclipses'))

    elif transits:
        ax1.set_title('{}'.format(plot_settings['TTV_PLOT']['title'] + ' Transits'))

    ax1.legend()
    legend1 = ax1.legend()
    for line in legend1.get_lines():
        line.set_linewidth(6)
        line.set_alpha(1)

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


def make_rv_plots(plot_settings, outfile, suffix='', model='constant'):
    """Radial velocity plot.

    Generates a 3-part plot of the RV models and data:

    1. The radial velocity measurements vs. time, with an option to plot the
       best-fit RV curve. The data are shifted by the zero velocity for each instrument.

    2. The residuals of the data after subtracting the signal of only the planet
       on the star. That is, subtracting the best-fit RV model without including the
       time-dependent radial velocity terms ('dvdt' and 'ddvdt').

    3. A phase-folded plot of the planet's orbit from 300 random samples from the model fit,
       along with the best-fit model.

    Parameters
    ----------
    plot_settings : dict
        A dictionary containing plot settings.
    outfile : str
        The filename to save the plot.
    suffix : str, optional
        An option to append a string to the end of the output files to differentiate fits.
    model : str, optional
        The chosen RV model ('constant', 'decay', or 'precession'), default is 'constant'.

    Returns
    -------
    None

    """
    CONSTANT = False
    DECAY = False
    PRECESSION = False

    TWOPI = 2 * np.pi

    print('-' * 100)
    print('Generating RV plot...')
    print('-' * 100)

    # load plot settings
    plt.rcParams.update(plot_settings['RV_PLOT']['rcParams'])

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

    # load full dataset
    data_file = plot_settings['RV_PLOT']['data_file'+suffix]
    data = utl.read_rv_data(data_file)

    # start subplots
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1], width_ratios=[4, 5], figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax3 = fig.add_subplot(gs[:, 1])

    if model == 'constant':
        CONSTANT = True
        results_file = plot_settings['RV_PLOT']['rv_constant_results_file'+suffix]
        samples_file = plot_settings['RV_PLOT']['rv_constant_samples_file'+suffix]

    elif model == 'decay':
        DECAY = True
        results_file = plot_settings['RV_PLOT']['rv_decay_results_file'+suffix]
        samples_file = plot_settings['RV_PLOT']['rv_decay_samples_file'+suffix]

    elif model == 'precession':
        PRECESSION = True
        results_file = plot_settings['RV_PLOT']['rv_precession_results_file'+suffix]
        samples_file = plot_settings['RV_PLOT']['rv_precession_samples_file'+suffix]

    else:
        raise ValueError('Invalid RV model, must be \'constant\', \'decay\', or \'precession\'.')

    # load fit results
    with open(results_file) as jf:
        rf = json.load(jf)
        res = rf['params']

    try:
        res['e0'] = res['e_derived']
        res['w0'] = res['w_derived']
    except KeyError:
        pass

    # load random samples
    s_orb, s_tdp, s_rv = read_random_samples(samples_file)

    # define time arrays
    times_observed = np.array(data['trv_all'])
    times_all = np.arange(min(times_observed), max(times_observed), 0.01 * res['P0'][0])
    times_fold_min = np.arange(min(times_observed), min(times_observed) + res['P0'][0], 0.01)
    times_fold_max = np.arange(max(times_observed) - res['P0'][0], max(times_observed), 0.01)

    # plot the best-fit radial velocity model
    if plot_settings["RV_PLOT"]['show_RV_curve'] == 'True':

        # compute the model with the best-fit parameters (not including systemic velocity 'v0')
        if CONSTANT:
            mod_all = rv.rv_constant(t0=res['t0'][0], P0=res['P0'][0], e0=res['e0'][0],
                                     w0=res['w0'][0], K=res['K'][0], v0=0.0, dvdt=res['dvdt'][0],
                                     ddvdt=res['ddvdt'][0], t=times_all)

        elif DECAY:
            mod_all = rv.rv_decay(t0=res['t0'][0], P0=res['P0'][0], e0=res['e0'][0],
                                  w0=res['w0'][0], K=res['K'][0], v0=0.0, dvdt=res['dvdt'][0],
                                  ddvdt=res['ddvdt'][0], PdE=res['PdE'][0], t=times_all)

        elif PRECESSION:
            mod_all = rv.rv_precession(t0=res['t0'][0], P0=res['P0'][0], e0=res['e0'][0],
                                       w0=res['w0'][0], K=res['K'][0], v0=0.0, dvdt=res['dvdt'][0],
                                       ddvdt=res['ddvdt'][0], wdE=res['wdE'][0], t=times_all)

        # plot the best-fit model for all times
        ax1.plot(times_all - res['t0'][0], mod_all, label='Best-Fit Model',
                 color='darkgrey', alpha=m_alpha, linewidth=m_lw, zorder=0)

    # plot data, separating instruments into different colours and labels
    for i in data['src_order']:

        # shift the RV measurement times by t0
        x = data['trv'][i] - res['t0'][0]

        # plot the RV data, adjusted by the best-fit systemic velocity 'v0'
        y = data['rvs'][i] - res['v0_' + data['src_tags'][i]][0]

        ax1.errorbar(x, y, yerr=data['err'][i], label=data['src_names'][i], color=data_colors[i],
                     fmt=dfmt, markersize=dms, elinewidth=delw, capsize=decap, zorder=1)

        # generate a best-fit model for each time step
        if CONSTANT:
            mod = rv.rv_constant(t0=res['t0'][0], P0=res['P0'][0], e0=res['e0'][0],
                                 w0=res['w0'][0], K=res['K'][0],
                                 v0=res['v0_' + data['src_tags'][i]][0],
                                 dvdt=res['dvdt'][0], ddvdt=res['ddvdt'][0], t=data['trv'][i])

        elif DECAY:
            mod = rv.rv_decay(t0=res['t0'][0], P0=res['P0'][0], e0=res['e0'][0], w0=res['w0'][0],
                              K=res['K'][0], v0=res['v0_' + data['src_tags'][i]][0],
                              dvdt=res['dvdt'][0], ddvdt=res['ddvdt'][0],
                              PdE=res['PdE'][0], t=data['trv'][i])

        elif PRECESSION:
            mod = rv.rv_precession(t0=res['t0'][0], P0=res['P0'][0],
                                   e0=res['e0'][0], w0=res['w0'][0], K=res['K'][0],
                                   v0=res['v0_' + data['src_tags'][i]][0], dvdt=res['dvdt'][0],
                                   ddvdt=res['ddvdt'][0], wdE=res['wdE'][0], t=data['trv'][i])

        # calculate the residual (data - model)
        residuals = data['rvs'][i] - mod

        # plot the residuals
        ax2.errorbar(x, residuals, yerr=data['err'][i], label=data['src_names'][i],
                     color=data_colors[i], fmt=dfmt, elinewidth=delw, markersize=dms,
                     capsize=decap, capthick=1)

        # phase-fold the data (subtracting long-term trends)
        y_fold = data['rvs'][i] - \
                 res['v0_' + data['src_tags'][i]][0] - \
                 res['dvdt'][0] * (data['trv'][i] - res['t0'][0]) - \
                 0.5 * res['ddvdt'][0] * (data['trv'][i] - res['t0'][0]) ** 2

        if CONSTANT:
            x_fold = []
            for t in data['trv'][i]:
                E = int((t - res['t0'][0]) / res['P0'][0])
                nu = TWOPI / res['P0'][0]
                f_tra = (np.pi / 2 - res['w0'][0]) % TWOPI
                E_tra = (2 * np.arctan(
                    np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_tra / 2))) % TWOPI
                M_tra = E_tra - res['e0'][0] * np.sin(E_tra)
                t_tra = res['t0'][0] + res['P0'][0] * E
                t_p = t_tra - (1 / nu) * M_tra
                MA = TWOPI / res['P0'][0] * (t - t_p)

                # save the phase
                x_fold.append(MA % TWOPI)

        elif DECAY:
            x_fold = []
            for t in data['trv'][i]:
                E = int((t - res['t0'][0]) / res['P0'][0])
                P_anom = res['P0'][0] + res['PdE'][0] * E
                nu = TWOPI / P_anom
                f_tra = (np.pi / 2 - res['w0'][0]) % TWOPI
                E_tra = (2 * np.arctan(
                    np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_tra / 2))) % TWOPI
                M_tra = E_tra - res['e0'][0] * np.sin(E_tra)
                t_tra = res['t0'][0] + res['P0'][0] * E + 0.5 * res['PdE'][0] * E ** 2
                t_p = t_tra - (1 / nu) * M_tra
                MA = TWOPI / P_anom * (t - t_p)

                # save the phase
                x_fold.append(MA % TWOPI)

        elif PRECESSION:

            x_fold = []
            for t in data['trv'][i]:

                E = int((t - res['t0'][0]) / res['P0'][0])
                P_anom = res['P0'][0] / (1 - res['wdE'][0] / (2 * np.pi))
                nu = TWOPI / P_anom
                w_p = (res['w0'][0] + res['wdE'][0] * E) % TWOPI

                f_tra = (np.pi / 2 - w_p) % TWOPI
                E_tra = (2 * np.arctan(
                    np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_tra / 2))) % TWOPI
                M_tra = E_tra - res['e0'][0] * np.sin(E_tra)
                t_tra = res['t0'][0] + res['P0'][0] * E - \
                        (res['e0'][0] * P_anom / np.pi) * np.cos(w_p)

                t_p = t_tra - (1 / nu) * M_tra
                MA = TWOPI / P_anom * (t - t_p)

                # save the phase
                x_fold.append(MA % TWOPI)

        # plot the phase-folded data
        ax3.errorbar(x_fold, y_fold, yerr=data['err'][i], label=data['src_names'][i], fmt='.',
                     color=data_colors[i], markersize=dms, elinewidth=delw, capsize=decap)

    # repeat the above for the 300 random posterior samples
    for i in range(np.shape(s_orb)[0]):
        times_s = np.arange(s_orb[i][0], s_orb[i][0] + s_orb[i][1], 0.01)

        if CONSTANT:
            x_fold_s = []
            for t in times_s:
                E = int((t - s_orb[i][0]) / s_orb[i][1])
                nu = TWOPI / s_orb[i][1]
                f_tra = (np.pi / 2 - s_orb[i][3]) % TWOPI
                E_tra = (2 * np.arctan(
                    np.sqrt((1 - s_orb[i][2]) / (1 + s_orb[i][2])) * np.tan(f_tra / 2))) % TWOPI
                M_tra = E_tra - s_orb[i][2] * np.sin(E_tra)
                t_tra = s_orb[i][0] + s_orb[i][1] * E
                t_p = t_tra - (1 / nu) * M_tra
                MA = TWOPI / s_orb[i][1] * (t - t_p)

                # save the phase
                x_fold_s.append(MA % TWOPI)

            y_fold_s = rv.rv_constant(t0=s_orb[i][0], P0=s_orb[i][1],
                                      e0=s_orb[i][2], w0=s_orb[i][3], K=s_rv[i][0],
                                      v0=0.0, dvdt=0.0, ddvdt=0.0, t=times_s)

            # plot the phase-folded model
            y_fold_s = [y for _, y in sorted(zip(x_fold_s, y_fold_s))]
            x_fold_s = sorted(x_fold_s)

            ax3.plot(x_fold_s, y_fold_s, label='_',
                     color='lightgrey', linestyle='-', linewidth=s_lw, alpha=s_alpha)

        elif DECAY:
            times_s = np.arange(s_orb[i][0], s_orb[i][0] + s_orb[i][1], 0.01)

            x_fold_s = []
            for t in times_s:
                E = int((t - s_orb[i][0]) / s_orb[i][1])
                P_anom = s_orb[i][1] + s_tdp[i][0] * E
                nu = TWOPI / P_anom
                f_tra = (np.pi / 2 - s_orb[i][3]) % TWOPI
                E_tra = (2 * np.arctan(
                    np.sqrt((1 - s_orb[i][2]) / (1 + s_orb[i][2])) * np.tan(f_tra / 2))) % TWOPI
                M_tra = E_tra - s_orb[i][2] * np.sin(E_tra)
                t_tra = s_orb[i][0] + s_orb[i][1] * E + 0.5 * s_tdp[i][0] * E ** 2
                t_p = t_tra - (1 / nu) * M_tra
                MA = TWOPI / P_anom * (t - t_p)

                # save the phase
                x_fold_s.append(MA % TWOPI)

            y_fold_s = rv.rv_decay(t0=s_orb[i][0], P0=s_orb[i][1], e0=s_orb[i][2],
                                   w0=s_orb[i][3], K=s_rv[i][0], v0=0.0, dvdt=0.0, ddvdt=0.0,
                                   PdE=s_tdp[i][0], t=times_s)

            # plot the phase-folded model
            y_fold_s = [y for _, y in sorted(zip(x_fold_s, y_fold_s))]
            x_fold_s = sorted(x_fold_s)

            ax3.plot(x_fold_s, y_fold_s, label='_',
                     color='lightgrey', linestyle='-', linewidth=s_lw, alpha=s_alpha)

    # phase-fold the best-fit model (for planet only)
    if CONSTANT:

        times_fold = np.arange(res['t0'][0], res['t0'][0] + res['P0'][0], 0.01)

        x_fold_m = []
        for t in times_fold:
            E = int((t - res['t0'][0]) / res['P0'][0])
            nu = TWOPI / res['P0'][0]
            f_tra = (np.pi / 2 - res['w0'][0]) % TWOPI
            E_tra = (2 * np.arctan(
                np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_tra / 2))) % TWOPI
            M_tra = E_tra - res['e0'][0] * np.sin(E_tra)
            t_tra = res['t0'][0] + res['P0'][0] * E
            t_p = t_tra - (1 / nu) * M_tra
            MA = TWOPI / res['P0'][0] * (t - t_p)

            # save the phase
            x_fold_m.append(MA % TWOPI)

        y_fold_m = rv.rv_constant(t0=res['t0'][0], P0=res['P0'][0],
                                  e0=res['e0'][0], w0=res['w0'][0], K=res['K'][0],
                                  v0=0.0, dvdt=0.0, ddvdt=0.0, t=times_fold)

        # plot the phase-folded model(s)
        y_fold_m = [y for _, y in sorted(zip(x_fold_m, y_fold_m))]
        x_fold_m = sorted(x_fold_m)

        ax3.plot(x_fold_m, y_fold_m, label='RV Model', color='dimgrey', linewidth=m_lw, alpha=0.75)

    elif DECAY:
        times_fold = np.arange(res['t0'][0], res['t0'][0] + res['P0'][0], 0.01)

        x_fold_m = []
        for t in times_fold:
            E = int((t - res['t0'][0]) / res['P0'][0])
            P_anom = res['P0'][0] + res['PdE'][0] * E
            nu = TWOPI / P_anom
            f_tra = (np.pi / 2 - res['w0'][0]) % TWOPI
            E_tra = (2 * np.arctan(
                np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_tra / 2))) % TWOPI
            M_tra = E_tra - res['e0'][0] * np.sin(E_tra)
            t_tra = res['t0'][0] + res['P0'][0] * E + 0.5 * res['PdE'][0] * E ** 2
            t_p = t_tra - (1 / nu) * M_tra
            MA = TWOPI / P_anom * (t - t_p)

            # save the phase
            x_fold_m.append(MA % TWOPI)

        # plot the phase-folded model
        y_fold_m = rv.rv_decay(t0=res['t0'][0], P0=res['P0'][0], e0=res['e0'][0],
                               w0=res['w0'][0], K=res['K'][0], v0=0.0, dvdt=0.0, ddvdt=0.0,
                               PdE=res['PdE'][0], t=times_fold)

        y_fold_m = [y for _, y in sorted(zip(x_fold_m, y_fold_m))]
        x_fold_m = sorted(x_fold_m)

        ax3.plot(x_fold_m, y_fold_m, label='RV Model',
                 color='dimgrey', linestyle='-', linewidth=1, alpha=0.75)

    elif PRECESSION:

        # calculate the anomalistic period
        P_anom = res['P0'][0] / (1 - res['wdE'][0] / (2 * np.pi))

        # plot the phase-folded model over all time
        x_fold_all = []
        for t in times_all:
            E = int((t - res['t0'][0]) / P_anom)
            nu = TWOPI / P_anom
            w_p = (res['w0'][0] + res['wdE'][0] * E) % TWOPI

            f_tra = (np.pi / 2 - w_p) % TWOPI
            E_tra = (2 * np.arctan(
                np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_tra / 2))) % TWOPI
            M_tra = E_tra - res['e0'][0] * np.sin(E_tra)
            t_tra = res['t0'][0] + res['P0'][0] * E  - \
                    (res['e0'][0] * P_anom / np.pi) * np.cos(w_p)

            t_p = t_tra - (1 / nu) * M_tra
            MA = TWOPI / P_anom * (t - t_p)

            # save the phase
            x_fold_all.append(MA % TWOPI)

        # plot the phase-folded model
        y_fold_all = rv.rv_precession(t0=res['t0'][0], P0=res['P0'][0],
                                      e0=res['e0'][0], w0=res['w0'][0], K=res['K'][0], v0=0.0,
                                      dvdt=0.0, ddvdt=0.0, wdE=res['wdE'][0], t=times_all)

        y_fold_all = [y for _, y in sorted(zip(x_fold_all, y_fold_all))]
        x_fold_all = sorted(x_fold_all)

        ax3.plot(x_fold_all, y_fold_all, label='Apsidal Precession Model',
                 color='dimgrey', linewidth=m_lw, alpha=0.06)

        # plot the phase-folded model at t_min
        x_fold_min = []
        for t in times_fold_min:
            E = int((t - res['t0'][0]) / P_anom)
            nu = TWOPI / P_anom
            w_p = (res['w0'][0] + res['wdE'][0] * E) % TWOPI

            f_tra = (np.pi / 2 - w_p) % TWOPI
            E_tra = (2 * np.arctan(
                np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_tra / 2))) % TWOPI
            M_tra = E_tra - res['e0'][0] * np.sin(E_tra)
            t_tra = res['t0'][0] + res['P0'][0] * E \
                    - (res['e0'][0] * P_anom / np.pi) * np.cos(w_p)

            t_p = t_tra - (1 / nu) * M_tra
            MA = TWOPI / P_anom * (t - t_p)

            # save the phase
            x_fold_min.append(MA % TWOPI)

        y_fold_min = rv.rv_precession(t0=res['t0'][0], P0=res['P0'][0],
                                      e0=res['e0'][0], w0=res['w0'][0], K=res['K'][0], v0=0.0,
                                      dvdt=0.0, ddvdt=0.0, wdE=res['wdE'][0], t=times_fold_min)

        # plot the phase-folded model
        y_fold_min = [y for _, y in sorted(zip(x_fold_min, y_fold_min))]
        x_fold_min = sorted(x_fold_min)

        tmin = Time([min(times_observed)], format='jd', scale='tdb')
        label_min = 'Orbit on ' + (tmin.to_value('iso', subfmt='date')[0]) + r' ($t_{\rm min}$)'

        ax3.plot(x_fold_min, y_fold_min, label=label_min,
                 color='dimgrey', linestyle='-.', linewidth=1)

        # plot the phase-folded model at t_max
        x_fold_max = []
        for t in times_fold_max:
            E = int((t - res['t0'][0]) / P_anom)
            nu = TWOPI / P_anom
            w_p = (res['w0'][0] + res['wdE'][0] * E) % TWOPI

            f_tra = (np.pi / 2 - w_p) % TWOPI
            E_tra = (2 * np.arctan(
                np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_tra / 2))) % TWOPI
            M_tra = E_tra - res['e0'][0] * np.sin(E_tra)
            t_tra = res['t0'][0] + res['P0'][0] * E \
                    - (res['e0'][0] * P_anom / np.pi) * np.cos(w_p)

            t_p = t_tra - (1 / nu) * M_tra
            MA = TWOPI / P_anom * (t - t_p)

            # save the phase
            x_fold_max.append(MA % TWOPI)

        y_fold_max = rv.rv_precession(t0=res['t0'][0], P0=res['P0'][0],
                                      e0=res['e0'][0], w0=res['w0'][0], K=res['K'][0], v0=0.0,
                                      dvdt=0.0, ddvdt=0.0, wdE=res['wdE'][0], t=times_fold_max)

        # plot the phase-folded model
        y_fold_max = [y for _, y in sorted(zip(x_fold_max, y_fold_max))]
        x_fold_max = sorted(x_fold_max)

        tmax = Time([max(times_observed)], format='jd', scale='tdb')
        label_max = 'Orbit on ' + (tmax.to_value('iso', subfmt='date')[0]) + r' ($t_{\rm max}$)'

        ax3.plot(x_fold_max, y_fold_max, label=label_max, color='k', linestyle='-',  linewidth=1)

    # plot zero lines for reference
    ax1.axhline(y=0, linestyle='-', color='grey', linewidth=1)
    ax2.axhline(y=0, linestyle='-', color='grey', linewidth=1)
    ax3.axhline(y=0, linestyle='-', color='grey', linewidth=0.5)
    ax1.axvline(x=0, linestyle='-', color='grey', linewidth=0.5)
    ax2.axvline(x=0, linestyle='-', color='grey', linewidth=0.5)

    # calculate the true anomaly, eccentric anomaly, and mean anomaly at t0
    f_t0 = (np.pi / 2 - res['w0'][0]) % (2 * np.pi)
    E_t0 = (2 * np.arctan(np.sqrt((1 - res['e0'][0]) / (1 + res['e0'][0])) * np.tan(f_t0 / 2))) % (2 * np.pi)
    M_t0 = E_t0 - res['e0'][0] * np.sin(E_t0)
    ax3.axvline(x=M_t0, linestyle='--', color='dimgrey', linewidth=1)
    ax3.text(M_t0-0.15, min(y_fold), 'transits', fontsize=11, rotation=90)

    # add text for t0 date
    t = Time([str(res['t0'][0])], format='jd', scale='tdb')
    ax2.text(.99, .01, 't$_0$ = ' + t.to_value('iso', subfmt='date')[0],
             ha='right', va='bottom', transform=ax2.transAxes, fontsize=11)

    # finish plots
    labels = ['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$',
              r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$', r'$2\pi$']
    ax3.set_xticks(np.arange(0, 2 * np.pi + 0.01, np.pi / 4))
    ax3.set_xticklabels(labels)
    ax3.set_xlim(0, 2*np.pi)

    plt.suptitle(plot_settings['RV_PLOT']['title'], fontsize=18)
    ax1.set_title('All Time')
    ax1.set_ylabel('Radial Velocity (m/s)')
    ax2.set_ylabel('Residuals (m/s)')
    ax2.set_xlabel('BJD $-$ t$_0$ (days)')
    ax3.set_ylabel('Radial Velocity (m/s)')
    ax3.set_xlabel('Mean Anomaly')
    ax3.set_title('Phase Folded')
    ax3.legend()
    plt.savefig(outfile, bbox_inches=bbox, dpi=dpi, pad_inches=pad_inches)
    plt.close()


def corner_plot(dic, samples, params, outfile):
    """Generate a corner plot from the weighted posterior samples.

    This method generates a corner plot from the weighted posterior samples using the
    corner.py package written by Daniel Foreman-Mackey [1]_.

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
    # specify number of dimensions
    ndim = len(params)

    # update plot settings
    plot_dic = {
          "figure.figsize": [10, 10], "font.family": "serif",
          "xtick.direction": "in", "ytick.direction": "in",
          "xtick.labelsize": 11, "ytick.labelsize": 11,
          "axes.labelsize": 14, "axes.titlesize": 14}
    plt.rcParams.update(plot_dic)

    # generate corner plot
    figure = corner.corner(samples, labels=params, color='k', top_ticks=False,
                           show_titles=False, plot_contours=True,
                           max_n_ticks=5, title_fmt='.2E', labelpad=0.2)

    axes = np.array(figure.axes).reshape((ndim, ndim))

    for ax in figure.get_axes():
        ax.tick_params(axis='both', labelsize=12)

    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(dic[params[i]][0], color='darkred')

    plt.savefig(outfile, bbox_inches='tight', pad_inches=0.5, dpi=300)
    plt.close()

    return


def read_random_samples(data_file, delim='\t'):
    """Read files with sets of random posterior samples into the format needed for plotting.

    Parameters
    ----------
    data_file : str
        Name of the file containing the posterior samples.
    delim : str, optional
        The data file delimiter, default is '\t' for a tab-separated file.

    Returns
    -------
    None

    """
    # data holders
    orbital_elements = []
    time_dependent = []
    radial_velocity = []

    # read samples from file
    with open(data_file, 'r') as file:
        reader = csv.reader(file, delimiter=delim)

        # skip header
        next(reader)

        # parse each row of the CSV file
        for row in reader:
            orb = [float(x) for x in row[0:6]]
            tdp = [float(x) for x in row[6:10]]
            rdv = [float(row[11])]

            v0 = list(map(str.strip, row[12].strip('][').replace('"', '').split(',')))
            v0 = [float(x) for x in v0]
            jit = list(map(str.strip, row[12].strip('][').replace('"', '').split(',')))
            jit = [float(x) for x in jit]

            rdv.append(v0)
            rdv.append(jit)
            rdv.append(float(row[14]))
            rdv.append(float(row[15]))

            orbital_elements.append(orb)   # t0 P0 e0 w0 i0 O0
            time_dependent.append(tdp)     # PdE wdE edE idE OdE
            radial_velocity.append(rdv)    # K v0 jit dvdt ddvdt

    # make arrays
    orbital_elements = np.array(orbital_elements)
    time_dependent = np.array(time_dependent)
    radial_velocity = np.array(radial_velocity, dtype=object)

    return orbital_elements, time_dependent, radial_velocity


def periodogram(data_file, outfile, min_period=3, max_period=100000.0, num=1000000):
    """Generate a periodogram from the residuals of an RV model fit.

    Notes
    -----
    The making of this plot is not currently implemented in the model fit routines, and this
    function should be considered to be under development. It is provided here for convenience.

    Parameters
    ----------
    data_file : str
        Name of the file containing the residuals.
    outfile : str
        File name of the saved plot.
    min_period : float
        The minimum period in days to include in the periodogram.
    max_period : float
        The maximum period in days to include in the periodogram.
    num : int
        The number of periods/frequencies in between 'min_period' and 'max_period'.

    Returns
    -------
    None

    """

    print('-' * 100)
    print('Generating periodogram of RV residuals...')
    print('-' * 100)

    # load full dataset
    data = utl.read_rv_data(data_file)

    times = data['trv_all']
    velocities = data['rvs_all']
    errors = data['err_all']

    frequency = np.linspace(1 / max_period, 1 / min_period, num)
    power = LombScargle(times, velocities, errors).power(frequency)
    mask = sci.find_peaks(power, distance=100000, prominence=0.4)[0]

    best_freq = np.array([f for _, f in sorted(zip(power[mask], frequency[mask]), reverse=True)])
    best_power = sorted(power[mask], reverse=True)
    best_period = 1 / best_freq

    figure_params = {'figure.figsize': [10, 7], 'font.family': 'serif',
                     'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.labelsize': 12, 'ytick.labelsize': 12,
                     'axes.labelsize': 13, 'axes.titlesize': 16,
                     'legend.fontsize': 12, 'legend.labelspacing': 0.9,
                     'legend.framealpha': 1., 'legend.borderpad': 0.7,
                     'legend.loc': 'upper right'}

    plt.rcParams.update(figure_params)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()

    ax1.plot(frequency, power, color='cornflowerblue')
    ax1.scatter(frequency[mask], power[mask], color='firebrick')

    print('best orbital periods:')
    for i, p in enumerate(best_period):
        print("     %.1f d, power = %.3f" % (p, best_power[i]))
        plt.text(best_freq[i] + 0.01 * best_freq[i],
                 best_power[i] + 0.03 * best_power[i],
                 '%.0f days' % best_period[i], fontsize=12)

    ax1.set_ylabel('Power')
    ax1.set_xlabel(r'Frequency [days$^{-1}$]')
    ax2.set_xlabel('Period [days]')
    ax1.set_ylim(0, best_power[0] + 0.1)
    new_tick_locations = np.linspace(1 / max_period, 1 / min_period, num=10)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(["%.0f" % z for z in 1 / new_tick_locations])
    plt.savefig(outfile, bbox_inches='tight', dpi=300, pad_inches=0.25)
    plt.close()
