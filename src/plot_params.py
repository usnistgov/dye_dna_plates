from .util import figure_name_to_abspath
from .parameter_extraction import Parameters
from .get_data import CombinedData
import numpy as np
import matplotlib.pyplot as plt


def plot_figure7(
    T_SS, f_SS, K_SS, f_std_SS, K_std_SS,
    T_DS, f_DS, K_DS, f_std_DS, K_std_DS):
    ds_color = "tab:purple"
    ss_color = "tab:olive"

    fig, (ax_K, ax_f) = plt.subplots(figsize=(3.25, 5.), ncols=1, nrows=2, sharex=True)
    ax_K.annotate("(a)", xy=(-0.18, 0.95), xycoords='axes fraction')
    ax_f.annotate("(b)", xy=(-0.18, 0.95), xycoords='axes fraction')

    ysps_kwargs = dict(mfc='None', markersize=2)
    ax_K.plot(T_SS, K_SS, '-', color=ss_color, label="SS", **ysps_kwargs)
    ax_K.fill_between(T_SS, K_SS - 3*K_std_SS, K_SS + 3*K_std_SS, color=ss_color, alpha=0.3, clip_on=False)
    ax_f.plot(T_SS, f_SS, '-', color=ss_color, label=r"SS", **ysps_kwargs)
    ax_f.fill_between(T_SS, f_SS - 3*f_std_SS, f_SS + 3*f_std_SS, color=ss_color, label=r"SS", alpha=0.3, clip_on=False)

    ax_K.plot(T_DS, K_DS, '-', color=ds_color, label="DS", **ysps_kwargs)
    ax_K.fill_between(T_DS, K_DS - 3*K_std_DS, K_DS + 3*K_std_DS, color=ds_color, alpha=0.3, clip_on=False)
    ax_f.plot(T_DS, f_DS, '-', color=ds_color, label=r"DS", **ysps_kwargs)
    ax_f.fill_between(T_DS, f_DS - 3*f_std_DS, f_DS + 3*f_std_DS, color=ds_color, alpha=0.3, clip_on=False)


    for ax in (ax_f, ax_K):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(which='both', direction="in")

    ax_f.set_ylabel("$\\mathbf{f}_j$", rotation=0., labelpad=12)
    ax_K.set_ylabel("$\\mathbf{K}_j$", rotation=0., labelpad=12)
    fig.subplots_adjust(right=0.99, top=0.99, bottom=0.09)
    leg_kwargs = dict(loc=(0.15, 0.85), framealpha=0.)
    leg_kwargs.pop('framealpha')

    ax_K.legend(framealpha=0., loc=(0.45, 0.8))
    # axes_yP[1].legend(framealpha=0., loc=(0.05, 0.05))
    ax_K.set_yticks([0., 1., 2., 3.])
    ax_K.set_yticks([
        0.2, 0.4, 0.6, 0.8,
        1.2, 1.4, 1.6, 1.8,
        2.2, 2.4, 2.6, 2.8
        ], minor=True)
    ax_K.set_ylim([0., 2.8])
    ax_f.set_yticks([1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.2])
    ax_f.set_yticks([0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1], minor=True)
    ax_f.set_ylim([0.9, 2.1])
    ax_f.set_xlabel(r"$T_j$ [K]")

    fig.subplots_adjust(left=0.17, hspace=0.06, bottom=0.08)
    fig.savefig(figure_name_to_abspath("figure7.png"), transparent=True, dpi=300)


def plot_figure8(SS: Parameters, DS: Parameters):
    fig, ax = plt.subplots(ncols=1, nrows=1, sharex=True, sharey=True, 
                                         figsize=(3.25, 3.25))
    dg_SS = SS.get_dg()
    dh_SS = SS.get_dh()
    dg_std_SS = SS.get_dg_std()
    dh_std_SS = SS.get_dh_std()
    ax.plot(
        SS.T, dg_SS, '-', color='C0', label="$\\Delta g_j$,SS", clip_on=False,
    )
    ax.fill_between(SS.T, dg_SS - 3.*dg_std_SS, dg_SS + 3.*dg_std_SS, color='C0', alpha=0.3, clip_on=False)

    ax.plot(
        SS.T[1:-1], dh_SS, 'x', color='black', markersize=4.,
        mfc='None', label="$\\Delta h_j$,SS", clip_on=True ## NB!
    )
    # axes.fill_between(
    #     SS.T[1:-1], dh_SS - dh_std_SS, dh_SS + dh_std_SS, color='C1', alpha=0.3, clip_on=True ## NB!
    # )


    dg_DS = DS.get_dg()
    dh_DS = DS.get_dh()
    dg_std_DS = DS.get_dg_std() #
    dh_std_DS = DS.get_dh_std()
    ax.plot(
        DS.T, dg_DS, '-', color='C2', label="$\\Delta g_j$,DS", clip_on=False,
    )
    ax.fill_between(DS.T, dg_DS - 3.*dg_std_DS, dg_DS + 3.*dg_std_DS, color='C2', alpha=0.3, clip_on=False)
    i_T = -24
    print("dg_SS at %3.2f K is %f +/- %f" % (SS.T[i_T], dg_SS[i_T], 3.*dg_std_SS[i_T]))
    print("dg_DS at %3.2f K is %f +/- %f" % (DS.T[i_T], dg_DS[i_T], 3.*dg_std_DS[i_T]))
    ax.plot(
        DS.T[1:-1], dh_DS, '.', color='C3', mfc='None', label="$\\Delta h_j$,DS", clip_on=False,
    )

    ax.legend(
        borderpad=0.1, handlelength=1., labelspacing=0.1,
        facecolor='None', edgecolor='black',
        loc=(0.33, 0.88), ncol=2)
    ax.set_ylabel("$\\Delta g_j$ or $\\Delta h_j$ [kJ/mol]", labelpad=0)
    ax.tick_params(which='both', direction='in', right=True)
    ax.set_xlabel("$T_j$ [K]")
    ax.set_yticks([-60., -50, -40., -30., -20, -10.])
    ax.set_yticks([
            -58, -56, -54, -52,
            -48, -46, -44, -42,
            -38, -36, -34, -32,
            -28, -26, -24, -22,
            -18, -16, -14, -12,
    ], minor=True)
    ax.set_ylim([-50., -14])
    ax.set_xticks([285., 300, 315])
    ax.set_xticks([277.5, 292.5, 307.5, 322.5], minor=True)
    ax.spines['top'].set_visible(False)
    ax.set_xlim([DS.T.min(), DS.T.max()])

    fig.subplots_adjust(right=0.99, left=0.16, bottom=0.13, top=0.98, hspace=0.0, wspace=0.0)
    fig.savefig(figure_name_to_abspath("figure8.png"), transparent=True, dpi=300)


def plot_figure_S3(SS_1: CombinedData, SS_2: CombinedData, DS_1: CombinedData, DS_2: CombinedData):
    # plot \|F - Fhat\| and \|C - Chat\| vs rho
    from .noise_removal import predictor_corrector
    fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True, figsize=(5., 5.))
    num_rhos = 10
    rhos = np.logspace(-2, 4, num_rhos)
    for (cls, ax) in [
        (SS_1, axes[0, 0]),
        (SS_2, axes[0, 1]),
        (DS_1, axes[1, 0]),
        (DS_2, axes[1, 1])
        ]:
        name = "%s,\n$\\mathbf{D}=%i$" % (cls.t, int(cls.D))
        norm_dF = np.zeros(num_rhos)
        norm_dC = np.zeros(num_rhos)
        for l in range(num_rhos):
            Mtls, Chat = predictor_corrector(cls.F, cls.C, rhos[l], maxiter=100000, print_iter=False)
            Fhat = np.outer(Mtls, Chat)
            dF = cls.F - Fhat
            dC = cls.C - Chat
            norm_dF[l] = np.sum(dF*dF)
            norm_dC[l] = np.sum(dC*dC)
    
        kwargs = dict(mfc='None', markersize=6)
        ax.plot(rhos, norm_dF, '-o', label=r"$\|\| \mathbf{F} - \widehat{\mathbf{F}}^\mathrm{TLS}\|\|_F^2$", **kwargs)
        ax.plot(rhos, norm_dC, '-o', label=r"$\|\|\mathbf{C} - \widehat{\mathbf{C}}\|\|_2^2$", **kwargs)
        ax.semilogx()
        ax.annotate(name, xy=(0.65, 0.3), xycoords='axes fraction')

    axes[1, 0].legend(framealpha=0., loc=(0.01, 0.55))
    axes[1, 0].set_xlabel(r"$\rho^2$")
    axes[1, 1].set_xlabel(r"$\rho^2$")
    fig.subplots_adjust(right=0.99, left=0.08, hspace=0.05, wspace=0.2, top=0.99, bottom=0.1)
    fig.savefig(figure_name_to_abspath("figureS3.png"), transparent=True, dpi=300)


def plot_figure_S4(SS: Parameters, DS: Parameters):
    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    def plot_bf(T, phi, dphi, color, label):
        ax.plot(T, phi, color=color, label=label)
        ax.fill_between(T, phi - 3*dphi, phi + 3*dphi, color=color, alpha=0.5)
    
    plot_bf(SS.T, SS.get_phi_1(), SS.get_std_phi_1(), "C0", "SS, $\\mathbf{D}=1$")
    plot_bf(SS.T, SS.get_phi_2(), SS.get_std_phi_2(), "C1", "SS, $\\mathbf{D}=2$")

    ax.errorbar(DS.T, DS.get_phi_1(), 
                                 fmt='d', yerr=3*DS.get_std_phi_1(), 
                                 color="C0", label="DS, $\\mathbf{D}=1$", mfc='None')
    ax.errorbar(DS.T, DS.get_phi_2(), fmt='d', 
                                 yerr=3*DS.get_std_phi_2(), color="C1", label="DS, $\\mathbf{D}=2$", mfc='None')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlabel("$T_j$ [K]")
    ax.set_ylabel("$\\varphi_{j\\mathbf{D}}$", rotation=0., labelpad=12)
    ax.legend(edgecolor='None', facecolor='None', ncol=2)
    ax.tick_params(direction='in', which='both')
    fig.subplots_adjust(left=0.17, right=0.99, top=0.97, bottom=0.14)
    ax.set_ylim([0., 1.])
    fig.savefig(figure_name_to_abspath("figureS4.png"), dpi=300, transparent=True)


def plot_figure_S5(SS_data: Parameters, DS_data: Parameters):
    # plot 
    fig, ax = plt.subplots(figsize=(5., 5.))
    ax.plot(SS_data.T, SS_data.get_theta_b_all_1(1.5), label="SS")
    ax.plot(DS_data.T, DS_data.get_theta_b_all_1(1.5), label="DS")
    ax.legend(framealpha=0.)
    ax.set_xlabel("$T_j$ [K]")
    ax.set_ylabel("$\\theta_{b,j,1}$", rotation=0, labelpad=12)

    fig.subplots_adjust(right=0.98, bottom=0.10, left=0.15, top=0.99)
    fig.savefig(figure_name_to_abspath("figureS5.png"), dpi=300, transparent=True)


def plot_figure_S6(T, rb, d_rb):
    """Plot figure S5, relative brightness

    Parameters
    ----------
    T : np.array
        array of temperatures
    rb : np.array
        Relative brightness, :math:`\\mathbf{f}^\\mathrm{DS}_j/\\mathbf{f}^\\mathrm{SS}_j`
    d_rb : np.array
        standard deviation in relative brightness
    """
    figbr, axbr = plt.subplots(figsize=(3.25, 3.25), dpi=300)
    axbr.plot(T, rb)
    axbr.fill_between(T, rb - d_rb,  rb + d_rb, 
                              color='C0', alpha=0.3)
    axbr.set_xlabel("$T$ [K]")
    axbr.set_ylabel("Relative Brightness, $f^\\mathrm{DS}/f^\\mathrm{SS}$")
    figbr.subplots_adjust(left=0.23, right=0.99, bottom=0.13, top=0.98)
    axbr.tick_params(which='both', direction='in')
    axbr.spines['top'].set_visible(False)
    axbr.spines['right'].set_visible(False)
    figbr.savefig(
        figure_name_to_abspath("figureS6.png"), 
        dpi=300, transparent=True
    )