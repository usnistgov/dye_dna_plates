import numpy as np
import matplotlib.pyplot as plt
from .util import figure_name_to_abspath


def plot_error_F(ax, F_k_tl, F_hat_k_tl, T):
    """Plot the predicted fluorescence :math:`\\widehat{\\mathbf{F}}_{ijk}^{t\ell}`
    against :math:`\\mathbf{F}_{ijk}^{t\ell}` 
    for :math:`i=1,\dots,n` and :math:`j=1,\dots,n`.

    Parameters
    ----------
    ax : matplotlib.axes
        axes to plot on
    F_kl : np.ndarray
        Experimentally measured fluorescence :math:`\\mathbf{F}_{k}^{t\ell}`.
    F_hat_kl : np.ndarray
        Predicted fluorescence :math:`\\widehat{\\mathbf{F}}_{k}^{t\ell}`.
    T : np.array
        Temperatures in K

    Returns
    -------
    img
        image for colorbar
    """
    m, n = F_k_tl.shape
    for i_T in range(m-1, 0, -1):
        img = ax.scatter(F_k_tl[i_T, :], F_hat_k_tl[i_T, :], s=2, marker='.', c=np.full(n, fill_value=T[i_T]), 
                         vmin = 283., vmax=330.,
                         cmap='coolwarm')
    
    E = F_hat_k_tl - F_k_tl
    
    ax.annotate("$R=%4.f$" % np.sum(E), xy=(0.8, 0.2), xycoords="data", color="purple", fontsize=8)

    ax.set_xticks([0., 1., 2., 3.])
    ax.set_yticks([0., 1., 2., 3.])
    ax.set_xticks([0.5, 1.5, 2.5], minor=True)
    ax.set_yticks([0.5, 1.5, 2.5], minor=True)
    ax.set_xlim([0., 2.5])
    ax.set_ylim([0., 2.5])
    ax.plot([0., 2.5], [0., 2.5], '--', color='black', linewidth=0.3)
    ax.tick_params(which='both', direction='in')
    return img


def plot_error_C(ax, X, x, x_std=None, **kwargs):
    # line, = ax.plot(X, x, '-o', mfc='None', markersize=2, **kwargs)
    line, = ax.plot(X, x, '-', linewidth=0.5, alpha=1, **kwargs)
    if 'color' not in kwargs.keys():
        kwargs['color'] = line.get_color() 
    if 'label' in kwargs.keys():
        kwargs.pop('label')
    if x_std is not None:
        ax.fill_between(X, x + 3*x_std, x - 3*x_std,
                                 color=kwargs['color'], alpha=0.5)
    ax.plot(X[::12], x[::12], '*', markersize=4, **kwargs)
    ax.set_xticks([0., 1., 2., 3., 4., 5.])
    ax.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], minor=True)
    ax.set_yticks([0., 1., 2., 3., 4., 5.])
    ax.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], minor=True)
    ax.set_xlim([0., 4.5])
    ax.set_ylim([0., 4.5])
    ax.plot([0., 4.5], [0., 4.5], '--', color='black', linewidth=0.6, zorder=-10)
    ax.tick_params(which='both', direction='in', right=True, top=True)


def plot_Chat_vs_C(Cs, C_hats, C_stds, fname):
    figls, axes = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True, figsize=(3.25, 3.25))
    axls = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]

    n = Cs[0].shape[0]
    C, Chat, C_std = Cs[0], C_hats[0], C_stds[0]
    ni = n // 3
    plot_error_C(axls[0], C[:ni], Chat[:ni], C_std[:ni], label="(SS,A,1)", color="C1")
    plot_error_C(axls[0], C[ni:2*ni], Chat[ni:2*ni], C_std[ni:2*ni], label="(SS,B,1)", color="C0")
    plot_error_C(axls[0], C[2*ni:], Chat[2*ni:], C_std[2*ni:], label="(SS,C,1)", color="C2")

    icol = 1
    n = Cs[icol].shape[0]
    C, Chat, C_std = Cs[icol], C_hats[icol], C_stds[icol]
    ni = n // 2
    plot_error_C(axls[icol], C[:ni], Chat[:ni], C_std[:ni], label="(SS,A,2)", color="C1")
    plot_error_C(axls[icol], C[ni:2*ni], Chat[ni:2*ni], C_std[ni:2*ni], label="(SS,B,2)", color="C2")
    icol += 1
    n = Cs[icol].shape[0]
    C, Chat, C_std = Cs[icol], C_hats[icol], C_stds[icol]
    ni = n // 2
    plot_error_C(axls[icol], C[:ni], Chat[:ni], C_std[:ni], label="(DS,A,1)", color="C2")
    plot_error_C(axls[icol], C[ni:2*ni], Chat[ni:2*ni], C_std[ni:2*ni],  label="(DS,B,1)", color="C3")
    icol += 1
    n = Cs[icol].shape[0]
    C, Chat, C_std = Cs[icol], C_hats[icol], C_stds[icol]
    ni = n // 2
    plot_error_C(axls[icol], C[:ni], Chat[:ni], C_std[:ni], label="(DS,A,2)", color="C0")
    plot_error_C(axls[icol], C[ni:2*ni], Chat[ni:2*ni], C_std[ni:2*ni],  label="(DS,B,2)", color="C3")

    axes[1, 1].set_xlabel(r"$\mathbf{C}_{i}$")
    kwargs = dict(borderpad=0., labelspacing=0.1, handlelength=1., handletextpad=0.1, fontsize=9)
    axes[1, 0].set_xlabel(r"$\mathbf{C}_{i}$")

    axes[0, 0].set_ylabel(r"$\widehat{\mathbf{C}}_i$", rotation=0.)
    axes[1, 0].set_ylabel(r"$\widehat{\mathbf{C}}_i$", rotation=0.)

    axes[0, 1].legend(framealpha=0., **kwargs)
    axes[1, 0].legend(framealpha=0., **kwargs)
    axes[1, 1].legend(framealpha=0., **kwargs)
    axes[0, 0].legend(framealpha=0., **kwargs)

    figls.subplots_adjust(left=0.08, right=0.99, top=0.99, hspace=0., bottom=0.12, wspace=0.)
    figls.savefig(figure_name_to_abspath(fname), transparent=True, dpi=300)


def plot_Fhat_vs_F(Fs, Fhats, Ts, fname, sname=r"$\widehat{\mathbf{F}}_{ji}^\mathrm{LS}$"):
    figls, axls = plt.subplots(ncols=3, nrows=4, sharex=True, sharey=True, figsize=(3.25, 5.))
    F_SS_1, F_SS_2, F_DS_1, F_DS_2 = Fs
    Fhat_SS_1, Fhat_SS_2, Fhat_DS_1, Fhat_DS_2 = Fhats
    T_SS_1, T_SS_2, T_DS_1, T_DS_2 = Ts

    m, n = F_SS_1.shape
    ni = n // 3
    plot_error_F(axls[0, 0], F_SS_1[:, :ni], Fhat_SS_1[:, :ni], T_SS_1)
    plot_error_F(axls[0, 1], F_SS_1[:, ni:2*ni], Fhat_SS_1[:, ni:2*ni], T_SS_1)
    img = plot_error_F(axls[0, 2], F_SS_1[:, 2*ni:], Fhat_SS_1[:, 2*ni:], T_SS_1)
    irow = 1
    for (F, Fhat, T) in [(F_SS_2, Fhat_SS_2, T_SS_2), 
                         (F_DS_1, Fhat_DS_1, T_DS_1), 
                         (F_DS_2, Fhat_DS_2, T_DS_2)]:
        m, n = F.shape
        ni = n // 2
        plot_error_F(axls[irow, 0], F[:, :ni], Fhat[:, :ni], T)
        plot_error_F(axls[irow, 1], F[:, ni:], Fhat[:, ni:], T)
        irow = irow + 1

    # l = 0.07; r = 0.98
    # figls.subplots_adjust(bottom=0.25, left=l, right=r, top=0.9,hspace=0.15, wspace=0.15)
    cax = figls.add_axes([0.75, 0.25, 0.04, 0.4])
    cbar = figls.colorbar(img, cax=cax, orientation="vertical", label="$T_j$ [K]")

    figls.text(
        0.84, 0.15,
        r"$R=\sum_{i,j}$ %s $ - \mathbf{F}_{j,i}$" % sname, 
        va='center', ha='center', color="purple"
        )

    for ax in [axls[1, 2], axls[2, 2], axls[3, 2]]:
        for k in ('top', 'bottom', 'right', 'left'):
            ax.spines[k].set_visible(False)
        ax.tick_params(which='both', length=0.)
        plt.setp(ax.get_xticklabels(), visible=False)

    kwargs = dict(xy=(0.05, 0.89), xycoords='axes fraction')
    axls[0, 0].annotate("(SS,A,1)", **kwargs)
    axls[0, 1].annotate("(SS,B,1)", **kwargs)
    axls[0, 2].annotate("(SS,C,1)", **kwargs)
    axls[3, 0].set_xlabel(r"$\mathbf{F}_{j,i}$")
    axls[3, 1].set_xlabel(r"$\mathbf{F}_{j,i}$")
    for ax in (axls[0, 2], ):
        ax.set_xlabel(r"$\mathbf{F}_{j,i}$")
        ax.tick_params(labelbottom=True)

    for irow in range(4):
        axls[irow, 0].set_ylabel(sname, rotation=0)
    axls[1, 0].annotate("(SS,A,2)", **kwargs)
    axls[1, 1].annotate("(SS,B,2)", **kwargs)
    axls[2, 0].annotate("(DS,A,1)", **kwargs)
    axls[2, 1].annotate("(DS,B,1)", **kwargs)
    axls[3, 0].annotate("(DS,A,2)", **kwargs)
    axls[3, 1].annotate("(DS,B,2)", **kwargs)

    for r in range(4):
        for c in range(3):
            axls[r, c].spines["top"].set_visible(False)
            axls[r, c].spines["right"].set_visible(False)

    figls.subplots_adjust(left=0.1, right=0.98, top=0.99, hspace=0.0, wspace=0., bottom=0.08)
    figls.savefig(figure_name_to_abspath(fname), transparent=True, dpi=300)


def plot_figure6(
    M_mean, M_std, T,
    color_r1="C0", color_r2="C1", 
    ):
    fig, axes = plt.subplots(ncols=2, nrows=1, sharex=True, sharey=True, figsize=(3.25, 2))

    def plot(T_k, M_k, M_std_k, ax, color, label=None):
        ax.plot(T_k, M_k, '-', lw=0.5, color=color, label=label)
        ax.fill_between(T_k, M_k+3*M_std_k, M_k - 3*M_std_k,
                                 color=color, alpha=0.5)
    
    M_SS_1, M_SS_2, M_DS_1, M_DS_2 = M_mean
    M_std_SS_1, M_std_SS_2, M_std_DS_1, M_std_DS_2 = M_std
    T_SS_1, T_SS_2, T_DS_1, T_DS_2 = T

    plot(T_SS_1, M_SS_1, M_std_SS_1, axes[0], color_r1, r"$\mathbf{D}=1$")
    plot(T_SS_2, M_SS_2, M_std_SS_2, axes[0], color_r2, r"$\mathbf{D}=2$")
    plot(T_SS_1, 2*M_SS_1, 2*M_std_SS_1, axes[0], "C2", r"$2\mathbf{M}_{j,1}^\mathrm{TLS}$")

    plot(T_DS_1, M_DS_1, M_std_DS_1, axes[1], color_r1)
    plot(T_DS_2, M_DS_2, M_std_DS_2, axes[1], color_r2)
    plot(T_DS_1, 2*M_DS_1, 2*M_std_DS_1, axes[1], "C2")

    axes[0].annotate("SS", xy=(0.5, 0.9), xycoords='axes fraction', ha='center')
    axes[1].annotate("DS", xy=(0.5, 0.9), xycoords='axes fraction', ha='center')

    kwargs = dict(borderpad=0., labelspacing=0.1, handlelength=1., handletextpad=0.1, fontsize=9)
    axes[0].legend(framealpha=0., loc=(0.55, 0.55), **kwargs)

    for ax in axes:
        ax.tick_params(which='both', direction='in')
        ax.set_xlabel(r"$T_j$ [K]")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    
    axes[0].set_ylabel(r"$\mathbf{M}_{j,\mathbf{D}}^\mathrm{TLS}$", rotation=0, labelpad=12)

    for ax in axes:
        ax.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9], minor=True)
        ax.set_ylim([0.1, 1.0])

    fig.subplots_adjust(left=0.18, bottom=0.20, right=0.99, top=0.97, hspace=0., wspace=0.03)
    fig.savefig(figure_name_to_abspath("figure6.png"), transparent=True, dpi=300)
    plt.close(fig)