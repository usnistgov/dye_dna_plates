import matplotlib.pyplot as plt
from .get_data import RawData
from .util import concentration_label, figure_name_to_abspath


def plot_linemap(
    cls: RawData,
    ax, colorbar=False, get_ticks=False, ordered_by_row=True
    ):
    """Plot :math:`F_d^t\\left(T_j, C_i\\right)` vs :math:`C_i`
    for various temperatures :math:`T_j` (colors)

    Parameters
    ----------
    cls : Data
        Instance of data (i.e., a dataset)
    ax : axis
        Matplotlib axis to plot on
    colorbar : bool, optional
        whether or not to plot colorbar, in which case the axis is colorbar axis, by default False
    get_ticks : bool, optional
        whether or not to return list of ticks, by default False
    ordered_by_row : bool, optional
        whether or not well concentrations are ordered by row, by default True

    Returns
    -------
    None or list
        Only returns list of :code:`get_ticks=True`.
    """
    from matplotlib import cm
    colors = cm.viridis((cls.T-cls.T[-1])/(cls.T[0]-cls.T[-1]))
    ticks = []
    if not ordered_by_row:
        ordered_wells = [i*12 + j for j in range(12) for i in range(8)]
    else:
        ordered_wells = list(range(96))
        
    for i_T in range(len(cls.T)):
        if ((cls.T[i_T] - 273.) % 10) < 1e-8:
            if colorbar:
                ax.plot([0, 1], [cls.T[i_T], cls.T[i_T]], '-', color=colors[i_T], lw=3, clip_on=False)
            elif get_ticks:
                ticks.append(cls.T[i_T])
            else:
                ax.plot(cls.C[ordered_wells] * 1e6, cls.F[i_T, ordered_wells] * 1e-6, '-', color=colors[i_T])
    
    if get_ticks:
        return ticks


def make_figure_S1():
    """Makes Figure S1."""

    pluto_plateA_rot = RawData(
        fluorescence_file_name="6-2022_intercalate-pluto-inter_PlateA_Rot_Pluto_6-14-22_data.xls",
        B_d=44e-6, t="SS", l="A", N=22,
        dye_conc_file_name="dye_conc_uM_rotated.csv"
    )
    fig, ax = plt.subplots(sharex=True, sharey=True, figsize=(3.25, 3.25))
    plot_linemap(pluto_plateA_rot, ax)

    # label_kwargs = dict(rotation=0, labelpad=12)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel(r"$F_{i,j,2}^{\mathrm{SS},\alpha}\times 10^{-6}$")  #, **label_kwargs)
    ax.tick_params(direction='in')
    ax.set_xlabel(concentration_label)
    ax.set_ylim([0., 0.3])

    cax = fig.add_axes([0.87, 0.59, 0.08, 0.2])
    for d in ('right', 'left', 'bottom', 'top'):
        cax.spines[d].set_visible(False)
    cax.tick_params(which="both", axis="both", length=0., labelbottom=False)
    plot_linemap(pluto_plateA_rot, cax, colorbar=True)
    cax.set_ylabel("$T_j$ [K]", rotation=0, labelpad=20.)
    ticks = plot_linemap(pluto_plateA_rot, cax, get_ticks=True)
    cax.set_yticks(ticks)
    cax.set_ylim([ticks[0], ticks[-1]])

    fig.subplots_adjust(bottom=0.18, left=0.2, right=0.99, top=0.98, wspace=0.01, hspace=0.01)
    fig.savefig(figure_name_to_abspath("figureS1.pdf"), transparent=True, dpi=300)


def make_figure_S3():
    """Makes Figure S3."""
    import os
    import numpy as np

    kwargs = dict(B_d=1e-6, t="SS", l="A", N=22)
    base_path = os.path.join("6-2022_intercalate", "Pluto")
    plateB_6_15 = RawData(
        fluorescence_file_name=os.path.join(base_path, "inter_PlateB_Pluto_6-15-22_data.xls"),  
        **kwargs
    )
    plateB_6_16 = RawData(
        fluorescence_file_name=os.path.join(base_path, "inter_PlateB_Pluto_6-16-22_data.xls"), 
        **kwargs
    )
    delta_f = RawData(
        fluorescence_file_name=os.path.join(base_path, "inter_PlateB_Pluto_6-15-22_data.xls"), 
        **kwargs
    )

    fig, ax = plt.subplots(ncols=1, nrows=1, sharex=True, sharey=True, figsize=(5., 3.25))
    # plot_linemap(plateB_6_15, axes[0])
    # plot_linemap(plateB_6_16, axes[1])
    delta_f.F[:, :] = plateB_6_15.F - plateB_6_16.F 
    plot_linemap(delta_f, ax)


    df = delta_f.F[:, :-32]
    print('Maximum change from 6/15 to 6/16: %i' % round(np.max(df)))
    print('Minimum change from 6/15 to 6/16: %i' % round(np.min(df)))
    print('Average change from 6/15 to 6/16: %i' % round(np.mean(df)))

    # for ax in (axes[0], axes[1]):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(direction='in')
    ax.set_xlabel(concentration_label)
        # ax.set_ylim([0., 0.3])
    
    # axes[0].annotate("(a)", xy=(0.05, 0.9), xycoords='axes fraction')
    # axes[1].annotate("(b)", xy=(0.05, 0.9), xycoords='axes fraction')

    # label_kwargs = dict(rotation=0, labelpad=10)
    ax.set_ylabel(r"$\Delta F_{i,j,2}^{\mathrm{SS}} \times 10^{-6}$") #, **label_kwargs)
    cax = fig.add_axes([0.94, 0.17, 0.05, 0.75])
    for d in ('right', 'left', 'bottom', 'top'):
        cax.spines[d].set_visible(False)
    cax.tick_params(which="both", axis="both", length=0., labelbottom=False)
    plot_linemap(plateB_6_15, cax, colorbar=True)
    cax.set_ylabel("$T_j$ [K]")
    ticks = plot_linemap(plateB_6_15, cax, get_ticks=True)
    cax.set_yticks(ticks)
    cax.set_ylim([ticks[0], ticks[-1]])

    fig.subplots_adjust(bottom=0.17, left=0.18, right=0.82, top=0.97, wspace=0.01, hspace=0.01)
    fig.savefig(figure_name_to_abspath("figureS3.pdf"), transparent=True, dpi=300)


def make_figure_2(
    SS_A_1: RawData, 
    SS_B_1: RawData, 
    SS_C_1: RawData, 
    SS_A_2: RawData, 
    SS_B_2: RawData, 
    DS_A_1: RawData, 
    DS_B_1: RawData, 
    DS_A_2: RawData, 
    DS_B_2: RawData, 
    A_1: RawData):
    """ Makes figure 2

    Parameters
    ----------
    SS_A_1 : RawData
        Data associated with :math:`(t,l,d)=(SS,A,1)` (see Table 1 of main text)
    SS_B_1 : RawData
        Data associated with :math:`(t,l,d)=(SS,B,1)` (see Table 1 of main text)
    SS_C_1 : RawData
        Data associated with :math:`(t,l,d)=(SS,C,1)` (see Table 1 of main text)
    SS_A_2 : RawData
        Data associated with :math:`(t,l,d)=(SS,A,2)` (see Table 1 of main text)
    SS_B_2 : RawData
        Data associated with :math:`(t,l,d)=(SS,B,2)` (see Table 1 of main text)
    DS_A_1 : RawData
        Data associated with :math:`(t,l,d)=(DS,A,1)` (see Table 1 of main text)
    DS_B_1 : RawData
        Data associated with :math:`(t,l,d)=(DS,B,1)` (see Table 1 of main text)
    DS_A_2 : RawData
        Data associated with :math:`(t,l,d)=(DS,A,2)` (see Table 1 of main text)
    DS_B_2 : RawData
        Data associated with :math:`(t,l,d)=(DS,B,2)` (see Table 1 of main text)
    A_1 : RawData
        Data associated with :math:`(l,d)=(A,1)` (without DNA, see Table 1 of main text)


    """

    # plot heatmaps
    figs, axes = plt.subplots(ncols=3, nrows=4, sharex=True, sharey=True, figsize=(3.25, 5.0))
    label_kwargs = dict(labelpad=0)

    axes[0, 0].annotate("$l = \\mathrm{A}$", xy=(0.05, 0.9), xycoords='axes fraction')
    axes[0, 0].set_ylabel(r"$F_{i,j,1}^{\mathrm{SS},l}\times 10^{-6}$", **label_kwargs)
    axes[0, 1].annotate("$l =\\mathrm{B}$", xy=(0.05, 0.9), xycoords='axes fraction')
    axes[1, 0].set_ylabel(r"$F_{i,j,2}^{\mathrm{SS},l}\times 10^{-6}$", **label_kwargs)
    axes[0, 2].annotate("$l = \\mathrm{C}$", xy=(0.05, 0.9), xycoords='axes fraction')
    axes[2, 0].set_ylabel(r"$F_{i,j,1}^{\mathrm{DS},l}\times 10^{-6}$", **label_kwargs)
    axes[3, 0].set_ylabel(r"$F_{i,j,2}^{\mathrm{DS},l}\times 10^{-6}$", **label_kwargs)

    for irow in range(axes.shape[0]):
        for icol in range(axes.shape[1]):
            axes[irow, icol].spines['top'].set_visible(False)
            axes[irow, icol].spines['right'].set_visible(False)
            axes[irow, icol].tick_params(axis='both', direction='in')
    
    plot_linemap(SS_A_1, axes[0, 0])
    plot_linemap(SS_B_1, axes[0, 1])
    plot_linemap(SS_C_1, axes[0, 2])
    plot_linemap(SS_A_2, axes[1, 0])
    axes[1, 0].annotate("$l = \\mathrm{A}$", xy=(0.05, 0.9), xycoords='axes fraction')
    plot_linemap(SS_B_2, axes[1, 1])
    axes[1, 1].annotate("$l = \\mathrm{B}$", xy=(0.05, 0.9), xycoords='axes fraction')

    for i in ((2, 2), (1, 2), (3, 2)):
        for kk in ('left', 'right', 'top', 'bottom'):
            axes[i].spines[kk].set_visible(False)
            axes[i].tick_params(which='both', length=0.)
            # axes[i].tick_params(axis='y', direction='in')

    axes[0, 0].set_xticks([0, 1, 2, 3, 4])
    plt.setp(axes[3, 2].get_xticklabels(), visible=False)
    plt.setp(axes[0, 2].get_xticklabels(), visible=True)
    axes[0, 2].tick_params(labelbottom=True)
    axes[3, 0].set_xlabel(concentration_label)
    axes[3, 1].set_xlabel(concentration_label)
    axes[0, 2].set_xlabel(concentration_label)

    plot_linemap(DS_A_1, axes[2, 0])
    axes[2, 0].annotate("$l = \\mathrm{A}$", xy=(0.05, 0.9), xycoords='axes fraction')
    plot_linemap(DS_B_1, axes[2, 1])
    axes[2, 1].annotate("$l = \\mathrm{B}$", xy=(0.05, 0.9), xycoords='axes fraction')
    plot_linemap(DS_A_2, axes[3, 0])
    axes[3, 0].annotate("$l = \\mathrm{A}$", xy=(0.05, 0.9), xycoords='axes fraction')
    plot_linemap(DS_B_2, axes[3, 1])
    axes[3, 1].annotate("$l = \\mathrm{B}$", xy=(0.05, 0.9), xycoords='axes fraction')

    ax_nodna = figs.add_axes([0.8, 0.2, 0.2, 0.2])
    c = 'blueviolet'
    # ax_nodna = axes[3, 2]
    ax_nodna.annotate("$l = \\mathrm{A}$", xy=(0.05, 1.02), xycoords='axes fraction', color=c)
    # ax_nodna.spines['left'].set_visible(False)
    # plt.setp(ax_nodna.get_xticklabels(), visible=True)
    ax_nodna.tick_params(labelbottom=True)
    # ax = ax_nodna.twinx()
    # ax.tick_params(axis='y', labelcolor=c, color=c)
    ax_nodna.tick_params(axis='both', labelcolor=c, color=c)
    ax_nodna.spines['right'].set_visible(False)
    ax_nodna.spines['top'].set_color(c)
    # ax.arrow(2., 11, 1., 0., color=c, length_includes_head=True, 
    #          width=1.e-4, head_width=1, head_length=0.2)
    ax_nodna.spines['bottom'].set_color(c)
    ax_nodna.spines['left'].set_color(c)
    A_1.F = A_1.F*1000.
    plot_linemap(A_1, ax_nodna)
    ax_nodna.set_ylabel(r"$F_{i,j,0}^\mathrm{A} \times 10^{-3}$", color=c, labelpad=-8)
                        

    ax_nodna.set_ylim([-5., 5.])
    ax_nodna.set_xlabel(concentration_label, color=c)
    ax_nodna.tick_params(which='both', direction='in', color=c)
    cax = figs.add_axes([0.88, 0.47, 0.05, 0.18])

    cax.tick_params(which="both", axis="x", length=0., labelbottom=False)
    plot_linemap(SS_A_1, cax, colorbar=True)
    # cax.set_ylabel("$T_j$", rotation=0, labelpad=10, fontsize=14)
    ticks = plot_linemap(SS_A_1, cax, get_ticks=True)
    cax.set_yticks(ticks)
    cax.set_ylim([ticks[0], ticks[-1]])
    # cax2 = cax.twinx()
    # ticks_Celcius = [i - 273 for i in ticks]
    # cax2.set_yticks(ticks_Celcius)
    # cax2.set_ylim([ticks_Celcius[0], ticks_Celcius[-1]])

    for d in ('right', 'left', 'top', 'bottom'):
        for a in (cax, 
                #   cax2
                  ):
            a.spines[d].set_visible(False)
            a.tick_params(axis='y', length=0.)
    
    cax.set_yticklabels(['343 K', '333 K', '323 K', '313 K', '303 K', '293 K'])
    # cax2.set_yticklabels(['$70^\\circ$C', '$60^\\circ$C', '$50^\\circ$C', 
    #                       '$40^\\circ$C', '$30^\\circ$C', '$20^\\circ$C'])

    figs.subplots_adjust(bottom=0.12, left=0.12, right=0.99, top=0.99, wspace=0.05, hspace=0.05)
    figs.savefig(
        figure_name_to_abspath("figure2.pdf"), 
        transparent=True, dpi=300
        )