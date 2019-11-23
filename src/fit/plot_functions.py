import numpy as np
from scipy import stats
import matplotlib

matplotlib.use('Agg')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import TransformedBbox, BboxPatch, \
    BboxConnector
import matplotlib.pyplot as plt
from utilities import em_algorithm as emCC
import seaborn as sns
import warnings

warnings.filterwarnings("ignore",
                        category=UserWarning)  # to suppress matplotlib warnings


bin_color = "#235789"
fit_line_color = "#CD534C99"
bin_line_width = 0
fit_line_width = 1
dots_color = "#CD534C99"
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))


def mark_inset(parent_axes, inset_axes, loc1a=1, loc1b=1, loc2a=2, loc2b=2,
               **kwargs):
    rect = TransformedBbox(inset_axes.viewLim, parent_axes.transData)
    pp = BboxPatch(rect, fill=False, linewidth=.5, **kwargs)
    parent_axes.add_patch(pp)

    p1 = BboxConnector(inset_axes.bbox, rect, loc1=loc1a, loc2=loc1b,
                       linewidth=.5, **kwargs)
    inset_axes.add_patch(p1)
    p1.set_clip_on(False)
    p2 = BboxConnector(inset_axes.bbox, rect, loc1=loc2a, loc2=loc2b,
                       linewidth=.5, **kwargs)
    inset_axes.add_patch(p2)
    p2.set_clip_on(False)

    return pp, p1, p2


def plot_EM_fitted_distribution(lambd,
                                phi,
                                half_genome_size,
                                GO_pairs_distances,
                                likelihood,
                                organism_name,
                                plot_file_name,
                                args):
    global bin_color
    global fit_line_color
    global bin_line_width
    global fit_line_width
    global dots_color
    global matplotlib_file

    with matplotlib.rc_context(fname=args["matplotlib_style_file"]):
        plt.rcParams.update({'axes.labelsize': 8,
                             'axes.titlesize': 8,
                             'xtick.labelsize': 7,
                             'ytick.labelsize': 7})

        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        points_x1 = np.linspace(0, half_genome_size, 1000)
        exp = emCC.expon_pdf(points_x1, lambd)
        uni = stats.uniform.pdf(points_x1, 0, half_genome_size)
        points_y1 = phi * exp + (1 - phi) * uni

        fig, ((sub1), (sub3)) = plt.subplots(figsize=(4.5, 6.0),
                                             nrows=2, ncols=1,
                                             facecolor="white",
                                             # constrained_layout=True,
                                             sharex=False,
                                             sharey=False)

        bins = np.arange(0, max(GO_pairs_distances), 900)

        sub1.plot(points_x1, points_y1,
                  linestyle="dashed",
                  color=fit_line_color,
                  linewidth=.8,
                  alpha=1)
        sub1.hist(GO_pairs_distances,
                  normed=1,
                  color=bin_color,
                  alpha=1,
                  bins=bins,
                  edgecolor=bin_color,
                  linewidth=0.8)
        sub1.set_xlim([0, half_genome_size])
        sub1.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda x, p: format(int(x / 1000), ',')))
        y1, y2 = sub1.get_ylim()

        # Fake inset to get the right shape of the marker box
        axins2 = inset_axes(sub1, height="40%", width="30%", loc=1, borderpad=1)
        axins2.set_xlim(-20, 0)
        axins2.set_ylim(0, y2 / 1000)
        plt.yticks(visible=False)
        plt.xticks(visible=False)
        mark_inset(sub1, axins2, loc1a=2, loc2a=3, loc1b=1, loc2b=4, fc="grey",
                   ec="grey", alpha=.6)

        # Correct inset
        x1 = 0
        axins = inset_axes(sub1, height="40%", width="30%",
                           loc=1, borderpad=1)
        axins.hist(GO_pairs_distances,
                   normed=1,
                   color=bin_color,
                   alpha=1,
                   bins=bins,
                   edgecolor=bin_color,
                   linewidth=0.8)
        axins.plot(points_x1, points_y1,
                   linestyle="dashed",
                   color=fit_line_color,
                   linewidth=.8,
                   alpha=1)
        axins.set_xlim(x1, 20000)
        axins.set_ylim(y1, y2)
        axins.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda x, p: format(int(x / 1000), ',')))
        axins.get_yaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: ""))
        axins.spines['top'].set_visible(True)
        axins.spines['right'].set_visible(True)
        axins.spines['left'].set_visible(True)
        axins.spines['bottom'].set_visible(True)

        for item in ([axins.title, axins.xaxis.label, axins.yaxis.label] +
                     axins.get_xticklabels() + axins.get_yticklabels()):
            item.set_fontsize(6)

        sub1.set_xlabel("Chromosomal domain pair distance / kbp")
        sub1.set_ylabel("Normalised domain pairs count")

        sub3.plot(range(0, len(likelihood)),
                  likelihood,
                  "o",
                  markersize=4,
                  color=dots_color,
                  markeredgewidth=0.0,
                  alpha=1)
        sub3.set_yticklabels(sub3.get_yticks())
        sub3.set_xlabel("EM Algorithm Iteration")
        sub3.set_ylabel("Marginal Log Likelihood")
        sub3.set_xlim([-2, len(likelihood)])

        mm_label = "lambda = {:.3g}, phi = {:.3g}, \n half chromosome size = {} bp" \
            .format(lambd, phi, half_genome_size)
        sub1.set_title(
            r"$\bf{Assembly}$" + " " + r"$\bf{{}" + organism_name.replace("_",
                                                                          "\_") + "}$"
            + "\n{}".format(mm_label))

        fig.set_tight_layout(True)
        fig.savefig(plot_file_name)
        plt.close()


def plot_histogram_bootstrap(parameters_data,
                             param_name,
                             args,
                             plot_file_name):
    global bin_color
    global bin_line_width
    global fit_line_color

    with matplotlib.rc_context(fname=args['matplotlib_style_file']):
        fig, ((ax)) = plt.subplots(1, 1, figsize=(8,5))

        sns.distplot(parameters_data,
                     norm_hist=True,
                     rug=False,
                     hist_kws={'alpha': .7},
                     ax=ax,
                     color="#235789")

        plt.ylabel("Density")
        name = "\overline{\phi^*}" if param_name == "phi" \
            else "\overline{\lambda^*}"
        plt.xlabel(r"$" + name + r"$")

        title = "B" if param_name == "phi" else "A"
        ax.set_title(title,
                     fontsize=16,
                     fontweight='bold',
                     x=-0.05,
                     y=1)
        ax.yaxis.set_ticklabels([])

        fig.set_tight_layout(True)
        fig.savefig(plot_file_name, bbox_inches='tight')
        plt.close()


def zeta_correlation_plots(true_zetas,
                           mean_values_zetas,
                           pearson_correlation,
                           lambd,
                           phi,
                           organism_name,
                           args,
                           plot_file_name):
    global dots_color

    with matplotlib.rc_context(fname=args['matplotlib_style_file']):
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        fig, ((sub1)) = plt.subplots(1, 1, figsize=(7, 6))
        sub1.plot(true_zetas,
                  mean_values_zetas,
                  "o",
                  markersize=3,
                  color=dots_color,
                  markeredgewidth=0.0,
                  alpha=1)
        sub1.set_xlabel("Zetas with organism-specific parameter values")
        sub1.set_ylabel("Zetas with mean parameter values")
        sub1.set_xlim([0, 1])
        sub1.set_ylim([0, 1])

        sub1.set_title(
            r"$\bf{Assembly}$" + " " + r"$\bf{{}" + organism_name.replace("_",
                                                                          "\_") + "}$")

        # place a text box in upper left in axes coords
        props = dict(facecolor='white', alpha=0.5)
        text_pearson = "Pearson correlation={:.3g}\np-value={:.3g}\n" \
                       "lambda={:.3g} \nphi={:.3g}".format(
            pearson_correlation[0],
            pearson_correlation[1],
            lambd,
            phi)
        sub1.text(0.05, 0.95, text_pearson, transform=sub1.transAxes,
                  fontsize=8,
                  verticalalignment='top', bbox=props)

        fig.set_tight_layout(True)
        fig.savefig(plot_file_name)
        plt.close()


def plot_pseudo_auc(tp_rates, fp_rates,
                    plot_file_name, args):
    global dots_color

    with matplotlib.rc_context(fname=args['matplotlib_style_file']):
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        fig, ((sub1)) = plt.subplots(1, 1, figsize=(5, 5))
        sub1.plot(fp_rates,
                  tp_rates,
                  "o",
                  markersize=3,
                  color=dots_color,
                  markeredgewidth=0.0,
                  alpha=1)
        sub1.set_xlabel("FPR")
        sub1.set_ylabel("TPR")
        sub1.set_xlim([0, 1])
        sub1.set_ylim([0, 1])

        fig.set_tight_layout(True)

        fig.savefig(plot_file_name)
        plt.close()
