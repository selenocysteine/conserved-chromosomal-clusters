import logging
import seaborn as sns
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.offsetbox
import matplotlib.transforms as mtransforms
import math
import numpy as np
import ete3
import os

logger = logging.getLogger(__name__)


def plot_perm(args):
    out_file = args['permutations_file']
    data = pd.read_csv(out_file, sep="\t")
    with matplotlib.rc_context(fname=args["matplotlib_style_file"]):
        sns.set(rc={'figure.figsize': (12, 16)})
        fig, ax = plt.subplots(1)
        data = data[['n_members', 'occurrences', 'pct_above_0.5']]
        data = data.groupby(['n_members', 'occurrences'], as_index=False
                            ).mean()
        data = data.pivot_table("pct_above_0.5", "n_members",
                                "occurrences")
        names = \
            np.transpose(np.array(pd.DataFrame(data.to_records())
                                  [['n_members']].values))[0]
        log_norm = LogNorm(vmin=0.0001, vmax=1.0001)
        cbar_ticks = [math.pow(10, i) for i in
                      range(math.floor(-4),
                            1)]
        ax = sns.heatmap(data,
                         cmap="OrRd_r",
                         cbar=True,
                         vmin=10 ** -4,
                         vmax=1,
                         cbar_kws={"ticks": cbar_ticks,
                                   "shrink": .4,
                                   "pad": 0.09},
                         ax=ax,
                         norm=log_norm,
                         yticklabels=names,
                         alpha=.9)

        ax.invert_yaxis()
        n_labels = len(ax.get_yticklabels()) - 1

        ax.tick_params(labeltop=True, labelright=True)
        ax.set_xticklabels(ax.get_xticklabels(),
                           rotation=0)
        ax.set_yticklabels(ax.get_yticklabels(),
                           rotation=0)
        for ind, label in enumerate(ax.get_yticklabels()):
            if ind % 8 == 0 or ind == n_labels:
                label.set_visible(True)
            else:
                label.set_visible(False)

        ax.set(ylabel='Organisms in subclade',
               xlabel='Pair instances in subclade')
        ax.collections[0].colorbar.set_label(
            r'Average $\hat p$ for $CCP$ > 0.5')
        ax.set_facecolor('white')
        plt.tight_layout()

        fig.savefig(os.path.join(args["permutations_dir"],
                                 "permutation_plot.png"))
        plt.close()




def count_distance(node, distances):
    if node.is_leaf():
        all_dist = 0
        node.name = node.name.split("__")[0]
    else:
        all_dist = 0
        for child in node.get_children():
            all_dist += max(0, child.dist) + count_distance(child, distances)

    distances[node.name] = all_dist
    return all_dist



def plot_occurrences_tree(data, occurrences, distances, args):
    data2 = \
        data[data['occurrences'] == occurrences]
    [['subclade_id', 'n_members', 'pct_above_0.5']]

    distances_array = np.zeros(data2['subclade_id'].values.shape[0])
    for index, subclade in enumerate(data2['subclade_id'].values):
        distances_array[index] = distances[str(subclade)]

    indexes = np.argsort(distances_array)
    distances_array = distances_array[indexes]
    distances_array = (distances_array - min(distances_array)) / (
            max(distances_array) - min(distances_array))
    values = np.array(data2['pct_above_0.5'].values)[indexes]
    n_organisms = np.array(data2['n_members'].values, dtype=int)[indexes]
    with matplotlib.rc_context(fname=args["matplotlib_style_file"]):
        n_bins = 50
        fig = plt.subplots(1, figsize=(6, 12))
        values = np.log10(values)
        g = sns.jointplot(distances_array,
                          values,
                          color='black',
                          space=0,
                          ylim=(-5.2, 0.1),
                          ratio=6)
        g.ax_joint.cla()
        # g.ax_joint.axvspan(-1e-6, 1e-6, color='red', alpha=0.1)
        g.ax_joint.scatter(distances_array,
                          values,
                             c=n_organisms,
                             s=15,
                             cmap="viridis_r",
                             linewidths=0,
                             alpha=.5,
                             norm=matplotlib.colors.LogNorm(vmin=1, vmax=14178))


        step = 0.05
        g.ax_marg_x.cla()
        g.ax_marg_y.cla()
        _ = g.ax_marg_x.hist(distances_array, color="grey", alpha=.4,
                             bins=np.arange(0, 1, 1e-5))
        _ = g.ax_marg_y.hist(values, color="grey", alpha=.4,
                             orientation="horizontal",
                             bins=np.arange(-5, 1.1,
                                            step),
                             normed=True)

        g.ax_joint.set_xscale('symlog', linthreshx=1e-5)
        g.ax_marg_x.set_xscale('symlog', linthreshx=1e-5)
        g.ax_joint.set(ylim=(-5.1, 0.1), xlim=(-1e-6,1.5))
        g.ax_marg_y.set_xscale('log')
        g.ax_marg_x.set_yscale('log')
        g.ax_joint.set_yticks([-5, -4, -3, -2, -1, 0])
        g.ax_joint.spines['top'].set_visible(True)
        g.ax_joint.spines['right'].set_visible(True)
        g.ax_joint.set_yticklabels([r"$10^{" + str(x) + "}$"
                                    for x in [-5, -4, -3, -2, -1, 0]])
        g.ax_marg_y.grid('off')
        g.ax_marg_y.axis('off')
        g.ax_marg_x.grid('off')
        g.ax_marg_x.axis('off')
        g.ax_joint.set(ylabel=r'$\hat \alpha$ for $CCP$ > 0.5',
               xlabel=r'Phylogenetic variation in subclade $s$')
        if str(occurrences).startswith('5'):
            g.ax_joint.set(ylabel='')
        plt.title(r'Simulation results for $M$ = {:,}'.format(occurrences),
                  y = 1.2, x=-3.2, fontsize = 16)
        # plt.tight_layout()

        plt.savefig(os.path.join(args["permutations_dir"],
                                 "permutation_plot_distances"
                                 "{}_scatterplot.png".format(occurrences)),
                                 figsize=(8,15),
        bbox_inches = 'tight')
        plt.close()

        # if occurrences == 100000:
        #     fig, ax = plt.subplots(figsize=(10, 1.2))
        #     cmap = 'cubehelix'
        #     norm = matplotlib.colors.LogNorm(vmin=1, vmax=np.max(occurrences))
        #
        #     # cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
        #     #                                        norm=norm,
        #     #                                        orientation='horizontal',
        #     #                                        alpha=.7,
        #     #                                        format=
        #     #                                        matplotlib.ticker.
        #     #                                        FuncFormatter(
        #     #                                            show_only_some))
        #     cb1.set_label(r'Number of organisms in subclade $s$')
        #
        #     plt.tight_layout()
        #     plt.savefig(os.path.join(args['permutations_dir'],
        #                                  "legend.png"))



def scatterplots(args):
    data = pd.read_csv(args['permutations_file'], sep="\t")
    n_organisms = np.array(data['n_members'].values)
    indexes = np.argsort(n_organisms)
    n_organisms = n_organisms[indexes]
    values = np.array(data['pct_above_0.5'].values)[indexes]
    n_occurrences = np.array(data['occurrences'].values)[indexes]
    sizes = (np.around(np.log10(n_occurrences), decimals=1) + 1)**3 * 10
    print(sizes)
    with matplotlib.rc_context(fname=args["matplotlib_style_file"]):
        fig = plt.subplots(1, figsize=(6, 12))
        values = np.log10(values)
        g = sns.jointplot(n_organisms,
                          values,
                          color='black',
                          space=0,
                          ylim=(-5.2, 0.1),
                          ratio=6)
        g.ax_joint.cla()
        scatter = g.ax_joint.scatter(n_organisms,
                          values,
                             s=sizes * 0.5,
                             linewidths=1,
                             alpha=.5,
                                     edgecolors='grey')

        g.ax_joint.legend().set_visible(False)

        step = 0.05
        g.ax_marg_x.cla()
        g.ax_marg_y.cla()
        _ = g.ax_marg_x.hist(n_organisms, color="grey", alpha=.4,
                             bins=np.arange(0, 14178, 100), normed=True)
        _ = g.ax_marg_y.hist(values, color="grey", alpha=.4,
                             orientation="horizontal",
                             bins=np.arange(-5, 1.1,
                                            step),
                             normed=True)

        g.ax_joint.set(ylim=(-5.1, 0.1), xlim=(-1000, 16000))
        # g.ax_marg_y.set_xscale('log')
        # g.ax_marg_x.set_yscale('log')
        g.ax_joint.set_yticks([-5, -4, -3, -2, -1, 0])
        g.ax_joint.spines['top'].set_visible(True)
        g.ax_joint.spines['right'].set_visible(True)
        g.ax_joint.set_yticklabels([r"$10^{" + str(x) + "}$"
                                    for x in [-5, -4, -3, -2, -1, 0]])
        g.ax_marg_y.grid('off')
        g.ax_marg_y.axis('off')
        g.ax_marg_x.grid('off')
        g.ax_marg_x.axis('off')
        g.ax_joint.set(ylabel=r'$\hat \alpha$ for $CCP$ > 0.5',
               xlabel=r'Number of organisms in subclade $s$')
        # plt.title(r'Simulation results for $M$ = {:,}'.format(occurrences),
        #           y = 1.2, x=-3.2, fontsize = 16)

        plt.savefig(os.path.join(args["permutations_dir"],
                                 "permutation_plot_distances"
                                 "_n_org.png"),
                                 figsize=(8,15),
        bbox_inches = 'tight')
        g.ax_joint.legend().set_visible(True)

        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)

        print(handles)
        labels = [5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                   ncol=2, mode="expand", borderaxespad=0.)
        legend = g.ax_joint.legend(handles, labels,
                                   labelspacing=1.2,
                            title=r"Total instances $M$",
                            borderpad=1,
                            frameon=True,
                            framealpha=0.6,
                            edgecolor="k", facecolor="w",
                            ncol=12,
                            fontsize=10)


        # Plot sizes legend
        # pws = np.unique(sizes)
        # pws.sort()
        # keys = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]
        # for index, pw in enumerate(pws):
        #     plt.scatter([], [], s=pw, label=str(keys[index]), color='C0')
        # legend = plt.legend(loc="upper right", labelspacing=1.2,
        #                     title=r"Total instances $M$",
        #                     borderpad=1,
        #                     frameon=True,
        #                     framealpha=0.6,
        #                     edgecolor="k", facecolor="w",
        #                     ncol=int(len(pws) / 2) + 1,
        #                     fontsize=10)
        # plt.cla()


        fig = legend.figure
        fig.canvas.draw()
        bbox = legend.get_window_extent()
        bbox = bbox.from_extents(*(bbox.extents + np.array([5, 5, 5, 5])))
        bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
        plt.tight_layout()
        plt.savefig(os.path.join(args['permutations_dir'],
                                 "legend_scatterplot.png"),
                    bbox_inches=bbox)

        # Plot n_organisms vs alpha
        fig = plt.subplots(1, figsize=(6, 12))
        g = sns.jointplot(n_occurrences,
                          values,
                          color='black',
                          space=0,
                          ylim=(-5.2, 0.1),
                          ratio=6)
        g.ax_joint.cla()
        g.ax_joint.scatter(n_occurrences,
                           values,
                           c=n_organisms,
                           s=20,
                           cmap="viridis_r",
                           linewidths=0,
                           alpha=.7,
                           norm=matplotlib.colors.LogNorm(vmin=1, vmax=14178))
        step = 0.05
        g.ax_marg_x.cla()
        g.ax_marg_y.cla()
        _ = g.ax_marg_x.hist(n_organisms, color="grey", alpha=.4,
                             bins=np.arange(1, 100000, 1), normed=True)
        _ = g.ax_marg_y.hist(values, color="grey", alpha=.4,
                             orientation="horizontal",
                             bins=np.arange(-5, 1.1,
                                            step),
                             normed=True)

        g.ax_joint.set_xscale('log')
        g.ax_joint.set(ylim=(-5.1, 0.1),
                       xlim=(0.9, 1000000))

        # g.ax_marg_y.set_xscale('log')
        # g.ax_marg_x.set_yscale('log')
        g.ax_joint.set_yticks([-5, -4, -3, -2, -1, 0])
        g.ax_joint.spines['top'].set_visible(True)
        g.ax_joint.spines['right'].set_visible(True)
        g.ax_joint.set_yticklabels([r"$10^{" + str(x) + "}$"
                                    for x in [-5, -4, -3, -2, -1, 0]])
        g.ax_marg_y.grid('off')
        g.ax_marg_y.axis('off')
        g.ax_marg_x.grid('off')
        g.ax_marg_x.axis('off')
        g.ax_joint.set(ylabel=r'$\hat \alpha$ for $CCP$ > 0.5',
                       xlabel=r'$M$')
        # plt.title(r'Simulation results for $M$ = {:,}'.format(occurrences),
        #           y = 1.2, x=-3.2, fontsize = 16)

        plt.savefig(os.path.join(args["permutations_dir"],
                                 "permutation_plot_distances"
                                 "_occurrences.png"),
                    figsize=(8, 15),
                    bbox_inches='tight')
        plt.close()
        sys.exit()


def show_only_some(x, pos):
    if x == int(10 ** (math.floor(math.log10(x)))):
        return '{:g}'.format(x)
    return ''


def plot_perm_tree(args):
    tree = ete3.Tree(args['labelled_tree_file'], format=1)
    tree.name = "0"
    distances = {}
    count_distance(tree, distances)
    out_file = args['permutations_file']
    data = pd.read_csv(out_file, sep="\t")
    for occurrences in [1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]:
        plot_occurrences_tree(data, occurrences, distances, args)
