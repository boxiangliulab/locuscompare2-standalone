import functools
import os
from math import sqrt
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from mayavi import mlab
from scipy.interpolate import RegularGridInterpolator

import sim_result_plot


def gwas_tissue(target_dir):
    for f1 in os.listdir(target_dir):
        if f1.startswith('.'):
            continue
        h1_coloc = None
        h1_predixcan = None
        h1_smr = None
        h1_fastenloc = None
        h1_twas = None
        h1_ecaviar = None
        tissue_full_path = os.path.join(target_dir, f1)
        if os.path.isfile(tissue_full_path):
            continue
        print(f'Processing files in dir {tissue_full_path}')
        for f2 in os.listdir(tissue_full_path):
            f2_full_path = os.path.join(tissue_full_path, f2)
            if f2.startswith('coloc_output_'):
                h1_coloc = f2_full_path
            elif f2.startswith('predixcan_output_'):
                h1_predixcan = f2_full_path
            elif f2.startswith('smr_output_'):
                h1_smr = f2_full_path
            elif f2.startswith('fastenloc_output_'):
                h1_fastenloc = f2_full_path
            elif f2.startswith('twas_output_'):
                h1_twas = f2_full_path
            elif f2.startswith('ecaviar_output_'):
                h1_ecaviar = f2_full_path
        fig_path = os.path.join(tissue_full_path, 'HEAT_raw.png')
        if h1_coloc is None and h1_smr is None and h1_fastenloc is None and h1_ecaviar is None and h1_twas is None:
            print(f'No reports found in dir {tissue_full_path}')
        try:
            sim_result_plot.plot_spearman_heatmap({'coloc': h1_coloc,
                                                   'smr': h1_smr,
                                                   'fastenloc': h1_fastenloc,
                                                   'predixcan': h1_predixcan,
                                                   'ecaviar': h1_ecaviar,
                                                   'twas': h1_twas}, genetic_model='',
                                                  output_figure_path=fig_path)
        except Exception as exc:
            print(
                f'plot corr failed for dir {tissue_full_path}, maybe too many genes are missing in reports of some tools')


def gwas_no_tissue(target_dir):
    for f1 in os.listdir(target_dir):
        tissue_full_path = os.path.join(target_dir, f1)
        if f1.startswith('.') or os.path.isfile(tissue_full_path) or not f1.startswith('GCST'):
            continue
        h1_coloc = None
        h1_predixcan = None
        h1_smr = None
        h1_fastenloc = None
        h1_twas = None
        h1_ecaviar = None
        if os.path.isfile(tissue_full_path):
            continue
        print(f'Processing files in dir {tissue_full_path}')
        for f2 in os.listdir(tissue_full_path):
            f2_full_path = os.path.join(tissue_full_path, f2)
            if f2.startswith('coloc_output_'):
                h1_coloc = f2_full_path
            elif f2.startswith('predixcan_output_'):
                h1_predixcan = f2_full_path
            elif f2.startswith('smr_output_'):
                h1_smr = f2_full_path
            elif f2.startswith('fastenloc_output_'):
                h1_fastenloc = f2_full_path
            elif f2.startswith('twas_output_'):
                h1_twas = f2_full_path
            elif f2.startswith('ecaviar_output_'):
                h1_ecaviar = f2_full_path
        if h1_coloc is None and h1_smr is None and h1_fastenloc is None and h1_ecaviar is None and h1_twas is None:
            print(f'No reports found in dir {tissue_full_path}')
        fig_path = os.path.join(tissue_full_path, 'HEAT_raw.png')
        try:
            sim_result_plot.plot_spearman_heatmap({'coloc': h1_coloc,
                                                   'smr': h1_smr,
                                                   'fastenloc': h1_fastenloc,
                                                   'predixcan': h1_predixcan,
                                                   'ecaviar': h1_ecaviar,
                                                   'twas': h1_twas}, genetic_model='',
                                                  output_figure_path=fig_path)
        except Exception as exc:
            print(
                f'plot corr failed for dir {tissue_full_path}, maybe too many genes are missing in reports of some tools')


def download_from_aws(file_list_file, dest_dir):
    with open(file_list_file) as f_list:
        for f in f_list:
            f = f.strip()
            splits = f.split(os.sep)
            tissue = splits[-5]
            file_name = splits[-1]
            tissue_dir = os.path.join(dest_dir, tissue)
            Path(tissue_dir).mkdir(exist_ok=True, parents=True)
            dest_file = os.path.join(tissue_dir, file_name)
            print(f'Downloading {f} to dest {dest_file}')
            os.system(f'rsync -e "ssh -i ~/.ssh/liulab-zhijiang.pem" -avhP linyi@65.1.12.42:{f} {dest_file}')
            if not os.path.exists(dest_file):
                raise ValueError(f'{f} synced failed')


def download_from_28(file_list_file, dest_dir):
    with open(file_list_file) as f_list:
        for f in f_list:
            f = f.strip()
            splits = f.split(os.sep)
            tissue = splits[-5]
            file_name = splits[-1]
            tissue_dir = os.path.join(dest_dir, tissue)
            Path(tissue_dir).mkdir(exist_ok=True, parents=True)
            dest_file = os.path.join(tissue_dir, file_name)
            print(f'Downloading {f} to dest {dest_file}')
            os.system(f'cp {f} {dest_file}')
            if not os.path.exists(dest_file):
                raise ValueError(f'{f} synced failed')


def calc_overview_data(target_dir, prefix=''):
    df_list = []
    for f1 in os.listdir(target_dir):
        tissue_full_path = os.path.join(target_dir, f1)
        if f1.startswith('.') or os.path.isfile(tissue_full_path):
            continue
        if prefix and not f1.startswith('GCST'):
            continue
        print(f'Processing files in dir {tissue_full_path}')
        # _pre_heatmap.tsv are output from gwas_tissue() and gwas_no_tissue()
        rank_file = os.path.join(tissue_full_path, '_pre_heatmap.tsv')
        if not os.path.exists(rank_file):
            print(f'Skipping {tissue_full_path} because _pre_heatmap.tsv does not exist')
            continue
        df = pd.read_table(rank_file)
        if df.shape[1] < 7:
            print(f'Skipping {tissue_full_path} because reports of some tools are missing')
            continue
        if (df.isna().sum() == df.shape[0]).values.any():
            print(f'Skipping {tissue_full_path} because ranking of some tools are all NA')
            continue
        corr = df.corr(method='spearman', numeric_only=True)
        if (corr.isna().sum() != 0).values.any():
            print(f'Skipping {tissue_full_path} because correlation has NAs')
            continue
        df_list.append(corr)
    if len(df_list) == 0:
        print(f'Target dir {target_dir} has no complete reports')
        return
    mean_df = functools.reduce(lambda a, b: a + b, df_list) / len(df_list)
    std_df = (functools.reduce(lambda a, b: a + b, map(lambda a: (a - mean_df).pow(2), df_list)) / len(df_list)).apply(
        np.sqrt)
    df_se = std_df / sqrt(len(df_list))
    # 95% CI,
    # for T-distribution, scipy.stats.t.ppf;
    # for normal distribution z=1.96
    z = scipy.stats.t.ppf(q=1 - 0.05 / 2, df=len(df_list) - 1)
    ci_lower = mean_df - z * df_se
    ci_upper = mean_df + z * df_se
    mean_out = os.path.join(target_dir, 'mean.tsv')
    mean_df.to_csv(mean_out, sep='\t', header=True, index=True, na_rep='NA')
    std_out = os.path.join(target_dir, 'std.tsv')
    std_df.to_csv(std_out, sep='\t', header=True, index=True, na_rep='NA')
    ci_lower_out = os.path.join(target_dir, 'ci_lower.tsv')
    ci_lower.to_csv(ci_lower_out, sep='\t', header=True, index=True, na_rep='NA')
    ci_upper_out = os.path.join(target_dir, 'ci_upper.tsv')
    ci_upper.to_csv(ci_upper_out, sep='\t', header=True, index=True, na_rep='NA')
    print(f'Mean\n: {mean_df}')
    print(f'Std\n: {mean_df}')
    print(f'95_ci_lower:\n {ci_lower}')
    print(f'95_ci_upper:\n {ci_upper}')
    return mean_out, std_out, ci_lower_out, ci_upper_out


def plot_surface_with_errorbar(mean_tsv, std_tsv, output_figure_path):
    mean_df = pd.read_table(mean_tsv, index_col=0)
    std_df = pd.read_table(std_tsv, index_col=0)
    # x,y 从1开始, 否则第一个工具名称label不显示
    x = range(1, mean_df.shape[0] + 1, 1)
    y = range(1, mean_df.shape[0] + 1, 1)
    X, Y = np.meshgrid(x, y)
    Z = mean_df.to_numpy()
    interp_step_size = 0.05
    # 保留原始的x和y
    total_num = int(1 + (mean_df.shape[0] - 1) / interp_step_size)
    x_new = np.linspace(1, mean_df.shape[0], total_num)
    y_new = np.linspace(1, mean_df.shape[0], total_num)
    X_new, Y_new = np.meshgrid(x_new, y_new)
    interp = RegularGridInterpolator((x, y), Z, method='linear', bounds_error=False, fill_value=None)
    Z_new = interp((X_new, Y_new))
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 8))
    ax.view_init(elev=45, azim=-45)
    cmap = 'vlag'
    surface = ax.plot_surface(X_new, Y_new, Z_new, cmap=cmap, rstride=1, cstride=1)
    ax.errorbar(X.flatten(), Y.flatten(), Z.flatten(), zerr=std_df.to_numpy().flatten(), fmt=',', ecolor='black',
                capsize=5, zorder=10)
    ax.set(xticks=[0] + [*x], yticks=[0] + [*y],
           # zticks=[0, .2, .4, .6, .8, 1.0],
           xticklabels=[''] + [*mean_df.columns],
           yticklabels=[''] + [*mean_df.columns])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('correlation')
    # # show colorbar from data
    # cbar = fig.colorbar(surface, shrink=0.2, aspect=10, location='right', pad=.1)
    # tick_locator = ticker.MaxNLocator(nbins=6)
    # cbar.locator = tick_locator
    # cbar.update_ticks()
    # cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # cbar.set_ticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1.0'])
    # # show colorbar from cmap
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap), ax=ax, shrink=0.2, aspect=10, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
                 location='top', pad=0.05)
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


# @mlab.show
def plot_mayavi_surface_with_errorbar(mean_tsv, std_tsv, output_figure_path):
    mean_df = pd.read_table(mean_tsv, index_col=0)
    std_df = pd.read_table(std_tsv, index_col=0)
    # x,y 从1开始, 否则第一个工具名称label不显示
    x = range(1, mean_df.shape[0] + 1, 1)
    y = range(1, mean_df.shape[0] + 1, 1)
    X, Y = np.meshgrid(x, y)
    Z = mean_df.to_numpy()
    interp_step_size = 0.05
    # 保留原始的x和y
    total_num = int(1 + (mean_df.shape[0] - 1) / interp_step_size)
    x_new = np.linspace(1, mean_df.shape[0], total_num)
    y_new = np.linspace(1, mean_df.shape[0], total_num)
    X_new, Y_new = np.meshgrid(x_new, y_new)
    interp = RegularGridInterpolator((x, y), Z, method='linear', bounds_error=False, fill_value=None)
    Z_new = interp((X_new, Y_new))
    cmap = matplotlib.colormaps['vlag']
    cmaplist = np.array([cmap(i) for i in range(cmap.N)]) * 255
    surface = mlab.surf(X_new.transpose(), Y_new.transpose(), Z_new, warp_scale='auto', colormap='jet')
    # surface = mlab.mesh(X_new, Y_new, Z_new, colormap='jet')
    surface.module_manager.scalar_lut_manager.lut.table = cmaplist
    # points3d error
    # mlab.points3d(X_new.transpose().flatten(), Y_new.transpose().flatten(), Z_new.flatten(), std_df.to_numpy().flatten(), scale_mode='none', scale_factor=0.1, color=(1, 0, 0))
    mlab.savefig(filename=output_figure_path)
    # mlab.show()
    mlab.close()
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 8))
    # ax.view_init(elev=45, azim=-45)
    # cmap = 'vlag'
    # surface = ax.plot_surface(X_new, Y_new, Z_new, cmap=cmap, rstride=1, cstride=1)
    # ax.errorbar(X.flatten(), Y.flatten(), Z.flatten(), zerr=std_df.to_numpy().flatten(), fmt=',', ecolor='black',
    #             capsize=5, zorder=10)
    # ax.set(xticks=[0] + [*x], yticks=[0] + [*y],
    #        # zticks=[0, .2, .4, .6, .8, 1.0],
    #        xticklabels=[''] + [*mean_df.columns],
    #        yticklabels=[''] + [*mean_df.columns])
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('correlation')
    # # # show colorbar from data
    # # cbar = fig.colorbar(surface, shrink=0.2, aspect=10, location='right', pad=.1)
    # # tick_locator = ticker.MaxNLocator(nbins=6)
    # # cbar.locator = tick_locator
    # # cbar.update_ticks()
    # # cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # # cbar.set_ticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1.0'])
    # # # show colorbar from cmap
    # fig.colorbar(plt.cm.ScalarMappable(cmap=cmap), ax=ax, shrink=0.2, aspect=10, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
    #              location='top', pad=0.05)
    # # if output_figure_path is not None:
    # #     plt.savefig(output_figure_path)
    # plt.show()
    # plt.close()


def plot_polar_with_errorbar(mean_tsv, std_tsv, output_figure_path):
    mean_df = pd.read_table(mean_tsv, index_col=0)
    std_df = pd.read_table(std_tsv, index_col=0)
    # !endpoint必须是false, 因为polar是圆型的, 否则第一个值和最后一个值重合
    x = np.linspace(0, 2 * np.pi, mean_df.shape[0], endpoint=False)  # theta: Angle values
    y = np.linspace(0, 1, mean_df.shape[0])  # r: Radius values
    X, Y = np.meshgrid(x, y)
    Z = mean_df.to_numpy()

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 8))
    cmap = 'vlag'
    mesh = ax.pcolormesh(X, Y, Z, cmap=cmap)
    plt.colorbar(mesh, ax=ax, shrink=0.6, aspect=20, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
                 location='right', pad=0.1)
    ax.errorbar(X.flatten(), Y.flatten(), yerr=std_df.to_numpy().flatten(), fmt=',', ecolor='black', capsize=5)
    ax.set(xticks=[*x], yticks=[*y],
           xticklabels=[*mean_df.columns],
           yticklabels=[*mean_df.columns])
    ax.tick_params(axis='x', which='major', pad=15)
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


def plot_clustermap_with_errorbar(mean_tsv, std_tsv, output_figure_path=None):
    spearman_corr = pd.read_table(mean_tsv, index_col=0)
    std_df = pd.read_table(std_tsv, index_col=0)
    sns.set(font_scale=1.2)
    cluster_grid = sns.clustermap(spearman_corr, vmin=0, vmax=1, annot=False, cmap="vlag", fmt='.3f',
                                  annot_kws={'ha': 'right', 'va': 'center_baseline'},
                                  cbar_kws=dict(ticks=[0, .2, .4, .6, .8, 1.0]))
    ax = cluster_grid.ax_heatmap
    # Plot error bars
    # clustermap 会进行cluster操作, 给定数据和图上数据顺序是不一样的! 所以需要通过tick查找CI
    xtool = [xlabel.get_text() for xlabel in ax.get_xticklabels()]
    ytool = [ylabel.get_text() for ylabel in ax.get_yticklabels()]
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    for xidx, x in enumerate(xticks):
        i = x - 0.5
        for yidx, y in enumerate(yticks):
            # cell center
            j = y - 0.5
            ax.text(i + 0.22, j + 0.5, '{:.3f}'.format(spearman_corr.loc[xtool[xidx], ytool[yidx]]), ha='left',
                    va='center_baseline')
            ax.errorbar(i + 0.7, j + 0.5, yerr=std_df.loc[xtool[xidx], ytool[yidx]], fmt=',', ecolor='black', capsize=5)
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


def plot_clustermap_with_ci(mean_tsv, ci_lower_tsv, ci_upper_tsv, output_figure_path=None):
    spearman_corr = pd.read_table(mean_tsv, index_col=0)
    ci_lower = pd.read_table(ci_lower_tsv, index_col=0)
    ci_upper = pd.read_table(ci_upper_tsv, index_col=0)
    sns.set(font_scale=1.2)
    # 这里面fmt, ha, va左右是反的
    cluster_grid = sns.clustermap(spearman_corr, vmin=0, vmax=1, annot=False, cmap='vlag', fmt='.3f',
                                  annot_kws={'ha': 'right', 'va': 'center_baseline'},
                                  cbar_kws=dict(ticks=[0, .2, .4, .6, .8, 1.0]))
    ax = cluster_grid.ax_heatmap
    # clustermap 会进行cluster操作, 给定数据和图上数据顺序是不一样的! 所以需要通过tick查找CI
    xtool = [xlabel.get_text() for xlabel in ax.get_xticklabels()]
    ytool = [ylabel.get_text() for ylabel in ax.get_yticklabels()]
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    cmap = matplotlib.colormaps['vlag']
    for xidx, x in enumerate(xticks):
        i = x - 0.5
        for yidx, y in enumerate(yticks):
            # cell center
            j = y - 0.5
            ax.text(i + 0.22, j + 0.5, '{:.3f}'.format(spearman_corr.loc[xtool[xidx], ytool[yidx]]), ha='left',
                    va='center_baseline')
            ax.scatter(i + 0.7, j + 0.4, color=cmap(ci_upper.loc[xtool[xidx], ytool[yidx]]))
            ax.scatter(i + 0.7, j + 0.6, color=cmap(ci_lower.loc[xtool[xidx], ytool[yidx]]))
    if output_figure_path is not None:
        plt.savefig(output_figure_path)
    # plt.show()
    plt.close()


if __name__ == '__main__':
    # gwas_no_tissue('/Users/haiyue.meng/Downloads/real_gwas/')
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/Chronotype')
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl')
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/SleepDuration')
    # download_from_aws('/Users/haiyue.meng/Downloads/idl.txt', '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl')
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl')
    # download_from_aws('/Users/haiyue.meng/Downloads/ldl.txt', '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl')
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl')

    # ----------------------------
    # gwas_no_tissue('/Users/haiyue.meng/Downloads/real_gwas/')
    # calc_overview_data('/Users/haiyue.meng/Downloads/real_gwas', 'GCST')
    plot_surface_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/mean.tsv',
                               '/Users/haiyue.meng/Downloads/real_gwas/std.tsv',
                               '/Users/haiyue.meng/Downloads/real_gwas/3d_heat_with_error.png')
    #
    # plot_polar_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/mean.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/std.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/polar_with_error.png')
    #
    # plot_clustermap_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/mean.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/std.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/cluster_with_err.png')
    #
    # plot_clustermap_with_ci('/Users/haiyue.meng/Downloads/real_gwas/mean.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/ci_lower.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/ci_upper.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/cluster_with_ci.png')
    #
    # plot_mayavi_surface_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/mean.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/std.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/mayavi_3d_heat_with_error.png')

    # ----------------------------
    # download_from_28('/Users/haiyue.meng/Downloads/real_gwas_height.txt',
    #                  '/Users/haiyue.meng/Downloads/real_gwas/height')
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/height')
    #
    # calc_overview_data('/Users/haiyue.meng/Downloads/real_gwas/height')
    #
    # plot_surface_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/height/mean.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/height/std.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/height/3d_heat_with_error.png')
    #
    # plot_polar_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/height/mean.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/height/std.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/height/polar_with_error.png')
    #
    # plot_clustermap_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/height/mean.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/height/std.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/height/cluster_with_err.png')
    #
    # plot_clustermap_with_ci('/Users/haiyue.meng/Downloads/real_gwas/height/mean.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/height/ci_lower.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/height/ci_upper.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/height/cluster_with_ci.png')

    # ----------------------------
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl')
    #
    # calc_overview_data('/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl')
    #
    # plot_surface_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/mean.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/std.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/3d_heat_with_error.png')
    #
    # plot_polar_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/mean.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/std.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/polar_with_error.png')
    #
    # plot_clustermap_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/mean.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/std.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/cluster_with_err.png')
    #
    # plot_clustermap_with_ci('/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/mean.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/ci_lower.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/ci_upper.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_ldl/cluster_with_ci.png')

    # ----------------------------
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl')
    #
    # calc_overview_data('/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl')
    #
    # plot_surface_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/mean.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/std.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/3d_heat_with_error.png')
    #
    # plot_polar_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/mean.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/std.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/polar_with_error.png')
    #
    # plot_clustermap_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/mean.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/std.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/cluster_with_err.png')
    #
    # plot_clustermap_with_ci('/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/mean.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/ci_lower.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/ci_upper.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_idl/cluster_with_ci.png')

    # # ----------------------------
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl')
    #
    # calc_overview_data('/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl')
    #
    # plot_surface_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/mean.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/std.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/3d_heat_with_error.png')
    #
    # plot_polar_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/mean.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/std.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/polar_with_error.png')
    #
    # plot_clustermap_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/mean.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/std.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/cluster_with_err.png')
    #
    # plot_clustermap_with_ci('/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/mean.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/ci_lower.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/ci_upper.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/magnetic_hdl/cluster_with_ci.png')

    # # ------------no complete reports----------------
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/SleepDuration')
    #
    # calc_overview_data('/Users/haiyue.meng/Downloads/real_gwas/SleepDuration')
    #
    # plot_surface_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/mean.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/std.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/3d_heat_with_error.png')
    #
    # plot_polar_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/mean.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/std.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/polar_with_error.png')
    #
    # plot_clustermap_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/mean.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/std.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/cluster_with_err.png')
    #
    # plot_clustermap_with_ci('/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/mean.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/ci_lower.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/ci_upper.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/SleepDuration/cluster_with_ci.png')

    # # ----------------------------
    # gwas_tissue('/Users/haiyue.meng/Downloads/real_gwas/Chronotype')
    #
    # calc_overview_data('/Users/haiyue.meng/Downloads/real_gwas/Chronotype')
    #
    # plot_surface_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/Chronotype/mean.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/std.tsv',
    #                            '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/3d_heat_with_error.png')
    #
    # plot_polar_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/Chronotype/mean.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/std.tsv',
    #                          '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/polar_with_error.png')
    #
    # plot_clustermap_with_errorbar('/Users/haiyue.meng/Downloads/real_gwas/Chronotype/mean.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/std.tsv',
    #                               '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/cluster_with_err.png')
    #
    # plot_clustermap_with_ci('/Users/haiyue.meng/Downloads/real_gwas/Chronotype/mean.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/ci_lower.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/ci_upper.tsv',
    #                         '/Users/haiyue.meng/Downloads/real_gwas/Chronotype/cluster_with_ci.png')
