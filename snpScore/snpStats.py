import sys
import delegator
import numpy as np
import pandas as pd
from io import StringIO
from functools import reduce
from pybedtools import BedTool
from pathlib import PurePath, Path

script_dir = Path(__file__).parent
SNP_SCORE_PLOT = script_dir / 'snpScorePlot.R'
OFFSET = 1e-05
GROUPS = ('mutant', 'wild', 'mutant_parent', 'wild_parent', 'background')


def load_hd5(hd5_file, index_cols=None):
    table_name = PurePath(hd5_file).stem
    table_df = pd.read_hdf(hd5_file, table_name)
    if index_cols is not None:
        table_df = table_df.set_index(index_cols)
    return table_df


def groups_from_params(params):
    groups_list = []
    for group_i in GROUPS:
        if params.get(group_i):
            groups_list.extend([group_i] * len(params[group_i]))
    return ','.join(groups_list)


def vcfs_from_params(params, data_dir):
    data_list = []
    for group_i in GROUPS:
        if params.get(group_i) is None:
            continue
        for each_sample in params[group_i]:
            sample_file = f'{data_dir}/{each_sample}.h5'
            data_list.append(sample_file)
    return ','.join(data_list)


def outdir_suffix_from_params(params):
    outdir_suffix_suffix = []
    for group_i in GROUPS:
        group_i_list = []
        if params.get(group_i) is None:
            continue
        sample_list = sorted(params[group_i])
        for sample_i in sample_list:
            sample_i_id = sample_i.split('.')[0]
            group_i_list.append(sample_i_id.replace('0', ''))
        outdir_suffix_suffix.append('_'.join(group_i_list))
    return '/'.join(outdir_suffix_suffix)


def filter_alt_stat(alt_freq_stat_df, cols, freq):
    if len(cols) == 2:
        if len(alt_freq_stat_df.columns.intersection(pd.Index(cols))) == 2:
            group_i, group_j = cols
            mask1 = alt_freq_stat_df.loc[:, group_i] <= freq[group_i][0]
            mask2 = alt_freq_stat_df.loc[:, group_j] >= freq[group_j][1]
            filter_df1 = alt_freq_stat_df[mask1 & mask2]
            mask3 = alt_freq_stat_df.loc[:, group_i] >= freq[group_i][1]
            mask4 = alt_freq_stat_df.loc[:, group_j] <= freq[group_j][0]
            filter_df2 = alt_freq_stat_df[mask3 & mask4]
            alt_freq_stat_df = pd.concat([filter_df1, filter_df2],
                                         sort=False)
            alt_freq_stat_df = alt_freq_stat_df.sort_index()
    else:
        if cols[0] in alt_freq_stat_df.columns:
            mask1 = alt_freq_stat_df.loc[:, cols[0]] <= freq[cols[0]][0]
            mask2 = alt_freq_stat_df.loc[:, cols[0]] >= freq[cols[0]][1]
            alt_freq_stat_df = alt_freq_stat_df[mask1 | mask2]
    return alt_freq_stat_df


def filter_snp(alt_freq_stat_df, freq_dict, filter_label, filter_freq_stats):
    alt_freq_stat_filter_df = alt_freq_stat_df.copy()
    for n, cols_i in enumerate([['mutant', 'wild'],
                                ['mutant_parent', 'wild_parent'],
                                ['background']]):
        alt_freq_stat_filter_df = filter_alt_stat(alt_freq_stat_filter_df,
                                                  cols_i,
                                                  freq_dict)
    alt_freq_stat_filter_df.to_csv(filter_freq_stats)
    return alt_freq_stat_filter_df


def slidewindow(obj, window, step):
    for i in range(0, len(obj), step):
        yield obj[i:i + window]


def make_snp_number_windows(stat_df, group_label, window, step, outdir):
    snp_num_window_file = outdir / \
        f'{group_label}.snp_num.window.w{window}.s{step}.bed'
    if not snp_num_window_file.exists():
        snp_num_window_list = []
        for slidewindow_i in slidewindow(stat_df.index, window, step):
            chrom = stat_df.Chr[slidewindow_i].unique()
            if len(chrom) == 1:
                score_chrom = chrom[0]
            else:
                continue
            start = stat_df.Pos[slidewindow_i[0]] - 1
            end = stat_df.Pos[slidewindow_i[-1]]
            snp_num_window_list.append([score_chrom, start, end])
        snp_num_window_df = pd.DataFrame(snp_num_window_list,
                                         columns=['Chrom', 'Start', 'End'])
        snp_num_window_df.to_csv(
            snp_num_window_file, sep='\t', index=False, header=False)
    return snp_num_window_file


def make_genome_windows(chr_size, window, step, outdir):
    window_file = outdir / f'genome.window.w{window}.s{step}.bed'
    if (not window_file.exists()):
        cmd = ('bedtools makewindows '
               '-g {chr_size} -w {window} -s {step}'
               ' > {window_file}'.format(**locals()))

        delegator.run(cmd)
        window_df = pd.read_csv(window_file, sep='\t', header=None,
                                names=['Chrom', 'Start', 'End'])
        window_df.to_csv(window_file, sep='\t', index=False,
                         header=False)
    return window_file


def snp_freq_by_window(stat_df,
                       group_label,
                       window_file,
                       outdir):
    groups = stat_df.columns[3:]
    if 'background' in groups:
        groups = groups.drop('background')
    alt_freq_stat_bed = outdir / f'{group_label}.snp.plot.bed'
    if not alt_freq_stat_bed.exists():
        alt_freq_stat_df = stat_df.copy()
        alt_freq_stat_df.loc[:, 'start'] = alt_freq_stat_df.Pos - 1
        bed_cols = ['Chr', 'start', 'Pos']
        bed_cols.extend(groups)
        alt_freq_stat_df.to_csv(alt_freq_stat_bed, sep='\t',
                                columns=bed_cols, header=None, index=False)
    window_bed = BedTool(str(window_file))
    snp_bed = BedTool(str(alt_freq_stat_bed))
    intersect_obj = window_bed.intersect(snp_bed, sorted=True, wo=True)
    intersect_obj_cols = ['Chrom', 'Start', 'End']
    intersect_obj_cols.extend(['snp_Chrom', 'snp_start', 'snp_end'])
    intersect_obj_cols.extend(groups)
    intersect_obj_cols.append('overlap')
    intersect_str = StringIO(str(intersect_obj))
    intersect_df = pd.read_csv(intersect_str, sep='\t', header=None,
                               names=intersect_obj_cols)
    intersect_df.drop(['snp_Chrom', 'snp_start', 'snp_end',
                       'overlap'], axis=1, inplace=True)
    return intersect_df


def log_varscore(row):
    return np.power(-np.log10(reduce(lambda a, b:  a * b, row)), 10)


def cal_score(intersect_df, method='var', min_snp_num=10):
    varscore_size_df = intersect_df.groupby(
        ['Chrom', 'Start', 'End']).size()
    mask = varscore_size_df >= min_snp_num
    if method == 'var':
        varscore_df = intersect_df.groupby(
            ['Chrom', 'Start', 'End']).agg(
                lambda x: np.var(x))
    # elif method == 'est':
    #     for n, group_i in enumerate(groups):
    #         intersect_df.loc[:, group_i] = intersect_df.loc[
    #             :, group_i] - est[n]

    #     varscore_df = intersect_df.groupby(
    #         ['Chrom', 'Start', 'End', 'Label']).agg(
    #             lambda x: np.average(np.power(x, 2))
    #     )
    elif method == 'snp_index':
        varscore_df = intersect_df.groupby(
            ['Chrom', 'Start', 'End']).agg('mean')
    else:
        sys.exit('Wrong analysis method.')
    varscore_df = varscore_df[mask]
    if method == 'snp_index':
        group0, group1 = varscore_df.columns[0:2]
        varscore_df.loc[:, 'snp_score'] = varscore_df.loc[:, group0] - \
            varscore_df.loc[:, group1]
    else:
        varscore_df = varscore_df.applymap(
            lambda x: x if x >= OFFSET else OFFSET)
        varscore_df.loc[:, 'snp_score'] = varscore_df.apply(
            log_varscore, axis=1)
    return varscore_df


def score_plot(score_file, method):
    out_prefix = score_file.with_suffix('.plot')
    plot_name = score_file.stem
    if method in ['var', 'est', 'density']:
        out_plot = score_file.with_suffix('.plot.jpg')
    elif method == 'snp_index':
        out_plot = score_file.with_suffix('')
        out_prefix = out_plot
    else:
        sys.exit(f'Wrong snp score method [{method}]!')
    cmd = ('Rscript {SNP_SCORE_PLOT} '
           '--input {score_file} '
           '--output {out_prefix} '
           '--plot_type {method}'.format(
               SNP_SCORE_PLOT=SNP_SCORE_PLOT,
               score_file=score_file,
               out_prefix=out_prefix,
               method=method
           ))
    if not out_plot.exists():
        return cmd
    else:
        None
