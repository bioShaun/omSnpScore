import asyncio
import shutil
import numpy as np
import pandas as pd
from io import StringIO
from loguru import logger
from decimal import Decimal, getcontext
from functools import reduce
from pybedtools import BedTool
from ._var import GROUPS, REF_FREQ, ALT_FREQ, OFFSET
from ._var import SnpGroup, SnpRep
from ._var import SNP_SCORE_PLOT

getcontext().prec = 3


class AppNotFound(Exception):
    pass


class SampleFileNotMatch(Exception):
    pass


class UnsupportedFormat(Exception):
    pass


class UnsupportedScoreMethod(Exception):
    pass


class UnsupportedPlot(Exception):
    pass


async def async_sh_job(cmd, sema):
    with (await sema):
        p = await asyncio.create_subprocess_shell(
            cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT)
        return (await p.communicate())[0].splitlines()


def async_batch_sh_jobs(cmd_list, thread=2):
    if cmd_list:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        loop = asyncio.get_event_loop()
        semaphore = asyncio.Semaphore(thread)
        coro_list = [async_sh_job(cmd, semaphore) for cmd in cmd_list]
        try:
            loop.run_until_complete(asyncio.wait(coro_list))
        finally:
            loop.close()


def check_app(app_name):
    if shutil.which(app_name) is None:
        raise AppNotFound(app_name)


def sample_and_group(*args):
    sample_list = []
    group_list = []
    for n, s_i in enumerate(args):
        if not s_i:
            continue
        s_i_list = s_i.split(',')
        g_i_list = [GROUPS[n]] * len(s_i_list)
        sample_list.extend(s_i_list)
        group_list.extend(g_i_list)
    return sample_list, group_list


def freq_accordance(freq_a, freq_b, message, equal=True):
    non_legal_freq = (freq_a != 0.5) and (freq_b != 0.5)
    assert non_legal_freq, 'alt freqency should not equal to 0.5.'
    if freq_a is None:
        return True
    elif freq_b is None:
        return True
    else:
        freq_direct = (freq_a - 0.5) * (freq_b - 0.5)
        freq_accord = freq_direct > 0
        flag_a = freq_accord and equal
        flag_b = not (freq_accord or equal)
        assert flag_a or flag_b, message


def alt_ref_cut(freq, is_ref=True):
    if freq is None:
        return -np.inf, np.inf
    if is_ref:
        ref_cut = freq
        alt_cut = float(Decimal(1) - Decimal(freq))
    else:
        if freq > 0.5:
            ref_cut = -np.inf
            alt_cut = freq
        else:
            ref_cut = freq
            alt_cut = np.inf
    return ref_cut, alt_cut


def snpfreq2rep(alt_freq, alt_cut, ref_cut):
    if alt_freq >= alt_cut:
        return SnpRep.alt.value
    elif alt_freq <= ref_cut:
        return SnpRep.ref.value
    elif pd.isna(alt_freq):
        return SnpRep.unkown.value
    else:
        return np.nan


def equal2parent(snp_rep_df, child, parent):
    if parent in snp_rep_df.columns:
        mask1 = snp_rep_df.loc[:, child] == snp_rep_df.loc[:, parent]
        mask2 = snp_rep_df.loc[:, parent] == SnpRep.unkown
        return snp_rep_df[mask1 | mask2]
    else:
        return snp_rep_df


def filter_snp(alt_freq_stat_df, freq_dict, filter_label, filter_freq_stats):
    alt_rep_df = alt_freq_stat_df.copy()
    for member in SnpGroup.__members__.values():
        if member.value not in alt_rep_df.columns:
            continue
        ref_cut, alt_cut = freq_dict[member.value]
        alt_rep_df.loc[:, member.value] = [
            snpfreq2rep(alt_freq_i, alt_cut, ref_cut)
            for alt_freq_i in alt_rep_df.loc[:, member.value]
        ]
    # step1 remove non ref/alt
    alt_rep_df.dropna(inplace=True)
    # step2 child equal to parent or parent unkown
    alt_rep_df = equal2parent(alt_rep_df, SnpGroup.mut.value,
                              SnpGroup.mut_pa.value)
    alt_rep_df = equal2parent(alt_rep_df, SnpGroup.wild.value,
                              SnpGroup.wild_pa.value)
    # step3 mutant not equal to wild
    mask = alt_rep_df.loc[:, SnpGroup.mut.value] \
        != alt_rep_df.loc[:, SnpGroup.wild.value]
    alt_rep_df = alt_rep_df[mask]
    alt_freq_stat_filter_df = alt_freq_stat_df.loc[alt_rep_df.index]
    alt_freq_stat_filter_df.to_csv(filter_freq_stats, index=False)
    return alt_freq_stat_filter_df


def slidewindow(obj, window, step):
    for i in range(0, len(obj), step):
        yield obj[i:i + window]


def make_snp_number_windows(stat_df, group_label, window, step, outdir):
    snp_num_window_file = outdir / \
        f'{group_label}.snp_num.window.w{window}.s{step}.bed'
    if not snp_num_window_file.is_file():
        logger.info(
            'Making snp number slidewindow bed file windows {w} step {s}...',
            w=window,
            s=step)
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
        snp_num_window_df.to_csv(snp_num_window_file,
                                 sep='\t',
                                 index=False,
                                 header=False)
    return snp_num_window_file


def snp_freq_by_window(stat_df, group_label, window_file, outdir):
    groups = [SnpGroup.mut.value, SnpGroup.wild.value]
    alt_freq_stat_bed = outdir / f'{group_label}.snp.plot.bed'
    if not alt_freq_stat_bed.is_file():
        alt_freq_stat_df = stat_df.copy()
        alt_freq_stat_df.loc[:, 'start'] = alt_freq_stat_df.Pos - 1
        bed_cols = ['Chr', 'start', 'Pos', 'Alt']
        bed_cols.extend(groups)
        alt_freq_stat_df.to_csv(alt_freq_stat_bed,
                                sep='\t',
                                columns=bed_cols,
                                header=None,
                                index=False)
    window_bed = BedTool(str(window_file))
    snp_bed = BedTool(str(alt_freq_stat_bed))
    intersect_obj = window_bed.intersect(snp_bed, sorted=True, wo=True)
    intersect_obj_cols = ['Chrom', 'Start', 'End']
    intersect_obj_cols.extend(['snp_Chr', 'snp_start', 'Pos', 'Alt'])
    intersect_obj_cols.extend(groups)
    intersect_obj_cols.append('overlap')
    intersect_str = StringIO(str(intersect_obj))
    intersect_df = pd.read_csv(intersect_str,
                               sep='\t',
                               header=None,
                               names=intersect_obj_cols)
    intersect_df.drop(['snp_Chr', 'snp_start', 'overlap'],
                      axis=1,
                      inplace=True)
    return intersect_df


def log_varscore(row):
    return np.power(-np.log10(reduce(lambda a, b: a * b, row)), 10)


def mut_wild_ext_freq(intersect_df, freq_dict, mut='alt'):
    if mut == 'alt':
        mut_freq = freq_dict[SnpGroup.mut.value][1]
        wild_freq = freq_dict[SnpGroup.wild.value][0]
    else:
        mut_freq = freq_dict[SnpGroup.mut.value][0]
        wild_freq = freq_dict[SnpGroup.wild.value][1]
    if np.isinf(mut_freq) or np.isinf(wild_freq):
        return None
    else:
        intersect_df.loc[:, SnpGroup.mut.value] = mut_freq
        intersect_df.loc[:, SnpGroup.wild.value] = wild_freq
        return intersect_df


def cal_score(intersect_df, freq_dict, method='var', min_snp_num=5):
    stats_cols = [
        'Chrom', 'Start', 'End', SnpGroup.mut.value, SnpGroup.wild.value
    ]
    stats_df = intersect_df.loc[:, stats_cols]
    varscore_size_df = stats_df.groupby(['Chrom', 'Start', 'End']).size()
    mask = varscore_size_df >= min_snp_num
    if method == 'var':
        varscore_df = stats_df.groupby(['Chrom', 'Start',
                                        'End']).agg(lambda x: np.var(x))
    elif 'est' in method:
        stats_df = stats_df.set_index(['Chrom', 'Start', 'End'])
        alt_freq_df = stats_df.copy()
        mut_stat = method.split('_')[-1]
        mut_wild_ext_df = mut_wild_ext_freq(alt_freq_df,
                                            freq_dict,
                                            mut=mut_stat)
        if mut_wild_ext_df is None:
            return None
        else:
            stats_df = stats_df - mut_wild_ext_df
        varscore_df = stats_df.groupby(
            ['Chrom', 'Start',
             'End']).agg(lambda x: np.average(np.power(x, 2)))
    elif method == 'snp_index':
        varscore_df = stats_df.groupby(['Chrom', 'Start', 'End']).agg('mean')
    else:
        raise UnsupportedScoreMethod(method)
    varscore_df = varscore_df[mask]
    if method == 'snp_index':
        group0, group1 = varscore_df.columns[0:2]
        varscore_df.loc[:, 'snp_score'] = varscore_df.loc[:, group0] - \
            varscore_df.loc[:, group1]
    else:
        varscore_df = varscore_df.applymap(lambda x: x
                                           if x >= OFFSET else OFFSET)
        varscore_df.loc[:, 'snp_score'] = varscore_df.apply(log_varscore,
                                                            axis=1)
    return varscore_df


def score_plot(score_file, method):
    out_prefix = score_file.with_suffix('.plot')
    if method in ['var', 'est_mut_alt', 'est_mut_ref', 'density']:
        out_plot = score_file.with_suffix('.plot.jpg')
    elif method == 'snp_index':
        out_plot = score_file.with_suffix('')
        out_prefix = out_plot
    else:
        raise UnsupportedPlot(method)
    cmd = (f'Rscript {SNP_SCORE_PLOT} '
           f'--input {score_file} '
           f'--output {out_prefix} '
           f'--plot_type {method}')
    if not out_plot.exists():
        return cmd
    else:
        return None


def extract_snpeff_anno(anno_line):
    anno_stats = []
    fileds = (1, 3, 6, 9, 10)
    gene_anno = anno_line.split(';')[0]
    anno_line_stats = gene_anno.split(",")
    for annStr in anno_line_stats:
        annDetailArray = annStr.split("|")
        filed_stats = []
        for filled_i in fileds:
            filed_stats.append(annDetailArray[filled_i])
        anno_stats.append(filed_stats)
    zip_anno_stats = list(map(lambda x: '|'.join(x), zip(*anno_stats)))
    return zip_anno_stats


def split_dataframe_rows(df, column_selectors, row_delimiter):
    # we need to keep track of the ordering of the columns
    def _split_list_to_rows(row, row_accumulator, column_selector,
                            row_delimiter):
        split_rows = {}
        max_split = 0
        for column_selector in column_selectors:
            split_row = row[column_selector].split(row_delimiter)
            split_rows[column_selector] = split_row
            if len(split_row) > max_split:
                max_split = len(split_row)

        for i in range(max_split):
            new_row = row.to_dict()
            for column_selector in column_selectors:
                try:
                    new_row[column_selector] = split_rows[column_selector].pop(
                        0)
                except IndexError:
                    new_row[column_selector] = ''
            row_accumulator.append(new_row)

    new_rows = []
    df.apply(_split_list_to_rows,
             axis=1,
             args=(new_rows, column_selectors, row_delimiter))
    new_df = pd.DataFrame(new_rows, columns=df.columns)
    return new_df
