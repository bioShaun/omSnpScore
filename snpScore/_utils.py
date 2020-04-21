import re
import asyncio
import shutil
import jinja2
import numpy as np
import pandas as pd
from io import StringIO
from loguru import logger
from decimal import Decimal, getcontext
from pathlib import Path, PurePath
from functools import reduce
from pybedtools import BedTool
from ._var import GROUPS, REF_FREQ, ALT_FREQ, OFFSET
from ._var import SnpGroup, SnpRep, SnpGroupFreq
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


class DuplicatedRecord(Exception):
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


def sample_and_group_for_web(parameters_obj):
    group_names = parameters_obj.get('group_names')
    if group_names:
        pass
    else:
        group_names = GROUPS
    sample_list = []
    group_list = []
    for group_i in group_names:
        if parameters_obj.get(group_i):
            group_list.extend([group_i] * len(parameters_obj[group_i]))
            sample_list.extend(parameters_obj[group_i])
    if len(group_names) == 1:
        if group_names[0] == 'not_a_group_id':
            group_list = sample_list[:]
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


def filter_snp(alt_freq_stat_df,
               freq_dict,
               filter_freq_stats,
               filter_method="nonsymmetrical"):
    alt_rep_df = alt_freq_stat_df.copy()
    for member in SnpGroupFreq.__members__.values():
        if member.value not in alt_rep_df.columns:
            continue
        ref_cut, alt_cut = freq_dict[getattr(SnpGroup, member.name).value]
        alt_rep_df.loc[:, member.value] = [
            snpfreq2rep(alt_freq_i, alt_cut, ref_cut)
            for alt_freq_i in alt_rep_df.loc[:, member.value]
        ]
    # step1 remove non ref/alt
    alt_rep_df.dropna(inplace=True)
    # step2 child equal to parent or parent unkown
    alt_rep_df = equal2parent(alt_rep_df, SnpGroupFreq.mut.value,
                              SnpGroupFreq.mut_pa.value)
    alt_rep_df = equal2parent(alt_rep_df, SnpGroupFreq.wild.value,
                              SnpGroupFreq.wild_pa.value)
    if filter_method == 'nonsymmetrical':
        # step3 mutant not equal to wild
        mask = alt_rep_df.loc[:, SnpGroupFreq.mut.value] \
            != alt_rep_df.loc[:, SnpGroupFreq.wild.value]
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
    groups = stat_df.columns[3:]
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
    score_s = np.power(-np.log10(reduce(lambda a, b: a * b, row)), 10)
    return score_s.astype('int')


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
        intersect_df.loc[:, SnpGroupFreq.mut.value] = mut_freq
        intersect_df.loc[:, SnpGroupFreq.wild.value] = wild_freq
        return intersect_df


def cal_score(intersect_df, freq_dict, method='var', min_snp_num=5):
    stats_cols = [
        'Chrom', 'Start', 'End', SnpGroupFreq.mut.value,
        SnpGroupFreq.wild.value
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
    varscore_df.drop([SnpGroupFreq.mut.value, SnpGroupFreq.wild.value],
                     axis=1,
                     inplace=True)
    return varscore_df


def score_plot(score_file,
               method,
               plot_title="",
               chr_size="",
               platform="local"):
    out_prefix = score_file.with_suffix('.plot')
    if method in ['var', 'est_mut_alt', 'est_mut_ref', 'density']:
        out_plot = score_file.with_suffix('.plot.jpg')
    elif method == 'snp_index':
        out_plot = score_file.with_suffix('')
        out_prefix = out_plot
    elif method in ['ED', 'snpIndex', 'Gprime']:
        out_prefix = score_file.parent / plot_title
        if method == 'Gprime':
            out_plot = score_file.with_suffix(f'.{method}.plot.jpg')
        else:
            out_plot = score_file.with_suffix(f'.{method}.plot.png')
    else:
        raise UnsupportedPlot(method)
    cmd = (f'Rscript {SNP_SCORE_PLOT} '
           f'--input {score_file} '
           f'--output {out_prefix} '
           f'--plot_type {method} '
           f'--title {plot_title} '
           f'--chr_size {chr_size}')
    if platform == 'web':
        cmd = f'{cmd} --web'
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


def valid_grp(grp_list):
    valid_grp_list = []
    for n, snp_group_i in enumerate(SnpGroup.__members__.items()):
        _, member = snp_group_i
        if member.value in grp_list:
            valid_grp_list.append(member.value)
    return valid_grp_list


def abbr_sample_id(sample_id):
    pattern = re.compile('(TC[A-Z])0+(\w+)')
    if pattern.match(sample_id):
        pre, suf = pattern.match(sample_id).groups()
        return f'{pre}{suf}'
    return sample_id


def outdir_suffix_from_params(params):
    outdir_suffix_suffix = []
    for group_i in GROUPS:
        group_i_list = []
        if params.get(group_i) is None:
            continue
        sample_list = sorted(params[group_i])
        for sample_i in sample_list:
            sample_i_id = sample_i.split('.')[0]
            group_i_list.append(abbr_sample_id(sample_i_id))
        outdir_suffix_suffix.append('_'.join(group_i_list))
    return '/'.join(outdir_suffix_suffix)


def replace_outdir(args, chrom):
    arg_list = []
    flag = False
    for arg_i in args:
        if flag:
            arg_i = f'{arg_i}/split/{chrom}'
            flag = False
        if arg_i == '-o' or arg_i == '--outdir':
            flag = True
        arg_list.append(arg_i)
    return arg_list


def wrap_param_arg(args):
    flag = False
    for arg_i in args:
        if flag:
            arg_i = f"'{arg_i}'"
            flag = False
        if arg_i == '-p' or arg_i == '--parameters':
            flag = True
        yield arg_i


def merge_split_file(file_dir, file_pattern, sortby=None):
    pattern_file = Path(file_dir).glob(f'split/*/{file_pattern}')
    df_list = []
    for file_i in pattern_file:
        # if file_i.name not in exist_file_name:
        df_list.append(pd.read_csv(file_i))
    outfile = Path(file_dir) / file_pattern
    df = pd.concat(df_list)
    if sortby:
        df.sort_values(sortby, inplace=True)
    if not outfile.is_file():
        df.to_csv(outfile, index=False)
    return outfile


def gene2pos(gene_bed, genes):
    gene_bed_df = pd.read_csv(gene_bed,
                              sep='\t',
                              header=None,
                              names=['chrom', 'start', 'end'],
                              index_col=3)
    for gene_i in genes:
        if gene_i in gene_bed_df.index:
            gene_i_pos = gene_bed_df.loc[gene_i]
            yield f'{gene_i_pos.chrom}:{gene_i_pos.start}-{gene_i_pos.end}'


def printdf(df):
    col_str = [str(col_i) for col_i in df.columns]
    print('\t'.join(col_str))
    for index_i in df.index:
        line_i = [str(col_i) for col_i in df.loc[index_i]]
        print('\t'.join(line_i))


def freq2qtlseqr(snpfreq):
    snpfreq = Path(snpfreq)
    snpfreq_df = pd.read_csv(snpfreq)

    def rename_col(col_i):
        pos_map = {'Chr': 'CHROM', 'Pos': 'POS', 'Alt': 'ALT'}
        if pos_map.get(col_i):
            return pos_map[col_i]
        else:
            return re.sub(r"(\w+).(\w+).AD", r"AD_\2.\1", col_i)

    snpfreq_df.columns = [rename_col(col_i) for col_i in snpfreq_df.columns]
    qtlseqr_table = snpfreq.with_suffix('.qtlseqr.csv')
    snpfreq_df.to_csv(qtlseqr_table, index=False)
    return qtlseqr_table


def circos_suffix(varscore_prefix, qtlseq_prefix):
    varscore_prefix = varscore_prefix.replace('mutant', 'mut')
    varscore_prefix = varscore_prefix.replace('wild', 'wd')
    varscore_prefix = varscore_prefix.replace('symmetrical', 'sym')
    varscore_prefix = varscore_prefix.replace('snp_num.window.', '')
    qtlseq_prefix = qtlseq_prefix.replace('window_', '')
    qtlseq_prefix = qtlseq_prefix.replace('popStru_', '')
    qtlseq_prefix = qtlseq_prefix.replace('refFreq_', '')
    qtlseq_prefix = qtlseq_prefix.replace('minDepth_', '')
    return f'{varscore_prefix}-{qtlseq_prefix}'


def circos_cfg(circos_prefix):
    circos_prefix.mkdir(parents=True, exist_ok=True)
    circos_path = circos_prefix.parent
    circos_file = f'{circos_prefix.name}.circos.png'
    # jinja2 load template
    PLOT_DIR = PurePath(__file__).parent / 'plot'
    JINJA_ENV = jinja2.Environment(loader=jinja2.FileSystemLoader(
        searchpath=f'{PLOT_DIR}'))
    CIRCOS_CONF = JINJA_ENV.get_template('circos.conf')
    cfgObj = CIRCOS_CONF.render({
        'circos_prefix': circos_prefix,
        'circos_path': circos_path,
        'circos_file': circos_file
    })
    cfgFile = circos_prefix / 'circos.conf'
    with open(cfgFile, 'w') as file_inf:
        file_inf.write(cfgObj)
    return circos_path / circos_file


def circos_plot(varScore_csv, qtlseqr_ed_csv, snp_freq_csv, out_prefix):
    circos_sh = PurePath(__file__).parent / 'plot' / 'data2circos.sh'
    cmd = f'sh {circos_sh} {varScore_csv} {qtlseqr_ed_csv} {snp_freq_csv} {out_prefix}'
    return cmd
