import fire
import gzip
import sys
import delegator
import gtfparse
import pandas as pd
import numpy as np
from pathlib import Path
from orderedset import OrderedSet
import datetime
from loguru import logger
from itertools import zip_longest
from functools import reduce
from pybedtools import BedTool
from io import StringIO


script_dir = Path(__file__).parent
SNP_SCORE_PLOT = script_dir / 'snpScorePlot.R'
SNP_INDEX_PLOT = script_dir / 'snpIndex.R'


def init_logger(outdir):
    logfile = outdir / 'log.txt'
    logger.add(logfile,
               format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}")


def extract_vcf_header(vcf):
    if vcf.suffix == '.gz':
        vcf_inf = gzip.open(vcf)
    else:
        vcf_inf = open(vcf)
    for eachline in vcf_inf:
        if vcf.suffix == '.gz':
            eachline = eachline.decode()
        if eachline[:6] == '#CHROM':
            return eachline.strip().split('\t')


def unique_snpeff_anno(snp_anno_obj):
    return '|'.join(OrderedSet(snp_anno_obj))


def extract_snpeff_anno(snpeff_col):
    anno_stats = []
    fileds = (1, 3, 9)
    anno_inf = sorted(snpeff_col.split(';'))
    anno_line = None
    for anno_i in anno_inf:
        if anno_i[:3] == 'ANN':
            anno_line = anno_i
            break
    if anno_line is None:
        sys.exit('No [ANN] record in vcf, check you use the right vcf file.')

    anno_line_stats = anno_line.split(",")
    for annStr in anno_line_stats:
        annDetailArray = annStr.split("|")
        filed_stats = []
        for filled_i in fileds:
            filed_stats.append(annDetailArray[filled_i])
        anno_stats.append(filed_stats)
    zip_anno_stats = list(map(unique_snpeff_anno, zip(*anno_stats)))
    return zip_anno_stats


def extractAlleFreq(reads_inf):
    if reads_inf.count('./.') > 0:
        return 0
    else:
        return reads_inf.split(":")[1]


def count_allele_num(x):
    return x.count(',') + 1


def allele_present(reads_inf):
    reads_inf_str_list = str(reads_inf).split(',')
    return '|'.join(reads_inf_str_list)


def allele_freq(reads_inf):
    reads_inf_str_list = str(reads_inf).split(',')
    reads_inf_int_list = [int(read_i) for read_i in reads_inf_str_list]
    total_reads = sum(reads_inf_int_list)
    if total_reads == 0:
        return 'NA'
    else:
        return min(reads_inf_int_list) / total_reads


def allele_direct(reads_inf):
    reads_inf_str_list = str(reads_inf).split(',')
    reads_inf_int_list = [int(read_i) for read_i in reads_inf_str_list]
    total_reads = sum(reads_inf_int_list)
    if total_reads == 0:
        return 'NA'
    else:
        return reads_inf_int_list.index(max(reads_inf_int_list))


def allele_depth(reads_inf):
    reads_inf_str_list = str(reads_inf).split(',')
    reads_inf_int_list = [int(read_i) for read_i in reads_inf_str_list]
    return sum(reads_inf_int_list)


def alt_allele_num(reads_inf):
    reads_inf_str_list = str(reads_inf).split('|')
    reads_inf_int_list = [int(read_i) for read_i in reads_inf_str_list]
    if len(reads_inf_int_list) == 1:
        return 0
    else:
        return reads_inf_int_list[1]


def alt_allele_freq(reads_inf):
    reads_inf_str_list = str(reads_inf).split('|')
    reads_inf_int_list = [int(read_i) for read_i in reads_inf_str_list]
    total_reads = sum(reads_inf_int_list)
    if total_reads == 0:
        return np.nan
    else:
        return 1 - (reads_inf_int_list[0] / total_reads)


def merge_allele_inf(row):
    row_list = []
    for col_i in row:
        col_i_reads = [int(each) for each in str(col_i).split(',')]
        row_list.append(col_i_reads)
    merge_row_list = zip_longest(*row_list)
    merge_row_str_list = [str(sum(filter(None, each)))
                          for each in merge_row_list]
    return ','.join(merge_row_str_list)


def slidewindow(obj, window, step):
    for i in range(0, len(obj), step):
        yield obj[i:i+window]


def make_snp_number_windows(stat_df, stat_label, window, step, outdir, force):
    snp_num_window_file = outdir / f'snp_num.{stat_label}.window.w{window}.s{step}.bed'
    if (not snp_num_window_file.exists()) or force:
        logger.info(
            'making {l} snp number slidewindow bed file windows {w} step {s}...',
            w=window, s=step, l=stat_label)
        snp_num_window_list = []
        for slidewindow_i in slidewindow(stat_df.index, window, step):
            chrom = stat_df.Chrom[slidewindow_i].unique()
            if len(chrom) == 1:
                score_chrom = chrom[0]
            else:
                continue
            start = stat_df.Pos[slidewindow_i[0]] - 1
            end = stat_df.Pos[slidewindow_i[-1]]
            snp_num_window_list.append([score_chrom, start, end])
        snp_num_window_df = pd.DataFrame(snp_num_window_list,
                                         columns=['Chrom', 'Start', 'End'])
        snp_num_window_df.loc[:, 'label'] = snp_num_window_df.index.map(lambda x: f'snp_group{x+1}')
        snp_num_window_df.to_csv(
            snp_num_window_file, sep='\t', index=False, header=False)
    return snp_num_window_file


def make_genome_windows(chr_file, window, step, outdir, force):
    window_file = outdir / f'genome.window.w{window}.s{step}.bed'
    if (not window_file.exists()) or force:
        logger.info(
            'making genome slidewindow bed file windows {w} step {s}...',
            w=window, s=step)
        cmd = f'bedtools makewindows -g {chr_file} -w {window} -s {step} > {window_file}'
        delegator.run(cmd)
        window_df = pd.read_csv(window_file, sep='\t', header=None,
                                names=['Chrom', 'Start', 'End'])
        window_df.loc[:, 'Label'] = window_df.apply(
            lambda x: '{0}:{1}-{2}'.format(
                x['Chrom'], x['Start'], x['End']), axis=1)
        window_df.to_csv(window_file, sep='\t', index=False,
                         header=False)
    return window_file


def make_gene_window(gtf, outdir, force, by='gene'):
    window_file = outdir / f'{by}.window.bed'
    if (not window_file.exists()) or force:
        logger.info(
            'making {by} bed file...',
            by=by)
        gtf_df = gtfparse.read_gtf(gtf)
        if by == 'gene':
            if 'gene' in gtf_df.feature.unique():
                pass
            else:
                if 'exon' not in gtf_df.feature.unique():
                    sys.exit('exon not in gtf file!')
                exon_gtf = gtf_df[gtf_df.feature == 'exon']
                gene_start = exon_gtf.groupby(['gene_id'])['start'].min() - 1
                gene_end = exon_gtf.groupby(['gene_id'])['end'].max()
                gene_chrom = exon_gtf.loc[
                    :, ['gene_id', 'seqname']
                ].drop_duplicates().set_index('gene_id')
                gene_df = pd.concat(
                    [gene_chrom, gene_start, gene_end],
                    axis=1, sort=False)
                gene_df = gene_df.reset_index()
                gene_df.sort_values(['seqname', 'start', 'end'], inplace=True)
                gene_df.to_csv(window_file, sep='\t',
                               index=False, header=False,
                               columns=['seqname', 'start', 'end', 'index'])
    return window_file


def varscore_by_window(stat_df,
                       window_file,
                       outdir,
                       stat_label,
                       by,
                       est=None, method='var',
                       offset=0.00001):
    logger.info('calculating score based on {stat_label} snp by {by} using {method}...',
                by=by, method=method, stat_label=stat_label)
    groups = stat_df.columns[3:-1]
    alt_freq_stat_bed = outdir / f'snp.{stat_label}.plot.bed'
    if not alt_freq_stat_bed.exists():
        alt_freq_stat_df = stat_df.copy()
        alt_freq_stat_df.loc[:, 'start'] = alt_freq_stat_df.Pos - 1
        bed_cols = ['Chrom', 'start', 'Pos']
        bed_cols.extend(groups)
        alt_freq_stat_df.to_csv(alt_freq_stat_bed, sep='\t',
                                columns=bed_cols, header=None, index=False)
    window_bed = BedTool(str(window_file))
    snp_bed = BedTool(str(alt_freq_stat_bed))
    intersect_obj = window_bed.intersect(snp_bed, sorted=True, wo=True)
    intersect_obj_cols = ['Chrom', 'Start', 'End', 'Label']
    intersect_obj_cols.extend(['snp_Chrom', 'snp_start', 'snp_end'])
    intersect_obj_cols.extend(groups)
    intersect_obj_cols.append('overlap')
    intersect_str = StringIO(str(intersect_obj))
    intersect_df = pd.read_csv(intersect_str, sep='\t', header=None,
                               names=intersect_obj_cols)
    intersect_df.drop(['snp_Chrom', 'snp_start', 'snp_end',
                       'overlap'], axis=1, inplace=True)
    varscore_size_df = intersect_df.groupby(
        ['Chrom', 'Start', 'End', 'Label']).size()
    if method == 'var':
        mask = varscore_size_df > 1
        varscore_df = intersect_df.groupby(
            ['Chrom', 'Start', 'End', 'Label']).agg(
                lambda x: np.var(x) + offset)
        varscore_df = varscore_df[mask]
    elif method == 'est':
        for n, group_i in enumerate(groups):
            intersect_df.loc[:, group_i] = intersect_df.loc[
                :, group_i] - est[n]
        varscore_df = intersect_df.groupby(
            ['Chrom', 'Start', 'End', 'Label']).agg(
                lambda x: np.average(np.power(x, 2))
        )
    elif method == 'snp_index':
        mask = varscore_size_df >= 10
        varscore_df = intersect_df.groupby(
            ['Chrom', 'Start', 'End', 'Label']).agg('mean')
        varscore_df = varscore_df[mask]
    else:
        sys.exit('Wrong analysis method.')
    return varscore_df.reset_index()


def log_varscore(row, offset=0.00001):
    row_vals = []
    for row_i in row:
        if row_i == 0:
            row_i = row_i + offset
        row_vals.append(row_i)
    return -np.log10(reduce(lambda a, b:  a * b, row_vals))


def score_plot(score_file, method, force):
    out_prefix = score_file.with_suffix('.plot')
    plot_name = score_file.stem
    if method in ['var', 'est', 'density']:
        out_plot = score_file.with_suffix('.plot.jpg')
        cmd = f'Rscript {SNP_SCORE_PLOT} --input {score_file} --output {out_prefix}'
    elif method == 'snp_index':
        out_plot = score_file.with_suffix('')
        cmd = f'Rscript {SNP_INDEX_PLOT} --snp_index_file {score_file} --plot_dir {out_plot}'
    else:
        sys.exit(f'Wrong snp score method [{method}]!')
    if not out_plot.exists() or force:
        logger.info('generating {} plot...', plot_name)
        r = delegator.run(cmd)
        return r.out
    else:
        logger.info('{} plot exists.', plot_name)


def vcf2df(vcf, outdir, force):
    out_table = outdir / 'snp.table.txt'
    if (not out_table.exists()) or force:
        vcf = Path(vcf)
        vcf_header = extract_vcf_header(vcf)
        logger.info('reading vcf data...')
        vcf_df = pd.read_csv(vcf,
                             compression='gzip', comment='#',
                             header=None, names=vcf_header,
                             sep='\t')
        out_table_df = vcf_df.loc[:, ['#CHROM', 'POS', 'REF', 'ALT']]
        out_table_df.columns = ["Chrom", "Pos", "Refer", "Alter"]
        # count allele number
        logger.info('counting allele number...')
        out_table_df.loc[:, 'AlleNum'] = out_table_df.Alter.map(
            count_allele_num)
        # snpeff_anno
        logger.info('extracting snpEff annotation...')
        snpeff_anno = list(vcf_df.INFO.map(extract_snpeff_anno))
        snpeff_anno_df = pd.DataFrame(snpeff_anno)
        snpeff_anno_df.columns = ['Feature', 'Gene', 'Alle']
        out_table_anno_df = pd.concat([out_table_df, snpeff_anno_df], axis=1)
        # reads number
        logger.info('extracting reads number...')
        samples = vcf_df.columns[9:]
        reads_number_df = vcf_df.loc[:, samples].applymap(extractAlleFreq)
        out_table_anno_num_df = pd.concat(
            [out_table_anno_df, reads_number_df], axis=1)
        # write snp table
        logger.info('writing snp table...')
        out_table_anno_num_df.to_csv(out_table, sep='\t', index=False)
    else:
        logger.info('snp table exists.')
        out_table_anno_num_df = pd.read_csv(out_table, sep='\t')
    return out_table_anno_num_df


def group_allele(out_table_anno_num_df, sample_group,
                 outdir, force):
    stats_table = outdir / 'snp.stats.txt'
    if (not stats_table.exists()) or force:
        out_table_anno_df = out_table_anno_num_df.iloc[:, 0:8]
        reads_number_df = out_table_anno_num_df.iloc[:, 8:]
        logger.info('merging allele info by group...')
        plot_header = ('mutant_bulk', 'wild_bulk',
                       'mutant_parent', 'wild_parent',
                       'background')
        analysis_groups = []
        group_allele_stats = []
        for n, sample_i in enumerate(sample_group):
            if sample_i is not None:
                group_df = pd.DataFrame(reads_number_df.loc[:, sample_i])
                group_allele_stats.append(
                    group_df.apply(merge_allele_inf, axis=1))
                analysis_groups.append(plot_header[n])
        group_allele_df = pd.concat(group_allele_stats, axis=1)
        group_allele_df.columns = analysis_groups
        logger.info('calculating group allele stats...')
        reads_present_df = group_allele_df.applymap(allele_present)
        reads_freq_df = group_allele_df.applymap(allele_freq)
        reads_freq_df.columns = [f'{group_i}_Freq' for group_i in analysis_groups]
        reads_direct_df = group_allele_df.applymap(allele_direct)
        reads_direct_df.columns = [f'{group_i}_Direction' for group_i in analysis_groups]
        reads_depth_df = group_allele_df.applymap(allele_depth)
        reads_depth_df.columns = [f'{group_i}_Depth' for group_i in analysis_groups]
        reads_stats_df = pd.concat([out_table_anno_df, reads_present_df, reads_freq_df,
                                    reads_direct_df, reads_depth_df], sort=False,
                                   axis=1)
        reads_stats_df.loc[:, 'TotalDepth'] = reads_depth_df.sum(1)
        logger.info('writing snp stats table...')
        reads_stats_df.to_csv(stats_table, sep='\t', index=False,
                              float_format='%.5f')
    else:
        logger.info('snp stats table exists.')
        reads_stats_df = pd.read_csv(stats_table, sep='\t')
    return reads_stats_df


def group_alt_freq(reads_stats_df, outdir, plot_min_depth, force):
    alt_freq_stat_file = outdir / 'snp.freq.txt'
    if (not alt_freq_stat_file.exists()) or force:
        group_num = int((len(reads_stats_df.columns) - 9) / 4)
        group_allele_df = reads_stats_df.iloc[:, 8:8+group_num]
        reads_depth_df = reads_stats_df.iloc[:, -group_num-1:-1]
        analysis_groups = group_allele_df.columns
        # generate plot data
        logger.info('generating snp alt freq data...')
        alt_freq_stats = []
        out_col = ['SNP', 'Chrom', 'Pos']
        out_col.extend(analysis_groups)
        alt_freq_stat_df = group_allele_df.applymap(alt_allele_freq)
        alt_freq_stat_pos_df = reads_stats_df.loc[:, [
            'Chrom', 'Pos', 'AlleNum']]
        alt_freq_stat_df = pd.concat(
            [alt_freq_stat_pos_df, alt_freq_stat_df], axis=1)
        alt_freq_stat_df.loc[:, 'SNP'] = alt_freq_stat_df.apply(
            lambda x: '{0}:{1}'.format(x['Chrom'], x['Pos']), axis=1)
        # filter freq stat
        passed_rows = reads_depth_df[
            reads_depth_df.min(1) >= plot_min_depth].index
        alt_freq_stat_df = alt_freq_stat_df.take(passed_rows)
        alt_freq_stat_df = alt_freq_stat_df[alt_freq_stat_df.AlleNum == 1]
        logger.info('writing plot data...')
        alt_freq_stat_df = alt_freq_stat_df.loc[:, out_col]
        alt_freq_stat_df.to_csv(alt_freq_stat_file, sep='\t', index=False,
                                float_format='%.5f')
    else:
        logger.info('snp alt freq data exists.')
        alt_freq_stat_df = pd.read_csv(alt_freq_stat_file, sep='\t')
    score_plot(alt_freq_stat_file, 'density', force)
    return alt_freq_stat_df


def filter_alt_stat(alt_freq_stat_df, cols, freq):
    if len(cols) == 2:
        if len(alt_freq_stat_df.columns.intersection(pd.Index(cols))) == 2:
            filter_dfs = []
            for cols_comb in [cols, cols[::-1]]:
                mask1 = alt_freq_stat_df.loc[:, cols_comb[0]] <= freq
                mask2 = alt_freq_stat_df.loc[:, cols_comb[1]] >= 1 - freq
                filter_dfs.append(alt_freq_stat_df[mask1 & mask2])
            alt_freq_stat_df = pd.concat(filter_dfs)
            alt_freq_stat_df = alt_freq_stat_df.sort_index()
    else:
        if cols[0] in alt_freq_stat_df.columns:
            mask1 = alt_freq_stat_df.loc[:, cols[0]] <= freq
            mask2 = alt_freq_stat_df.loc[:, cols[0]] >= 1 - freq
            alt_freq_stat_df = alt_freq_stat_df[mask1 | mask2]
    return alt_freq_stat_df


def filter_snp(alt_freq_stat_df, freq_order, outdir, force):
    filter_freq_stas = outdir / 'snp.freq.filter.txt'
    if (not filter_freq_stas.exists()) or force:
        logger.info('filtering snp ...')
        alt_freq_stat_filter_df = alt_freq_stat_df.copy()
        for n, cols_i in enumerate([['mutant_bulk', 'wild_bulk'],
                                    ['mutant_parent', 'wild_parent'],
                                    ['background']]):
            alt_freq_stat_filter_df = filter_alt_stat(alt_freq_stat_filter_df,
                                                      cols_i,
                                                      freq_order[n])
        alt_freq_stat_filter_df.to_csv(filter_freq_stas, sep='\t', index=False)
    else:
        logger.info('snp filtered.')
        alt_freq_stat_filter_df = pd.read_csv(filter_freq_stas, sep='\t')
    score_plot(filter_freq_stas, 'density', force)
    return alt_freq_stat_filter_df


def cal_varscore(alt_freq_stat_df, freq_est,
                 outdir, force,
                 by='snp_num', method='var',
                 snp_number_window=10,
                 snp_number_step=5,
                 genome_window=2000000,
                 genome_step=1000000,
                 stat_label='all',
                 chr_size=None,
                 gtf=None):
    varscore_file = outdir / f'snp.{stat_label}.{by}.{method}.score.txt'
    if (not varscore_file.exists()) or force:
        alt_freq_stat_pos_df = alt_freq_stat_df.iloc[:, 0:3]
        analysis_groups = alt_freq_stat_df.columns[3:]
        if by in 'snp_num':
            snp_number_window_file = make_snp_number_windows(
                alt_freq_stat_df, stat_label, snp_number_window,
                snp_number_step, outdir, force)
            varscore_df = varscore_by_window(alt_freq_stat_df,
                                             snp_number_window_file,
                                             outdir,
                                             stat_label,
                                             by,
                                             est=freq_est,
                                             method=method)
        elif by == 'genome_window':
            if chr_size is None:
                return 0
            window_file = make_genome_windows(chr_size, genome_window,
                                              genome_step, outdir, force)
            varscore_df = varscore_by_window(alt_freq_stat_df,
                                             window_file,
                                             outdir,
                                             stat_label,
                                             by,
                                             est=freq_est,
                                             method=method)
        elif by == 'gene':
            if gtf is None:
                return 0
            window_file = make_gene_window(gtf, outdir, force)
            varscore_df = varscore_by_window(alt_freq_stat_df,
                                             window_file,
                                             outdir,
                                             stat_label,
                                             by,
                                             est=freq_est,
                                             method=method)
        if method in ['var', 'est']:
            varscore_df.loc[:, 'varscore'] = varscore_df.loc[
                :, analysis_groups[:-1]].apply(log_varscore, axis=1)
        elif method == 'snp_index':
            group0, group1 = analysis_groups[0:2]
            varscore_df.loc[:, 'varscore'] = varscore_df.loc[:, group0] - \
                varscore_df.loc[:, group1]
        logger.info('writing varscore by snp number data...')
        varscore_df.to_csv(varscore_file, sep='\t', index=False)
    else:
        logger.info('snp score based on {stat_label} snp by {by} using {method} exist.',
                    by=by, method=method, stat_label=stat_label)
        varscore_df = pd.read_csv(varscore_file, sep='\t')
    score_plot(varscore_file, method, force)
    return varscore_df


def varscore(vcf, outdir, mutant_bulk, mutant_bulk_freq,
             wild_bulk, wild_bulk_freq, chr_size=None,
             mutant_parent=None, mutant_parent_freq=None,
             wild_parent=None, wild_parent_freq=None,
             background=None, plot_min_depth=5,
             force=False, f2_freq=0.3, f1_freq=0.1,
             bg_freq=0.2, gtf=None):
    # init log
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    init_logger(outdir)
    # read vcf file
    out_table_anno_num_df = vcf2df(vcf, outdir, force)
    # reads stats by group
    reads_stats_df = group_allele(out_table_anno_num_df,
                                  [mutant_bulk, wild_bulk, mutant_parent,
                                   wild_parent, background],
                                  outdir, force)
    # generate plot data
    alt_freq_stat_df = group_alt_freq(
        reads_stats_df, outdir, plot_min_depth, force)
    # filter data
    freq_order = (f2_freq, f1_freq, bg_freq)
    alt_freq_stat_filter_df = filter_snp(alt_freq_stat_df, freq_order,
                                         outdir, force)

    # varscore
    freq_est = [mutant_bulk_freq, wild_bulk_freq,
                mutant_parent_freq, wild_parent_freq]
    by = ['snp_num', 'genome_window', 'gene']
    method = ['var', 'est', 'snp_index']
    freq_df_dict = {'all': alt_freq_stat_df,
                    'filter': alt_freq_stat_filter_df}
    for by_i in by:
        for method_i in method:
            for stat_label in freq_df_dict:
                varsore_df = cal_varscore(freq_df_dict[stat_label],
                                          freq_est,
                                          outdir, force, by=by_i,
                                          method=method_i,
                                          stat_label=stat_label,
                                          chr_size=chr_size, gtf=gtf)


if __name__ == '__main__':
    fire.Fire(varscore)
