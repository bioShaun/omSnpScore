import sys
import attr
import asyncio
import pkg_resources
import numpy as np
import pandas as pd
from loguru import logger
from pathlib import Path, PurePath
from functools import reduce
from collections import OrderedDict
from asyncio.subprocess import PIPE, STDOUT
from . import snpStats, snpAnn

DATA_DIR = PurePath(pkg_resources.resource_filename('snpScore', 'data'))
CHR_SIZE = DATA_DIR / 'chr.size'
SEMA = asyncio.Semaphore(20)
SNP_FREQ_BIAS = 0.1
ALT_FREQ = np.round(2 / 3 - SNP_FREQ_BIAS, 2)
REF_FREQ = np.round(1 / 3 + SNP_FREQ_BIAS, 2)

MUT_NAME = snpStats.SnpGroup.mut.value
WILD_NAME = snpStats.SnpGroup.wild.value


async def launch_cmd(cmd):
    with (await SEMA):
        p = await asyncio.create_subprocess_shell(cmd,
                                                  stdout=PIPE,
                                                  stderr=STDOUT)
        return (await p.communicate())[0].splitlines()


def default_value(value, default_value):
    if not value:
        return default_value
    return value


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


def alt_ref_cut(freq):
    if freq is None:
        ref_cut = REF_FREQ
        alt_cut = ALT_FREQ
    else:
        if freq > 0.5:
            ref_cut = -np.inf
            alt_cut = freq
        else:
            ref_cut = freq
            alt_cut = np.inf
    return ref_cut, alt_cut


@attr.s(kw_only=True)
class SNPscore:
    vcf_table_files = attr.ib()
    vcf_ann_file = attr.ib()
    group_labels = attr.ib()
    outdir = attr.ib(converter=Path)
    min_depth = attr.ib(default=5)
    snp_number_window = attr.ib(default=20)
    snp_number_step = attr.ib(default=5)
    mutant_alt_exp = attr.ib(default=None)
    wild_alt_exp = attr.ib(default=None)
    mutant_parent_alt_exp = attr.ib(default=None)
    wild_parent_alt_exp = attr.ib(default=None)
    background_alt_exp = attr.ib(default=None)
    methods = attr.ib(default='var,snp_index')
    qtlseqr_window = attr.ib(default=1e7)
    qtlseqr_ref_freq = attr.ib(default=REF_FREQ)
    pop_stru = attr.ib(default='RIL')

    def __attrs_post_init__(self):
        self.vcf_table_file_list = self.vcf_table_files.split(',')
        self.methods_list = self.methods.split(',')
        self.group_label_list = self.group_labels.split(',')
        self.group_order = []
        self.freq_dict = OrderedDict()
        self.plot_cmds = []
        self.group_sample_dict = dict()
        self.grp_alt_df = None
        self.grp_dep_df = None
        self.snp_ann_df = None
        self.mutant_wild_group = None
        self.group_label = None
        self.snp_stats_df = None
        self.snp_group_stats_df = None
        self.snp_alt_freq_df = None

    def init_logger(self):
        logfile = self.outdir / 'log.txt'
        logger.add(
            logfile,
            format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}")

    def check_freq(self):
        freq_accordance(self.mutant_alt_exp,
                        self.wild_alt_exp,
                        message=('Mutant and Wild alt freqency direction '
                                 'should not consistent.'),
                        equal=False)
        freq_accordance(
            self.mutant_alt_exp,
            self.mutant_parent_alt_exp,
            message='Mutant and parent alt freqency direction not consistent.')
        freq_accordance(
            self.wild_alt_exp,
            self.wild_parent_alt_exp,
            message='Wild and parent alt freqency direction not consistent.')

    def make_freq_dict(self):
        alt_freq_list = [
            self.mutant_alt_exp, self.wild_alt_exp, self.mutant_parent_alt_exp,
            self.wild_parent_alt_exp, self.background_alt_exp
        ]
        for n, snp_group_i in enumerate(snpStats.SnpGroup.__members__.items()):
            name, member = snp_group_i
            ref_cut, alt_cut = alt_ref_cut(alt_freq_list[n])
            self.freq_dict.update({member.value: [ref_cut, alt_cut]})
        if self.mutant_alt_exp and self.wild_alt_exp:
            self.methods_list.extend(['est_mut_alt', 'est_mut_ref'])

    def check_groups(self):
        group_pairs = [[
            snpStats.SnpGroup.mut.value, snpStats.SnpGroup.wild.value
        ], [snpStats.SnpGroup.mut_pa.value], [snpStats.SnpGroup.wild_pa.value],
            [snpStats.SnpGroup.bg.value]]
        group_label_set = set(self.group_label_list)
        group_out_label = []
        for n, group_pair_i in enumerate(group_pairs):
            union_num = len(set(group_pair_i).intersection(group_label_set))
            if n == 0:
                if union_num != 2:
                    sys.exit('Mutant and Wild is needed!')
                else:
                    self.group_order.extend(group_pair_i)
            else:
                if union_num:
                    self.group_order.extend(group_pair_i)
                else:
                    pass
        for group_i in self.group_order:
            if (self.freq_dict[group_i][0] is None
                    or self.freq_dict[group_i][1] is None):
                sys.exit(f'{group_i} frequency not spcified!')
            else:
                label_group = [
                    str(each) for each in self.freq_dict[group_i]
                    if not np.isinf(each)
                ]
                group_out_label.append(group_i)
                group_out_label.extend(label_group)
        self.group_label = '_'.join(group_out_label)

    def load_stats(self):
        logger.info('Loading tables...')
        self.snp_stats_dfs = [
            snpStats.load_hd5(table_i) for table_i in self.vcf_table_file_list
        ]
        logger.info('Concatinating tables...')
        self.snp_stats_df = reduce(
            lambda x, y: pd.merge(x, y, on=['Chr', 'Pos', 'Alt']),
            self.snp_stats_dfs)
        self.snp_stats_df.fillna(0, inplace=True)

    def group_stats(self):
        logger.info('Merging group reads...')
        self.alt_df = self.snp_stats_df.loc[:, 'alt_count'].copy()
        self.dep_df = self.snp_stats_df.loc[:, 'dep_count'].copy()
        self.alt_df.columns = self.group_label_list
        self.dep_df.columns = self.group_label_list
        self.grp_alt_df = self.alt_df.T.groupby(level=0).agg('sum').T
        self.grp_dep_df = self.dep_df.T.groupby(level=0).agg('sum').T

    def alt_freq(self):
        logger.info('Filtering allele depth...')
        self.mutant_wild_group = [MUT_NAME, WILD_NAME]
        dep_passed_snp = self.grp_dep_df.loc[:, self.mutant_wild_group].min(
            1) >= self.min_depth
        self.passed_grp_dep_df = self.grp_dep_df[dep_passed_snp]
        self.passed_grp_alt_df = self.grp_alt_df[dep_passed_snp]
        self.passed_grp_dep_df.applymap(lambda x: x
                                        if x >= self.min_depth else np.nan)
        logger.info('Calculating alt allele freq...')
        self.grp_alt_freq_df = self.passed_grp_alt_df / self.passed_grp_dep_df
        self.grp_alt_freq_df = self.grp_alt_freq_df.loc[:, self.group_order]
        self.grp_alt_freq_df = self.grp_alt_freq_df.reset_index()
        self.grp_alt_freq_df.to_csv(self.grp_alt_freq_file, index=False)

    def snp_filter(self):
        logger.info('Filtering snp by freq...')
        self.snp_alt_filter_df = snpStats.filter_snp(self.grp_alt_freq_df,
                                                     self.freq_dict,
                                                     self.group_label,
                                                     self.snp_alt_filter_file)

    def make_windows(self):
        logger.info(
            'Making snp number slidewindow bed file windows {w} step {s}...',
            w=self.snp_number_window,
            s=self.snp_number_step)
        self.snp_num_window_file = snpStats.make_snp_number_windows(
            self.snp_alt_filter_df, self.group_label, self.snp_number_window,
            self.snp_number_step, self.outdir)
        self.windows_files = [self.snp_num_window_file]

    def load_snp_ann(self):
        logger.info('Loading snp annotation...')
        if self.snp_ann_df is None:
            self.snp_ann_df = pd.read_pickle(self.vcf_ann_file)

    def annotate_snp_window(self):
        # add snp annotation to snp score region
        self.load_snp_ann()
        self.freq_dis_ann_df = self.freq_dis_df.merge(
            self.snp_ann_df,
            left_on=['Chrom', 'Pos', 'Alt'],
            right_on=['#CHROM', 'POS', 'ALT'],
            how='left')
        self.freq_dis_ann_df.drop(['#CHROM', 'POS', 'Alt'],
                                  inplace=True,
                                  axis=1)
        self.freq_dis_ann_df.rename(columns={
            MUT_NAME: f'{MUT_NAME}_alt_freq',
            WILD_NAME: f'{WILD_NAME}_alt_freq',
        },
            inplace=True)
        return self.freq_dis_ann_df

    def annotate_snp_score(self):
        # add snp annotation to snp score table and flat
        self.score_ann_df = self.score_df.merge(
            self.freq_dis_ann_df,
            left_on=['Chrom', 'Start', 'End'],
            right_on=['Chrom', 'Start', 'End'],
            how='left')
        snpeff_anno = list(
            self.score_ann_df.INFO.map(snpAnn.extract_snpeff_anno))
        snpeff_anno_df = pd.DataFrame(snpeff_anno)
        snpeff_anno_df.columns = [
            'Feature', 'Gene', 'Transcript', 'Variant_DNA_Level',
            'Variant_Protein_Level'
        ]
        self.score_ann_df = pd.concat([self.score_ann_df, snpeff_anno_df],
                                      axis=1)
        self.score_ann_df.drop('INFO', axis=1, inplace=True)
        self.score_ann_df = snpAnn.split_dataframe_rows(
            self.score_ann_df,
            column_selectors=[
                'Feature', 'Gene', 'Transcript', 'Variant_DNA_Level',
                'Variant_Protein_Level'
            ],
            row_delimiter='|')
        return self.score_ann_df

    def snp_score(self):
        for n, window_file in enumerate(self.windows_files):
            # group snp by window
            self.freq_dis_df = snpStats.snp_freq_by_window(
                self.snp_alt_filter_df, self.group_label, window_file,
                self.outdir)
            self.freq_dis_ann_df = self.annotate_snp_window()

            # calculating snp score using different methods
            for method in self.methods_list:
                window_file_name = window_file.name
                logger.info(
                    'Calculating snp score using {window} by {method}...',
                    window=window_file_name,
                    method=method)
                if self.group_label in window_file.stem:
                    score_name = window_file.stem
                else:
                    score_name = f'{self.group_label}.{window_file.stem}'
                self.score_file = self.outdir / \
                    f'{score_name}.{method}.score.csv'
                if snpStats.is_valid_file(self.score_file):
                    self.score_df = pd.read_csv(self.score_file)
                else:
                    self.score_df = snpStats.cal_score(self.freq_dis_df,
                                                       self.freq_dict,
                                                       method=method)

                    if self.score_df is None:
                        continue
                    self.score_df.to_csv(self.score_file)
                self.score_ann_file = self.outdir / \
                    f'{score_name}.{method}.score.ann.csv'
                if not snpStats.is_valid_file(self.score_ann_file):
                    self.score_ann_df = self.annotate_snp_score()
                    self.score_ann_df.to_csv(self.score_ann_file, index=False)

                self.plot_cmds.append(
                    snpStats.score_plot(self.score_file, method))

    def run_qtlseqr(self):
        logger.info('Running QTLseqr...')
        if not snpStats.is_valid_file(self.qtlseqr_input):
            if (not self.grp_dep_df) or (not self.grp_alt_df):
                self.load_stats()
                self.group_stats()
            self.grp_ref_df = self.grp_dep_df - self.grp_alt_df
            self.grp_ref_df.columns = [
                f'AD_REF.{sp_i}' for sp_i in self.grp_ref_df.columns
            ]
            self.grp_alt_df.columns = [
                f'AD_ALT.{sp_i}' for sp_i in self.grp_alt_df.columns
            ]
            self.grp_ref_df = self.grp_ref_df.astype('int')
            self.qtlseqr_df = self.grp_ref_df.merge(self.grp_alt_df,
                                                    on=['Chr', 'Pos', 'Alt'])
            self.qtlseqr_df.index.names = ['CHROM', 'POS', 'ALT']
            self.qtlseqr_df = self.qtlseqr_df[self.qtlseqr_df.sum(1) > 0]
            self.qtlseqr_df = self.qtlseqr_df.reset_index()
            self.qtlseqr_df.to_csv(self.qtlseqr_input, index=False)
        out_prefix = self.outdir / 'QTLseqr'
        cmd = snpStats.run_qtlseqr_cmd(self.qtlseqr_input,
                                       h_bulk=MUT_NAME,
                                       l_bulk=WILD_NAME,
                                       out_prefix=out_prefix,
                                       window=self.qtlseqr_window,
                                       ref_freq=self.qtlseqr_ref_freq,
                                       min_sample_dp=self.min_depth,
                                       pop_stru=self.pop_stru)
        self.plot_cmds.append(cmd)

    def plot(self):
        logger.info('Ploting ...')
        self.grp_alt_freq_file = self.outdir / 'snp.freq.csv'
        self.plot_cmds.append(
            snpStats.score_plot(self.grp_alt_freq_file, 'density'))
        self.plot_cmds.append(
            snpStats.score_plot(self.snp_alt_filter_file, 'density'))
        self.plot_cmds = list(filter(None, self.plot_cmds))
        if self.plot_cmds:
            loop = asyncio.get_event_loop()
            f = asyncio.wait([launch_cmd(each) for each in self.plot_cmds])
            loop.run_until_complete(f)
            loop.close()

    @property
    def run(self):
        self.init_logger()
        self.check_freq()
        self.make_freq_dict()
        self.check_groups()
        self.qtlseqr_input = self.outdir / 'qtlseqr.csv'
        self.snp_alt_filter_file = self.outdir / \
            f'{self.group_label}.snp.freq.csv'
        if snpStats.is_valid_file(self.snp_alt_filter_file):
            self.snp_alt_filter_df = pd.read_csv(self.snp_alt_filter_file)
        else:
            self.grp_alt_freq_file = self.outdir / 'snp.freq.csv'
            if snpStats.is_valid_file(self.grp_alt_freq_file):
                self.grp_alt_freq_df = pd.read_csv(self.grp_alt_freq_file)
            else:
                self.load_stats()
                self.group_stats()
                self.alt_freq()
            self.snp_filter()
        self.make_windows()
        self.snp_score()
        self.run_qtlseqr()
        self.plot()
        logger.info('The End.')
