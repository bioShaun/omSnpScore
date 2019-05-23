import sys
import asyncio
import pkg_resources
import numpy as np
import pandas as pd
from loguru import logger
from pathlib import Path, PurePath
from collections import OrderedDict
from asyncio.subprocess import PIPE, STDOUT
from . import snpStats


DATA_DIR = PurePath(pkg_resources.resource_filename('snpScore', 'data'))
CHR_SIZE = DATA_DIR / 'chr.size'
SEMA = asyncio.Semaphore(20)
SNP_FREQ_BIAS = 0.1
ALT_FREQ = np.round(2 / 3 - SNP_FREQ_BIAS, 2)
REF_FREQ = np.round(1 / 3 + SNP_FREQ_BIAS, 2)
DOMINANT_ALT_FREQ = 0.9
DOMINANT_REF_freq = 0.1


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


def alt_ref_cut(freq, dominant):
    if freq is None:
        if dominant == 'yes':
            ref_cut = DOMINANT_REF_freq
            alt_cut = DOMINANT_ALT_FREQ
        else:
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


class SNPscore:

    def __init__(self,
                 vcf_table_files,
                 group_labels,
                 outdir,
                 mutant_alt_exp=None,
                 wild_alt_exp=None,
                 mutant_parent_alt_exp=None,
                 wild_parent_alt_exp=None,
                 background_alt_exp=None,
                 mutant_dominant='no',
                 wild_dominant='no',
                 background_dominant='no',
                 min_depth=5,
                 snp_number_window=10,
                 snp_number_step=5,
                 genome_window=1000000,
                 genome_step=500000,
                 methods='var,snp_index,est_mut_alt,est_mut_ref',
                 tow_side=False):
        self.outdir = Path(outdir)
        self.vcf_table_file_list = vcf_table_files.split(',')
        self.group_label_list = group_labels.split(',')
        self.group_order = []
        self.mutant_wild_group = None
        self.group_label = None
        self.snp_stats_df = None
        self.snp_group_stats_df = None
        self.snp_alt_freq_df = None
        self.min_depth = default_value(min_depth, 5)
        self.snp_number_window = default_value(snp_number_window, 20)
        self.snp_number_step = default_value(snp_number_step, 5)
        self.genome_window = default_value(genome_window, 1000000)
        self.genome_step = default_value(genome_step, 500000)
        methods = default_value(
            methods, 'var,snp_index,ext_mut_alt,ext_mut_ref')
        self.methods_list = methods.split(',')
        self.mutant_alt_exp = mutant_alt_exp
        self.wild_alt_exp = wild_alt_exp
        self.mutant_parent_alt_exp = mutant_parent_alt_exp
        self.wild_parent_alt_exp = wild_parent_alt_exp
        self.background_alt_exp = background_alt_exp
        self.mutant_dominant = mutant_dominant
        self.wild_dominant = wild_dominant
        self.background_dominant = background_dominant
        self.freq_dict = OrderedDict()
        self.plot_cmds = []

    def init_logger(self):
        logfile = self.outdir / 'log.txt'
        logger.add(
            logfile,
            format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}")

    def check_dominant(self):
        check_mut_wild = (self.mutant_dominant,
                          self.wild_dominant) != ('yes', 'yes')
        assert check_mut_wild, 'Only one of mutant/wild can be dominant'

    def check_freq(self):
        freq_accordance(
            self.mutant_alt_exp,
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
            self.mutant_alt_exp,
            self.wild_alt_exp,
            self.mutant_parent_alt_exp,
            self.wild_parent_alt_exp,
            self.background_alt_exp]
        dominant_list = [
            self.mutant_dominant,
            self.wild_dominant,
            self.mutant_dominant,
            self.wild_dominant,
            self.background_dominant]
        freq_dominant_list = list(zip(alt_freq_list, dominant_list))
        for n, snp_group_i in enumerate(snpStats.SnpGroup.__members__.items()):
            name, member = snp_group_i
            freq, dominant = freq_dominant_list[n]
            ref_cut, alt_cut = alt_ref_cut(freq, dominant)
            self.freq_dict.update({
                member.value: [ref_cut, alt_cut]
            })

    def check_groups(self):
        group_pairs = [[snpStats.SnpGroup.mut.value,
                        snpStats.SnpGroup.wild.value],
                       [snpStats.SnpGroup.mut_pa.value],
                       [snpStats.SnpGroup.wild_pa.value],
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
            if (self.freq_dict[group_i][0] is None or
                    self.freq_dict[group_i][1] is None):
                sys.exit(f'{group_i} frequency not spcified!')
            else:
                label_group = [str(each) for each in self.freq_dict[group_i]
                               if not np.isinf(each)]
                group_out_label.append(group_i)
                group_out_label.extend(label_group)
        self.group_label = '_'.join(group_out_label)

    def load_stats(self):
        logger.info('Loading tables...')
        self.snp_stats_dfs = [
            snpStats.load_hd5(
                table_i)
            for table_i in self.vcf_table_file_list
        ]
        logger.info('Concatinating tables...')
        self.snp_stats_df = pd.concat(self.snp_stats_dfs, sort=False, axis=1)
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
        self.mutant_wild_group = [
            snpStats.SnpGroup.mut.value, snpStats.SnpGroup.wild.value]
        dep_passed_snp = self.grp_dep_df.loc[
            :, self.mutant_wild_group].min(1) >= self.min_depth
        self.passed_grp_dep_df = self.grp_dep_df[dep_passed_snp]
        self.passed_grp_alt_df = self.grp_alt_df[dep_passed_snp]
        self.passed_grp_dep_df.applymap(
            lambda x: x if x >= self.min_depth else np.nan)
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
            'making snp number slidewindow bed file windows {w} step {s}...',
            w=self.snp_number_window, s=self.snp_number_step)
        self.snp_num_window_file = snpStats.make_snp_number_windows(
            self.snp_alt_filter_df,
            self.group_label,
            self.snp_number_window,
            self.snp_number_step,
            self.outdir)
        # logger.info(
        #     'making genome slidewindow bed file windows {w} step {s}...',
        #     w=self.genome_window, s=self.genome_step)
        # self.snp_genome_window_file = snpStats.make_genome_windows(
        #     CHR_SIZE, self.genome_window,
        #     self.genome_step, self.outdir)
        self.windows_files = [self.snp_num_window_file]

    def snp_score(self):
        for n, window_file in enumerate(self.windows_files):
            self.freq_dis_df = snpStats.snp_freq_by_window(
                self.snp_alt_filter_df,
                self.group_label,
                window_file,
                self.outdir
            )
            for method in self.methods_list:
                window_file_name = window_file.name
                logger.info(
                    'calculating snp score using {window} by {method}...',
                    window=window_file_name,
                    method=method)
                if self.group_label in window_file.stem:
                    score_name = window_file.stem
                else:
                    score_name = f'{self.group_label}.{window_file.stem}'
                self.score_file = self.outdir / \
                    f'{score_name}.{method}.score.csv'
                self.score_df = snpStats.cal_score(self.freq_dis_df,
                                                   self.freq_dict,
                                                   method=method)
                if self.score_df is None:
                    continue
                self.score_df.to_csv(self.score_file)
                self.plot_cmds.append(
                    snpStats.score_plot(self.score_file, method))

    def plot(self):
        logger.info('ploting ...')
        self.grp_alt_freq_file = self.outdir / 'snp.freq.csv'
        self.plot_cmds.append(snpStats.score_plot(
            self.grp_alt_freq_file, 'density'))
        self.plot_cmds.append(snpStats.score_plot(
            self.snp_alt_filter_file, 'density'))
        self.plot_cmds = list(filter(None, self.plot_cmds))
        if self.plot_cmds:
            loop = asyncio.get_event_loop()
            f = asyncio.wait([launch_cmd(each) for each in self.plot_cmds])
            loop.run_until_complete(f)
            loop.close()

    @property
    def run(self):
        self.init_logger()
        self.check_dominant()
        self.check_freq()
        self.make_freq_dict()
        self.check_groups()
        self.snp_alt_filter_file = self.outdir / \
            f'{self.group_label}.snp.freq.csv'
        if snpStats.is_valid_file(self.snp_alt_filter_file):
            self.snp_alt_filter_df = pd.read_csv(
                self.snp_alt_filter_file)
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
        self.plot()
        logger.info('all done ...')
