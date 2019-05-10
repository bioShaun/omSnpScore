import sys
import asyncio
import pkg_resources
from asyncio.subprocess import PIPE, STDOUT
import pandas as pd
from loguru import logger
from pathlib import Path, PurePath
from collections import OrderedDict
from . import snpStats


DATA_DIR = PurePath(pkg_resources.resource_filename('snpScore', 'data'))
CHR_SIZE = DATA_DIR / 'chr.size'

SEMA = asyncio.Semaphore(20)


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


class SNPscore:

    def __init__(self,
                 vcf_table_files,
                 group_labels,
                 outdir,
                 mutant_freq_low=0.3,
                 wild_freq_low=0.3,
                 mutant_parent_freq_low=None,
                 wild_parent_freq_low=None,
                 background_freq_low=None,
                 mutant_freq_up=0.7,
                 wild_freq_up=0.7,
                 mutant_parent_freq_up=None,
                 wild_parent_freq_up=None,
                 background_freq_up=None,
                 min_depth=5,
                 snp_number_window=10,
                 snp_number_step=5,
                 genome_window=1000000,
                 genome_step=500000,
                 methods='var,snp_index'):
        self.outdir = Path(outdir)
        self.vcf_table_file_list = vcf_table_files.split(',')
        self.group_label_list = group_labels.split(',')
        self.group_order = []
        self.group_label = None
        self.snp_stats_df = None
        self.snp_group_stats_df = None
        self.snp_alt_freq_df = None
        self.min_depth = default_value(min_depth, 5)
        self.snp_number_window = default_value(snp_number_window, 10)
        self.snp_number_step = default_value(snp_number_step, 5)
        self.genome_window = default_value(genome_window, 1000000)
        self.genome_step = default_value(genome_step, 500000)
        methods = default_value(methods, 'var,snp_index')
        self.methods_list = methods.split(',')
        mutant_freq_low = default_value(mutant_freq_low, 0.3)
        mutant_freq_up = default_value(mutant_freq_up, 0.7)
        wild_freq_low = default_value(wild_freq_low, 0.3)
        wild_freq_up = default_value(wild_freq_up, 0.7)
        mutant_parent_freq_low = default_value(mutant_parent_freq_low, None)
        mutant_parent_freq_up = default_value(mutant_parent_freq_up, None)
        wild_parent_freq_low = default_value(wild_parent_freq_low, None)
        wild_parent_freq_up = default_value(wild_parent_freq_up, None)
        background_freq_low = default_value(background_freq_low, None)
        background_freq_up = default_value(background_freq_up, None)
        self.plot_cmds = []
        self.freq_dict = OrderedDict({
            'mutant': [mutant_freq_low, mutant_freq_up],
            'wild': [wild_freq_low, wild_freq_up],
            'mutant_parent': [mutant_parent_freq_low, mutant_parent_freq_up],
            'wild_parent': [wild_parent_freq_low, wild_parent_freq_up],
            'background': [background_freq_low, background_freq_up],
        })

    def init_logger(self):
        logfile = self.outdir / 'log.txt'
        logger.add(
            logfile,
            format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}")

    def check_groups(self):
        group_pairs = [['mutant', 'wild'],
                       ['mutant_parent', 'wild_parent'],
                       ['background']]
        group_label_set = set(self.group_label_list)
        group_out_label = []
        for n, group_pair_i in enumerate(group_pairs):
            union_num = len(set(group_pair_i).intersection(group_label_set))
            if n == 0:
                if union_num != 2:
                    sys.exit('Mutant and Wild is needed!')
                else:
                    self.group_order.extend(group_pair_i)
            elif n == 1:
                if union_num == 1:
                    sys.exit('Mutant and Wild Parent should be paired!')
                elif union_num == 2:
                    self.group_order.extend(group_pair_i)
                else:
                    pass
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
                group_out_label.extend(
                    [group_i,
                     str(self.freq_dict[group_i][0]),
                     str(self.freq_dict[group_i][1])
                     ])
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
        self.snp_stats_df.dropna(inplace=True)

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
        dep_passed_snp = self.grp_dep_df.min(1) >= self.min_depth
        self.passed_grp_dep_df = self.grp_dep_df[dep_passed_snp]
        self.passed_grp_alt_df = self.grp_alt_df[dep_passed_snp]
        logger.info('Calculating alt allele freq...')
        self.grp_alt_freq_df = self.passed_grp_alt_df / self.passed_grp_dep_df
        self.grp_alt_freq_df = self.grp_alt_freq_df.loc[:, self.group_order]

    def snp_filter(self):
        logger.info('Filtering snp by freq...')
        self.snp_alt_filter_df = snpStats.filter_snp(self.grp_alt_freq_df,
                                                     self.freq_dict,
                                                     self.group_label,
                                                     self.snp_alt_filter_file)
        self.snp_alt_filter_df = self.snp_alt_filter_df.reset_index()

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
        logger.info(
            'making genome slidewindow bed file windows {w} step {s}...',
            w=self.genome_window, s=self.genome_step)
        self.snp_genome_window_file = snpStats.make_genome_windows(
            CHR_SIZE, self.genome_window,
            self.genome_step, self.outdir)
        self.windows_files = [self.snp_num_window_file,
                              self.snp_genome_window_file]

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
                if not self.score_file.exists():
                    self.score_df = snpStats.cal_score(self.freq_dis_df,
                                                       method=method)
                    self.score_df.to_csv(self.score_file)
                self.plot_cmds.append(
                    snpStats.score_plot(self.score_file, method))

    def plot(self):
        logger.info('ploting ...')
        self.grp_alt_freq_file = self.outdir / 'snp.freq.csv'
        if not self.grp_alt_freq_file.exists():
            out_grp_alt_freq_df = self.grp_alt_freq_df.reset_index()
            out_grp_alt_freq_df.to_csv(self.grp_alt_freq_file,
                                       index=False)
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
        self.check_groups()
        self.init_logger()
        self.snp_alt_filter_file = self.outdir / \
            f'{self.group_label}.snp.freq.csv'
        if not self.snp_alt_filter_file.exists():
            self.load_stats()
            self.group_stats()
            self.alt_freq()
            self.snp_filter()
        else:
            self.snp_alt_filter_df = pd.read_csv(
                self.snp_alt_filter_file)
        self.make_windows()
        self.snp_score()
        self.plot()
        logger.info('all done ...')
