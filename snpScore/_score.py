import attr
import numpy as np
import pandas as pd
from loguru import logger
from pathlib import Path
from collections import OrderedDict
from ._var import REF_FREQ, QTLSEQR_PLOT
from ._var import MUT_NAME, WILD_NAME
from ._utils import freq_accordance
from ._utils import alt_ref_cut
from ._utils import filter_snp
from ._utils import make_snp_number_windows
from ._utils import snp_freq_by_window
from ._utils import cal_score, score_plot
from ._utils import extract_snpeff_anno
from ._utils import split_dataframe_rows


@attr.s
class snpScoreBox:
    alt_freq_df = attr.ib()
    outdir = attr.ib(converter=Path)
    grp_list = attr.ib(factory=list)
    method_list = attr.ib(factory=list)
    min_depth = attr.ib(default=5)
    snp_number_window = attr.ib(default=20)
    snp_number_step = attr.ib(default=5)
    ref_freq = attr.ib(default=REF_FREQ)
    mutant_alt_exp = attr.ib(default=None)
    wild_alt_exp = attr.ib(default=None)
    mutant_parent_alt_exp = attr.ib(default=None)
    wild_parent_alt_exp = attr.ib(default=None)
    background_alt_exp = attr.ib(default=None)
    vcf_ann_file = attr.ib(default=None)

    def __attrs_post_init__(self):
        self._freq_dict = OrderedDict()
        self.plot_cmds = []
        self._alt_filter_freq_df = None
        self._group_label = None
        self._alt_freq_dis_df = None
        self._snp_ann_df = None
        self._snp_window_ann_df = None
        if self.ref_freq > 0.5:
            self.ref_freq = 1 - self.ref_freq

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

    @property
    def freq_dict(self):
        if not self._freq_dict:
            alt_freq_list = [
                self.mutant_alt_exp, self.wild_alt_exp,
                self.mutant_parent_alt_exp, self.wild_parent_alt_exp,
                self.background_alt_exp
            ]
            for n, member in enumerate(self.grp_list):
                ref_cut, alt_cut = alt_ref_cut(alt_freq_list[n], self.ref_freq)
                self._freq_dict.update({member: [ref_cut, alt_cut]})
        return self._freq_dict

    @property
    def group_label(self):
        if self._group_label is None:
            group_out_label = []
            for group_i in self.grp_list:
                label_group = [
                    str(each) for each in self.freq_dict[group_i]
                    if not np.isinf(each)
                ]
                group_out_label.append(group_i)
                group_out_label.extend(label_group)
            self._group_label = '_'.join(group_out_label)
        return self._group_label

    @property
    def alt_filter_freq_file(self):
        return self.outdir / f'{self.group_label}.snp.freq.csv'

    @property
    def alt_filter_freq_df(self):
        if self._alt_filter_freq_df is None:
            if self.alt_filter_freq_file.is_file():
                self._alt_filter_freq_df = pd.read_csv(
                    self.alt_filter_freq_file)
            else:
                logger.info('Filtering snp by freq...')
                self._alt_filter_freq_df = filter_snp(
                    self.alt_freq_df, self.freq_dict, self.group_label,
                    self.alt_filter_freq_file)
        return self._alt_filter_freq_df

    @property
    def snp_number_window_file(self):
        return make_snp_number_windows(self.alt_filter_freq_df,
                                       self.group_label,
                                       self.snp_number_window,
                                       self.snp_number_step, self.outdir)

    @property
    def alt_freq_dis_df(self):
        if self._alt_freq_dis_df is None:
            self._alt_freq_dis_df = snp_freq_by_window(
                self.alt_filter_freq_df, self.group_label,
                self.snp_number_window_file, self.outdir)
        return self._alt_freq_dis_df

    @property
    def snp_ann_df(self):
        if self._snp_ann_df is None:
            logger.info('Loading snp annotation...')
            self._snp_ann_df = pd.read_pickle(self.vcf_ann_file)
        return self._snp_ann_df

    @property
    def snp_window_ann_df(self):
        if self._snp_window_ann_df is None:
            # add snp annotation to snp score region
            logger.info('Generating snp window annotation...')
            self._snp_window_ann_df = self.alt_freq_dis_df.merge(
                self.snp_ann_df,
                left_on=['Chrom', 'Pos', 'Alt'],
                right_on=['#CHROM', 'POS', 'ALT'],
                how='left')
            self._snp_window_ann_df.drop(['#CHROM', 'POS', 'Alt'],
                                         inplace=True,
                                         axis=1)
            self._snp_window_ann_df.rename(columns={
                MUT_NAME: f'{MUT_NAME}_alt_freq',
                WILD_NAME: f'{WILD_NAME}_alt_freq',
            },
                                           inplace=True)
        return self._snp_window_ann_df

    @property
    def score_ann_df(self):
        # add snp annotation to snp score table and flat
        logger.info('Annotating snp score...')
        self._score_ann_df = self.score_df.merge(
            self.snp_window_ann_df,
            left_on=['Chrom', 'Start', 'End'],
            right_on=['Chrom', 'Start', 'End'],
            how='left')
        snpeff_anno = list(self._score_ann_df.INFO.map(extract_snpeff_anno))
        snpeff_anno_df = pd.DataFrame(snpeff_anno)
        snpeff_anno_df.columns = [
            'Feature', 'Gene', 'Transcript', 'Variant_DNA_Level',
            'Variant_Protein_Level'
        ]
        self._score_ann_df = pd.concat([self._score_ann_df, snpeff_anno_df],
                                       axis=1)
        self._score_ann_df.drop('INFO', axis=1, inplace=True)
        self._score_ann_df = split_dataframe_rows(self._score_ann_df,
                                                  column_selectors=[
                                                      'Feature', 'Gene',
                                                      'Transcript',
                                                      'Variant_DNA_Level',
                                                      'Variant_Protein_Level'
                                                  ],
                                                  row_delimiter='|')
        return self._score_ann_df

    @property
    def score_jobs(self):
        # calculating snp score using different methods
        for method in self.method_list:
            window_file_name = self.snp_number_window_file.name
            score_name = self.snp_number_window_file.stem
            self.score_file = self.outdir / \
                f'{score_name}.{method}.score.csv'
            if not self.score_file.is_file():
                logger.info(
                    'Calculating snp score using {window} by {method}...',
                    window=window_file_name,
                    method=method)
                self.score_df = cal_score(self.alt_freq_dis_df,
                                          self.freq_dict,
                                          method=method)
                if self.score_df is None:
                    continue
                self.score_df.to_csv(self.score_file)
            else:
                self.score_df = pd.read_csv(self.score_file)
            self.plot_cmds.append(score_plot(self.score_file, method))
            self.score_ann_file = self.outdir / \
                f'{score_name}.{method}.score.ann.csv'
            if not self.score_ann_file.is_file():
                if self.vcf_ann_file:
                    self.score_ann_df.to_csv(self.score_ann_file, index=False)
        self.grp_alt_freq_file = self.outdir / 'snp.freq.csv'
        self.plot_cmds.append(score_plot(self.grp_alt_freq_file, 'density'))
        self.plot_cmds.append(score_plot(self.alt_filter_freq_file, 'density'))
        self.plot_cmds = list(filter(None, self.plot_cmds))
        return self.plot_cmds


@attr.s
class qtlSeqr:

    input_table = attr.ib(converter=Path)
    out_dir = attr.ib(converter=Path)
    run_qtlseqr = attr.ib(factory=bool)
    run_ed = attr.ib(factory=bool)
    window = attr.ib(default=1e7)
    ref_freq = attr.ib(default=REF_FREQ)
    pop_stru = attr.ib(default='RIL')
    min_sample_dp = attr.ib(default=5)

    @property
    def qtlseqr_job(self):
        cmd_flag = ''
        if self.run_qtlseqr:
            cmd_flag += '--qtlseqr '
        if self.run_ed:
            cmd_flag += '--ed '
        logger.info('Generating QTLseqr command...')
        cmd_line = (f'Rscript {QTLSEQR_PLOT} '
                    f'--input {self.input_table} '
                    f'--high_bulk {MUT_NAME} '
                    f'--low_bulk {WILD_NAME} '
                    f'--out_dir {self.out_dir} '
                    f'--window {self.window} '
                    f'--ref_freq {self.ref_freq} '
                    f'--min_sample_dp {self.min_sample_dp} '
                    f'--pop_stru {self.pop_stru} '
                    f'{cmd_flag}')
        return cmd_line
