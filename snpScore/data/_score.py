import attr


@attr.s
class qtlSeqr:
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