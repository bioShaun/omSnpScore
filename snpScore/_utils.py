import re
import sys
import json
import shutil
import jinja2
import asyncio
import delegator
import functools
import numpy as np
import pandas as pd
from io import StringIO
from loguru import logger
from decimal import Decimal, getcontext
from typing import Callable, Iterable, List, Optional, Union
from pathlib import Path, PurePath
from functools import reduce
from pybedtools import BedTool
from datetime import datetime
from ._var import GROUPS, OFFSET
from ._var import SnpGroup, SnpRep, SnpGroupFreq, VarScoreParams
from ._var import SNP_SCORE_PLOT, COLUMN_NAME_MAP, SNP_BASIC_COL, QTLSERQ_BASIC_COL
from ._var import (
    VAR_SCORE_OUT_COL,
    ANN_POS_COLS,
    SNP_DENSITY_POS_COLS,
    QTLSEQR_POS_COLS,
)
from ._var import SCIENTIFIC_NUMBER_COLS, ED_SPECIFIC_COLS, QTLSEQR_SPECIFIC_COLS
from ._var import (
    ANN_OUT_COLS,
    QTLSEQR_CHROM_NAME,
    SNP_DENSITY_SORT_COL,
    VAR_SCORE_SORT_COL,
    VAR_SCORE_ANN_SORT_COL,
)
from ._var import (
    SnpDensityStatsTable,
    SnpTableConstants,
    VarFilterParams,
    ANN_TABLE_COL,
    ANN_TABLE_NAME,
)

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


class SnpDensityWindowFailed(Exception):
    pass


async def async_sh_job(cmd, sema):
    with (await sema):
        p = await asyncio.create_subprocess_shell(
            cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.STDOUT
        )
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
        s_i_list = s_i.split(",")
        g_i_list = [GROUPS[n]] * len(s_i_list)
        sample_list.extend(s_i_list)
        group_list.extend(g_i_list)
    return sample_list, group_list


def sample_and_group_for_web(parameters_obj):
    group_names = parameters_obj.get("group_names")
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
        if group_names[0] == "not_a_group_id":
            group_list = sample_list[:]
    return sample_list, group_list


def freq_accordance(freq_a, freq_b, message, equal=True):
    non_legal_freq = (freq_a != 0.5) and (freq_b != 0.5)
    assert non_legal_freq, "alt freqency should not equal to 0.5."
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


def filter_snp(
    alt_freq_stat_df, freq_dict, filter_freq_stats, filter_method="nonsymmetrical"
):
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
    alt_rep_df = equal2parent(
        alt_rep_df, SnpGroupFreq.mut.value, SnpGroupFreq.mut_pa.value
    )
    alt_rep_df = equal2parent(
        alt_rep_df, SnpGroupFreq.wild.value, SnpGroupFreq.wild_pa.value
    )
    if filter_method == "nonsymmetrical":
        # step3 mutant not equal to wild
        mask = (
            alt_rep_df.loc[:, SnpGroupFreq.mut.value]
            != alt_rep_df.loc[:, SnpGroupFreq.wild.value]
        )
        alt_rep_df = alt_rep_df[mask]
    alt_freq_stat_filter_df = alt_freq_stat_df.loc[alt_rep_df.index]
    alt_freq_stat_filter_df.to_csv(filter_freq_stats, index=False)
    return alt_freq_stat_filter_df


def slidewindow(obj, window, step):
    for i in range(0, len(obj), step):
        yield obj[i : i + window]


def make_snp_number_windows(stat_df, group_label, window, step, outdir):
    snp_num_window_file = outdir / f"{group_label}.snp_num.window.w{window}.s{step}.bed"
    if not snp_num_window_file.is_file():
        logger.info(
            "Making snp number slidewindow bed file windows {w} step {s}...",
            w=window,
            s=step,
        )
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
        snp_num_window_df = pd.DataFrame(
            snp_num_window_list, columns=["Chrom", "Start", "End"]
        )
        snp_num_window_df.to_csv(
            snp_num_window_file, sep="\t", index=False, header=False
        )
    return snp_num_window_file


def snp_freq_by_window(stat_df, group_label, window_file, outdir):
    groups = stat_df.columns[4:]
    alt_freq_stat_bed = outdir / f"{group_label}.snp.plot.bed"
    if not alt_freq_stat_bed.is_file():
        alt_freq_stat_df = stat_df.copy()
        alt_freq_stat_df.loc[:, "start"] = alt_freq_stat_df.Pos - 1
        bed_cols = ["Chr", "start", "Pos", "Alt"]
        bed_cols.extend(groups)
        alt_freq_stat_df.to_csv(
            alt_freq_stat_bed, sep="\t", columns=bed_cols, header=None, index=False
        )
    window_bed = BedTool(str(window_file))
    snp_bed = BedTool(str(alt_freq_stat_bed))
    intersect_obj = window_bed.intersect(snp_bed, sorted=True, wo=True)
    intersect_obj_cols = ["Chrom", "Start", "End"]
    intersect_obj_cols.extend(["snp_Chr", "snp_start", "Pos", "Alt"])
    intersect_obj_cols.extend(groups)
    intersect_obj_cols.append("overlap")
    intersect_str = StringIO(str(intersect_obj))
    intersect_df = pd.read_csv(
        intersect_str, sep="\t", header=None, names=intersect_obj_cols
    )
    intersect_df.drop(["snp_Chr", "snp_start", "overlap"], axis=1, inplace=True)
    return intersect_df


def log_varscore(row):
    score_s = np.power(-np.log10(reduce(lambda a, b: a * b, row)), 10)
    return score_s.astype("int")


def mut_wild_ext_freq(intersect_df, freq_dict, mut="alt"):
    if mut == "alt":
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


def cal_score(intersect_df, freq_dict, method="var", min_snp_num=3):
    stats_cols = [
        "Chrom",
        "Start",
        "End",
        SnpGroupFreq.mut.value,
        SnpGroupFreq.wild.value,
    ]
    stats_df = intersect_df.loc[:, stats_cols]
    varscore_size_df = stats_df.groupby(["Chrom", "Start", "End"]).size()
    mask = varscore_size_df >= min_snp_num
    if method == "var":
        varscore_df = stats_df.groupby(["Chrom", "Start", "End"]).agg(
            lambda x: np.var(x)
        )

    elif "est" in method:
        stats_df = stats_df.set_index(["Chrom", "Start", "End"])
        alt_freq_df = stats_df.copy()
        mut_stat = method.split("_")[-1]
        mut_wild_ext_df = mut_wild_ext_freq(alt_freq_df, freq_dict, mut=mut_stat)
        if mut_wild_ext_df is None:
            return None
        else:
            stats_df = stats_df - mut_wild_ext_df
        varscore_df = stats_df.groupby(["Chrom", "Start", "End"]).agg(
            lambda x: np.average(np.power(x, 2))
        )
    elif method == "snp_index":
        varscore_df = stats_df.groupby(["Chrom", "Start", "End"]).agg("mean")
    else:
        raise UnsupportedScoreMethod(method)
    varscore_df = varscore_df[mask]
    if method == "snp_index":
        group0, group1 = varscore_df.columns[0:2]
        varscore_df.loc[:, "snp_score"] = (
            varscore_df.loc[:, group0] - varscore_df.loc[:, group1]
        )
    else:
        varscore_df = varscore_df.applymap(lambda x: x if x >= OFFSET else OFFSET)
        varscore_df.loc[:, "snp_score"] = varscore_df.apply(log_varscore, axis=1)
    varscore_df.drop(
        [SnpGroupFreq.mut.value, SnpGroupFreq.wild.value], axis=1, inplace=True
    )
    return varscore_df


def score_plot(score_file, method, plot_title="", chr_size="", platform="local"):
    out_prefix = score_file.with_suffix(".plot")
    if method in ["var", "est_mut_alt", "est_mut_ref", "density", "ED", "density-new"]:
        out_plot = score_file.with_suffix(".plot.jpg")
    elif method == "snp_index":
        out_plot = score_file.with_suffix("")
        out_prefix = out_plot
    elif method in ["snpIndex", "Gprime"]:
        out_prefix = score_file.parent / plot_title
        out_plot = score_file.with_suffix(f".{method}.plot.jpg")
    else:
        raise UnsupportedPlot(method)
    cmd = (
        f"Rscript {SNP_SCORE_PLOT} "
        f"--input {score_file} "
        f"--output {out_prefix} "
        f"--plot_type {method} "
        f"--title {plot_title} "
        f"--chr_size {chr_size}"
    )
    if platform == "web":
        cmd = f"{cmd} --web"
    if not out_plot.exists():
        return cmd
    else:
        return None


def extract_snpeff_anno(anno_line):
    anno_stats = []
    fileds = (1, 3, 6, 9, 10)
    try:
        gene_anno = anno_line.split(";")[0]
        anno_line_stats = gene_anno.split(",")
    except Exception:
        print(anno_line)
        sys.exit("snp annotation error!")
    for annStr in anno_line_stats:
        annDetailArray = annStr.split("|")
        filed_stats = []
        for filled_i in fileds:
            filed_stats.append(annDetailArray[filled_i])
        anno_stats.append(filed_stats)
    zip_anno_stats = list(map(lambda x: "|".join(x), zip(*anno_stats)))
    return zip_anno_stats


def split_dataframe_rows(df, column_selectors, row_delimiter):
    # we need to keep track of the ordering of the columns
    def _split_list_to_rows(row, row_accumulator, column_selector, row_delimiter):
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
                    new_row[column_selector] = split_rows[column_selector].pop(0)
                except IndexError:
                    new_row[column_selector] = ""
            row_accumulator.append(new_row)

    new_rows = []
    df.apply(
        _split_list_to_rows, axis=1, args=(new_rows, column_selectors, row_delimiter)
    )
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
    pattern = re.compile("(TC[A-Z])0+(\w+)")
    if pattern.match(sample_id):
        pre, suf = pattern.match(sample_id).groups()
        return f"{pre}{suf}"
    return sample_id


def outdir_suffix_from_params(params):
    outdir_suffix_suffix = []
    for group_i in GROUPS:
        group_i_list = []
        if params.get(group_i) is None:
            continue
        sample_list = sorted(params[group_i])
        for sample_i in sample_list:
            sample_i_id = sample_i.split(".")[0]
            group_i_list.append(abbr_sample_id(sample_i_id))
        outdir_suffix_suffix.append("_".join(group_i_list))
    return "/".join(outdir_suffix_suffix)


def replace_outdir(args, chrom):
    arg_list = []
    flag = False
    for arg_i in args:
        if flag:
            arg_i = f"{arg_i}/split/{chrom}"
            flag = False
        if arg_i == "-o" or arg_i == "--outdir":
            flag = True
        arg_list.append(arg_i)
    return arg_list


def wrap_param_arg(args):
    flag = False
    for arg_i in args:
        if flag:
            arg_i = f"'{arg_i}'"
            flag = False
        if arg_i == "-p" or arg_i == "--parameters":
            flag = True
        yield arg_i


def valid_output_cols(output_cols: List[str], df: pd.DataFrame) -> List[str]:
    return [col_i for col_i in output_cols if col_i in df.columns]


def add_qtlserq_like_cols(df: pd.DataFrame, out_column: List[str]) -> pd.DataFrame:
    df.loc[:, "REF_FRQ"] = (df["mutant.REF.AD"] + df["wild.REF.AD"]) / (
        df["mutant.REF.AD"]
        + df["wild.REF.AD"]
        + df["mutant.ALT.AD"]
        + df["wild.ALT.AD"]
    )
    for grp in [each.value for each in SnpGroup.__members__.values()]:
        if f"{grp}.REF.AD" in df.columns:
            df.loc[:, f"{grp}.DP"] = df[f"{grp}.REF.AD"] + df[f"{grp}.ALT.AD"]
    df.loc[:, "AFD(deltaSNP)"] = df["mutant.FREQ"] - df["wild.FREQ"]

    df.rename(columns=COLUMN_NAME_MAP, inplace=True)
    df = df.loc[:, ~df.columns.duplicated()]
    real_out_col = valid_output_cols(out_column, df)
    return df[real_out_col]


def reformat_df(
    df: pd.DataFrame, file_name: str, chr_list: Optional[Iterable[str]] = None
) -> pd.DataFrame:
    sortby = SNP_DENSITY_SORT_COL
    if ".snp.freq.csv" in file_name:
        df = add_qtlserq_like_cols(df, SNP_BASIC_COL)
    elif "var.filter.csv" in file_name:
        out_cols = SNP_BASIC_COL + ["INFO"]
        df = add_qtlserq_like_cols(df, out_cols)
        df = flat_snpeff_ann(df, out_cols)
    elif ".var.score.ann.csv" in file_name:
        df = add_qtlserq_like_cols(df, VAR_SCORE_OUT_COL)
        sortby = VAR_SCORE_ANN_SORT_COL
    elif ".var.score.csv" in file_name:
        df.rename(columns=COLUMN_NAME_MAP, inplace=True)
        sortby = VAR_SCORE_SORT_COL
    elif ".snp.plot.bed" in file_name:
        df.loc[:, 0] = df[0].astype("str")
        sortby = [0, 1]
    else:
        pass
    if chr_list is not None:
        df = df[df[QTLSEQR_CHROM_NAME].isin(chr_list)]
    df.sort_values(sortby, inplace=True)
    return df


def merge_split_file(
    file_dir,
    file_pattern,
    chr_list: Optional[Iterable[str]] = None,
    top_rate: float = 0.05,
    out_dir=None,
    input_header="infer",
    input_sep=",",
    out_header=True,
    out_sep=",",
):
    if (
        "qtlseqr" in file_pattern
        or "snp.freq.csv" in file_pattern
        or "var.filter.csv" in file_pattern
    ):
        pattern_file = Path(file_dir).glob(f"split/*/fmt/{file_pattern}")
    else:
        pattern_file = Path(file_dir).glob(f"split/*/{file_pattern}")
    df_list = []
    for file_i in pattern_file:
        try:
            df_list.append(pd.read_csv(file_i, header=input_header, sep=input_sep))
        except pd.errors.EmptyDataError:
            logger.warning(f"File is Empty: {file_i}")
    if out_dir is None:
        out_dir = Path(file_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    outfile = Path(out_dir) / file_pattern
    df = pd.concat(df_list)
    if not (
        "qtlseqr" in file_pattern
        or "snp.freq.csv" in file_pattern
        or "var.filter.csv" in file_pattern
    ):
        df = reformat_df(df, file_pattern, chr_list=chr_list)
    else:
        df.sort_values(SNP_DENSITY_SORT_COL, inplace=True)
    if not outfile.is_file():
        if "qtlseqr" in file_pattern:
            df.to_csv(outfile, index=False)
        elif ".var.score.ann.csv" in file_pattern:
            outFilePath = outfile.with_suffix(f".top{top_rate}.csv")
            out_cols = valid_output_cols(VAR_SCORE_OUT_COL[:-1], df)
            top_ranked_results(
                df, "varBScore", out_cols, outFilePath, select_rate=top_rate
            )
        else:
            df.to_csv(
                outfile,
                index=False,
                float_format="%.3f",
                header=out_header,
                sep=out_sep,
            )
    return outfile


def format_outfile(
    filePath: Path,
    outDir: Path,
    top_rate: float = 0.05,
    merge_cols: List[str] = SNP_DENSITY_POS_COLS,
    float_format: Optional[str] = "%.3f",
    chr_list: Optional[Iterable[str]] = None,
    ann_df: Optional[pd.DataFrame] = None,
) -> Path:
    df = pd.read_csv(filePath)
    outFilePath = outDir / filePath.name
    if not outFilePath.is_file():
        if ann_df is not None:
            df = add_snp_ann(df, merge_cols=merge_cols, ann_df=ann_df)
        df = reformat_df(df, filePath.name, chr_list=chr_list)
        if df.empty:
            return outFilePath
        if ".var.score.ann.csv" in filePath.name:
            outFilePath = outFilePath.with_suffix(f".top{top_rate}.csv")
            out_cols = valid_output_cols(VAR_SCORE_OUT_COL[:-1], df)
            top_ranked_results(
                df, "varBScore", out_cols, outFilePath, select_rate=top_rate
            )
        else:
            df.to_csv(outFilePath, index=False, float_format=float_format)
    return outFilePath


def gene2pos(gene_bed, genes):
    gene_bed_df = pd.read_csv(
        gene_bed, sep="\t", header=None, names=["chrom", "start", "end"], index_col=3
    )
    for gene_i in genes:
        if gene_i in gene_bed_df.index:
            gene_i_pos = gene_bed_df.loc[gene_i]
            yield f"{gene_i_pos.chrom}:{gene_i_pos.start}-{gene_i_pos.end}"


def printdf(df):
    col_str = [str(col_i) for col_i in df.columns]
    print("\t".join(col_str))
    for index_i in df.index:
        line_i = [str(col_i) for col_i in df.loc[index_i]]
        print("\t".join(line_i))


def freq2qtlseqr(snpfreq):
    snpfreq = Path(snpfreq)
    qtlseqr_table = snpfreq.with_suffix(".qtlseqr.csv")
    if not qtlseqr_table.is_file():
        snpfreq_df = pd.read_csv(snpfreq)

        def rename_col(col_i):
            pos_map = {"Chr": "CHROM", "Pos": "POS", "Alt": "ALT"}
            if pos_map.get(col_i):
                return pos_map[col_i]
            else:
                return re.sub(r"(\w+).(\w+).AD", r"AD_\2.\1", col_i)

        snpfreq_df.columns = [rename_col(col_i) for col_i in snpfreq_df.columns]
        snpfreq_df.to_csv(qtlseqr_table, index=False)
    return qtlseqr_table


def circos_suffix(varscore_prefix, qtlseq_prefix):
    varscore_prefix = varscore_prefix.replace("mutant", "mut")
    varscore_prefix = varscore_prefix.replace("wild", "wd")
    varscore_prefix = varscore_prefix.replace("symmetrical", "sym")
    varscore_prefix = varscore_prefix.replace("snp_num.window.", "")
    qtlseq_prefix = qtlseq_prefix.replace("window_", "")
    qtlseq_prefix = qtlseq_prefix.replace("popStru_", "")
    qtlseq_prefix = qtlseq_prefix.replace("refFreq_", "")
    qtlseq_prefix = qtlseq_prefix.replace("minDepth_", "")
    return f"{varscore_prefix}-{qtlseq_prefix}"


def circos_cfg(circos_prefix, circos_path: Path = None) -> Path:
    circos_prefix.mkdir(parents=True, exist_ok=True)
    if circos_path is None:
        circos_path = circos_prefix.parent.parent
    else:
        circos_path.mkdir(parents=True, exist_ok=True)
    circos_file = f"{circos_prefix.name}.circos.png"
    # jinja2 load template
    PLOT_DIR = PurePath(__file__).parent / "plot"
    JINJA_ENV = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath=f"{PLOT_DIR}")
    )
    CIRCOS_CONF = JINJA_ENV.get_template("circos.conf")
    cfgObj = CIRCOS_CONF.render(
        {
            "circos_prefix": circos_prefix,
            "circos_path": circos_path,
            "circos_file": circos_file,
        }
    )
    cfgFile = circos_prefix / "circos.conf"
    with open(cfgFile, "w") as file_inf:
        file_inf.write(cfgObj)
    return circos_path / circos_file


def add_default_params(param_obj: dict) -> dict:
    for name, member in VarScoreParams.__members__.items():
        if not param_obj.get(name):
            param_obj[name] = member.value
    return param_obj


def add_filter_default_params(param_obj: dict) -> dict:
    for name, member in VarFilterParams.__members__.items():
        if not param_obj.get(name):
            param_obj[name] = member.value
    return param_obj


def is_same_param(param_obj: dict, param_json: Path) -> bool:
    with open(param_json) as json_inf:
        old_param = json.load(json_inf)
        return param_obj == old_param


def is_new_cmd(param_obj: dict, cmd_history_dir: Path) -> bool:
    cmd_history_dir.mkdir(parents=True, exist_ok=True)
    old_params = cmd_history_dir.glob("*.json")
    for file_i in old_params:
        if is_same_param(param_obj, file_i):
            return False
    return True


def now_str() -> str:
    return datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S")


def save_params(param_obj: dict, cmd_history_dir: Path) -> None:
    new_param_file = cmd_history_dir / f"{now_str()}.params.json"
    with open(str(new_param_file), "w") as json_inf:
        json.dump(param_obj, json_inf)


def params_cfg(
    cfg_file: Path, cfg_value: dict, cmd_history_dir: Path, cfg_type="varBscore"
) -> None:
    save_params(cfg_value, cmd_history_dir)
    cfg_value["report_time"] = now_str()
    cfg_dir = PurePath(__file__).parent / "config"
    jinja_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath=f"{cfg_dir}")
    )
    if cfg_type == "varBscore":
        template_name = "params.cfg"
    else:
        template_name = "varFilter.params.cfg"
    cfg_temp = jinja_env.get_template(template_name)
    cfg_obj = cfg_temp.render(cfg_value)
    with open(str(cfg_file), "a") as file_inf:
        file_inf.write(cfg_obj)


def circos_plot(
    varScore_csv, qtlseqr_ed_csv, snp_freq_csv, out_prefix, window_file=None
):
    if window_file:
        circos_sh = PurePath(__file__).parent / "plot" / "data2circos-plant.sh"
        cmd = f"sh {circos_sh} {varScore_csv} {qtlseqr_ed_csv} {snp_freq_csv} {window_file} {out_prefix}"
        print(cmd)
    else:
        circos_sh = PurePath(__file__).parent / "plot" / "data2circos.sh"
        cmd = f"sh {circos_sh} {varScore_csv} {qtlseqr_ed_csv} {snp_freq_csv} {out_prefix}"
    return cmd


def extract_qtlseqr_result(
    df: pd.DataFrame, selected_cols: List[str], outFile: Path
) -> None:
    if not outFile.is_file():
        out_cols = QTLSERQ_BASIC_COL + selected_cols
        df = df[out_cols].copy()
        df = sientific_number_col_to_str(df)
        df.to_csv(outFile, index=False, float_format="%.3f")


def flat_snpeff_ann(df: pd.DataFrame, out_cols: List[str]) -> pd.DataFrame:
    snpeff_anno = list(df.INFO.map(extract_snpeff_anno))
    snpeff_anno_df = pd.DataFrame(snpeff_anno)
    snpeff_anno_df.columns = [
        "Feature",
        "Gene",
        "Transcript",
        "Variant_DNA_Level",
        "Variant_Protein_Level",
    ]
    ann_df = pd.concat([df, snpeff_anno_df], axis=1)
    ann_df.drop("INFO", axis=1, inplace=True)
    ann_df = split_dataframe_rows(
        ann_df,
        column_selectors=[
            "Feature",
            "Gene",
            "Transcript",
            "Variant_DNA_Level",
            "Variant_Protein_Level",
        ],
        row_delimiter="|",
    )
    out_cols = out_cols + ANN_OUT_COLS
    real_out_cols = valid_output_cols(out_cols, ann_df)
    out_df = ann_df[real_out_cols]
    return out_df


def top_ranked_results(
    df: pd.DataFrame,
    rank_col: str,
    out_cols: List[str],
    out_file: Path,
    select_rate: float = 0.05,
) -> None:
    if not out_file.is_file():
        top_cutoff = np.quantile(sorted(df[rank_col].unique()), 1 - select_rate)
        df = df[df[rank_col] >= top_cutoff]
        df.reset_index(drop=True, inplace=True)
        out_df = flat_snpeff_ann(df, out_cols)
        out_df = sientific_number_col_to_str(out_df)
        out_df.to_csv(out_file, index=False)


def sientific_number_col_to_str(df):
    for col_i in SCIENTIFIC_NUMBER_COLS:
        if col_i in df.columns:
            df.loc[:, col_i] = df[col_i].astype("str")
    return df


def split_qtlseqr_results(
    qtlseqrFile: Path,
    qtlseqrAloneFile: Path,
    edFile: Path,
    chr_list: Optional[Iterable[str]] = None,
    top_rate: float = 0.05,
    ann_df: Optional[pd.DataFrame] = None,
) -> None:
    df = pd.read_csv(qtlseqrFile)
    df.loc[:, QTLSEQR_CHROM_NAME] = df[QTLSEQR_CHROM_NAME].astype("str")
    if chr_list is not None:
        df = df[df[QTLSEQR_CHROM_NAME].isin(chr_list)]
    ad_cols = [each for each in df.columns if "AD" in each]
    df.loc[:, "AFD(deltaSNP)"] = df["deltaSNP"]
    df.rename(columns=COLUMN_NAME_MAP, inplace=True)
    df.drop(ad_cols, axis=1, inplace=True)
    if ann_df is not None:
        df = add_snp_ann(df, merge_cols=QTLSEQR_POS_COLS, ann_df=ann_df)

    if "ED" in df.columns:
        extract_qtlseqr_result(df, ED_SPECIFIC_COLS, edFile)
        edTopFile = edFile.with_suffix(f".top{top_rate}.csv")
        out_cols = QTLSERQ_BASIC_COL + ED_SPECIFIC_COLS
        top_ranked_results(df, "fitted", out_cols, edTopFile, top_rate)

    if "Gprime" in df.columns:
        extract_qtlseqr_result(df, QTLSEQR_SPECIFIC_COLS, qtlseqrAloneFile)
        qtlseqrTopFile = qtlseqrAloneFile.with_suffix(f".top{top_rate}.csv")
        out_cols = QTLSERQ_BASIC_COL + QTLSEQR_SPECIFIC_COLS
        top_ranked_results(df, "Gprime", out_cols, qtlseqrTopFile, top_rate)


def snp_density_stats(
    window_bed: PurePath, snp_density_bed: Path, density_stats_file: Path
) -> None:
    if not density_stats_file.is_file():
        window_bed = BedTool(str(window_bed))
        snp_bed = BedTool(str(snp_density_bed))
        cov_res = window_bed.coverage(snp_bed, counts=True, sorted=True)
        cov_str = StringIO(str(cov_res))
        cov_df = pd.read_csv(
            cov_str,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "variantCount"],
        )
        cov_df.to_csv(density_stats_file, index=False)


def cp_if_not_exist(fileItem: Path, outPath: Path) -> None:
    outFile = outPath / fileItem.name
    if not outFile.exists():
        shutil.copy(fileItem, outPath)


def cp_files(fileList: List[Path], outPath: Path) -> None:
    for file_i in fileList:
        cp_if_not_exist(file_i, outPath)


def window_number_format(number):
    megabase = 1000 * 1000
    kilobase = 1000
    if number >= megabase:
        return f"{number / megabase}M"
    elif number >= kilobase:
        return f"{int(number / kilobase)}K"
    else:
        return str(number)


def make_chr_window(chr_size: Path, window: int, step: int, outPath: Path) -> Path:
    if step:
        fileName = f"w{window_number_format(window)}.s{window_number_format(step)}.bed"
    else:
        fileName = f"w{window_number_format(window)}.bed"
    window_file = outPath / fileName
    if not window_file.is_file():
        step_param = f"-s {step}" if step else ""
        cmd = f"bedtools makewindows -g {chr_size} -w {window} {step_param} | bedtools sort -i - > {window_file}"
        delegator.run(cmd)
        if not window_file.is_file():
            raise SnpDensityWindowFailed
    return window_file


def add_snp_ann(
    df: pd.DataFrame, merge_cols: List[str], ann_df: pd.DataFrame
) -> pd.DataFrame:
    chrom_col = merge_cols[0]
    df.loc[:, chrom_col] = df[chrom_col].astype("str")
    return df.merge(ann_df, left_on=merge_cols, right_on=ANN_POS_COLS)


class CheckOutPutFunc:
    def __init__(self, output: Path, func: Callable):
        self.output = output
        self.func = func

    def __call__(self, *args, **kwargs):
        if not self.output.exists():
            return self.func(*args, **kwargs)


def check_output(outFile: Path):
    return functools.partial(CheckOutPutFunc, outFile)


def filter_by_freq(df, column, freq):
    return (df[column] <= freq) | (df[column] >= 1 - freq)


def filter_by_afd(df, column, afd, deviation):
    return (df[column].abs() >= afd - deviation) & (df[column].abs() <= afd + deviation)


def non_score_filter_snp(
    alt_freq_df, mut_freq, wild_freq, afd, afd_deviation, isParent=False
):
    if not isParent:
        mut_name = SnpGroupFreq.mut.value
        wild_name = SnpGroupFreq.wild.value
        afd_name = "AFD"
    else:
        mut_name = SnpGroupFreq.mut_pa.value
        wild_name = SnpGroupFreq.wild_pa.value
        afd_name = "P_AFD"
    # filter mut freq
    mut_freq_mask = filter_by_freq(alt_freq_df, mut_name, mut_freq)
    # filter wild freq
    wild_freq_mask = filter_by_freq(alt_freq_df, wild_name, wild_freq)
    # filter by afd
    afd_mask = filter_by_afd(alt_freq_df, afd_name, afd, afd_deviation)
    return mut_freq_mask & wild_freq_mask & afd_mask


def chrom_bin_snp_number_df(
    start: int, window: int, chr_len: int, chrom: str, df: pd.DataFrame
) -> pd.DataFrame:
    cut_range = range(start, chr_len + window, window)
    range_count_df = pd.DataFrame(
        pd.cut(df[SnpTableConstants.ANN_TABLE_COL[1]], cut_range)
        .value_counts()
        .sort_index()
    )
    range_count_df.columns = [SnpDensityStatsTable.COUNT]
    range_count_df.loc[:, SnpDensityStatsTable.CHROM] = chrom
    range_count_df.loc[:, SnpDensityStatsTable.START] = [
        each.left for each in range_count_df.index
    ]
    range_count_df.loc[:, SnpDensityStatsTable.END] = [
        each.right for each in range_count_df.index
    ]
    return range_count_df[
        [
            SnpDensityStatsTable.CHROM,
            SnpDensityStatsTable.START,
            SnpDensityStatsTable.END,
            SnpDensityStatsTable.COUNT,
        ]
    ].reset_index(drop=True)


def var_density_stats(
    chr_size_file: Path,
    ann_alt_freq_df: pd.DataFrame,
    window: int = 1000 * 1000,
    step: Optional[int] = None,
) -> pd.DataFrame:
    chr_size_df = pd.read_csv(chr_size_file, sep="\t", index_col=0, names=["chr_len"])
    stats_df_list = []
    ann_alt_freq_df.drop_duplicates(
        subset=SnpTableConstants.ANN_TABLE_COL[:2], inplace=True
    )
    for chrom, chrom_df in ann_alt_freq_df.groupby(
        [SnpTableConstants.ANN_TABLE_COL[0]]
    ):
        chr_len = chr_size_df.loc[chrom, "chr_len"]
        if step is None:
            step = window
        for start in range(0, window, step):
            stats_df_i = chrom_bin_snp_number_df(
                start=start,
                chr_len=chr_len,
                window=window,
                chrom=str(chrom),
                df=chrom_df,
            )
            stats_df_list.append(stats_df_i)
    return pd.concat(stats_df_list)


def var_density_file_suffix(window: int, step: Optional[int]):
    window_str = window_number_format(window)
    if step:
        step_str = window_number_format(step)
        return f".window{window_str}_step{step_str}.density.csv"
    else:
        return f".window{window_str}.density.csv"


def has_parent(group: List[str]) -> bool:
    return SnpGroup.mut_pa.value in group


def table2vcf(input_df: pd.DataFrame, vcf_file: Path, column_map: dict) -> None:
    vcf_header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    df = input_df.copy()
    df.loc[:, "ID"] = "."
    df.loc[:, "QUAL"] = "0"
    df.loc[:, "FILTER"] = "PASS"
    df.loc[:, "INFO"] = "."
    df.rename(columns=column_map, inplace=True)
    df.to_csv(vcf_file, sep="\t", index=False, columns=vcf_header)


def get_info_ann(info: str) -> str:
    return info.split("ANN=")[1]


def load_snpeff(vcf_file: str) -> pd.DataFrame:
    ann_df = pd.read_csv(vcf_file, sep="\t", header=None, comment="#", dtype={0: str})
    ann_df = ann_df.loc[:, [0, 1, 3, 4, 7]]
    ann_df.columns = ANN_TABLE_COL
    ann_df.loc[:, "INFO"] = ann_df.INFO.map(get_info_ann)
    return ann_df


def extract_snpeff_annotation(vcf_file: str, annotation_file: Path) -> pd.DataFrame:
    ann_df = load_snpeff(vcf_file)
    ann_df.to_pickle(annotation_file)
    return ann_df


def table2annotation_df(
    input_table: str,
    column_map: dict,
    snpeff_cfg: str,
    snpeff_db: str,
    annotation_dir: Union[str, Path],
    prefix: str,
) -> pd.DataFrame:
    annotation_dir = Path(annotation_dir).absolute()
    annotation_file = annotation_dir / f"{prefix}.{ANN_TABLE_NAME}"
    if not annotation_file.is_file():
        input_table = Path(input_table)
        input_table_df = pd.read_csv(input_table)
        alt_freq_vcf = input_table.with_suffix(".vcf")
        alt_ann_vcf = alt_freq_vcf.with_suffix(".ann.vcf.gz")
        table2vcf(input_table_df, alt_freq_vcf, column_map)
        annotate_cmd = (
            f"snpEff -Xmx4g -c {snpeff_cfg} "
            f"-v {snpeff_db} -quiet {alt_freq_vcf} "
            f"| bgzip > {alt_ann_vcf}"
        )
        delegator.run(annotate_cmd)
        return extract_snpeff_annotation(alt_ann_vcf, annotation_file)
    annotation_df = pd.read_pickle(annotation_file)
    annotation_df.loc[:, "#CHROM"] = annotation_df["#CHROM"].astype("str")
    return annotation_df
