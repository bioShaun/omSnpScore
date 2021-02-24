import pkg_resources
import numpy as np
from enum import Enum, IntEnum
from pathlib import PurePath, Path


class SnpGroup(Enum):
    mut = "mutant"
    wild = "wild"
    mut_pa = "mutant_parent"
    wild_pa = "wild_parent"
    bg = "background"


class SnpGroupFreq(Enum):
    mut = "mutant.FREQ"
    wild = "wild.FREQ"
    mut_pa = "mutant_parent.FREQ"
    wild_pa = "wild_parent.FREQ"
    bg = "background.FREQ"


class SnpRep(IntEnum):
    alt = 0
    ref = 1
    unkown = 2


class VarScoreParams(Enum):
    min_depth = 5
    snp_number_window = 20
    snp_number_step = 5
    ref_freq = 0.3
    p_ref_freq = 0.3
    background_ref_freq = 0.3
    qtlseqr_window = 10_000_000
    qtlseqr_ref_freq = 0.3
    pop_stru = "RIL"
    qtlseqr_min_depth = 5
    snp_density_window = 1_000_000
    snp_density_step = 0


class VarFilterParams(Enum):
    min_depth = 5
    mut_freq = 0.4
    wild_freq = 0.4
    p_mut_freq = 0
    p_wild_freq = 0
    afd_deviation = 0.05
    p_afd = 1
    p_afd_deviation = 0.05


script_dir = Path(__file__).parent

SNP_SCORE_PLOT = script_dir / "plot/snpScorePlot.R"
QTLSEQR_PLOT = script_dir / "plot/run_qtlseqr.R"
QTLSEQR_PLOT_WEB = script_dir / "plot/run_qtlseqr_web.R"
VAR_DENSITY_PLOT = script_dir / "plot/var_density_compare.R"

DATA_DIR = PurePath(pkg_resources.resource_filename("snpScore", "data"))
CHR_SIZE = DATA_DIR / "chr.size"
CHR_WINDOW = DATA_DIR / "chr.1m.window.bed"

OFFSET = 1e-05
SNP_FREQ_BIAS = 0.1
ALT_FREQ = np.round(2 / 3 - SNP_FREQ_BIAS, 2)
REF_FREQ = np.round(1 / 3 + SNP_FREQ_BIAS, 2)
VCF_SAMPLE_INDEX = 9

GROUPS = tuple([grp_i.value for grp_i in SnpGroup.__members__.values()])
MUT_NAME = SnpGroup.mut.value
WILD_NAME = SnpGroup.wild.value

# position cols
SNP_DENSITY_POS_COLS = ["Chr", "Pos", "Alt"]
ANN_POS_COLS = ["#CHROM", "POS", "ALT"]
QTLSEQR_POS_COLS = ["CHROM", "POS", "ALT"]
SNP_DENSITY_SORT_COL = ["CHROM", "POS"]
# OUTPUT TABLE

COLUMN_NAME_MAP = {
    "Chr": "CHROM",
    "Chrom": "CHROM",
    "Pos": "POS",
    "Alt": "ALT",
    "P_AFD": "PARENT_AFD",
    "snp_score": "varBScore",
    "mutant.FREQ": "mutant.AF",
    "wild.FREQ": "wild.AF",
    "mutant_parent.FREQ": "mutant_parent.AF",
    "wild_parent.FREQ": "wild_parent.AF",
    "background.FREQ": "background.AF",
    "LOW.FREQ": "wild.AF",
    "DP.LOW": "wild.DP",
    "SNPindex.LOW": "wild.SNPindex",
    "HIGH.FREQ": "mutant.AF",
    "DP.HIGH": "mutant.DP",
    "SNPindex.HIGH": "mutant.SNPindex",
    "REF_FRQ": "REF_FRQ",
    "euc": "ED",
}

SNP_BASIC_COL = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "wild.AF",
    "mutant.AF",
    "wild_parent.AF",
    "mutant_parent.AF",
    "background.AF",
    "AFD(deltaSNP)",
    "PARENT_AFD",
    "REF_FRQ",
    "wild.DP",
    "mutant.DP",
    "wild_parent.DP",
    "mutant_parent.DP",
    "background.DP",
]

QTLSERQ_BASIC_COL = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "wild.AF",
    "mutant.AF",
    "AFD(deltaSNP)",
    "REF_FRQ",
    "wild.DP",
    "mutant.DP",
]

# varscore
VAR_SCORE_OUT_COL = [
    "CHROM",
    "Start",
    "End",
    "varBScore",
    "wild.AF",
    "mutant.AF",
    "wild_parent.AF",
    "mutant_parent.AF",
    "background.AF",
    "AFD(deltaSNP)",
    "REF_FRQ",
    "wild.DP",
    "mutant.DP",
    "wild_parent.DP",
    "mutant_parent.DP",
    "background.DP",
    "POS",
    "REF",
    "ALT",
    "INFO",
]

VAR_SCORE_ANN_SORT_COL = ["CHROM", "Start", "POS"]
VAR_SCORE_SORT_COL = ["CHROM", "Start"]
## qtlseqr
QTLSEQR_CHROM_NAME = "CHROM"

SCIENTIFIC_NUMBER_COLS = [
    "pvalue",
    "negLog10Pval",
    "qvalue",
    "ED",
    "fitted",
    "unfitted",
    "dis2edcutoff",
]

ED_SPECIFIC_COLS = ["ED", "fitted", "unfitted", "dis2edcutoff"]

QTLSEQR_SPECIFIC_COLS = [
    "nSNPs",
    "tricubeDeltaSNP",
    "Gprime",
    "negLog10Pval",
    "qvalue",
    "minDP",
    "CI_95",
    "CI_99",
]

# ANNOTATION OUT
ANN_OUT_COLS = ["Feature", "Gene", "Transcript"]


# OUTPUT DIR
class VarScoreOutDirName(Enum):
    snp_density = "1.SNP_Density"
    var_score = "2.varBScore"
    ed = "3.ed"
    qtlseqr = "4.QTLseqr"
    circos = "5.circosPlot"


# result readme
DOC_DIR = script_dir / "doc"


class VarScoreDocName(Enum):
    snp_density = DOC_DIR / "SNP_Density.readme.txt"
    var_score = DOC_DIR / "varBScore.readme.txt"
    ed = DOC_DIR / "ed.readme.txt"
    qtlseqr = DOC_DIR / "qtlseqr.readme.txt"
    circos = DOC_DIR / "circos.readme.txt"


class SnpTableConstants:
    ANN_TABLE_COL = ["CHROM", "POS", "REF", "ALT", "INFO"]


class SnpDensityStatsTable:
    CHROM = "chrom"
    START = "start"
    END = "end"
    COUNT = "variantCount"


ANN_TABLE_COL = ("#CHROM", "POS", "REF", "ALT", "INFO")

ANN_TABLE_NAME = "snp.ann.table.pkl"

VAR_TO_VCF_COLUMN_MAP = {"Chr": "#CHROM", "Pos": "POS", "Ref": "REF", "Alt": "ALT"}
QTLSEQR_TO_VCF_COLUMN_MAP = {"CHROM": "#CHROM", "Ref": "REF"}
