import platform
from Bio.Seq import Seq
import os
from hyperopt import hp

# General Configuration
RUNNING_ON_WINDOWS = (platform.system() == "Windows")
MAX_WORKERS_MULTIPLIER = 2.5
MAX_WORKERS = int(os.cpu_count() * MAX_WORKERS_MULTIPLIER)


# File/Directory Paths
DIR_START = "learning_midgam_final" if not RUNNING_ON_WINDOWS else "midgam"
RESULTS = "filtering_results"
READS = "raw_reads"
SEQUENCE_OUTPUT = "sequence_output"
SEQUENCE_OUTPUT_WITH_PRIMERS = "with_primers"
SEQUENCE_OUTPUT_WO_PRIMERS = "without_primers"
DEBUG = "debug"
DATA_SET_FILE = "data_after_indices_two_files.txt"
LEARNING = "learning_info"

DIRECTORIES = {
    "input_dir": os.path.join(DIR_START, READS),
    "output_with_primers_dir": os.path.join(DIR_START, RESULTS, SEQUENCE_OUTPUT, SEQUENCE_OUTPUT_WITH_PRIMERS),
    "output_wo_primers_dir": os.path.join(DIR_START, RESULTS, SEQUENCE_OUTPUT, SEQUENCE_OUTPUT_WO_PRIMERS),
    "debug_dir": os.path.join(DIR_START, RESULTS, DEBUG),
    "learning_dir": os.path.join(DIR_START, RESULTS, LEARNING),
}


# Alignment Parameters
INDEX_LEN = 12
LIB_LENGTH = 140

# Primers
FRONT_PRIMER = Seq("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
BACK_PRIMER = Seq("CTGTCTCTTATACACATCTCCGAGCCCACGAGAC")
PRIMERS = {"front": FRONT_PRIMER, "back": BACK_PRIMER}

# Default Scores
DEFAULT_SCORES = {
    "match_score": 4.254276492,
    "mismatch_score": -3.515462288,
    "open_gap_score": -4.032420994,
    "extend_gap_score": -2.999836807
}


# Machine Learning Settings
MACHINE_LEARNING_MODE = True
EVAL_NUM = 50
DATA_SET_PATH = os.path.join(DIR_START, DATA_SET_FILE)

DEFAULT_SEARCH_SPACE = {
 'match_score': hp.uniform('match_score', 1, 5),
 'mismatch_score': hp.uniform('mismatch_score', -5, -1),
 'open_gap_score': hp.uniform('open_gap_score', -5, -1),
 'extend_gap_score': hp.uniform('extend_gap_score', -5, -1)
}

# Debug Settings
DEBUG_MODE = True


# Creating Important Directories
if not MACHINE_LEARNING_MODE:
    os.makedirs(DIRECTORIES["output_with_primers_dir"], exist_ok=True)
    os.makedirs(DIRECTORIES["output_wo_primers_dir"], exist_ok=True)

if DEBUG_MODE:
    os.makedirs(DIRECTORIES["debug_dir"], exist_ok=True)