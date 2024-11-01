import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from time import time
from Bio.Seq import Seq
from helper_functions import read_FASTQ_file, set_aligner_parameters
from learning_algorithm import run_learning_algorithm
from search_functions import search_seq_and_write, MatchResult
from config import MAX_WORKERS, DIRECTORIES, PRIMERS, INDEX_LEN, DEFAULT_SCORES, MACHINE_LEARNING_MODE, DEBUG_MODE, \
                    RUNNING_ON_WINDOWS

INVALID_FIELD = -1


class Statistics:
    def __init__(self):
        self.seq_count = 0
        self.filtered_out_count = 0
        self.valid_seq_count = 0
        self.invalid_seq_count = 0
        self.percent_filtered = 0.0
        self.loss = INVALID_FIELD
        self.file_name = ""


def process_file(file_index, fastq_file, aligner, min_valid_score, compare_aligner,
                 compare_min_valid_score):
    """
    Process a single FASTQ file, applying sequence filtering and comparing results to naive filtering.

    :param file_index: Index of the file being processed.
    :param fastq_file: Name of the FASTQ file.
    :param aligner: Aligner for local alignment.
    :param min_valid_score: Minimum valid score for alignment.
    :param compare_aligner: Aligner for global alignment comparison.
    :param compare_min_valid_score: Minimum valid score for comparison.
    :return: Statistics object for the file.
    """
    print(f"---Processing file #{file_index}: {fastq_file}")

    file_stats = Statistics()
    file_start_time = time()

    input_dir = DIRECTORIES['input_dir']
    output_with_primers_dir = DIRECTORIES['output_with_primers_dir']
    output_wo_primers_dir = DIRECTORIES['output_wo_primers_dir']
    debug_dir = DIRECTORIES['debug_dir']

    input_path = os.path.join(input_dir, fastq_file)

    base_filename = os.path.splitext(fastq_file)[0]
    output_with_primers_file_path = os.path.join(output_with_primers_dir, f"{base_filename}_with_primers.fastq")
    output_wo_primers_file_path = os.path.join(output_wo_primers_dir, f"{base_filename}_wo_primers.fastq")
    debug_path = os.path.join(debug_dir, f"{base_filename}_debug.txt")
    compare_params = {"aligner": compare_aligner, "min_valid_score": compare_min_valid_score}

    sequences_and_scores_dict, seq_count = read_FASTQ_file(input_path)
    filtered_out_count = 0

    # for learning algorithm
    valid_seq_count = 0
    invalid_seq_count = 0

    output_with_primers_file = None
    output_wo_primers_file = None
    debug_file = None
    if not MACHINE_LEARNING_MODE:
        output_with_primers_file = open(output_with_primers_file_path, "w")
        output_wo_primers_file = open(output_wo_primers_file_path, "w")

    if DEBUG_MODE:
        debug_file = open(debug_path, "w")
        debug_file.write(f"-------- Start of file #{file_index} --------\n")

    for seq_input in sequences_and_scores_dict.keys():
        sequence = Seq(seq_input)
        output_files = {"wo_primers": output_wo_primers_file, "with_primers": output_with_primers_file}

        search_res = search_seq_and_write(aligner, min_valid_score, compare_params, sequence, output_files, debug_file)
        if search_res == MatchResult.FOUND_VALID_MATCH:
            valid_seq_count += 1
        elif search_res == MatchResult.FOUND_INVALID_MATCH:
            invalid_seq_count += 1
        else:
            filtered_out_count += 1

    execution_time = time() - file_start_time

    if DEBUG_MODE:
        debug_file.write(f"Filtered out (we want it to be small){filtered_out_count} out of {seq_count} sequences\n")
        debug_file.write(f"Execution time: {execution_time: .2f} sec\n")
        debug_file.write(f"-------- End of file number {file_index}--------\n")
        debug_file.close()

    print(f"----Finished processing file #{file_index}, Execution time: {execution_time: .2f} sec.\n"
          f"    Filtered total: {filtered_out_count} out of {seq_count} sequences")

    if MACHINE_LEARNING_MODE:
        file_stats.seq_count = seq_count
        file_stats.filtered_out_count = filtered_out_count
        file_stats.valid_seq_count = valid_seq_count
        file_stats.invalid_seq_count = invalid_seq_count
        file_stats.file_name = base_filename
        if seq_count != 0:
            file_stats.percent_filtered = (filtered_out_count/seq_count)*100
    else:
        output_with_primers_file.close()
        output_wo_primers_file.close()

    return file_stats


def run_filtering(scores=DEFAULT_SCORES):
    """
    Run the filtering process with given parameters across multiple files using parallel processing.

    :param scores: Scores for the alignment.
    :return: eval_stats, file_stats_list - Statistics for the evaluation and for each file.
    """
    start_time = time()

    input_dir = DIRECTORIES['input_dir']
    aligner, min_valid_score = set_aligner_parameters(scores=scores, mode="local",
                                                      wanted_match_percentage=0.85, length=len(PRIMERS['front']))
    compare_aligner, compare_min_valid_score = set_aligner_parameters(scores=scores, mode="global",
                                                                      wanted_match_percentage=0.85, length=INDEX_LEN)

    eval_stats = Statistics()
    file_stats_list = []

    # find the sequences in the files
    print(f"--Starting to iterate over all files in {input_dir}")
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [
            executor.submit(process_file, file_index + 1, fastq_file, aligner, min_valid_score, compare_aligner,
                            compare_min_valid_score)
            for file_index, fastq_file in enumerate(os.listdir(input_dir)) if fastq_file.endswith(".fastq")
        ]

        if MACHINE_LEARNING_MODE:
            for future in as_completed(futures):
                try:
                    file_stats = future.result()
                    file_stats_list.append(file_stats)
                    eval_stats.seq_count += file_stats.seq_count
                    eval_stats.filtered_out_count += file_stats.filtered_out_count
                    eval_stats.valid_seq_count += file_stats.valid_seq_count
                    eval_stats.invalid_seq_count += file_stats.invalid_seq_count
                    if eval_stats.seq_count != 0:  # if it is 0, it indicates an error
                        eval_stats.percent_filtered = (eval_stats.filtered_out_count / eval_stats.seq_count) * 100
                except Exception as exc:
                    print(f"---Error processing file: {exc}")

    print("-----Finished iterating over all files")
    print(f"-----Total execution time: {time() - start_time: .2f} sec")

    return eval_stats, file_stats_list


def main():
    if MACHINE_LEARNING_MODE:
        os.makedirs(DIRECTORIES["learning_dir"], exist_ok=True)

        if RUNNING_ON_WINDOWS:
            print("Running on Windows. Adjusting for Windows environment...")
            from multiprocessing import freeze_support  # For Windows compatibility
            freeze_support()  # Necessary for Windows multiprocessing

        run_learning_algorithm()
    else:
        run_filtering()


if __name__ == "__main__":
    main()
