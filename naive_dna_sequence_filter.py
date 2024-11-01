from time import time
from helper_functions import (read_FASTQ_file, copy_reverse_complement)
from config import LIB_LENGTH, PRIMERS
import os


# ******************* Helper Functions *******************

def naive_search_and_write_seq(sequence, front_primer, back_primer):
    """
    Search for valid sequences containing front and back primers within a given sequence.

    :param sequence: DNA sequence to search.
    :param front_primer: The front primer sequence.
    :param back_primer: The back primer sequence.
    :return: Boolean indicating whether a valid match was found.
    """
    front_index = sequence.find(front_primer)
    back_index = sequence.find(back_primer)

    if front_index != -1 and back_index != -1:
        seq_without_primers = sequence[front_index + len(front_primer):back_index]
        if abs(len(seq_without_primers) - LIB_LENGTH) <= 5:
            return True

    elif front_index != -1:
        seq_without_primers = sequence[front_index + len(front_primer):min(len(sequence), front_index + len(front_primer) + LIB_LENGTH)]
        if abs(len(seq_without_primers) - LIB_LENGTH) <= 5:
            return True

    elif back_index != -1:
        seq_without_primers = sequence[max(0, back_index - LIB_LENGTH):back_index]
        if abs(len(seq_without_primers) - LIB_LENGTH) <= 5:
            return True

    return False

# ******************* Filtering Functions *******************


def naive_filtering_of_file(input_file_name, primers):
    """
    Perform naive filtering on a FASTQ file to count the number of filtered sequences.

    :param input_file_name: Path to the input FASTQ file.
    :param primers: Dictionary containing the front and back primers.
    :return: Tuple containing the number of filtered sequences and the total number of sequences.
    """
    sequences = read_FASTQ_file(input_file_name)[0].keys()
    start_time = time()
    count_filtered = 0

    front_primer = primers['front']
    back_primer = primers['back']
    front_primer_rc = primers['front_rc']
    back_primer_rc = primers['back_rc']

    for sequence in sequences:
        if not naive_search_and_write_seq(sequence, front_primer, back_primer):
            if not naive_search_and_write_seq(sequence, front_primer_rc, back_primer_rc):
                count_filtered += 1

    # Uncomment for timing and debug info
    # print(f"Execution time: {(time() - start_time): .2f} sec")
    # print("Filtered sequences count:", count_filtered)

    return count_filtered, len(sequences)


def naive_filtering_percent_of_dir(input_dir):
    """
    Perform naive filtering on all FASTQ files in a directory and calculate the percentage of filtered sequences.

    :param input_dir: Path to the input directory containing FASTQ files.
    :return: Percentage of naively filtered sequences.
    """
    total_filtered = 0
    total_sequences = 0

    primers = {
        'front': str(PRIMERS['front']),
        'back': str(PRIMERS['back']),
        'front_rc': str(copy_reverse_complement(PRIMERS['back'])),
        'back_rc': str(copy_reverse_complement(PRIMERS['front']))
    }

    for file in os.listdir(input_dir):
        if file.endswith(".fastq"):
            filtered_count, sequence_count = naive_filtering_of_file(os.path.join(input_dir, file), primers)
            total_filtered += filtered_count
            total_sequences += sequence_count

    assert(total_sequences > 0)
    return (total_filtered / total_sequences) * 100
