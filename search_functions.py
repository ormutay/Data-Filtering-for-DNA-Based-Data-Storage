from helper_functions import (get_alignment_positions, read_sequences_from_data_set, percentage_to_min_score,
                              copy_reverse_complement, get_valid_alignments)
from Bio.Seq import Seq
from enum import Enum, auto
from config import INDEX_LEN, LIB_LENGTH, PRIMERS, MACHINE_LEARNING_MODE, DEBUG_MODE


class MatchResult(Enum):
    """Enumeration for different algorithm types used for hyperparameter optimization."""
    FOUND_VALID_MATCH = auto()
    FOUND_INVALID_MATCH = auto()
    NO_MATCH_FOUND = auto()

# ******************* DataSet *******************


def compare_seq_to_dataset(compare_params, sequence, debug_file):
    """
    Compare a given sequence to a dataset and find the closest match using the provided aligner.

    :param compare_params: Dictionary containing compare_aligner and compare_min_valid_score.
    :param sequence: Sequence to compare against the dataset.
    :param debug_file: File object to write debug information.
    :return: MatchResult indicating whether a valid sequence is found that meets the minimum alignment score.
    """
    data_set_sequences = read_sequences_from_data_set()
    min_valid_score = percentage_to_min_score(compare_params['aligner'], 0.85, len(sequence))
    sequence_index = sequence[:INDEX_LEN]  # Take the first 12 bases of the sequence

    for data_set_seq in data_set_sequences:
        if len(data_set_seq) == 0:
            continue

        data_set_seq_index_as_string = data_set_seq[:INDEX_LEN]
        data_set_seq_index = Seq(data_set_seq_index_as_string)

        # Align the sequence index with the dataset sequence index
        alignment_idx = compare_params['aligner'].align(sequence_index, data_set_seq_index)[0]

        if alignment_idx.score >= compare_params['min_valid_score']:
            alignment_seq = compare_params['aligner'].align(sequence, data_set_seq)[0]

            if DEBUG_MODE:
                # Streamlined debug prints for consistency
                debug_msg = (
                    f"Index: {data_set_seq_index}\n"
                    f"Index alignment score: {alignment_idx.score}, minimum valid score is {compare_params['min_valid_score']}\n"
                    f"Best index alignment:\n{alignment_idx}\n"
                    "------------------\n"
                )
                debug_file.write(debug_msg)

            if min_valid_score <= alignment_seq.score:
                if DEBUG_MODE:
                    debug_full_msg = (
                        "---------Full sequence details:---------\n"
                        f"The current sequence (full): {sequence}\n"
                        f"Matching sequence out of DB (full): {data_set_seq}\n"
                        f"Sequence alignment score (full): {alignment_seq.score}, minimum valid score is {min_valid_score}\n"
                        f"Alignment:\n{alignment_seq}\n"
                        "\n" + "="*64 + "\n\n"
                    )
                    debug_file.write(debug_full_msg)
                return MatchResult.FOUND_VALID_MATCH

    if DEBUG_MODE:
        debug_file.write(
            "!!!!!!!!!!!!!! No valid sequence found !!!!!!!!\n"
            "\n" + "="*64 + "\n\n"
        )
    return MatchResult.FOUND_INVALID_MATCH

# ******************* Both primers *******************


def search_both_primers_helper(valid_alignments, sequence, output_files, compare_params, debug_file, is_reversed):
    """
    Search for sequences containing both front and back primers within a valid distance, considering reverse complement
    if needed.

    :param valid_alignments: Valid alignments for the regular/rc orientation.
    :param sequence: The full DNA sequence to search.
    :param output_files: Dictionary containing output_wo_primers and output_with_primers file objects.
    :param compare_params: Dictionary containing aligner and min_valid_score for comparison.
    :param debug_file: File object for writing debug logs.
    :param is_reversed: Boolean indicating if the sequence is reverse complemented.
    :return: MatchResult indicating whether a match was found.
    """

    valid_alignments_front, valid_alignments_back = valid_alignments.values()

    for f_align in valid_alignments_front:
        for b_align in valid_alignments_back:
            f_start, f_end = get_alignment_positions(f_align)
            b_start, b_end = get_alignment_positions(b_align)
            distance = b_start - f_end

            if LIB_LENGTH - 5 <= distance <= LIB_LENGTH + 5:
                seq_without_primers = sequence[f_end:b_start]
                seq_with_primers = sequence[f_start:b_end]

                if is_reversed:
                    seq_without_primers = copy_reverse_complement(seq_without_primers)
                    seq_with_primers = copy_reverse_complement(seq_with_primers)

                if not MACHINE_LEARNING_MODE:
                    output_files['wo_primers'].write(str(seq_without_primers) + "\n")
                    output_files['with_primers'].write(str(seq_with_primers) + "\n")

                if MACHINE_LEARNING_MODE or DEBUG_MODE:
                    return compare_seq_to_dataset(compare_params, seq_without_primers, debug_file)
                return MatchResult.FOUND_VALID_MATCH

    return MatchResult.NO_MATCH_FOUND


def search_both_primers(valid_alignments, valid_alignments_rc, sequence, output_files, compare_params, debug_file):
    """
    Searches for valid matches with both primers (front and back) in regular and reverse orientations.

    :param valid_alignments: Valid alignments for the regular orientation.
    :param valid_alignments_rc: Valid alignments for the reverse complement orientation.
    :param sequence: The sequence being analyzed.
    :param output_files: File handles for writing sequences.
    :param compare_params: Comparison parameters for alignment.
    :param debug_file: File handle for logging debug information.
    :return: MatchResult indicating the result of the search.
    """
    # Search for both primers - regular
    search_res = search_both_primers_helper(valid_alignments, sequence, output_files, compare_params,
                                            debug_file, False)
    if search_res != MatchResult.NO_MATCH_FOUND:
        return search_res

    # Search for both primers - reversed
    search_res_rc = search_both_primers_helper(valid_alignments_rc, sequence, output_files, compare_params,
                                               debug_file, True)
    if search_res_rc != MatchResult.NO_MATCH_FOUND:
        return search_res_rc

    return MatchResult.NO_MATCH_FOUND

# ******************* Single primer *******************


# f_end/_b_end are the index after the primer
# front: start_pos = f_end, end_pos = f_end + LIB_LENGTH
def search_single_front_primer_helper(valid_alignments_front, sequence, output_files, compare_params,
                                      debug_file, is_reversed):
    """
    Search for sequences with a valid front primer match.

    :param valid_alignments_front: List of valid front primer alignments.
    :param sequence: Full DNA sequence to search.
    :param output_files: Dictionary containing output_wo_primers and output_with_primers file objects.
    :param compare_params: Dictionary containing aligner and min_valid_score for comparison.
    :param debug_file: File for logging debug information.
    :param is_reversed: Boolean indicating if the sequence is reverse complemented.
    :return: MatchResult indicating whether a valid match was found.
    """
    for f_align in valid_alignments_front:
        f_start, f_end = get_alignment_positions(f_align)
        start_pos = f_end
        end_pos = start_pos + LIB_LENGTH
        seq_with_primers = sequence[f_start:end_pos]

        # Skip if out of bounds
        if end_pos >= len(sequence):
            continue

        seq_without_primers = sequence[start_pos:end_pos]

        if is_reversed:
            seq_without_primers = copy_reverse_complement(seq_without_primers)
            seq_with_primers = copy_reverse_complement(seq_with_primers)

        if not MACHINE_LEARNING_MODE:
            output_files['wo_primers'].write(str(seq_without_primers) + "\n")
            output_files['with_primers'].write(str(seq_with_primers) + "\n")

        if MACHINE_LEARNING_MODE or DEBUG_MODE:
            return compare_seq_to_dataset(compare_params, seq_without_primers, debug_file)
        return MatchResult.FOUND_VALID_MATCH

    return MatchResult.NO_MATCH_FOUND


# back: start_pos = b_start - LIB_LENGTH, end_pos = b_start
def search_single_back_primer_helper(valid_alignments_back, sequence, output_files, compare_params,
                                     debug_file, is_reversed):
    """
    Search for sequences with a valid back primer match.

    :param valid_alignments_back: List of valid back primer alignments.
    :param sequence: Full DNA sequence to search.
    :param output_files: Dictionary containing output_wo_primers and output_with_primers file objects.
    :param compare_params: Dictionary containing aligner and min_valid_score for comparison.
    :param debug_file: File for logging debug information.
    :param is_reversed: Boolean indicating if the sequence is reverse complemented.
    :return: MatchResult indicating whether a valid match was found.
    """
    for b_align in valid_alignments_back:
        b_start, b_end = get_alignment_positions(b_align)
        end_pos = b_start
        start_pos = end_pos - LIB_LENGTH
        seq_with_primers = sequence[start_pos:b_end]

        # Skip if out of bounds
        if start_pos < 0:
            continue

        seq_without_primers = sequence[start_pos:end_pos]

        # Apply ternary condition for reversed sequences
        if is_reversed:
            seq_without_primers = copy_reverse_complement(seq_without_primers)
            seq_with_primers = copy_reverse_complement(seq_with_primers)

        if not MACHINE_LEARNING_MODE:
            output_files['wo_primers'].write(str(seq_without_primers) + "\n")
            output_files['with_primers'].write(str(seq_with_primers) + "\n")

        if MACHINE_LEARNING_MODE or DEBUG_MODE:
            return compare_seq_to_dataset(compare_params, seq_without_primers, debug_file)
        return MatchResult.FOUND_VALID_MATCH

    return MatchResult.NO_MATCH_FOUND


def search_single_primer(valid_alignments, valid_alignments_rc, sequence, output_files, compare_params, debug_file):
    """
    Searches for a valid match with either the front or back primer, in both regular and reverse orientations.

    :param valid_alignments: Dictionary of valid alignments for the regular orientation (keys: 'front', 'back').
    :param valid_alignments_rc: Dictionary of valid alignments for the reverse complement orientation.
    :param sequence: The sequence being analyzed.
    :param output_files: File handles for writing sequences (with and without primers).
    :param compare_params: Comparison parameters (aligner, minimum score).
    :param debug_file: File handle for logging debug information.
    :return: MatchResult indicating the result
    """
    # Search for single front primer - regular
    search_res_front = search_single_front_primer_helper(valid_alignments['front'], sequence, output_files,
                                                         compare_params, debug_file, False)
    if search_res_front != MatchResult.NO_MATCH_FOUND:
        return search_res_front

    # Search for single front primer - reversed
    search_res_front_rc = search_single_front_primer_helper(valid_alignments_rc['front'], sequence, output_files,
                                                            compare_params, debug_file, True)
    if search_res_front_rc != MatchResult.NO_MATCH_FOUND:
        return search_res_front_rc

    # Search for single back primer - regular
    search_res_back = search_single_back_primer_helper(valid_alignments['back'], sequence, output_files, compare_params,
                                                       debug_file, False)
    if search_res_back != MatchResult.NO_MATCH_FOUND:
        return search_res_back

    # Search for single back primer - reversed
    search_res_back_rc = search_single_back_primer_helper(valid_alignments_rc['back'], sequence, output_files,
                                                          compare_params, debug_file, True)
    if search_res_back_rc != MatchResult.NO_MATCH_FOUND:
        return search_res_back_rc

    return MatchResult.NO_MATCH_FOUND

# ******************* Exact match *******************


def search_exact_match(sequence, primers, output_files):
    """
    Search for exact matches of both front and back primers within the sequence.

    :param sequence: The full DNA sequence to search.
    :param primers: The primers to search for. Could be reverse complement.
    :param output_files: Dictionary containing output_wo_primers and output_with_primers file objects.
    :return: Boolean indicating whether a valid match was found.
    """
    front_primer, back_primer = primers.values()
    f_end = sequence.find(front_primer)
    b_start = sequence.find(back_primer)
    f_start = f_end - len(front_primer)  # has to be in bounds, because we found exact matches for the primers
    b_end = b_start + len(back_primer)

    if f_end != -1 and b_start != -1:
        seq_without_primers = sequence[f_end:b_start]
        if LIB_LENGTH - 5 <= len(seq_without_primers) <= LIB_LENGTH + 5:
            seq_with_primers = sequence[f_start:b_end]
            if not MACHINE_LEARNING_MODE:
                output_files['wo_primers'].write(f"{seq_without_primers}\n")
                output_files['with_primers'].write(f"{seq_with_primers}\n")
            return True
    return False

# ******************* Main search function *******************


def search_seq_and_write(aligner, min_valid_score, compare_params, sequence, output_files, debug_file):
    """
    Search for sequences and write the results to the output files.

    :param aligner: PairwiseAligner object for alignment.
    :param min_valid_score: Minimum valid alignment score.
    :param compare_params: Dictionary containing aligner and min_valid_score for dataset comparison.
    :param sequence: The full DNA sequence to search.
    :param output_files: Dictionary containing output_wo_primers and output_with_primers file objects.
    :param debug_file: File object for writing debug information.
    :return: MatchResult indicating whether a valid match was found.
    """
    # Reverse complement the primers
    PRIMERS_RC = {
        'front': Seq(copy_reverse_complement(PRIMERS['back'])),
        'back': Seq(copy_reverse_complement(PRIMERS['front'])),
    }

    # Do an initial exact match search - without using an aligner(and without learning hyperparameters)
    if (search_exact_match(sequence, PRIMERS, output_files) or
            search_exact_match(sequence, PRIMERS_RC, output_files)):
        return MatchResult.FOUND_VALID_MATCH

    # Get valid alignments for both regular and reversed primers
    valid_alignments = get_valid_alignments(aligner, min_valid_score, sequence, PRIMERS)
    valid_alignments_rc = get_valid_alignments(aligner, min_valid_score, sequence, PRIMERS_RC)

    # Search for both primers
    search_res = search_both_primers(valid_alignments, valid_alignments_rc, sequence, output_files,
                                     compare_params, debug_file)
    if search_res != MatchResult.NO_MATCH_FOUND:
        return search_res

    # Search for single primer
    search_res = search_single_primer(valid_alignments, valid_alignments_rc, sequence, output_files,
                                      compare_params, debug_file)
    if search_res != MatchResult.NO_MATCH_FOUND:
        return search_res

    return MatchResult.NO_MATCH_FOUND
