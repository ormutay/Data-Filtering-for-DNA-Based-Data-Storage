from Bio.Align import PairwiseAligner
from config import DATA_SET_PATH


def read_FASTQ_file(fastq_file):
    """
     Reads a FASTQ file and extracts sequences and their corresponding ASCII scores.

     :param fastq_file: Path to the FASTQ file.
     :return: A dictionary with sequences as keys and ASCII average scores as values,
              and the length of the sequence dictionary.
     """
    seq_data = []
    seq, ascii_scores = [], []

    # Read the FASTQ file and filter out header lines (the first line) and separators
    with open(fastq_file, 'r') as FASTQ:
        lines = FASTQ.readlines()[1:]
        for i, line in enumerate(lines):
            if i % 2 == 0:
                seq_data.append(line.strip())

    # Parse the sequences and ASCII scores
    empty_seq = False
    for i, data in enumerate(seq_data):
        if i % 2 == 0:
            if len(data) == 0:
                empty_seq = True
                continue
            seq.append(data)
        else:
            if empty_seq:
                empty_seq = False
                continue
            scores = [ord(char) - 33 for char in data]
            ascii_scores.append(sum(scores) / len(scores))

    # Create a dictionary mapping sequence to score
    seq_dict = {s: score for s, score in zip(seq, ascii_scores)}
    return seq_dict, len(seq_dict)


def copy_reverse_complement(sequence):
    """
    Generates the reverse complement of a given DNA sequence.

    :param sequence: DNA sequence (e.g., 'ATCG')
    :return: The reverse complement of the input DNA sequence.
    """
    sequence = sequence.rstrip()  # Remove any trailing whitespace
    complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}  # DNA complement rules
    complement_sequence = ''.join([complement_map[base] for base in sequence])  # Replace bases with complements
    return complement_sequence[::-1]  # Reverse the complemented sequence to get the final result


def get_alignment_positions(alignment):
    """
    Extracts the start and end positions from a given alignment.

    :param alignment: PairwiseAlignment object from BioPython.
    :return: Tuple (start_pos, end_pos) representing the alignment start and end positions.
    """
    start_pos = alignment.aligned[0][0][0]
    end_pos = alignment.aligned[0][-1][1]
    return start_pos, end_pos


def is_overlap(start, end, regions):
    """
    Checks if a given range (start, end) overlaps with any regions in a list of regions.

    :param start: Start position of the range.
    :param end: End position of the range.
    :param regions: List of regions, where each region is represented as (start, end).
    :return: True if overlap is found, otherwise False.
    """
    for region in regions:
        if not (end < region[0] or start > region[1]):
            return True
    return False


def percentage_to_min_score(aligner, wanted_match_percentage, length):
    """
    Calculates the minimum valid alignment score based on a desired match percentage.

    :param aligner: PairwiseAligner object with match/mismatch/gap scores set.
    :param wanted_match_percentage: Desired match percentage (e.g., 0.8 for 80% match).
    :param length: The length of the sequence.
    :return: The minimum valid score for alignment.
    """
    max_score = aligner.match_score * length
    avg_penalty = (aligner.mismatch_score + aligner.open_gap_score + aligner.extend_gap_score) / 3
    penalty_score = avg_penalty * (1 - wanted_match_percentage) * length
    match_score = max_score * wanted_match_percentage
    return match_score + penalty_score


def set_aligner_parameters(scores, mode, wanted_match_percentage, length):
    """
    Configures the pairwise aligner with the provided scoring parameters.

    :param scores: Dictionary of match, mismatch, open gap, and extend gap scores.
    :param mode: Alignment mode ('global' or 'local').
    :param wanted_match_percentage: Desired match percentage.
    :param length: Length of the sequence to align.
    :return: Configured PairwiseAligner object and the minimum valid score.
    """
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = scores['match_score']
    aligner.mismatch_score = scores['mismatch_score']
    aligner.open_gap_score = scores['open_gap_score']
    aligner.extend_gap_score = scores['extend_gap_score']

    min_valid_score = percentage_to_min_score(aligner, wanted_match_percentage, length)

    return aligner, min_valid_score


def read_sequences_from_data_set():
    """
    Reads sequences from a text file, stripping whitespace from each line.

    :return: A list of sequences.
    """
    with open(DATA_SET_PATH, 'r') as file:
        sequences = [line.strip() for line in file.readlines()]
    return sequences


def get_valid_alignments(aligner, min_valid_score, sequence, primers):
    """
    Get valid alignments for both front and back primers within a sequence.

    :param aligner: PairwiseAligner object to perform the alignment.
    :param min_valid_score: Minimum score for valid alignments.
    :param sequence: The full DNA sequence to search.
    :param primers: The primers. Could be reverse complement.
    :return: dictionary of valid front and back alignments.
    """
    valid_alignments_front = []
    valid_alignments_back = []

    alignments_front = aligner.align(sequence, primers['front'])
    alignments_back = aligner.align(sequence, primers['back'])

    for f_align in alignments_front:
        if f_align.score >= min_valid_score:
            valid_alignments_front.append(f_align)

    for b_align in alignments_back:
        if b_align.score >= min_valid_score:
            valid_alignments_back.append(b_align)

    valid_alignments = {'front': valid_alignments_front, 'back': valid_alignments_back}
    return valid_alignments
