import matplotlib.pyplot as plt
import numpy as np
import os
from config import DIRECTORIES
import pandas as pd


# ============================= Helper Functions =============================

def save_hyperparameter_eval_mapping_to_excel(params, losses):
    """
    Save hyperparameter eval mapping to an Excel file.

    :param params: Hyperparameter values over evals.
    :param losses: Loss values over evals.
    """
    match_scores, mismatch_scores, open_gap_scores, extend_gap_scores = get_scores_lists(params)
    formatted_losses = [round(loss, 2) for loss in losses]

    eval_data = {
        'Eval': range(len(match_scores)),
        'match_score': match_scores,
        'mismatch_score': mismatch_scores,
        'open_gap_score': open_gap_scores,
        'extend_gap_score': extend_gap_scores,
        'Loss': formatted_losses
    }

    df = pd.DataFrame(eval_data)
    df.to_excel(os.path.join(DIRECTORIES["learning_dir"], "eval_hyperparameter_mapping.xlsx"), index=False)


def get_trials_data(trials):
    """
    Extract losses, percent filtered, and parameters from trials.

    :param trials: Trials object containing trial data.
    :return: Tuple (losses, percent_filtered, params).
    """
    losses = [x['result']['loss'] for x in trials.trials]
    percent_filtered = [x['result']['percent_filtered'] for x in trials.trials]
    params = [x['misc']['vals'] for x in trials.trials]

    return losses, percent_filtered, params


def get_scores_lists(params):
    """
    Extract lists of match, mismatch, open gap, and extend gap scores from parameters.

    :param params: Hyperparameter values over evals.
    :return: Tuple of match_scores, mismatch_scores, open_gap_scores, extend_gap_scores.
    """
    match_scores = [p['match_score'][0] for p in params]
    mismatch_scores = [p['mismatch_score'][0] for p in params]
    open_gap_scores = [p['open_gap_score'][0] for p in params]
    extend_gap_scores = [p['extend_gap_score'][0] for p in params]

    return match_scores, mismatch_scores, open_gap_scores, extend_gap_scores


def calculate_pie_sizes(filtered_out, correct_classified, incorrect_classified):
    """
    Calculate the sizes for pie chart sections based on classification and filtering results.

    :param filtered_out: Number of sequences filtered out.
    :param correct_classified: Number of correctly classified sequences.
    :param incorrect_classified: Number of incorrectly classified sequences.
    :return: List of percentages for the pie chart.
    """
    sizes = [filtered_out, correct_classified, incorrect_classified]
    total = sum(sizes)

    if total == 0:
        return []

    return [(s / total) * 100 for s in sizes]


# ============================= Plot Functions =============================

def plot_parameter_distributions(match_scores, mismatch_scores, open_gap_scores, extend_gap_scores, output_dir):
    """
    Plot the distribution of hyperparameter values for match, mismatch, open gap, and extend gap scores.

    :param match_scores: List of match scores.
    :param mismatch_scores: List of mismatch scores.
    :param open_gap_scores: List of open gap scores.
    :param extend_gap_scores: List of extend gap scores.
    :param output_dir: Directory where the plot will be saved.
    """
    # Create a 2x2 grid of subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Define score data, titles, and colors for each plot
    score_data = [
        (match_scores, 'Distribution of Match Scores', 'blue'),
        (mismatch_scores, 'Distribution of Mismatch Scores', 'green'),
        (open_gap_scores, 'Distribution of Open Gap Scores', 'orange'),
        (extend_gap_scores, 'Distribution of Extend Gap Scores', 'red')
    ]

    # Plot histograms for each set of scores
    for ax, (scores, title, color) in zip(axes.flat, score_data):
        ax.hist(scores, bins=20, color=color, edgecolor='black')
        ax.set_title(title, fontsize=18)
        ax.set_xlabel('Value', fontsize=16)
        ax.set_ylabel('Frequency', fontsize=16)
        ax.grid(True)

    # Adjust layout and add a super title
    plt.tight_layout()
    plt.suptitle('Parameter Distributions', fontsize=24, y=0.98)
    plt.subplots_adjust(top=0.9)  # Adjust top margin to fit the title
    plt.savefig(os.path.join(output_dir, "parameter_distributions.png"))
    plt.close()


def plot_percent_filtered_and_loss_over_evals(max_evals, percent_filtered, naive_percent_filtered, losses, output_dir):
    """
    Plot the percent of sequences filtered and loss over evals on two y-axes.

    :param max_evals: Maximum number of evals.
    :param percent_filtered: List of percent filtered values.
    :param naive_percent_filtered: Percent of sequences filtered by the naive algorithm.
    :param losses: List of loss values.
    :param output_dir: Directory where the plot will be saved.
    """
    ax1_color = 'tab:blue'
    ax2_color = 'tab:red'

    fig, ax1 = plt.subplots(figsize=(18, 8))

    ax1.set_xlabel('Eval', fontsize=20, labelpad=20)
    ax1.set_ylim(0, 100)
    ax1.set_ylabel('Percent\nFiltered', color=ax1_color, labelpad=50, rotation=0, fontsize=20)
    ax1.tick_params(axis='y', labelcolor=ax1_color)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    percent_filtered_line = ax1.plot(range(max_evals), percent_filtered, color=ax1_color, marker='o', linestyle='-',
                                     label='Percent Filtered')

    naive_filter_line = ax1.axhline(y=naive_percent_filtered, color='green', linestyle='--',
                                    label=f'Naive Algorithm: {naive_percent_filtered:.2f}%')

    ax2 = ax1.twinx()
    ax2.set_ylabel('Loss', color=ax2_color, fontsize=20, labelpad=30, rotation=0)
    ax2.tick_params(axis='y', labelcolor=ax2_color)
    ax2.tick_params(axis='both', which='major', labelsize=16)
    loss_line = ax2.plot(range(max_evals), losses, color=ax2_color, marker='x', linestyle='-', label='Loss')

    ax1.set_xticks(np.arange(0, max_evals, step=5))
    ax1.set_xticklabels(np.arange(0, max_evals, step=5), fontsize=16)

    # Create a legend for the plot
    lines = [percent_filtered_line[0], loss_line[0], naive_filter_line]
    labels = ['Percent Filtered', 'Loss', f'Naive Algorithm\'s\n'
                                          f'Percent Filtered: {naive_percent_filtered:.2f}%']
    fig.legend(lines, labels, loc="upper right", bbox_to_anchor=(0.98, 0.96), fontsize=18)

    plt.title('Percent Filtered and Loss over Evals', fontsize=26, pad=20)
    ax1.grid(True, which='both', axis='both')
    plt.subplots_adjust(top=0.85, right=0.7)
    plt.savefig(os.path.join(output_dir, "percent_filtered_and_loss_over_evals.png"))
    plt.close()


# ============================= Useful Functions =============================

def plot_results(trials, naive_percent_filtered, output_dir="plots"):
    """
    Plot results including loss over evals, parameter distributions, and more.

    :param trials: Trials object containing data on each eval.
    :param naive_percent_filtered: Percent of sequences filtered by the naive algorithm.
    :param output_dir: Directory to save the plots.
    """
    os.makedirs(output_dir, exist_ok=True)

    losses, percent_filtered, params = get_trials_data(trials)
    match_scores, mismatch_scores, open_gap_scores, extend_gap_scores = get_scores_lists(params)
    max_evals = len(losses)

    # Save hyperparameter-eval mapping
    save_hyperparameter_eval_mapping_to_excel(params, losses)

    # Generate plots
    plot_parameter_distributions(match_scores, mismatch_scores, open_gap_scores, extend_gap_scores, output_dir)
    plot_percent_filtered_and_loss_over_evals(max_evals, percent_filtered, naive_percent_filtered, losses, output_dir)


def plot_pie_chart(stats, output_dir):
    """
    Plot a pie chart for file-specific filtering and classification results.

    :param stats: Stats object containing filtering and classification results.
    :param output_dir: Directory to save the plot.
    """
    file_name = stats.file_name
    sizes = calculate_pie_sizes(stats.filtered_out_count, stats.valid_seq_count, stats.invalid_seq_count)

    plt.figure(figsize=(8, 8))
    wedges, _, autotexts = plt.pie(sizes, autopct='%1.1f%%', startangle=140, colors=['gold', 'skyblue', 'lightcoral'])
    plt.axis('equal')

    labels = ['Filtered Out', 'Classified as Correct', 'Classified as Incorrect']
    plt.legend(wedges, labels, loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

    plt.figtext(0.5, 0.87, "File-Specific Filtering and Classification Distribution", ha="center", fontsize=17)
    plt.figtext(0.5, 0.12, file_name, ha="center", fontsize=12)

    output_file = os.path.join(output_dir, f"file_specific_pie_chart_{file_name}.png")
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()


def plot_eval_pie_chart(eval_stats, eval_dir, eval_counter):
    """
    Plot a eval pie chart for filtering and classification results across all files.

    :param eval_stats: EvalStats object containing filtering and classification results.
    :param eval_dir: Directory to save the plot.
    :param eval_counter: Counter for the evaluation number.
    """
    sizes = calculate_pie_sizes(eval_stats.filtered_out_count, eval_stats.valid_seq_count, eval_stats.invalid_seq_count)

    plt.figure(figsize=(8, 8))
    wedges, _, autotexts = plt.pie(sizes, autopct='%1.1f%%', startangle=140, colors=['gold', 'skyblue', 'lightcoral'])
    plt.axis('equal')

    labels = ['Filtered Out', 'Classified as Correct', 'Classified as Incorrect']
    plt.legend(wedges, labels, loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

    plt.figtext(0.5, 0.87, f"Eval {eval_counter} Filtering and Classification Distribution", ha="center", fontsize=17)

    output_file = os.path.join(eval_dir, f"eval_{eval_counter}_pie_chart.png")
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()


# ============================= Utility Functions =============================

def score_comparison(naive_percent_filtered, eval_stats):
    """
    Compare the score of filtered sequences and apply a penalty if the filtering
    performed worse than the naive algorithm.

    :param naive_percent_filtered: Percent of sequences filtered by the naive algorithm.
    :param eval_stats: Dictionary containing statistics from the filtering process.
    :return: Computed score based on invalid sequence counts and penalty.
    """
    score = 0

    # Apply a penalty if the filtering performed worse than naive filtering
    if naive_percent_filtered < eval_stats.percent_filtered:
        score += 100

    assert(eval_stats.seq_count > 0)
    score += (eval_stats.invalid_seq_count / eval_stats.seq_count) * 100

    return score
