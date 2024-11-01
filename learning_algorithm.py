from hyperopt import fmin, tpe, Trials
import os
from config import DEFAULT_SEARCH_SPACE, EVAL_NUM, DIRECTORIES
from learning_helper_functions import score_comparison, plot_results, \
    plot_pie_chart, plot_eval_pie_chart
from naive_dna_sequence_filter import naive_filtering_percent_of_dir
from functools import partial


# ============================= Helper Functions =============================

def objective(params, naive_percent_filtered, eval_counter):
    """
    Objective function for the optimization algorithm. Runs the sequence filtering algorithm with the given
    hyperparameters

    :param params: Dictionary of hyperparameters (match_score, mismatch_score, open_gap_score, extend_gap_score)
    :param naive_percent_filtered: Percent of sequences filtered out by the naive algorithm
    :param eval_counter: Counter for the current evaluation eval. Starts its count from 1.
    :return: Dictionary containing the loss, percent_filtered, and status of the evaluation
    """
    scores = {
        "match_score": params['match_score'],
        "mismatch_score": params['mismatch_score'],
        "open_gap_score": params['open_gap_score'],
        "extend_gap_score": params['extend_gap_score']
    }

    eval_dir = os.path.join(DIRECTORIES["learning_dir"], "eval_plots", f"eval_{eval_counter-1}")
    os.makedirs(eval_dir, exist_ok=True)

    # Filter and get statistics
    from main_sequence_filtering import run_filtering
    eval_stats, file_stats_list = run_filtering(scores=scores)
    eval_stats.loss = score_comparison(naive_percent_filtered, eval_stats)

    # Plot individual file statistics pie charts
    for i, file_stats in enumerate(file_stats_list):
        os.makedirs(eval_dir, exist_ok=True)
        plot_pie_chart(file_stats, eval_dir)

    # Plot overall statistics pie chart for the current eval
    plot_eval_pie_chart(eval_stats, eval_dir, eval_counter-1)

    return {'loss': eval_stats.loss, 'percent_filtered': eval_stats.percent_filtered, 'status': 'ok'}


def plot_and_generate_results(trials_list, naive_percent_filtered, alg_dir_name):
    """
    Plots results from trials and generates pie charts for algorithm analysis.

    :param trials_list: List of trial dictionaries containing parameter values and results
    :param naive_percent_filtered: Naive filtering percentage of the input directory
    :param alg_dir_name: Subdirectory for specific algorithm plots
    """
    class MockTrials:
        def __init__(self, trials):
            self.trials = trials

    mock_trials_obj = MockTrials(trials_list)
    output_dir = os.path.join(DIRECTORIES["learning_dir"], alg_dir_name)
    plot_results(mock_trials_obj, naive_percent_filtered, output_dir=output_dir)


# ============================= Optimization Algorithm Functions =============================


def run_tpe(naive_percent_filtered):
    """
    Runs Tree-structured Parzen Estimator (TPE) optimization using Hyperopt to find the best hyperparameters.

    :param naive_percent_filtered: Naive filtering percentage of the input directory
    :return: Tuple of best hyperparameters and corresponding score
    """
    print("\n" + "="*64)
    print(f"=============== Running TPE with Hyperopt... ===============")
    print("="*64 + "\n")

    # Optimize using TPE
    def objective_with_args(params, naive_percent_filtered):
        eval_counter = len(trials.trials)
        return objective(params, naive_percent_filtered, eval_counter)

    objective_partial = partial(objective_with_args, naive_percent_filtered=naive_percent_filtered)

    trials = Trials()
    best = fmin(fn=objective_partial, space=DEFAULT_SEARCH_SPACE, algo=tpe.suggest, max_evals=EVAL_NUM, trials=trials)

    # Plot and generate charts
    plot_and_generate_results(trials.trials, naive_percent_filtered, "tpe_plots")

    # Extract the best score from the trials
    best_score = min(trial['result']['loss'] for trial in trials.trials)
    return best, best_score


# ============================= Useful Functions =============================

def run_learning_algorithm():
    """
    Runs the learning algorithm to find the best hyperparameters for the sequence filtering algorithm.

    :return: Dictionary containing the best parameters and scores for each algorithm
    """
    results = {}

    naive_percent_filtered = naive_filtering_percent_of_dir(DIRECTORIES['input_dir'])

    # Run TPE with Hyperopt
    best_params_tpe, best_score_tpe = run_tpe(naive_percent_filtered)
    results['TPE (Hyperopt)'] = {'params': best_params_tpe, 'score': best_score_tpe}

    # Display Results
    print("\n--- Results ---")
    for algo, result in results.items():
        print(f"{algo}: Best Score = {result['score']}, Best Params = {result['params']}")

    return results
