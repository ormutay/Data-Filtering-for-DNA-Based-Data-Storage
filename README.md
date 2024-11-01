# DNA Sequence Filtering and Classification Project

This project applies sequence alignment algorithms to filter DNA sequences. It includes **Filtering**, **Debug**, and **Learning** modes, with configurations accessible in a centralized configuration file for easy adjustments.

## Configuration (`config.py`)

All key parameters for this project are in `config.py`. Adjustments to these constants control the filtering behavior, alignment parameters, file paths, and hyperparameter search spaces used during execution.

- **File Paths**: Set the directories for input sequences, output results and more.
- **Library Length and Primers**: Defines expected library lengths and primer sequences to align against.
- **Alignment Scores**: Customize match, mismatch, and gap scores, impacting how sequence alignments are scored.
- **Learning Settings**: Defines the range and search space for hyperparameters used in Learning mode. Set the learning mode on/off.
- **Debug Settings**: Set the debug mode on/off.

## Files Overview

The main functions in this project are contained across several fiels, each targeting specific operations within the filtering pipeline:

1. **main_sequence_filtering.py**: Orchestrates the sequence filtering process by calling appropriate functions based on the chosen mode.
2. **naive_dna_sequence_filter.py**: Provides basic filtering operations based on primer matching.
3. **helper_functions.py**: Handles file reading, sequence manipulation, and helper functions for alignment.
4. **search_functions.py**: Implements functions to search for primers and filter sequences based on specified criteria.
5. **learning_helper_functions.py**: Supports Learning mode by managing evaluation statistics, parameter extraction, and plots creation.
6. **learning_algorithm.py**: Runs hyperparameter optimization using Tree-structured Parzen Estimator (TPE) to refine scoring parameters.

## Usage Modes

- **Dependencies**: Ensure you have the needed dependencies installed to run the project: `pip install -r requirements.txt`
- **Command**:
  - Make sure you set the config file according to your needs and the details in this README beforehand.
  - To run this project you can use the command: `python main_sequence_filtering.py`
  - You can also use the following command to ensure your run will not stop (due to a server disconnection, ect.): `nohup python3 main_sequence_filtering.py &; disown;`
> [!WARNING]
>  After using the `disown` command, the run cannot be stopped from the terminal. Ensure you are ready for the process to complete fully before initiating


### 1. **Filtering Mode**

This mode executes the core filtering process using the parameters specified in `config.py`. It filters sequences based on primer matching and alignment scores.

- **How to activate**: Under the config file the flag `MACHINE_LEARNING_MODE` needs to be set to False.
- **Outcome**: Generates output files with filtered sequences, distinguishing those that match the criteria from those that don’t.

### 2. **Debug Mode**

In Debug mode, the program provides detailed logs and intermediate results for each sequence processed. This mode helps identify errors in filtering and alignment configurations.

- **How to activate**: Under the config file the flag `DEBUG_MODE` needs to be set to True.
- **Outcome**: Creates detailed logs for each sequence, indicating filtering decisions, alignment scores, and any mismatches. When running with `MACHINE_LEARNING_MODE` on, only the files from the last eval will remain.

### 3. **Learning Mode**

This mode optimizes the filtering process by exploring different hyperparameter configurations to find the best match, mismatch, and gap scores. It uses a Tree-structured Parzen Estimator (TPE) to search the hyperparameter space defined in `config.py`.

- **How to activate**: Under the config file the flag `MACHINE_LEARNING_MODE` needs to be set to True.
- **Outcome**:
  - **TPE Plots**: Saves charts of precent filtered and loss over iterations as well as parameter distributions over the evaluations.
  - **Eval Plots**: Saves pie charts for filtering/classification percentages for each file and evaluation.
  - **Excel**: Saves an Excel file with the hyperparameters and loss for each evaluation.
