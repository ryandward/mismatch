import pandas as pd
import numpy as np
import argparse
import logging
import sys


def gc_content(seq):
    """Calculate GC content of a DNA sequence."""
    return (seq.count('G') + seq.count('C')) / len(seq)


def calculate_y_pred(original, variant, gc_weight, params):
    """Calculate y_pred based on original and variant sequences."""
    y_pred = params['intercept']

    for pos, (orig_base, var_base) in enumerate(zip(original, variant)):
        if orig_base != var_base:
            y_pred += params[f'{pos}']
            y_pred += params[f'{orig_base}{var_base}']

    y_pred += gc_weight * gc_content(original)
    return y_pred


def read_parameters(file_path):
    """Read parameters from a CSV file."""
    try:
        df = pd.read_csv(file_path)
        params = df.set_index('feature')['weight'].to_dict()
        return params
    except Exception as e:
        logging.error(f'Error reading parameters file: {e}')
        sys.exit(1)


def generate_header():
    """Generate and print header row for output."""
    header = ['original', 'variant', 'change_description', 'y_pred']
    print('\t'.join(header))


def find_closest_mismatch(score, mismatches, mismatch_list):
    """Find the closest mismatch based on the desired score."""
    closest_score = None
    closest_mismatch = None
    for mismatch, mismatch_score in mismatches:
        if closest_score is None or abs(mismatch_score - score) < abs(closest_score - score):
            if mismatch not in [x[0] for x in mismatch_list]:
                closest_score = mismatch_score
                closest_mismatch = mismatch
    return closest_mismatch, closest_score


def print_mismatches(mismatch_list, spacer):
    """Print the mismatches in a formatted manner."""
    for i, (mismatch, score) in enumerate(mismatch_list):
        if mismatch is not None:
            target_mismatch = spacer[:mismatch[0]] + mismatch[1] + spacer[mismatch[0]+1:]
            change_description = f"{spacer[mismatch[0]]}{mismatch[0]+1}{mismatch[1]}"
            row = [spacer, target_mismatch, change_description, f"{score:.4f}"]
            print('\t'.join(row))


def generate_mismatches(spacers, min_score, max_score, step, parameters):
    """Generate mismatches for a list of spacers based on desired scores."""
    nucleotides = ['A', 'C', 'G', 'T']

    for spacer in spacers:
        desired_scores = np.arange(min_score, max_score + step, step)

        mismatches = []
        for pos in range(len(spacer)):
            for nt in nucleotides:
                if nt == spacer[pos]:
                    continue
                target_mismatch = spacer[:pos] + nt + spacer[pos+1:]
                y_pred = calculate_y_pred(spacer, target_mismatch, parameters['GC_content'], parameters)
                mismatches.append(((pos, nt), y_pred))

        mismatch_list = []
        for score in desired_scores:
            closest_mismatch, closest_score = find_closest_mismatch(score, mismatches, mismatch_list)
            if closest_mismatch is not None:
                mismatch_list.append((closest_mismatch, closest_score))

        print_mismatches(mismatch_list, spacer)


def main(args):
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=getattr(logging, args.verbosity.upper()))

    logging.info(f'Reading parameters from {args.parameters_file}...')
    params = read_parameters(args.parameters_file)

    if args.mode == 'mismatches':
        with open(args.spacers_file, 'r') as f:
            spacers = [line.strip() for line in f.readlines() if line.strip()]

        generate_header()
        generate_mismatches(spacers, args.min, args.max, args.step, params)
    elif args.mode == 'recalculate':
        try:
            data = pd.read_csv(args.existing_mismatches, sep='\t')
        except Exception as e:
            logging.error(f'Error reading input data file: {e}')
            sys.exit(1)

        if not {'original', 'variant'}.issubset(data.columns):
            logging.error("Input data file must have 'original' and 'variant' columns.")
            sys.exit(1)

        logging.info('Calculating y_pred for each row...')
        new_y_pred_column_name = 'y_pred_new' if 'y_pred' in data.columns else 'y_pred'
        data[new_y_pred_column_name] = data.apply(lambda row: f"{calculate_y_pred(row['original'], row['variant'], params['GC_content'], params):.4f}", axis=1)

        logging.info('Displaying results:')
        print(data.to_csv(sep='\t', index=False))

        logging.info('Done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate mismatches for a list of spacers and/or recalculate y_pred.')
    parser.add_argument('mode', choices=['mismatches', 'recalculate'], help="Choose the mode of operation: 'mismatches' to generate mismatches, 'recalculate' to recalculate y_pred.")
    parser.add_argument('--spacers_file', help='Path to the file containing original spacers (one per line) (required for mismatches mode).')
    parser.add_argument('--existing_mismatches', help='Path to the input data file (TSV format) (required for recalculate mode).')
    parser.add_argument('--parameters_file', required=True, help='Path to the parameters file (CSV format).')
    parser.add_argument('--verbosity', choices=['debug', 'info', 'warning', 'error', 'critical'], default='info', help='Set the logging verbosity level (default: info).')
    parser.add_argument('--min', type=float, default=0, help='Minimum desired efficacy (default: 0) (required for mismatches mode).')
    parser.add_argument('--max', type=float, default=1, help='Maximum desired efficacy (default: 1) (required for mismatches mode).')
    parser.add_argument('--step', type=float, default=0.1, help='Step between desired efficacies (default: 0.1) (required for mismatches mode).')
    args = parser.parse_args()

    if args.mode == 'mismatches' and args.spacers_file is None:
        parser.error('The --spacers_file option is required for mismatches mode.')
    elif args.mode == 'recalculate' and args.existing_mismatches is None:
        parser.error('The --existing_mismatches option is required for recalculate mode.')

    main(args)
