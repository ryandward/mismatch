import pandas as pd
import numpy as np
import argparse
import logging
import sys

def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

def calculate_y_pred(original, variant, gc_weight, params):
    y_pred = params['intercept']

    for pos, (orig_base, var_base) in enumerate(zip(original, variant)):
        if orig_base != var_base:
            y_pred += params[f'{pos}']
            y_pred += params[f'{orig_base}{var_base}']

    y_pred += gc_weight * gc_content(original)
    return y_pred

def read_parameters(file_path):
    try:
        df = pd.read_csv(file_path)
        params = df.set_index("feature")["weight"].to_dict()
        return params
    except Exception as e:
        logging.error(f"Error reading parameters file: {e}")
        sys.exit(1)

def generate_header():
    header = ["original", "variant", "change_description", "y_pred"]
    print("\t".join(header))

def generate_mismatches(spacers, min_score, max_score, step, parameters):
    nucleotides = ["A", "C", "G", "T"]

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
            closest_score = None
            closest_mismatch = None
            for mismatch, mismatch_score in mismatches:
                if closest_score is None or abs(mismatch_score - score) < abs(closest_score - score):
                    if mismatch not in [x[0] for x in mismatch_list]:
                        closest_score = mismatch_score
                        closest_mismatch = mismatch
            mismatch_list.append((closest_mismatch, closest_score))

        for i, (mismatch, score) in enumerate(mismatch_list):
            target_mismatch = spacer[:mismatch[0]] + mismatch[1] + spacer[mismatch[0]+1:]
            change_description = f"{spacer[mismatch[0]]}{mismatch[0]+1}{mismatch[1]}"
            print(f"{spacer}\t{target_mismatch}\t{change_description}\t{score:.4f}")


def main(args):
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=getattr(logging, args.verbosity.upper()))

    logging.info(f"Reading parameters from {args.parameters_file}...")
    params = read_parameters(args.parameters_file)

    if args.mode == "mismatches":
        with open(args.spacers_file, "r") as f:
            spacers = [line.strip() for line in f.readlines() if line.strip()]

        generate_header()
        generate_mismatches(spacers, args.min, args.max, args.step, params)
    elif args.mode == "recalculate":
        try:
            data = pd.read_csv(args.input_data_file, sep="\t")
        except Exception as e:
            logging.error(f"Error reading input data file: {e}")
            sys.exit(1)

        if not {'original', 'variant'}.issubset(data.columns):
            logging.error("Input data file must have 'original' and 'variant' columns.")
            sys.exit(1)

        logging.info("Calculating y_pred for each row...")
        new_y_pred_column_name = 'y_pred_new' if 'y_pred' in data.columns else 'y_pred'
        data[new_y_pred_column_name] = data.apply(lambda row: calculate_y_pred(row['original'], row['variant'], params['GC_content'], params), axis=1)

        logging.info("Displaying results:")
        print(data.to_string(index=False, col_space=0, justify='left'))

        logging.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mismatches for a list of spacers and/or recalculate y_pred.")
    parser.add_argument("mode", choices=["mismatches", "recalculate"], help="Choose the mode of operation: 'mismatches' to generate mismatches, 'recalculate' to recalculate y_pred.")
    parser.add_argument("--spacers_file", help="Path to the file containing original spacers (one per line) (required for mismatches mode).")
    parser.add_argument("--input_data_file", help="Path to the input data file (TSV format) (required for recalculate mode).")
    parser.add_argument("--parameters_file", required=True, help="Path to the parameters file (CSV format).")
    parser.add_argument("--verbosity", choices=["debug", "info", "warning", "error", "critical"], default="info", help="Set the logging verbosity level (default: info).")
    args = parser.parse_args()

    if args.mode == "mismatches" and args.spacers_file is None:
        parser.error("The --spacers_file option is required for mismatches mode.")
    elif args.mode == "recalculate" and args.input_data_file is None:
        parser.error("The --input_data_file option is required for recalculate mode.")

    main(args)
